#!/usr/bin/env julia

include(joinpath(@__DIR__, "generate_oracle.jl"))

const STATIONARITY_PRECISIONS_BITS = Int[256, 512]
const STATIONARITY_SETTINGS = [
    (label = "baseline",   tolerance = 1.0e-16, max_iters = 4096),
    (label = "tighter_18", tolerance = 1.0e-18, max_iters = 8192),
    (label = "tighter_20", tolerance = 1.0e-20, max_iters = 8192),
    (label = "tighter_24", tolerance = 1.0e-24, max_iters = 16384),
    (label = "tighter_28", tolerance = 1.0e-28, max_iters = 32768),
    (label = "tighter_32", tolerance = 1.0e-32, max_iters = 65536),
]
const DERIVATIVE_RELATIVE_H = "1e-6"
const TARGET_LABEL = "state=4 phase=liquid"

function big_string(value)
    return string(value)
end

function state4_liquid(source, frozen)
    target = nothing
    for (flat_index, oracle_state) in enumerate(frozen["states"])
        String(oracle_state["label"]) == TARGET_LABEL || continue
        source_index = div(flat_index - 1, 2) + 1
        phase = isodd(flat_index) ? "liquid" : "vapor"
        require(phase == "liquid", "state 4 audit matched a non-liquid record")
        source_state = source["states"][source_index]
        target = (
            temperature_k = Float64(source_state["temperature_k"]),
            target_pressure_pa = Float64(source_state["pressure_pa"]),
            rho = Float64(
                source_state["liquid"]["molar_density_mol_per_m3"]),
            x_meoh = Float64(
                source_state["liquid"]["composition"]["MEOH"]),
            frozen_clapeyron_pressure_pa =
                Float64(oracle_state["pressure_pa"]),
            frozen_f_res =
                Float64(oracle_state["f_res_reduced_extensive"]),
        )
        break
    end
    target === nothing &&
        error("state 4 liquid not found in frozen oracle")
    return target
end

function association_center(model, volume, temperature, z)
    x = Clapeyron.assoc_fractions(model, volume, temperature, z)
    x_values = collect(x.v)
    if !(eltype(x_values) == BigFloat && all(isfinite, x_values))
        return Dict(
            "status" => "nonfinite_association_solution",
            "x" => [big_string(value) for value in x_values],
        )
    end

    matrix = Clapeyron.assoc_site_matrix(
        model, volume, temperature, z)
    residual_vector =
        x_values .* (1 .+ matrix * x_values) .- 1
    if !all(isfinite, residual_vector)
        return Dict(
            "status" => "nonfinite_mass_action_residual",
            "x" => [big_string(value) for value in x_values],
        )
    end

    a_res = Clapeyron.eos_res(model, volume, temperature, z)
    if !(a_res isa BigFloat && isfinite(a_res))
        return Dict(
            "status" => "nonfinite_scalar_helmholtz",
            "x" => [big_string(value) for value in x_values],
            "mass_action_max_residual" =>
                big_string(maximum(abs, residual_vector)),
        )
    end

    rt = BigFloat(Clapeyron.Rgas(model)) * temperature
    f_res = a_res / rt

    return Dict(
        "status" => "ok",
        "x" => [big_string(value) for value in x_values],
        "mass_action_max_residual" =>
            big_string(maximum(abs, residual_vector)),
        "mass_action_residual_vector" =>
            [big_string(value) for value in residual_vector],
        "a_res_joule" => big_string(a_res),
        "f_res_reduced_extensive" => big_string(f_res),
    )
end

function derivative_7(values, h)
    return (
        -values[-3] + 9 * values[-2] - 45 * values[-1] +
        45 * values[1] - 9 * values[2] + values[3]
    ) / (60 * h)
end

function derivative_probe(
    model,
    volume,
    temperature,
    z,
    target_pressure,
    internal_pressure)

    relative_h = parse(BigFloat, DERIVATIVE_RELATIVE_H)
    h = relative_h * volume
    require(volume - 3h > 0,
            "stationarity derivative stencil crossed V=0")

    values = Dict{Int,BigFloat}()
    for multiplier in (-3, -2, -1, 1, 2, 3)
        shifted_volume = volume + BigFloat(multiplier) * h
        scalar = Clapeyron.eos_res(
            model, shifted_volume, temperature, z)
        if !(scalar isa BigFloat && isfinite(scalar))
            return Dict(
                "status" => "nonfinite_shifted_scalar",
                "failed_multiplier" => multiplier,
            )
        end
        values[multiplier] = scalar
    end

    derivative = derivative_7(values, h)
    total_moles = sum(z)
    rt = BigFloat(Clapeyron.Rgas(model)) * temperature
    pressure = total_moles * rt / volume - derivative

    if !(isfinite(derivative) && isfinite(pressure) && pressure > 0)
        return Dict("status" => "nonfinite_derivative")
    end

    return Dict(
        "status" => "ok",
        "relative_h" => DERIVATIVE_RELATIVE_H,
        "stencil" => "central_7_o6",
        "dAres_dV_pa" => big_string(derivative),
        "pressure_pa" => big_string(pressure),
        "abs_pressure_delta_vs_target_pa" =>
            Float64(abs(pressure - target_pressure)),
        "abs_pressure_delta_vs_clapeyron_internal_pa" =>
            Float64(abs(pressure - internal_pressure)),
    )
end

function run_setting(state, precision_bits, setting)
    return setprecision(BigFloat, precision_bits) do
        model = build_model(
            assoc_rtol = setting.tolerance,
            assoc_atol = setting.tolerance,
            assoc_max_iters = setting.max_iters,
        )

        require(
            model.assoc_options.rtol == setting.tolerance &&
            model.assoc_options.atol == setting.tolerance &&
            model.assoc_options.max_iters == setting.max_iters,
            "stationarity solver options did not round-trip")

        temperature = BigFloat(state.temperature_k)
        z = BigFloat[state.x_meoh, 1.0 - state.x_meoh]
        total_moles = sum(z)
        volume = total_moles / BigFloat(state.rho)
        target_pressure = BigFloat(state.target_pressure_pa)
        internal_pressure =
            BigFloat(state.frozen_clapeyron_pressure_pa)

        center = association_center(
            model, volume, temperature, z)
        derivative = derivative_probe(
            model,
            volume,
            temperature,
            z,
            target_pressure,
            internal_pressure,
        )

        return Dict(
            "precision_bits" => precision_bits,
            "solver_label" => setting.label,
            "rtol" => setting.tolerance,
            "atol" => setting.tolerance,
            "max_iters" => setting.max_iters,
            "dampingfactor" => model.assoc_options.dampingfactor,
            "implicit_ad" => model.assoc_options.implicit_ad,
            "center" => center,
            "derivative" => derivative,
        )
    end
end

function annotate_against_baseline!(records)
    for precision_bits in STATIONARITY_PRECISIONS_BITS
        selected = [
            record for record in records
            if record["precision_bits"] == precision_bits
        ]
        baseline_index = findfirst(
            record -> record["solver_label"] == "baseline",
            selected)
        baseline_index === nothing &&
            error("missing baseline stationarity row")
        baseline = selected[baseline_index]

        require(
            baseline["center"]["status"] == "ok" &&
            baseline["derivative"]["status"] == "ok",
            "baseline stationarity row failed")

        baseline_x = parse.(BigFloat, baseline["center"]["x"])
        baseline_f = parse(
            BigFloat,
            baseline["center"]["f_res_reduced_extensive"])
        baseline_dv = parse(
            BigFloat,
            baseline["derivative"]["dAres_dV_pa"])

        for record in selected
            center = record["center"]
            derivative = record["derivative"]

            if center["status"] == "ok"
                current_x = parse.(BigFloat, center["x"])
                current_f = parse(
                    BigFloat,
                    center["f_res_reduced_extensive"])
                record["max_abs_x_delta_vs_baseline"] =
                    big_string(maximum(abs.(current_x .- baseline_x)))
                record["abs_f_res_delta_vs_baseline"] =
                    big_string(abs(current_f - baseline_f))
            else
                record["max_abs_x_delta_vs_baseline"] = nothing
                record["abs_f_res_delta_vs_baseline"] = nothing
            end

            if derivative["status"] == "ok"
                current_dv = parse(
                    BigFloat, derivative["dAres_dV_pa"])
                record["abs_dAres_dV_delta_vs_baseline_pa"] =
                    big_string(abs(current_dv - baseline_dv))
            else
                record["abs_dAres_dV_delta_vs_baseline_pa"] =
                    nothing
            end
        end
    end
end

function main_stationarity_floor()
    states_path = arg_value("--states")
    oracle_path = arg_value("--oracle")
    output_path = arg_value("--out")

    source = JSON.parsefile(states_path)
    frozen = JSON.parsefile(oracle_path)

    require(
        source["schema"] ==
            "MPMC_HNU/CPA/ThermoPack-phase-kernel-reference/v2",
        "unexpected stationarity coordinate-source schema")
    require(
        frozen["schema"] ==
            "MPMC_HNU/CPA/Clapeyron-phase-kernel-reference/v1",
        "unexpected stationarity frozen-oracle schema")
    require(
        frozen["source"]["commit"] == PINNED_CLAPEYRON_COMMIT,
        "stationarity audit Clapeyron revision mismatch")

    state = state4_liquid(source, frozen)

    baseline_model = build_model()
    model_audit = audit_model(baseline_model)

    records = Any[]
    for precision_bits in STATIONARITY_PRECISIONS_BITS
        for setting in STATIONARITY_SETTINGS
            push!(records, run_setting(
                state, precision_bits, setting))
        end
    end

    annotate_against_baseline!(records)

    result = Dict(
        "schema" =>
            "MPMC_HNU/CPA/Clapeyron-association-stationarity-floor-audit/v1",
        "source" => Dict(
            "software" => "ClapeyronThermo/Clapeyron.jl",
            "commit" => PINNED_CLAPEYRON_COMMIT,
            "julia_version" => string(VERSION),
            "scalar_api" => "Clapeyron.eos_res(model,V,T,z)",
        ),
        "model_audit" => model_audit,
        "state" => Dict(
            "label" => TARGET_LABEL,
            "temperature_k" => state.temperature_k,
            "molar_density_mol_per_m3" => state.rho,
            "composition_methanol" => state.x_meoh,
            "target_pressure_pa" => state.target_pressure_pa,
            "clapeyron_internal_pressure_pa" =>
                state.frozen_clapeyron_pressure_pa,
            "frozen_f_res" => state.frozen_f_res,
        ),
        "precision_bits" => STATIONARITY_PRECISIONS_BITS,
        "solver_matrix" => [
            Dict(
                "label" => setting.label,
                "rtol" => setting.tolerance,
                "atol" => setting.tolerance,
                "max_iters" => setting.max_iters,
            )
            for setting in STATIONARITY_SETTINGS
        ],
        "derivative_probe" => Dict(
            "stencil" => "central_7_o6",
            "relative_h" => DERIVATIVE_RELATIVE_H,
        ),
        "records" => records,
    )

    open(output_path, "w") do io
        JSON.print(io, result, 2)
        println(io)
    end

    println(
        "CLAPEYRON_ASSOC_STATIONARITY_AUDIT_OK",
        " state=", replace(TARGET_LABEL, " " => "_"),
        " precisions=", join(STATIONARITY_PRECISIONS_BITS, ","),
        " settings=", length(STATIONARITY_SETTINGS),
        " derivative_r=", DERIVATIVE_RELATIVE_H,
    )

    for record in records
        center = record["center"]
        derivative = record["derivative"]
        println(
            "CLAPEYRON_ASSOC_STATIONARITY_ROW",
            " precision_bits=", record["precision_bits"],
            " solver=", record["solver_label"],
            " rtol=", record["rtol"],
            " max_iters=", record["max_iters"],
            " center_status=", center["status"],
            " mass_action_residual=",
            get(center, "mass_action_max_residual", "NA"),
            " max_dX_baseline=",
            get(record, "max_abs_x_delta_vs_baseline", "NA"),
            " dF_baseline=",
            get(record, "abs_f_res_delta_vs_baseline", "NA"),
            " derivative_status=", derivative["status"],
            " dDV_baseline_pa=",
            get(
                record,
                "abs_dAres_dV_delta_vs_baseline_pa",
                "NA"),
            " dP_target_pa=",
            get(
                derivative,
                "abs_pressure_delta_vs_target_pa",
                "NA"),
            " dP_internal_pa=",
            get(
                derivative,
                "abs_pressure_delta_vs_clapeyron_internal_pa",
                "NA"),
        )
    end
end

if abspath(PROGRAM_FILE) == @__FILE__
    main_stationarity_floor()
end
