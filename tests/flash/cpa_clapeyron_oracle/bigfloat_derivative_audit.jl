#!/usr/bin/env julia

include(joinpath(@__DIR__, "generate_oracle.jl"))

const BIGFLOAT_PRECISIONS_BITS = Int[128, 256, 512]
const RELATIVE_VOLUME_STEP_STRINGS = String[
    "2e-3", "1e-3", "5e-4", "2e-4", "1e-4",
    "5e-5", "2e-5", "1e-5", "5e-6", "2e-6", "1e-6",
    "5e-7", "2e-7", "1e-7", "5e-8", "2e-8", "1e-8",
    "2e-9", "1e-10", "1e-11", "1e-12",
]
const PROBLEM_LABELS = Set([
    "state=0 phase=liquid",
    "state=2 phase=liquid",
    "state=3 phase=liquid",
    "state=4 phase=liquid",
])
const BIGFLOAT_STENCILS = ("central_5_o4", "central_7_o6")

function derivative_5_big(values, h)
    return (
        values[-2] - 8 * values[-1] +
        8 * values[1] - values[2]
    ) / (12 * h)
end

function derivative_7_big(values, h)
    return (
        -values[-3] + 9 * values[-2] - 45 * values[-1] +
        45 * values[1] - 9 * values[2] + values[3]
    ) / (60 * h)
end

function bigfloat_string(x)
    return string(x)
end

function frozen_problem_states(source, frozen)
    selected = Any[]
    for (flat_index, oracle_state) in enumerate(frozen["states"])
        label = String(oracle_state["label"])
        label in PROBLEM_LABELS || continue

        source_index = div(flat_index - 1, 2) + 1
        phase = isodd(flat_index) ? "liquid" : "vapor"
        require(phase == "liquid", "BigFloat audit selected non-liquid state")
        source_state = source["states"][source_index]
        require(
            label == "state=$(source_index-1) phase=liquid",
            "BigFloat audit state ordering mismatch for $label")

        push!(selected, (
            label = label,
            temperature_k = Float64(source_state["temperature_k"]),
            target_pressure_pa = Float64(source_state["pressure_pa"]),
            rho = Float64(
                source_state["liquid"]["molar_density_mol_per_m3"]),
            x_meoh = Float64(
                source_state["liquid"]["composition"]["MEOH"]),
            frozen_pressure_pa = Float64(oracle_state["pressure_pa"]),
            frozen_mu_meoh = Float64(
                oracle_state["mu_residual_over_rt"]["MEOH"]),
            frozen_mu_h2o = Float64(
                oracle_state["mu_residual_over_rt"]["H2O"]),
            frozen_lnphi_meoh = Float64(
                oracle_state["ln_phi"]["MEOH"]),
            frozen_lnphi_h2o = Float64(
                oracle_state["ln_phi"]["H2O"]),
        ))
    end
    require(length(selected) == 4, "BigFloat audit must select exactly four liquid states")
    return selected
end

function association_residual_big(model, volume, temperature, z)
    x = Clapeyron.assoc_fractions(model, volume, temperature, z)
    require(eltype(x.v) == BigFloat,
            "BigFloat association fractions did not preserve BigFloat")
    require(all(isfinite, x.v),
            "BigFloat association fractions became nonfinite")
    matrix = Clapeyron.assoc_site_matrix(model, volume, temperature, z)
    residual = x.v .* (1 .+ matrix * x.v) .- 1
    require(all(isfinite, residual),
            "BigFloat association equation residual became nonfinite")
    return maximum(abs, residual), x
end

function capability_probe(model, state)
    result = Dict{String,Any}()

    setprecision(BigFloat, 256) do
        volume = BigFloat(1.0 / state.rho)
        temperature = BigFloat(state.temperature_k)
        z = BigFloat[state.x_meoh, 1.0 - state.x_meoh]

        scalar = Clapeyron.eos_res(model, volume, temperature, z)
        require(scalar isa BigFloat,
                "Clapeyron eos_res did not return BigFloat for BigFloat state inputs")
        require(isfinite(scalar),
                "Clapeyron BigFloat eos_res returned nonfinite value")

        association_residual, x =
            association_residual_big(model, volume, temperature, z)

        repeated = Clapeyron.eos_res(model, volume, temperature, z)
        require(repeated isa BigFloat,
                "repeated BigFloat eos_res lost BigFloat type")

        result["direct_state_bigfloat_supported"] = true
        result["scalar_type"] = string(typeof(scalar))
        result["association_fraction_eltype"] = string(eltype(x.v))
        result["scalar_repeat_delta"] =
            bigfloat_string(abs(scalar - repeated))
        result["association_equation_max_residual"] =
            bigfloat_string(association_residual)
    end

    try
        promoted = Clapeyron.promote_model(BigFloat, model)
        result["promote_model_bigfloat_supported"] = true
        result["promoted_model_type"] = string(typeof(promoted))
        result["promoted_model_eltype"] = string(eltype(promoted))
    catch error
        result["promote_model_bigfloat_supported"] = false
        result["promote_model_error_type"] = string(typeof(error))
        result["promote_model_error"] = sprint(showerror, error)
    end

    return result
end

function audit_precision_state(model, state, precision_bits::Int)
    return setprecision(BigFloat, precision_bits) do
        temperature = BigFloat(state.temperature_k)
        z = BigFloat[state.x_meoh, 1.0 - state.x_meoh]
        total_moles = sum(z)
        volume = total_moles / BigFloat(state.rho)
        rgas = BigFloat(Clapeyron.Rgas(model))
        rt = rgas * temperature
        target_pressure = BigFloat(state.target_pressure_pa)
        frozen_pressure = BigFloat(state.frozen_pressure_pa)
        frozen_mu = BigFloat[state.frozen_mu_meoh, state.frozen_mu_h2o]

        center_a = Clapeyron.eos_res(model, volume, temperature, z)
        require(center_a isa BigFloat && isfinite(center_a),
                "BigFloat center eos_res is invalid")
        center_a_repeat = Clapeyron.eos_res(model, volume, temperature, z)
        center_repeat_delta = abs(center_a - center_a_repeat)

        association_residual, _ =
            association_residual_big(model, volume, temperature, z)

        target_z = target_pressure * volume / (total_moles * rt)
        frozen_z = frozen_pressure * volume / (total_moles * rt)
        require(
            target_z > 0 && frozen_z > 0 &&
            isfinite(target_z) && isfinite(frozen_z),
            "BigFloat audit encountered invalid target/internal Z")

        target_lnphi = frozen_mu .- log(target_z)
        frozen_pressure_lnphi = frozen_mu .- log(frozen_z)

        step_records = Any[]
        for step_string in RELATIVE_VOLUME_STEP_STRINGS
            relative_h = parse(BigFloat, step_string)
            h = relative_h * volume
            require(volume - 3h > 0,
                    "BigFloat finite-difference stencil crossed V=0")

            values = Dict{Int,BigFloat}()
            for multiplier in (-3, -2, -1, 1, 2, 3)
                shifted_volume = volume + BigFloat(multiplier) * h
                scalar = Clapeyron.eos_res(
                    model, shifted_volume, temperature, z)
                require(scalar isa BigFloat && isfinite(scalar),
                        "BigFloat shifted eos_res became invalid")
                values[multiplier] = scalar
            end

            derivatives = Dict(
                "central_5_o4" => derivative_5_big(values, h),
                "central_7_o6" => derivative_7_big(values, h),
            )

            stencil_records = Dict{String,Any}()
            for stencil in BIGFLOAT_STENCILS
                derivative = derivatives[stencil]
                pressure =
                    total_moles * rt / volume - derivative
                require(isfinite(pressure) && pressure > 0,
                        "BigFloat reconstructed pressure became invalid")
                z_fd = pressure * volume / (total_moles * rt)
                require(isfinite(z_fd) && z_fd > 0,
                        "BigFloat reconstructed Z became invalid")
                lnphi_fd = frozen_mu .- log(z_fd)

                stencil_records[stencil] = Dict(
                    "pressure_pa" => bigfloat_string(pressure),
                    "abs_pressure_delta_vs_target_pa" =>
                        Float64(abs(pressure - target_pressure)),
                    "abs_pressure_delta_vs_clapeyron_internal_pa" =>
                        Float64(abs(pressure - frozen_pressure)),
                    "max_abs_ln_phi_delta_vs_target_z" =>
                        Float64(maximum(abs.(lnphi_fd .- target_lnphi))),
                    "max_abs_ln_phi_delta_vs_clapeyron_internal_z" =>
                        Float64(maximum(
                            abs.(lnphi_fd .- frozen_pressure_lnphi))),
                )
            end

            push!(step_records, Dict(
                "relative_h" => step_string,
                "h_m3" => bigfloat_string(h),
                "stencils" => stencil_records,
            ))
        end

        return Dict(
            "label" => state.label,
            "precision_bits" => precision_bits,
            "temperature_k" => state.temperature_k,
            "molar_density_mol_per_m3" => state.rho,
            "molar_volume_m3_per_mol" => Float64(volume),
            "composition_methanol" => state.x_meoh,
            "target_pressure_pa" => state.target_pressure_pa,
            "clapeyron_internal_pressure_pa" => state.frozen_pressure_pa,
            "center_scalar_Ares_joule" => bigfloat_string(center_a),
            "center_scalar_repeat_delta_joule" =>
                bigfloat_string(center_repeat_delta),
            "association_equation_max_residual" =>
                bigfloat_string(association_residual),
            "steps" => step_records,
        )
    end
end

function precision_summary(records, precision_bits, stencil)
    selected = [
        record for record in records
        if record["precision_bits"] == precision_bits
    ]
    require(length(selected) == 4,
            "BigFloat precision summary expected four states")

    rows = Any[]
    for (step_index, step_string) in enumerate(RELATIVE_VOLUME_STEP_STRINGS)
        max_p_target = 0.0
        max_p_internal = 0.0
        max_lnphi_target = 0.0
        max_lnphi_internal = 0.0
        for state in selected
            metric = state["steps"][step_index]["stencils"][stencil]
            max_p_target = max(
                max_p_target,
                metric["abs_pressure_delta_vs_target_pa"])
            max_p_internal = max(
                max_p_internal,
                metric["abs_pressure_delta_vs_clapeyron_internal_pa"])
            max_lnphi_target = max(
                max_lnphi_target,
                metric["max_abs_ln_phi_delta_vs_target_z"])
            max_lnphi_internal = max(
                max_lnphi_internal,
                metric["max_abs_ln_phi_delta_vs_clapeyron_internal_z"])
        end

        push!(rows, Dict(
            "relative_h" => step_string,
            "max_abs_pressure_delta_vs_target_pa" => max_p_target,
            "max_abs_pressure_delta_vs_clapeyron_internal_pa" =>
                max_p_internal,
            "max_abs_ln_phi_delta_vs_target_z" => max_lnphi_target,
            "max_abs_ln_phi_delta_vs_clapeyron_internal_z" =>
                max_lnphi_internal,
        ))
    end

    best_p_target = rows[argmin([
        row["max_abs_pressure_delta_vs_target_pa"] for row in rows
    ])]
    best_p_internal = rows[argmin([
        row["max_abs_pressure_delta_vs_clapeyron_internal_pa"] for row in rows
    ])]
    best_lnphi_target = rows[argmin([
        row["max_abs_ln_phi_delta_vs_target_z"] for row in rows
    ])]

    return Dict(
        "precision_bits" => precision_bits,
        "stencil" => stencil,
        "steps" => rows,
        "best_global_pressure_vs_target" => best_p_target,
        "best_global_pressure_vs_clapeyron_internal" => best_p_internal,
        "best_global_ln_phi_vs_target_z" => best_lnphi_target,
    )
end

function state_best_summary(records, precision_bits, stencil)
    selected = [
        record for record in records
        if record["precision_bits"] == precision_bits
    ]
    result = Any[]
    for state in selected
        rows = Any[]
        for (step_index, step_string) in enumerate(RELATIVE_VOLUME_STEP_STRINGS)
            metric = state["steps"][step_index]["stencils"][stencil]
            push!(rows, Dict(
                "relative_h" => step_string,
                "abs_pressure_delta_vs_target_pa" =>
                    metric["abs_pressure_delta_vs_target_pa"],
                "max_abs_ln_phi_delta_vs_target_z" =>
                    metric["max_abs_ln_phi_delta_vs_target_z"],
            ))
        end
        best_p = rows[argmin([
            row["abs_pressure_delta_vs_target_pa"] for row in rows
        ])]
        push!(result, Dict(
            "label" => state["label"],
            "best_pressure_vs_target" => best_p,
        ))
    end
    return result
end

function main_bigfloat_derivative_audit()
    states_path = arg_value("--states")
    oracle_path = arg_value("--oracle")
    output_path = arg_value("--out")

    source = JSON.parsefile(states_path)
    frozen = JSON.parsefile(oracle_path)
    require(
        source["schema"] ==
            "MPMC_HNU/CPA/ThermoPack-phase-kernel-reference/v2",
        "unexpected BigFloat coordinate-source schema")
    require(
        frozen["schema"] ==
            "MPMC_HNU/CPA/Clapeyron-phase-kernel-reference/v1",
        "unexpected BigFloat frozen-oracle schema")
    require(
        frozen["source"]["commit"] == PINNED_CLAPEYRON_COMMIT,
        "BigFloat frozen Clapeyron revision mismatch")

    model = build_model()
    model_audit = audit_model(model)
    problem_states = frozen_problem_states(source, frozen)

    capability = capability_probe(model, first(problem_states))
    require(
        capability["direct_state_bigfloat_supported"] == true,
        "pinned Clapeyron direct BigFloat scalar path is unsupported")

    records = Any[]
    for precision_bits in BIGFLOAT_PRECISIONS_BITS
        for state in problem_states
            push!(records, audit_precision_state(
                model, state, precision_bits))
        end
    end

    summaries = Any[]
    state_best = Dict{String,Any}()
    for precision_bits in BIGFLOAT_PRECISIONS_BITS
        for stencil in BIGFLOAT_STENCILS
            push!(summaries, precision_summary(
                records, precision_bits, stencil))
            state_best["$(precision_bits)-$(stencil)"] =
                state_best_summary(records, precision_bits, stencil)
        end
    end

    result = Dict(
        "schema" =>
            "MPMC_HNU/CPA/Clapeyron-BigFloat-derivative-audit/v1",
        "source" => Dict(
            "software" => "ClapeyronThermo/Clapeyron.jl",
            "commit" => PINNED_CLAPEYRON_COMMIT,
            "scalar_api" => "Clapeyron.eos_res(model,V,T,z)",
            "julia_version" => string(VERSION),
        ),
        "capability" => capability,
        "model_audit" => model_audit,
        "problem_states" => collect(PROBLEM_LABELS),
        "precision_bits" => BIGFLOAT_PRECISIONS_BITS,
        "relative_volume_steps" => RELATIVE_VOLUME_STEP_STRINGS,
        "stencils" => collect(BIGFLOAT_STENCILS),
        "association_numerics_unchanged" => Dict(
            "rtol" => model.assoc_options.rtol,
            "atol" => model.assoc_options.atol,
            "max_iters" => model.assoc_options.max_iters,
            "dampingfactor" => model.assoc_options.dampingfactor,
            "implicit_ad" => model.assoc_options.implicit_ad,
        ),
        "records" => records,
        "precision_summaries" => summaries,
        "state_best_summaries" => state_best,
    )

    open(output_path, "w") do io
        JSON.print(io, result, 2)
        println(io)
    end

    println(
        "CLAPEYRON_BIGFLOAT_CAPABILITY_OK",
        " direct_state_bigfloat=",
        capability["direct_state_bigfloat_supported"],
        " scalar_type=", capability["scalar_type"],
        " association_eltype=",
        capability["association_fraction_eltype"],
        " promote_model_supported=",
        capability["promote_model_bigfloat_supported"],
    )

    for summary in summaries
        p_target = summary["best_global_pressure_vs_target"]
        p_internal =
            summary["best_global_pressure_vs_clapeyron_internal"]
        ln_target = summary["best_global_ln_phi_vs_target_z"]
        println(
            "CLAPEYRON_BIGFLOAT_DERIVATIVE_SUMMARY",
            " precision_bits=", summary["precision_bits"],
            " stencil=", summary["stencil"],
            " best_r_pressure_target=", p_target["relative_h"],
            " max_dP_target_pa=",
            p_target["max_abs_pressure_delta_vs_target_pa"],
            " best_r_pressure_internal=", p_internal["relative_h"],
            " max_dP_internal_pa=",
            p_internal[
                "max_abs_pressure_delta_vs_clapeyron_internal_pa"],
            " best_r_lnphi_target=", ln_target["relative_h"],
            " max_dlnphi_target=",
            ln_target["max_abs_ln_phi_delta_vs_target_z"],
        )
    end
end

if abspath(PROGRAM_FILE) == @__FILE__
    main_bigfloat_derivative_audit()
end
