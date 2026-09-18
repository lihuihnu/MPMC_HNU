#!/usr/bin/env julia

include(joinpath(@__DIR__, "generate_oracle.jl"))

const RELATIVE_VOLUME_STEPS = Float64[
    1.0e-2, 5.0e-3, 2.0e-3, 1.0e-3,
    5.0e-4, 2.0e-4, 1.0e-4, 5.0e-5,
    2.0e-5, 1.0e-5, 5.0e-6, 2.0e-6,
    1.0e-6, 5.0e-7, 2.0e-7, 1.0e-7,
    5.0e-8, 2.0e-8, 1.0e-8,
]

const STENCILS = ("central_3_o2", "central_5_o4", "central_7_o6")

function derivative_3(values, h)
    return (values[1] - values[-1]) / (2.0 * h)
end

function derivative_5(values, h)
    return (
        values[-2] - 8.0 * values[-1] +
        8.0 * values[1] - values[2]
    ) / (12.0 * h)
end

function derivative_7(values, h)
    return (
        -values[-3] + 9.0 * values[-2] - 45.0 * values[-1] +
        45.0 * values[1] - 9.0 * values[2] + values[3]
    ) / (60.0 * h)
end

function state_coordinate(source, flat_index::Int)
    source_index = div(flat_index - 1, 2) + 1
    phase = isodd(flat_index) ? "liquid" : "vapor"
    source_state = source["states"][source_index]
    phase_record = source_state[phase]
    return (
        source_index = source_index - 1,
        phase = phase,
        temperature_k = Float64(source_state["temperature_k"]),
        target_pressure_pa = Float64(source_state["pressure_pa"]),
        rho = Float64(phase_record["molar_density_mol_per_m3"]),
        x_meoh = Float64(phase_record["composition"]["MEOH"]),
    )
end

function audit_state(model, source, oracle_state, flat_index::Int)
    coord = state_coordinate(source, flat_index)
    require(
        oracle_state["label"] ==
            "state=$(coord.source_index) phase=$(coord.phase)",
        "oracle/source state ordering mismatch at flat index $flat_index")
    require_close(
        Float64(oracle_state["temperature_k"]),
        coord.temperature_k,
        "oracle/source temperature mismatch")
    require_close(
        Float64(oracle_state["molar_density_mol_per_m3"]),
        coord.rho,
        "oracle/source density mismatch")
    require_close(
        Float64(oracle_state["composition"]["MEOH"]),
        coord.x_meoh,
        "oracle/source composition mismatch")

    z = [coord.x_meoh, 1.0 - coord.x_meoh]
    n = sum(z)
    volume = n / coord.rho
    rgas = Clapeyron.Rgas(model)
    rt = rgas * coord.temperature_k

    center_a = Float64(Clapeyron.eos_res(
        model, volume, coord.temperature_k, z))
    center_a_repeat = Float64(Clapeyron.eos_res(
        model, volume, coord.temperature_k, z))
    center_repeat_delta = abs(center_a - center_a_repeat)

    oracle_f = Float64(oracle_state["f_res_reduced_extensive"])
    center_f = center_a / rt
    center_f_delta = abs(center_f - oracle_f)

    internal_pressure = Float64(oracle_state["pressure_pa"])
    internal_mu = Float64[
        oracle_state["mu_residual_over_rt"]["MEOH"],
        oracle_state["mu_residual_over_rt"]["H2O"],
    ]
    internal_lnphi = Float64[
        oracle_state["ln_phi"]["MEOH"],
        oracle_state["ln_phi"]["H2O"],
    ]

    z_target = coord.target_pressure_pa * volume / (n * rt)
    require(z_target > 0.0 && isfinite(z_target),
            "target state has invalid Z")
    logz_target = log(z_target)
    target_z_lnphi = internal_mu .- logz_target

    steps = Any[]
    for relative_h in RELATIVE_VOLUME_STEPS
        h = relative_h * volume
        require(volume - 3.0 * h > 0.0,
                "finite-difference stencil crossed V=0")

        scalar = Dict{Int,Float64}()
        for multiplier in (-3, -2, -1, 1, 2, 3)
            shifted_v = volume + multiplier * h
            scalar[multiplier] = Float64(Clapeyron.eos_res(
                model, shifted_v, coord.temperature_k, z))
        end

        derivatives = Dict(
            "central_3_o2" => derivative_3(scalar, h),
            "central_5_o4" => derivative_5(scalar, h),
            "central_7_o6" => derivative_7(scalar, h),
        )

        stencil_results = Dict{String,Any}()
        for stencil in STENCILS
            d_a_d_v = derivatives[stencil]
            pressure_fd = n * rt / volume - d_a_d_v
            require(isfinite(pressure_fd) && pressure_fd > 0.0,
                    "finite-difference pressure became invalid")

            z_fd = pressure_fd * volume / (n * rt)
            require(isfinite(z_fd) && z_fd > 0.0,
                    "finite-difference Z became invalid")
            logz_fd = log(z_fd)
            lnphi_fd_z = internal_mu .- logz_fd

            stencil_results[stencil] = Dict(
                "dAres_dV_pa" => d_a_d_v,
                "pressure_pa" => pressure_fd,
                "abs_pressure_delta_vs_clapeyron_internal_pa" =>
                    abs(pressure_fd - internal_pressure),
                "abs_pressure_delta_vs_target_pa" =>
                    abs(pressure_fd - coord.target_pressure_pa),
                "log_z" => logz_fd,
                "minus_log_z" => -logz_fd,
                "abs_log_z_delta_vs_target" =>
                    abs(logz_fd - logz_target),
                "ln_phi_from_fixed_mu_and_fd_z" => Dict(
                    "MEOH" => lnphi_fd_z[1],
                    "H2O" => lnphi_fd_z[2],
                ),
                "max_abs_ln_phi_delta_vs_clapeyron_internal" =>
                    maximum(abs.(lnphi_fd_z .- internal_lnphi)),
                "max_abs_ln_phi_delta_vs_target_z" =>
                    maximum(abs.(lnphi_fd_z .- target_z_lnphi)),
            )
        end

        push!(steps, Dict(
            "relative_h" => relative_h,
            "h_m3" => h,
            "scalar_Ares_joule" => Dict(
                string(multiplier) => scalar[multiplier]
                for multiplier in (-3, -2, -1, 1, 2, 3)
            ),
            "stencils" => stencil_results,
        ))
    end

    return Dict(
        "label" => oracle_state["label"],
        "temperature_k" => coord.temperature_k,
        "molar_density_mol_per_m3" => coord.rho,
        "molar_volume_m3_per_mol" => volume,
        "composition" => Dict("MEOH" => z[1], "H2O" => z[2]),
        "target_pressure_pa" => coord.target_pressure_pa,
        "clapeyron_internal_pressure_pa" => internal_pressure,
        "clapeyron_internal_mu_res_over_rt" => Dict(
            "MEOH" => internal_mu[1],
            "H2O" => internal_mu[2],
        ),
        "clapeyron_internal_ln_phi" => Dict(
            "MEOH" => internal_lnphi[1],
            "H2O" => internal_lnphi[2],
        ),
        "target_z_ln_phi_from_same_mu" => Dict(
            "MEOH" => target_z_lnphi[1],
            "H2O" => target_z_lnphi[2],
        ),
        "center_scalar_f_res" => center_f,
        "center_scalar_f_res_delta_vs_frozen" => center_f_delta,
        "center_scalar_repeat_delta_joule" => center_repeat_delta,
        "steps" => steps,
    )
end

function global_step_summary(states)
    result = Dict{String,Any}()
    for stencil in STENCILS
        rows = Any[]
        for (step_index, relative_h) in enumerate(RELATIVE_VOLUME_STEPS)
            max_p_internal = 0.0
            max_p_target = 0.0
            max_lnphi_internal = 0.0
            max_lnphi_target_z = 0.0
            max_logz_target = 0.0

            for state in states
                record = state["steps"][step_index]["stencils"][stencil]
                max_p_internal = max(
                    max_p_internal,
                    record["abs_pressure_delta_vs_clapeyron_internal_pa"])
                max_p_target = max(
                    max_p_target,
                    record["abs_pressure_delta_vs_target_pa"])
                max_lnphi_internal = max(
                    max_lnphi_internal,
                    record["max_abs_ln_phi_delta_vs_clapeyron_internal"])
                max_lnphi_target_z = max(
                    max_lnphi_target_z,
                    record["max_abs_ln_phi_delta_vs_target_z"])
                max_logz_target = max(
                    max_logz_target,
                    record["abs_log_z_delta_vs_target"])
            end

            push!(rows, Dict(
                "relative_h" => relative_h,
                "max_abs_pressure_delta_vs_clapeyron_internal_pa" =>
                    max_p_internal,
                "max_abs_pressure_delta_vs_target_pa" => max_p_target,
                "max_abs_ln_phi_delta_vs_clapeyron_internal" =>
                    max_lnphi_internal,
                "max_abs_ln_phi_delta_vs_target_z" => max_lnphi_target_z,
                "max_abs_log_z_delta_vs_target" => max_logz_target,
            ))
        end

        best_pressure_target = rows[argmin(
            [row["max_abs_pressure_delta_vs_target_pa"] for row in rows])]
        best_pressure_internal = rows[argmin(
            [row["max_abs_pressure_delta_vs_clapeyron_internal_pa"] for row in rows])]
        best_lnphi_target = rows[argmin(
            [row["max_abs_ln_phi_delta_vs_target_z"] for row in rows])]

        result[stencil] = Dict(
            "steps" => rows,
            "best_global_pressure_vs_target" => best_pressure_target,
            "best_global_pressure_vs_clapeyron_internal" =>
                best_pressure_internal,
            "best_global_ln_phi_vs_target_z" => best_lnphi_target,
        )
    end
    return result
end

function main_derivative_floor()
    states_path = arg_value("--states")
    oracle_path = arg_value("--oracle")
    output_path = arg_value("--out")

    source = JSON.parsefile(states_path)
    frozen = JSON.parsefile(oracle_path)
    require(
        source["schema"] ==
            "MPMC_HNU/CPA/ThermoPack-phase-kernel-reference/v2",
        "unexpected phase-coordinate source schema")
    require(
        frozen["schema"] ==
            "MPMC_HNU/CPA/Clapeyron-phase-kernel-reference/v1",
        "unexpected frozen Clapeyron oracle schema")
    require(
        frozen["source"]["commit"] == PINNED_CLAPEYRON_COMMIT,
        "frozen Clapeyron revision mismatch")
    require(length(frozen["states"]) == 10,
            "derivative-floor audit requires ten frozen states")

    model = build_model()
    model_audit = audit_model(model)

    states = Any[]
    for flat_index in 1:10
        push!(states, audit_state(
            model, source, frozen["states"][flat_index], flat_index))
    end

    max_center_f_delta = maximum(
        state["center_scalar_f_res_delta_vs_frozen"] for state in states)
    max_center_repeat_delta = maximum(
        state["center_scalar_repeat_delta_joule"] for state in states)

    result = Dict(
        "schema" =>
            "MPMC_HNU/CPA/Clapeyron-derivative-floor-audit/v1",
        "source" => Dict(
            "software" => "ClapeyronThermo/Clapeyron.jl",
            "commit" => PINNED_CLAPEYRON_COMMIT,
            "scalar_api" => "Clapeyron.eos_res(model,V,T,z)",
            "julia_version" => string(VERSION),
        ),
        "model_audit" => model_audit,
        "frozen_relative_volume_steps" => RELATIVE_VOLUME_STEPS,
        "stencils" => Dict(
            "central_3_o2" =>
                "(A(V+h)-A(V-h))/(2h)",
            "central_5_o4" =>
                "(A(V-2h)-8A(V-h)+8A(V+h)-A(V+2h))/(12h)",
            "central_7_o6" =>
                "(-A(V-3h)+9A(V-2h)-45A(V-h)+45A(V+h)-9A(V+2h)+A(V+3h))/(60h)",
        ),
        "max_center_scalar_f_res_delta_vs_frozen" =>
            max_center_f_delta,
        "max_center_scalar_repeat_delta_joule" =>
            max_center_repeat_delta,
        "states" => states,
        "global_step_summary" => global_step_summary(states),
    )

    open(output_path, "w") do io
        JSON.print(io, result, 2)
        println(io)
    end

    println(
        "CLAPEYRON_DERIVATIVE_FLOOR_AUDIT_OK",
        " states=", length(states),
        " steps=", length(RELATIVE_VOLUME_STEPS),
        " stencils=", length(STENCILS),
        " max_center_scalar_f_res_delta_vs_frozen=",
        max_center_f_delta,
        " max_center_scalar_repeat_delta_joule=",
        max_center_repeat_delta,
    )

    for stencil in STENCILS
        summary = result["global_step_summary"][stencil]
        p_target = summary["best_global_pressure_vs_target"]
        p_internal = summary["best_global_pressure_vs_clapeyron_internal"]
        ln_target = summary["best_global_ln_phi_vs_target_z"]
        println(
            "CLAPEYRON_DERIVATIVE_FLOOR_SUMMARY",
            " stencil=", stencil,
            " best_r_pressure_target=", p_target["relative_h"],
            " max_dP_target_pa=",
            p_target["max_abs_pressure_delta_vs_target_pa"],
            " best_r_pressure_internal=", p_internal["relative_h"],
            " max_dP_internal_pa=",
            p_internal["max_abs_pressure_delta_vs_clapeyron_internal_pa"],
            " best_r_lnphi_target=", ln_target["relative_h"],
            " max_dlnphi_target=",
            ln_target["max_abs_ln_phi_delta_vs_target_z"],
        )
    end
end

if abspath(PROGRAM_FILE) == @__FILE__
    main_derivative_floor()
end
