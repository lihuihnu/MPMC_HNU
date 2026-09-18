#!/usr/bin/env julia

include(joinpath(@__DIR__, "independent_exact2_audit.jl"))

const TEN_STATE_PRECISIONS_BITS = Int[256, 512]
const TEN_STATE_PRESSURE_THRESHOLD_PA = 5.0e-6
const TEN_STATE_LNPHI_THRESHOLD = 1.0e-10
const TEN_STATE_DERIVATIVE_RELATIVE_H = "1e-6"

function ten_state_records(source, frozen)
    records = Any[]
    require(length(frozen["states"]) == 10,
            "ten-state diagnostic requires exactly ten frozen states")

    for (flat_index, oracle_state) in enumerate(frozen["states"])
        source_index = div(flat_index - 1, 2) + 1
        phase = isodd(flat_index) ? "liquid" : "vapor"
        source_state = source["states"][source_index]
        source_phase = source_state[phase]
        expected_label =
            "state=$(source_index-1) phase=$phase"
        require(
            String(oracle_state["label"]) == expected_label,
            "ten-state source/oracle ordering mismatch")

        push!(records, (
            label = expected_label,
            temperature_k = Float64(source_state["temperature_k"]),
            target_pressure_pa = Float64(source_state["pressure_pa"]),
            rho = Float64(source_phase["molar_density_mol_per_m3"]),
            x_meoh = Float64(source_phase["composition"]["MEOH"]),
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
    return records
end

function exact2_path_probe(model, volume, temperature, z)
    delta = Clapeyron.delta_assoc(
        model, volume, temperature, z, nothing)
    original = Clapeyron.assoc_site_matrix(
        model, volume, temperature, z, nothing, delta)

    compress = Clapeyron.__maybe_compress(original)
    if !compress
        return (
            path = "uncorrected_raw_path",
            original = original,
            reduced = original,
            mapping = collect(1:size(original, 1)),
            actual_initializer_success = false,
            pseudodiagonal_initializer = false,
        )
    end

    reduced, mapping =
        Clapeyron.compress_assoc_matrix(original)
    mapping_vector = collect(mapping)

    if size(reduced) != (2, 2)
        return (
            path = "uncorrected_raw_path",
            original = original,
            reduced = reduced,
            mapping = mapping_vector,
            actual_initializer_success = false,
            pseudodiagonal_initializer = false,
        )
    end

    pseudo_buffer = Vector{eltype(reduced)}(
        undef, size(reduced, 1))
    pseudo_buffer, pseudo_init, _ =
        let
            result = Clapeyron.X_maybe_exact_pseudodiag!(
                reduced, pseudo_buffer)
            (pseudo_buffer, result[1], result[2])
        end

    actual_buffer = Vector{eltype(reduced)}(
        undef, size(reduced, 1))
    actual_buffer, actual_success =
        Clapeyron.assoc_matrix_x0!(
            reduced, actual_buffer)

    exact2_path =
        !pseudo_init && actual_success

    return (
        path = exact2_path ?
            "compressed_exact2_corrected" :
            "uncorrected_raw_path",
        original = original,
        reduced = reduced,
        mapping = mapping_vector,
        actual_initializer_success = actual_success,
        pseudodiagonal_initializer = pseudo_init,
    )
end

function corrected_scalar_state(
    model,
    volume,
    temperature,
    z)

    path = exact2_path_probe(
        model, volume, temperature, z)

    if path.path != "compressed_exact2_corrected"
        raw_A = Clapeyron.eos_res(
            model, volume, temperature, z)
        return Dict(
            "path" => path.path,
            "corrected" => false,
            "original_dimension" => size(path.original, 1),
            "reduced_dimension" => size(path.reduced, 1),
            "mapping" => path.mapping,
            "pseudodiagonal_initializer" =>
                path.pseudodiagonal_initializer,
            "actual_initializer_success" =>
                path.actual_initializer_success,
            "A_res_joule" => big_string(raw_A),
        )
    end

    require(
        size(path.original) == (4, 4) &&
        size(path.reduced) == (2, 2),
        "corrected exact2 path has unexpected dimensions")

    solution = independent_newton(path.reduced)
    require(
        solution.max_abs_residual <=
            parse(BigFloat, NEWTON_RESIDUAL_TARGET),
        "ten-state independent Newton missed residual target")

    expanded = BigFloat[
        solution.x[index] for index in path.mapping
    ]
    expanded_residual =
        residual_vector(path.original, expanded)
    expanded_norm = maximum(abs, expanded_residual)

    q_independent = independent_q_from_original_matrix(
        model, path.original, z, expanded)

    total_moles = sum(z)
    rt = BigFloat(Clapeyron.Rgas(model)) * temperature
    A_assoc = rt * q_independent.q_stationary
    A_cubic = Clapeyron.eos_res(
        model.cubicmodel, volume, temperature, z)
    A_total = A_cubic + A_assoc

    return Dict(
        "path" => path.path,
        "corrected" => true,
        "original_dimension" => size(path.original, 1),
        "reduced_dimension" => size(path.reduced, 1),
        "mapping" => path.mapping,
        "pseudodiagonal_initializer" =>
            path.pseudodiagonal_initializer,
        "actual_initializer_success" =>
            path.actual_initializer_success,
        "newton_iterations" => solution.iterations,
        "reduced_max_abs_residual" =>
            big_string(solution.max_abs_residual),
        "expanded_max_abs_residual" =>
            big_string(expanded_norm),
        "A_assoc_joule" => big_string(A_assoc),
        "A_cubic_joule" => big_string(A_cubic),
        "A_res_joule" => big_string(A_total),
    )
end

function corrected_pressure(
    model,
    volume,
    temperature,
    z,
    frozen_pressure)

    center = corrected_scalar_state(
        model, volume, temperature, z)
    corrected = Bool(center["corrected"])

    if !corrected
        return (
            center = center,
            pressure = BigFloat(frozen_pressure),
            derivative = nothing,
            shifted_paths = String[],
        )
    end

    relative_h =
        parse(BigFloat, TEN_STATE_DERIVATIVE_RELATIVE_H)
    h = relative_h * volume
    scalars = Dict{Int,BigFloat}()
    shifted_paths = String[]

    for multiplier in (-3, -2, -1, 1, 2, 3)
        shifted_volume =
            volume + BigFloat(multiplier) * h
        shifted = corrected_scalar_state(
            model, shifted_volume, temperature, z)

        require(
            Bool(shifted["corrected"]),
            "exact2 correction path changed inside derivative stencil")
        require(
            parse(
                BigFloat,
                shifted["reduced_max_abs_residual"]) <=
                parse(BigFloat, NEWTON_RESIDUAL_TARGET),
            "shifted independent Newton missed residual target")

        scalars[multiplier] =
            parse(BigFloat, shifted["A_res_joule"])
        push!(shifted_paths, String(shifted["path"]))
    end

    derivative = seven_point(scalars, h)
    total_moles = sum(z)
    rt = BigFloat(Clapeyron.Rgas(model)) * temperature
    pressure =
        total_moles * rt / volume - derivative

    require(
        isfinite(pressure) && pressure > 0,
        "ten-state corrected pressure became invalid")

    return (
        center = center,
        pressure = pressure,
        derivative = derivative,
        shifted_paths = shifted_paths,
    )
end

function run_ten_state_precision(
    states,
    precision_bits)

    return setprecision(BigFloat, precision_bits) do
        model = build_model()
        records = Any[]

        max_pressure_delta = 0.0
        max_lnphi_delta = 0.0
        corrected_count = 0
        raw_count = 0

        for state in states
            temperature = BigFloat(state.temperature_k)
            z = BigFloat[
                state.x_meoh,
                1.0 - state.x_meoh,
            ]
            total_moles = sum(z)
            volume = total_moles / BigFloat(state.rho)
            rt =
                BigFloat(Clapeyron.Rgas(model)) * temperature

            pressure_result = corrected_pressure(
                model,
                volume,
                temperature,
                z,
                state.frozen_pressure_pa)

            pressure = pressure_result.pressure
            corrected =
                Bool(pressure_result.center["corrected"])
            corrected_count += corrected ? 1 : 0
            raw_count += corrected ? 0 : 1

            target_pressure =
                BigFloat(state.target_pressure_pa)
            pressure_delta =
                Float64(abs(pressure - target_pressure))

            mu = BigFloat[
                state.frozen_mu_meoh,
                state.frozen_mu_h2o,
            ]
            z_corrected =
                pressure * volume / (total_moles * rt)
            z_target =
                target_pressure * volume /
                (total_moles * rt)

            require(
                z_corrected > 0 && z_target > 0 &&
                isfinite(z_corrected) &&
                isfinite(z_target),
                "ten-state diagnostic encountered invalid Z")

            lnphi_corrected =
                mu .- log(z_corrected)
            lnphi_target_z =
                mu .- log(z_target)
            lnphi_delta =
                Float64(maximum(abs.(
                    lnphi_corrected .-
                    lnphi_target_z)))

            max_pressure_delta =
                max(max_pressure_delta, pressure_delta)
            max_lnphi_delta =
                max(max_lnphi_delta, lnphi_delta)

            push!(records, Dict(
                "label" => state.label,
                "precision_bits" => precision_bits,
                "path" =>
                    pressure_result.center["path"],
                "corrected" => corrected,
                "original_dimension" =>
                    pressure_result.center[
                        "original_dimension"],
                "reduced_dimension" =>
                    pressure_result.center[
                        "reduced_dimension"],
                "mapping" =>
                    pressure_result.center["mapping"],
                "pseudodiagonal_initializer" =>
                    pressure_result.center[
                        "pseudodiagonal_initializer"],
                "actual_initializer_success" =>
                    pressure_result.center[
                        "actual_initializer_success"],
                "reduced_max_abs_residual" =>
                    get(
                        pressure_result.center,
                        "reduced_max_abs_residual",
                        nothing),
                "expanded_max_abs_residual" =>
                    get(
                        pressure_result.center,
                        "expanded_max_abs_residual",
                        nothing),
                "corrected_pressure_pa" =>
                    big_string(pressure),
                "target_pressure_pa" =>
                    state.target_pressure_pa,
                "raw_frozen_pressure_pa" =>
                    state.frozen_pressure_pa,
                "abs_pressure_delta_vs_target_pa" =>
                    pressure_delta,
                "corrected_z" =>
                    big_string(z_corrected),
                "target_z" =>
                    big_string(z_target),
                "corrected_z_only_ln_phi" => Dict(
                    "MEOH" =>
                        big_string(lnphi_corrected[1]),
                    "H2O" =>
                        big_string(lnphi_corrected[2]),
                ),
                "target_z_only_ln_phi" => Dict(
                    "MEOH" =>
                        big_string(lnphi_target_z[1]),
                    "H2O" =>
                        big_string(lnphi_target_z[2]),
                ),
                "max_abs_z_only_ln_phi_delta_vs_target" =>
                    lnphi_delta,
                "pressure_threshold_pass" =>
                    pressure_delta <=
                        TEN_STATE_PRESSURE_THRESHOLD_PA,
                "ln_phi_threshold_pass" =>
                    lnphi_delta <=
                        TEN_STATE_LNPHI_THRESHOLD,
            ))
        end

        diagnostic_pass =
            max_pressure_delta <=
                TEN_STATE_PRESSURE_THRESHOLD_PA &&
            max_lnphi_delta <=
                TEN_STATE_LNPHI_THRESHOLD

        return Dict(
            "precision_bits" => precision_bits,
            "corrected_state_count" => corrected_count,
            "uncorrected_raw_state_count" => raw_count,
            "max_abs_pressure_delta_vs_target_pa" =>
                max_pressure_delta,
            "max_abs_z_only_ln_phi_delta_vs_target" =>
                max_lnphi_delta,
            "pressure_threshold_pa" =>
                TEN_STATE_PRESSURE_THRESHOLD_PA,
            "ln_phi_threshold" =>
                TEN_STATE_LNPHI_THRESHOLD,
            "diagnostic_pass" => diagnostic_pass,
            "records" => records,
        )
    end
end

function main_ten_state_diagnostic()
    states_path = arg_value("--states")
    oracle_path = arg_value("--oracle")
    output_path = arg_value("--out")

    source = JSON.parsefile(states_path)
    frozen = JSON.parsefile(oracle_path)

    require(
        source["schema"] ==
            "MPMC_HNU/CPA/ThermoPack-phase-kernel-reference/v2",
        "unexpected ten-state coordinate-source schema")
    require(
        frozen["schema"] ==
            "MPMC_HNU/CPA/Clapeyron-phase-kernel-reference/v1",
        "unexpected ten-state frozen-oracle schema")
    require(
        frozen["source"]["commit"] ==
            PINNED_CLAPEYRON_COMMIT,
        "ten-state diagnostic Clapeyron revision mismatch")

    states = ten_state_records(source, frozen)
    model_audit = audit_model(build_model())

    summaries = Any[]
    for precision_bits in TEN_STATE_PRECISIONS_BITS
        push!(summaries, run_ten_state_precision(
            states, precision_bits))
    end

    result = Dict(
        "schema" =>
            "MPMC_HNU/CPA/Clapeyron-ten-state-independent-compressed-diagnostic/v1",
        "source" => Dict(
            "software" =>
                "ClapeyronThermo/Clapeyron.jl",
            "commit" => PINNED_CLAPEYRON_COMMIT,
            "julia_version" => string(VERSION),
        ),
        "model_audit" => model_audit,
        "precision_bits" =>
            TEN_STATE_PRECISIONS_BITS,
        "independent_newton_contract" => Dict(
            "residual_target" =>
                NEWTON_RESIDUAL_TARGET,
            "initial_x" => ["0.5", "0.5"],
            "max_iterations" =>
                NEWTON_MAX_ITERS,
            "max_backtracks" =>
                NEWTON_MAX_BACKTRACKS,
        ),
        "derivative_probe" => Dict(
            "stencil" => "central_7_o6",
            "relative_h" =>
                TEN_STATE_DERIVATIVE_RELATIVE_H,
        ),
        "frozen_acceptance_envelope" => Dict(
            "pressure_pa" =>
                TEN_STATE_PRESSURE_THRESHOLD_PA,
            "z_only_ln_phi" =>
                TEN_STATE_LNPHI_THRESHOLD,
        ),
        "summaries" => summaries,
    )

    open(output_path, "w") do io
        JSON.print(io, result, 2)
        println(io)
    end

    println(
        "CLAPEYRON_TEN_STATE_COMPRESSED_DIAGNOSTIC_OK",
        " states=", length(states),
        " precisions=",
        join(TEN_STATE_PRECISIONS_BITS, ","),
        " residual_target=",
        NEWTON_RESIDUAL_TARGET,
        " derivative_r=",
        TEN_STATE_DERIVATIVE_RELATIVE_H,
    )

    for summary in summaries
        println(
            "CLAPEYRON_TEN_STATE_COMPRESSED_SUMMARY",
            " precision_bits=",
            summary["precision_bits"],
            " corrected_states=",
            summary["corrected_state_count"],
            " raw_states=",
            summary["uncorrected_raw_state_count"],
            " max_dP_target_pa=",
            summary[
                "max_abs_pressure_delta_vs_target_pa"],
            " max_dlnphi_target=",
            summary[
                "max_abs_z_only_ln_phi_delta_vs_target"],
            " result=",
            summary["diagnostic_pass"] ?
                "PASS" : "FAIL",
        )

        for record in summary["records"]
            println(
                "CLAPEYRON_TEN_STATE_COMPRESSED_ROW",
                " precision_bits=",
                summary["precision_bits"],
                " label=",
                replace(
                    record["label"],
                    " " => "_"),
                " path=", record["path"],
                " dP_target_pa=",
                record[
                    "abs_pressure_delta_vs_target_pa"],
                " dlnphi_target=",
                record[
                    "max_abs_z_only_ln_phi_delta_vs_target"],
                " p_pass=",
                record["pressure_threshold_pass"],
                " lnphi_pass=",
                record["ln_phi_threshold_pass"],
            )
        end
    end
end

if abspath(PROGRAM_FILE) == @__FILE__
    main_ten_state_diagnostic()
end
