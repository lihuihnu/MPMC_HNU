#!/usr/bin/env julia

using LinearAlgebra

include(joinpath(@__DIR__, "association_stationarity_floor_audit.jl"))

const EXACT2_PRECISIONS_BITS = Int[256, 512]
const NEWTON_RESIDUAL_TARGET = "1e-30"
const NEWTON_MAX_ITERS = 100
const NEWTON_MAX_BACKTRACKS = 40
const EXACT2_DERIVATIVE_RELATIVE_H = "1e-6"
const EXPECTED_MAPPING = Int[1, 1, 2, 2]

function residual_vector(matrix, x)
    return x .* (1 .+ matrix * x) .- 1
end

function residual_norm(matrix, x)
    return maximum(abs, residual_vector(matrix, x))
end

function analytic_jacobian(matrix, x)
    kx = matrix * x
    jacobian = Matrix{eltype(x)}(undef, length(x), length(x))
    for i in eachindex(x), j in eachindex(x)
        value = x[i] * matrix[i, j]
        if i == j
            value += 1 + kx[i]
        end
        jacobian[i, j] = value
    end
    return jacobian
end

function independent_newton(matrix)
    require(size(matrix) == (2, 2),
            "independent exact2 audit requires a 2x2 reduced matrix")

    target = parse(BigFloat, NEWTON_RESIDUAL_TARGET)
    x = BigFloat[0.5, 0.5]
    history = Any[]

    for iteration in 0:NEWTON_MAX_ITERS
        residual = residual_vector(matrix, x)
        norm_value = maximum(abs, residual)
        push!(history, Dict(
            "iteration" => iteration,
            "residual" => [big_string(value) for value in residual],
            "max_abs_residual" => big_string(norm_value),
            "x" => [big_string(value) for value in x],
        ))

        if norm_value <= target
            return (
                x = x,
                residual = residual,
                max_abs_residual = norm_value,
                iterations = iteration,
                history = history,
            )
        end

        iteration == NEWTON_MAX_ITERS &&
            error("independent Newton did not hit residual target")

        jacobian = analytic_jacobian(matrix, x)
        step = jacobian \ (-residual)

        accepted = false
        lambda = one(BigFloat)
        for backtrack in 0:NEWTON_MAX_BACKTRACKS
            candidate = x .+ lambda .* step
            if all(value -> isfinite(value) && 0 < value <= 1, candidate)
                candidate_norm = residual_norm(matrix, candidate)
                if isfinite(candidate_norm) && candidate_norm < norm_value
                    x = candidate
                    accepted = true
                    push!(history, Dict(
                        "accepted_backtrack" => backtrack,
                        "accepted_lambda" => big_string(lambda),
                    ))
                    break
                end
            end
            lambda *= BigFloat(0.5)
        end
        accepted ||
            error("independent Newton line search failed")
    end

    error("unreachable independent Newton exit")
end

function extract_association_problem(model, volume, temperature, z)
    delta = Clapeyron.delta_assoc(
        model, volume, temperature, z, nothing)
    original = Clapeyron.assoc_site_matrix(
        model, volume, temperature, z, nothing, delta)

    require(size(original) == (4, 4),
            "state-4 original association matrix is not 4x4")
    require(Clapeyron.__maybe_compress(original),
            "state-4 association matrix is no longer compressible")

    reduced, mapping =
        Clapeyron.compress_assoc_matrix(original)
    mapping_vector = collect(mapping)

    require(size(reduced) == (2, 2),
            "state-4 reduced association matrix is not 2x2")
    require(mapping_vector == EXPECTED_MAPPING,
            "state-4 compression mapping drifted")

    return (
        delta = delta,
        original = original,
        reduced = reduced,
        mapping = mapping_vector,
    )
end

function site_weights(model, z)
    sites = Clapeyron.getsites(model)
    offsets = sites.n_sites.p
    multiplicities = sites.n_sites.v

    weights = Vector{BigFloat}(undef, length(multiplicities))
    for component in eachindex(z)
        for flat_index in offsets[component]:(offsets[component + 1] - 1)
            weights[flat_index] =
                z[component] * BigFloat(multiplicities[flat_index])
        end
    end
    return weights
end

function independent_q_from_original_matrix(model, original, z, x)
    require(all(iszero, diag(original)),
            "independent Q identity assumes no same-site diagonal association")
    weights = site_weights(model, z)
    kx = original * x

    q2 = zero(BigFloat)
    for i in eachindex(x)
        q2 += weights[i] * (log(x[i]) + 1 - x[i])
    end

    q1 = -BigFloat(0.5) * sum(weights .* x .* kx)
    q_stationary = q1 + q2

    q_direct = zero(BigFloat)
    for i in eachindex(x)
        q_direct += weights[i] *
            (log(x[i]) - BigFloat(0.5) * x[i] + BigFloat(0.5))
    end

    return (
        q1 = q1,
        q2 = q2,
        q_stationary = q_stationary,
        q_direct = q_direct,
        q_direct_delta = abs(q_stationary - q_direct),
    )
end

function solve_independent_state(model, volume, temperature, z)
    problem = extract_association_problem(
        model, volume, temperature, z)

    solution = independent_newton(problem.reduced)
    expanded = BigFloat[
        solution.x[index] for index in problem.mapping
    ]
    expanded_residual =
        residual_vector(problem.original, expanded)
    expanded_norm = maximum(abs, expanded_residual)

    pinned = Clapeyron.assoc_fractions(
        model, volume, temperature, z)
    pinned_x = collect(pinned.v)
    require(eltype(pinned_x) == BigFloat && all(isfinite, pinned_x),
            "pinned exact2 site fractions are not finite BigFloat values")

    q_independent = independent_q_from_original_matrix(
        model, problem.original, z, expanded)
    q_pinned_formula = independent_q_from_original_matrix(
        model, problem.original, z, pinned_x)

    total_moles = sum(z)
    rgas = BigFloat(Clapeyron.Rgas(model))
    rt = rgas * temperature

    a_assoc_independent =
        q_independent.q_stationary / total_moles
    A_assoc_independent =
        rt * q_independent.q_stationary

    pinned_a_assoc =
        Clapeyron.a_assoc(model, volume, temperature, z)
    pinned_q = pinned_a_assoc * total_moles
    pinned_A_assoc = rt * pinned_q

    cubic_A = Clapeyron.eos_res(
        model.cubicmodel, volume, temperature, z)
    total_A_independent = cubic_A + A_assoc_independent
    total_A_pinned = Clapeyron.eos_res(
        model, volume, temperature, z)

    return Dict(
        "reduced_matrix" => [
            [big_string(problem.reduced[i, j]) for j in 1:2]
            for i in 1:2
        ],
        "mapping" => problem.mapping,
        "newton_iterations" => solution.iterations,
        "reduced_x" => [big_string(value) for value in solution.x],
        "reduced_residual" =>
            [big_string(value) for value in solution.residual],
        "reduced_max_abs_residual" =>
            big_string(solution.max_abs_residual),
        "expanded_x" => [big_string(value) for value in expanded],
        "expanded_residual" =>
            [big_string(value) for value in expanded_residual],
        "expanded_max_abs_residual" => big_string(expanded_norm),
        "pinned_x" => [big_string(value) for value in pinned_x],
        "max_abs_x_delta_vs_pinned" =>
            big_string(maximum(abs.(expanded .- pinned_x))),
        "q_independent" => Dict(
            "q1" => big_string(q_independent.q1),
            "q2" => big_string(q_independent.q2),
            "q_stationary" =>
                big_string(q_independent.q_stationary),
            "q_direct" => big_string(q_independent.q_direct),
            "abs_q_stationary_minus_direct" =>
                big_string(q_independent.q_direct_delta),
        ),
        "q_pinned_x_via_independent_formula" => Dict(
            "q_stationary" =>
                big_string(q_pinned_formula.q_stationary),
            "q_direct" =>
                big_string(q_pinned_formula.q_direct),
            "abs_q_stationary_minus_direct" =>
                big_string(q_pinned_formula.q_direct_delta),
        ),
        "pinned_q_from_clapeyron_a_assoc" =>
            big_string(pinned_q),
        "abs_independent_formula_pinned_q_delta" =>
            big_string(abs(
                q_pinned_formula.q_stationary - pinned_q)),
        "a_assoc_independent" =>
            big_string(a_assoc_independent),
        "A_assoc_independent_joule" =>
            big_string(A_assoc_independent),
        "A_assoc_pinned_joule" =>
            big_string(pinned_A_assoc),
        "abs_A_assoc_delta_vs_pinned_joule" =>
            big_string(abs(
                A_assoc_independent - pinned_A_assoc)),
        "A_cubic_joule" => big_string(cubic_A),
        "A_res_independent_joule" =>
            big_string(total_A_independent),
        "A_res_pinned_joule" =>
            big_string(total_A_pinned),
        "abs_A_res_delta_vs_pinned_joule" =>
            big_string(abs(
                total_A_independent - total_A_pinned)),
    )
end

function seven_point(values, h)
    return (
        -values[-3] + 9 * values[-2] - 45 * values[-1] +
        45 * values[1] - 9 * values[2] + values[3]
    ) / (60 * h)
end

function derivative_audit(model, volume, temperature, z,
                          target_pressure, internal_pressure)
    relative_h = parse(BigFloat, EXACT2_DERIVATIVE_RELATIVE_H)
    h = relative_h * volume

    independent_assoc = Dict{Int,BigFloat}()
    pinned_assoc = Dict{Int,BigFloat}()
    independent_total = Dict{Int,BigFloat}()
    pinned_total = Dict{Int,BigFloat}()
    shifted_residuals = Dict{String,String}()

    for multiplier in (-3, -2, -1, 1, 2, 3)
        shifted_volume =
            volume + BigFloat(multiplier) * h
        state = solve_independent_state(
            model, shifted_volume, temperature, z)

        independent_assoc[multiplier] =
            parse(BigFloat, state["A_assoc_independent_joule"])
        pinned_assoc[multiplier] =
            parse(BigFloat, state["A_assoc_pinned_joule"])
        independent_total[multiplier] =
            parse(BigFloat, state["A_res_independent_joule"])
        pinned_total[multiplier] =
            parse(BigFloat, state["A_res_pinned_joule"])
        shifted_residuals[string(multiplier)] =
            state["reduced_max_abs_residual"]
    end

    d_assoc_independent =
        seven_point(independent_assoc, h)
    d_assoc_pinned =
        seven_point(pinned_assoc, h)
    d_total_independent =
        seven_point(independent_total, h)
    d_total_pinned =
        seven_point(pinned_total, h)

    total_moles = sum(z)
    rt = BigFloat(Clapeyron.Rgas(model)) * temperature
    pressure_independent =
        total_moles * rt / volume - d_total_independent
    pressure_pinned_fd =
        total_moles * rt / volume - d_total_pinned

    return Dict(
        "stencil" => "central_7_o6",
        "relative_h" => EXACT2_DERIVATIVE_RELATIVE_H,
        "shifted_reduced_residuals" => shifted_residuals,
        "dA_assoc_independent_dV_pa" =>
            big_string(d_assoc_independent),
        "dA_assoc_pinned_dV_pa" =>
            big_string(d_assoc_pinned),
        "abs_dA_assoc_delta_pa" =>
            big_string(abs(
                d_assoc_independent - d_assoc_pinned)),
        "dA_res_independent_dV_pa" =>
            big_string(d_total_independent),
        "dA_res_pinned_dV_pa" =>
            big_string(d_total_pinned),
        "abs_dA_res_delta_pa" =>
            big_string(abs(
                d_total_independent - d_total_pinned)),
        "pressure_independent_pa" =>
            big_string(pressure_independent),
        "pressure_pinned_fd_pa" =>
            big_string(pressure_pinned_fd),
        "abs_pressure_independent_vs_target_pa" =>
            Float64(abs(
                pressure_independent - target_pressure)),
        "abs_pressure_independent_vs_internal_pa" =>
            Float64(abs(
                pressure_independent - internal_pressure)),
        "abs_pressure_pinned_fd_vs_target_pa" =>
            Float64(abs(
                pressure_pinned_fd - target_pressure)),
        "abs_pressure_independent_minus_pinned_fd_pa" =>
            big_string(abs(
                pressure_independent - pressure_pinned_fd)),
    )
end

function run_precision(state, precision_bits)
    return setprecision(BigFloat, precision_bits) do
        model = build_model()
        temperature = BigFloat(state.temperature_k)
        z = BigFloat[state.x_meoh, 1.0 - state.x_meoh]
        volume = sum(z) / BigFloat(state.rho)
        target_pressure = BigFloat(state.target_pressure_pa)
        internal_pressure =
            BigFloat(state.frozen_clapeyron_pressure_pa)

        center = solve_independent_state(
            model, volume, temperature, z)

        require(
            parse(
                BigFloat,
                center["reduced_max_abs_residual"]) <=
                parse(BigFloat, NEWTON_RESIDUAL_TARGET),
            "independent Newton missed frozen residual target")

        derivative = derivative_audit(
            model,
            volume,
            temperature,
            z,
            target_pressure,
            internal_pressure,
        )

        return Dict(
            "precision_bits" => precision_bits,
            "center" => center,
            "derivative" => derivative,
        )
    end
end

function main_exact2_independent()
    states_path = arg_value("--states")
    oracle_path = arg_value("--oracle")
    output_path = arg_value("--out")

    source = JSON.parsefile(states_path)
    frozen = JSON.parsefile(oracle_path)
    require(
        source["schema"] ==
            "MPMC_HNU/CPA/ThermoPack-phase-kernel-reference/v2",
        "unexpected exact2 coordinate-source schema")
    require(
        frozen["schema"] ==
            "MPMC_HNU/CPA/Clapeyron-phase-kernel-reference/v1",
        "unexpected exact2 frozen-oracle schema")
    require(
        frozen["source"]["commit"] == PINNED_CLAPEYRON_COMMIT,
        "exact2 audit Clapeyron revision mismatch")

    state = state4_liquid(source, frozen)
    model_audit = audit_model(build_model())

    records = Any[]
    for precision_bits in EXACT2_PRECISIONS_BITS
        push!(records, run_precision(
            state, precision_bits))
    end

    result = Dict(
        "schema" =>
            "MPMC_HNU/CPA/Clapeyron-independent-exact2-audit/v1",
        "source" => Dict(
            "software" =>
                "ClapeyronThermo/Clapeyron.jl",
            "commit" => PINNED_CLAPEYRON_COMMIT,
            "julia_version" => string(VERSION),
        ),
        "state" => Dict(
            "label" => TARGET_LABEL,
            "temperature_k" => state.temperature_k,
            "molar_density_mol_per_m3" => state.rho,
            "composition_methanol" => state.x_meoh,
            "target_pressure_pa" => state.target_pressure_pa,
            "clapeyron_internal_pressure_pa" =>
                state.frozen_clapeyron_pressure_pa,
        ),
        "model_audit" => model_audit,
        "newton_contract" => Dict(
            "initial_x" => ["0.5", "0.5"],
            "residual_target" => NEWTON_RESIDUAL_TARGET,
            "max_iterations" => NEWTON_MAX_ITERS,
            "max_backtracks" => NEWTON_MAX_BACKTRACKS,
            "backtrack_factor" => "0.5",
            "allowed_interval" => "0 < x_i <= 1",
        ),
        "precision_bits" => EXACT2_PRECISIONS_BITS,
        "derivative_probe" => Dict(
            "stencil" => "central_7_o6",
            "relative_h" => EXACT2_DERIVATIVE_RELATIVE_H,
        ),
        "records" => records,
    )

    open(output_path, "w") do io
        JSON.print(io, result, 2)
        println(io)
    end

    println(
        "CLAPEYRON_INDEPENDENT_EXACT2_AUDIT_OK",
        " state=", replace(TARGET_LABEL, " " => "_"),
        " precisions=", join(EXACT2_PRECISIONS_BITS, ","),
        " residual_target=", NEWTON_RESIDUAL_TARGET,
        " derivative_r=", EXACT2_DERIVATIVE_RELATIVE_H,
    )

    for record in records
        center = record["center"]
        derivative = record["derivative"]
        println(
            "CLAPEYRON_INDEPENDENT_EXACT2_ROW",
            " precision_bits=", record["precision_bits"],
            " reduced_residual=",
            center["reduced_max_abs_residual"],
            " expanded_residual=",
            center["expanded_max_abs_residual"],
            " max_dX_pinned=",
            center["max_abs_x_delta_vs_pinned"],
            " dQ_direct=",
            center["q_independent"][
                "abs_q_stationary_minus_direct"],
            " dAassoc_pinned_j=",
            center["abs_A_assoc_delta_vs_pinned_joule"],
            " dAres_pinned_j=",
            center["abs_A_res_delta_vs_pinned_joule"],
            " dDassoc_pa=",
            derivative["abs_dA_assoc_delta_pa"],
            " dDres_pa=",
            derivative["abs_dA_res_delta_pa"],
            " dP_target_pa=",
            derivative[
                "abs_pressure_independent_vs_target_pa"],
            " dP_internal_pa=",
            derivative[
                "abs_pressure_independent_vs_internal_pa"],
            " dP_pinnedfd_pa=",
            derivative[
                "abs_pressure_independent_minus_pinned_fd_pa"],
        )
    end
end

if abspath(PROGRAM_FILE) == @__FILE__
    main_exact2_independent()
end
