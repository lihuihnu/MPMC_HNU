#!/usr/bin/env julia

using Clapeyron
using JSON

const PINNED_CLAPEYRON_COMMIT = "229b09452f36c2f812486150df0bb43b197bb4e5"
const COMPONENTS = ("MEOH", "H2O")
const NAMES = ("methanol", "water")
const R_MPMC = 8.31446261815324
const TEMPERATURE_K = 333.15

function require(condition::Bool, message::AbstractString)
    condition || error(message)
end

function require_close(actual::Real, expected::Real, message::AbstractString)
    a = Float64(actual)
    e = Float64(expected)
    scale = max(1.0, abs(a), abs(e))
    tolerance = 4096.0 * eps(Float64) * scale
    abs(a - e) <= tolerance ||
        error("$message: actual=$a expected=$e delta=$(abs(a-e)) tolerance=$tolerance")
    return abs(a - e)
end

function arg_value(flag::String)
    idx = findfirst(==(flag), ARGS)
    idx === nothing && error("missing required argument $flag")
    idx == length(ARGS) && error("missing value after $flag")
    return ARGS[idx + 1]
end

function assoc_value(param, model, i::Int, j::Int, site_i::String, site_j::String)
    ai = findfirst(==(site_i), model.sites.sites[i])
    aj = findfirst(==(site_j), model.sites.sites[j])
    ai === nothing && error("missing site $site_i on $(model.components[i])")
    aj === nothing && error("missing site $site_j on $(model.components[j])")
    return param.values[i, j][ai, aj]
end

function site_multiplicity(model, component::Int, site::String)
    index = findfirst(==(site), model.sites.sites[component])
    index === nothing && return 0
    return model.sites.n_sites[component][index]
end

function build_model()
    epsilon_methanol_k = 24591.0 / R_MPMC
    epsilon_water_k = 16655.0 / R_MPMC

    parameters = (;
        # Constructor metadata from the pinned Clapeyron molar-mass database.
        # Mw is not used by the CPA residual Helmholtz observables in this oracle.
        Mw = [32.042, 18.015],
        Tc = [512.6, 647.3],
        a = [0.40531, 0.12277],
        b = [3.0978e-5, 1.4515e-5],
        c1 = [0.43102, 0.67359],
        k = [0.0 -0.055; -0.055 0.0],
        n_H = [1, 2],
        n_e = [1, 2],
        epsilon_assoc = Dict(
            (("methanol", "e"), ("methanol", "H")) => epsilon_methanol_k,
            (("water", "e"), ("water", "H")) => epsilon_water_k,
            # Explicit CR-1 cross pairs. The pinned revision's native :cr1
            # recombination path fails to retain the mixed epsilon in the model,
            # so the source-complete CR-1 values are injected without modifying
            # Clapeyron itself.
            (("methanol", "H"), ("water", "e")) => 20623.0 / R_MPMC,
            (("methanol", "e"), ("water", "H")) => 20623.0 / R_MPMC,
        ),
        bondvol = Dict(
            (("methanol", "e"), ("methanol", "H")) => 0.0161,
            (("water", "e"), ("water", "H")) => 0.0692,
            (("methanol", "H"), ("water", "e")) => 0.03337843615270194,
            (("methanol", "e"), ("water", "H")) => 0.03337843615270194,
        ),
    )

    options = Clapeyron.AssocOptions(
        rtol = 1.0e-12,
        atol = 1.0e-12,
        max_iters = 512,
        dampingfactor = 0.5,
        # Explicit CR-1 pairs are already supplied above. :nocombining ensures
        # the pinned implementation cannot silently rewrite them.
        combining = :nocombining,
        implicit_ad = false,
    )

    return Clapeyron.sCPA(
        collect(NAMES);
        userlocations = parameters,
        assoc_options = options,
        verbose = false,
    )
end

function audit_model(model)
    max_delta = 0.0
    update(actual, expected, label) =
        (max_delta = max(max_delta, require_close(actual, expected, label)))

    require(model.components == collect(NAMES),
            "component order mismatch: $(model.components)")
    require(model.radial_dist == :KG,
            "Clapeyron model is not simplified/Kontogeorgis CPA")
    require(nameof(typeof(model.cubicmodel)) == :RK,
            "Clapeyron cubic model is not RK/SRK: $(typeof(model.cubicmodel))")
    require(nameof(typeof(model.cubicmodel.alpha)) in (:CPAAlpha, :sCPAAlpha),
            "Clapeyron alpha model is not CPA alpha: $(typeof(model.cubicmodel.alpha))")
    require(nameof(typeof(model.cubicmodel.mixing)) == :vdW1fRule,
            "Clapeyron mixing rule is not vdW1f: $(typeof(model.cubicmodel.mixing))")
    require(nameof(typeof(model.cubicmodel.translation)) == :NoTranslation,
            "Clapeyron volume translation is not disabled: $(typeof(model.cubicmodel.translation))")
    require(model.assoc_options.combining == :nocombining,
            "Clapeyron explicit CR-1 oracle must disable runtime combining")
    require(model.assoc_options.rtol == 1.0e-12 &&
            model.assoc_options.atol == 1.0e-12 &&
            model.assoc_options.dampingfactor == 0.5 &&
            model.assoc_options.max_iters == 512 &&
            !model.assoc_options.implicit_ad,
            "Clapeyron association numerical options drifted")

    rgas = Clapeyron.Rgas(model)
    require(rgas == R_MPMC,
            "gas constant mismatch: Clapeyron=$rgas MPMC=$R_MPMC")

    expected_mw = (32.042, 18.015)
    expected_tc = (512.6, 647.3)
    expected_a = (0.40531, 0.12277)
    expected_b = (3.0978e-5, 1.4515e-5)
    expected_c1 = (0.43102, 0.67359)
    for i in 1:2
        update(model.params.Mw.values[i], expected_mw[i], "Mw[$i]")
        update(model.params.Tc.values[i], expected_tc[i], "Tc[$i]")
        update(model.params.a.values[i, i], expected_a[i], "a0[$i]")
        update(model.params.b.values[i, i], expected_b[i], "b[$i]")
        update(model.params.c1.values[i], expected_c1[i], "c1[$i]")
    end

    expected_a12 = sqrt(expected_a[1] * expected_a[2]) * (1.0 - (-0.055))
    expected_b12 = 0.5 * (expected_b[1] + expected_b[2])
    update(model.params.a.values[1, 2], expected_a12, "mixed a12")
    update(model.params.b.values[1, 2], expected_b12, "mixed b12")
    derived_kij =
        1.0 - model.params.a.values[1, 2] /
        sqrt(model.params.a.values[1, 1] * model.params.a.values[2, 2])
    update(derived_kij, -0.055, "derived kij")

    require(site_multiplicity(model, 1, "H") == 1 &&
            site_multiplicity(model, 1, "e") == 1,
            "methanol site topology is not 2B (H=1,e=1)")
    require(site_multiplicity(model, 2, "H") == 2 &&
            site_multiplicity(model, 2, "e") == 2,
            "water site topology is not 4C (H=2,e=2)")
    require(length(model.sites.sites[1]) == 2 &&
            length(model.sites.sites[2]) == 2,
            "unexpected additional association site classes")

    epsilon = model.params.epsilon_assoc
    beta = model.params.bondvol
    expected_eps_j = [24591.0, 16655.0]
    expected_beta = [0.0161, 0.0692]
    for i in 1:2
        eps_k = assoc_value(epsilon, model, i, i, "e", "H")
        bet = assoc_value(beta, model, i, i, "e", "H")
        update(eps_k * rgas, expected_eps_j[i], "self epsilon J/mol[$i]")
        update(bet, expected_beta[i], "self beta[$i]")
    end

    expected_cross_eps_j = 20623.0
    expected_cross_beta = 0.03337843615270194
    for (i, j) in ((1, 2), (2, 1))
        eps_k = assoc_value(epsilon, model, i, j, "e", "H")
        bet = assoc_value(beta, model, i, j, "e", "H")
        update(eps_k * rgas, expected_cross_eps_j,
               "CR-1 cross epsilon J/mol[$i,$j]")
        update(bet, expected_cross_beta, "CR-1 cross beta[$i,$j]")
    end

    # Donor-donor and acceptor-acceptor interactions must remain absent.
    for i in 1:2, j in 1:2, sites in (("H", "H"), ("e", "e"))
        value_eps = assoc_value(epsilon, model, i, j, sites[1], sites[2])
        value_beta = assoc_value(beta, model, i, j, sites[1], sites[2])
        require(iszero(value_eps) && iszero(value_beta),
                "unexpected same-type association interaction at ($i,$j,$sites)")
    end

    return Dict(
        "status" => "exactly_representable",
        "max_abs_parameter_readback_delta" => max_delta,
        "component_order" => collect(NAMES),
        "radial_distribution" => "KG: g=1/(1-1.9eta)",
        "cubic" => string(nameof(typeof(model.cubicmodel))),
        "alpha" => string(nameof(typeof(model.cubicmodel.alpha))),
        "mixing" => string(nameof(typeof(model.cubicmodel.mixing))),
        "translation" => string(nameof(typeof(model.cubicmodel.translation))),
        "association_combining" =>
            "explicit CR-1 cross pairs; runtime combining=:nocombining",
        "gas_constant_j_per_mol_k" => rgas,
        "molar_mass_g_per_mol_constructor_metadata" => collect(expected_mw),
        "kij" => derived_kij,
        "methanol_sites" => Dict("H" => 1, "e" => 1),
        "water_sites" => Dict("H" => 2, "e" => 2),
        "cr1_cross_epsilon_j_per_mol" => expected_cross_eps_j,
        "cr1_cross_beta" => expected_cross_beta,
        "pinned_native_cr1_note" =>
            "Pinned assoc_mix! does not retain epsilon_assoc_mix return in recombine_assoc!; explicit source-complete CR-1 pairs are therefore injected.",
    )
end

function state_result(model, temperature_k::Float64, rho::Float64,
                      methanol_fraction::Float64, label::String)
    z = [methanol_fraction, 1.0 - methanol_fraction]
    require(abs(sum(z) - 1.0) <= 16eps(Float64),
            "$label composition does not sum to one")
    volume = 1.0 / rho
    rgas = Clapeyron.Rgas(model)
    rt = rgas * temperature_k

    a_res_j = Clapeyron.eos_res(model, volume, temperature_k, z)
    f_res = a_res_j / rt
    f_res_from_reduced = sum(z) * Clapeyron.a_res(model, volume, temperature_k, z)
    require_close(f_res, f_res_from_reduced,
                  "$label residual Helmholtz API self-consistency")

    pressure_pa = Clapeyron.pressure(model, volume, temperature_k, z)
    require(isfinite(pressure_pa) && pressure_pa > 0.0,
            "$label Clapeyron pressure must be positive and finite")

    mu_j_per_mol =
        Clapeyron.VT_chemical_potential_res(model, volume, temperature_k, z)
    require(length(mu_j_per_mol) == 2 &&
            all(isfinite, mu_j_per_mol),
            "$label invalid residual chemical potential")
    mu_over_rt = [Float64(value / rt) for value in mu_j_per_mol]

    zfactor = pressure_pa * volume / (rt * sum(z))
    require(isfinite(zfactor) && zfactor > 0.0,
            "$label invalid compressibility factor")
    log_z = log(zfactor)
    ln_phi = [value - log_z for value in mu_over_rt]

    phi = Clapeyron.VT_fugacity_coefficient(model, volume, temperature_k, z)
    require(length(phi) == 2 && all(v -> isfinite(v) && v > 0.0, phi),
            "$label invalid Clapeyron fugacity coefficient")
    for i in 1:2
        require_close(log(phi[i]), ln_phi[i],
                      "$label fugacity API self-consistency component $i")
    end

    x = Clapeyron.assoc_fractions(model, volume, temperature_k, z)
    require(all(v -> isfinite(v) && 0.0 < v <= 1.0, x.v),
            "$label invalid association site fraction")

    return Dict(
        "label" => label,
        "temperature_k" => temperature_k,
        "molar_density_mol_per_m3" => rho,
        "molar_volume_m3_per_mol" => volume,
        "composition" => Dict("MEOH" => z[1], "H2O" => z[2]),
        "f_res_reduced_extensive" => Float64(f_res),
        "pressure_pa" => Float64(pressure_pa),
        "compressibility_factor" => Float64(zfactor),
        "mu_residual_over_rt" => Dict(
            "MEOH" => mu_over_rt[1],
            "H2O" => mu_over_rt[2],
        ),
        "ln_phi" => Dict(
            "MEOH" => Float64(ln_phi[1]),
            "H2O" => Float64(ln_phi[2]),
        ),
        "association_site_fractions" => [Float64(v) for v in x.v],
    )
end

function generate(states_path::String)
    source = JSON.parsefile(states_path)
    require(source["schema"] ==
            "MPMC_HNU/CPA/ThermoPack-phase-kernel-reference/v2",
            "unexpected frozen phase-state source schema")

    model = build_model()
    audit = audit_model(model)

    states = Any[]
    for (state_index, source_state) in enumerate(source["states"])
        temperature_k = Float64(source_state["temperature_k"])
        require_close(temperature_k, TEMPERATURE_K,
                      "frozen state temperature $state_index")
        for phase in ("liquid", "vapor")
            record = source_state[phase]
            rho = Float64(record["molar_density_mol_per_m3"])
            methanol = Float64(record["composition"]["MEOH"])
            push!(states, state_result(
                model, temperature_k, rho, methanol,
                "state=$(state_index-1) phase=$phase"))
        end
    end
    require(length(states) == 10, "Gate-E oracle must contain exactly ten states")

    commit_from_env = get(ENV, "CLAPEYRON_COMMIT", "")
    require(commit_from_env == PINNED_CLAPEYRON_COMMIT,
            "workflow Clapeyron revision mismatch: $commit_from_env")

    return Dict(
        "schema" => "MPMC_HNU/CPA/Clapeyron-phase-kernel-reference/v1",
        "source" => Dict(
            "software" => "ClapeyronThermo/Clapeyron.jl",
            "commit" => PINNED_CLAPEYRON_COMMIT,
            "unmodified_source_required" => true,
            "julia_version" => string(VERSION),
            "state_coordinates_from" =>
                "ThermoPack frozen phase-kernel oracle; only T,V,n coordinates reused",
        ),
        "model" => Dict(
            "profile" =>
                "CPA/SRK-physical/simplified-rdf-1.9eta/explicit-site-pairs/v1",
            "parameter_snapshot" =>
                "thermopack-parity Tc + MPMC literature a0/b/c1/association + Folas kij=-0.055",
            "audit" => audit,
        ),
        "gas_constant_j_per_mol_k" => Clapeyron.Rgas(model),
        "states" => states,
    )
end

function main()
    states_path = arg_value("--states")
    output_path = arg_value("--out")
    result = generate(states_path)

    open(output_path, "w") do io
        JSON.print(io, result, 2)
        println(io)
    end

    audit = result["model"]["audit"]
    println(
        "CLAPEYRON_CPA_FORMULATION_AUDIT_OK",
        " commit=", PINNED_CLAPEYRON_COMMIT,
        " profile=", result["model"]["profile"],
        " R=", result["gas_constant_j_per_mol_k"],
        " kij=", audit["kij"],
        " max_abs_parameter_readback_delta=",
        audit["max_abs_parameter_readback_delta"],
    )
    println("CLAPEYRON_CPA_ORACLE_OK states=", length(result["states"]))
end

main()
