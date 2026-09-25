#include <mpmc/flow/thermodynamics_absent_phase_extension.hpp>
#include <mpmc/flow/single_phase_natural_variable.hpp>

#include <test_support.hpp>

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <iostream>
#include <optional>
#include <span>
#include <stdexcept>
#include <string>
#include <string_view>
#include <utility>
#include <vector>

namespace {

namespace flow = mpmc::flow;
namespace th = mpmc::thermodynamics;
namespace swt = sw92_test;

void require_provider(
    bool condition,
    std::string_view message) {
    if (!condition) {
        throw std::runtime_error(
            std::string{message});
    }
}

th::Provenance synthetic_source(
    std::string locator) {
    return {
        th::SourceKind::synthetic_test,
        "MPMC_HNU absent-phase thermodynamic provider regression",
        "v1",
        std::move(locator),
        "Manufactured parameters exercise branch binding and analytic density propagation only",
        "tests/flow/core/thermodynamics_absent_phase_extension_test.cpp",
        "Repository structural regression"};
}

th::SourcedScalar scalar(
    double value,
    th::Unit unit,
    std::string locator) {
    return {
        value,
        unit,
        synthetic_source(
            std::move(locator)),
        "SI",
        "identity"};
}

th::Component component(
    std::string id,
    double molar_mass_kg_per_mol) {
    const auto definition =
        synthetic_source(
            "component-" + id);
    return {
        id,
        id,
        th::ComponentKind::pure,
        definition,
        scalar(
            molar_mass_kg_per_mol,
            th::Unit::kilogram_per_mole,
            "molar-mass-" + id)};
}

th::PrParameterSet pr_parameters() {
    std::vector<th::Component> catalog{
        component("A", 0.020),
        component("B", 0.040)};
    const std::vector<std::string>
        order{"A", "B"};
    th::PrParameterInput input;
    input.model_id =
        std::string{th::pr76_profile};
    input.dataset_id =
        "synthetic-absent-phase-pr76";
    input.revision = "v1";
    input.applicability = {
        std::nullopt,
        std::nullopt,
        synthetic_source(
            "pr-applicability")};
    input.pure.push_back({
        "A",
        scalar(
            400.0,
            th::Unit::kelvin,
            "A-Tc"),
        scalar(
            4.0e6,
            th::Unit::pascal,
            "A-Pc"),
        scalar(
            0.2,
            th::Unit::dimensionless,
            "A-omega")});
    input.pure.push_back({
        "B",
        scalar(
            500.0,
            th::Unit::kelvin,
            "B-Tc"),
        scalar(
            5.0e6,
            th::Unit::pascal,
            "B-Pc"),
        scalar(
            0.1,
            th::Unit::dimensionless,
            "B-omega")});
    input.binary.push_back({
        "A",
        "B",
        scalar(
            0.0,
            th::Unit::dimensionless,
            "A-B-kij")});
    return th::PrParameterSet::create(
        catalog,
        order,
        input,
        th::DataPolicy::
            allow_synthetic_tests);
}

th::CpaParameterSet cpa_parameters() {
    std::vector<th::Component> catalog{
        component("A", 0.044),
        component("B", 0.030)};
    const std::vector<std::string>
        order{"A", "B"};
    th::CpaParameterInput input;
    input.dataset_id =
        "synthetic-absent-phase-cpa";
    input.revision = "v1";
    input.applicability = {
        std::nullopt,
        std::nullopt,
        synthetic_source(
            "cpa-applicability")};
    input.pure.push_back({
        "A",
        {500.0, synthetic_source("CPA-Tc-A")},
        {0.25, synthetic_source("CPA-a0-A")},
        {4.0e-5, synthetic_source("CPA-b-A")},
        {0.7, synthetic_source("CPA-c1-A")},
        {}});
    input.pure.push_back({
        "B",
        {400.0, synthetic_source("CPA-Tc-B")},
        {0.18, synthetic_source("CPA-a0-B")},
        {5.0e-5, synthetic_source("CPA-b-B")},
        {0.5, synthetic_source("CPA-c1-B")},
        {}});
    input.binary.push_back({
        "A",
        "B",
        {0.04, synthetic_source("CPA-kij-A-B")}});
    return th::CpaParameterSet::create(
        catalog,
        order,
        input,
        th::DataPolicy::
            allow_synthetic_tests);
}

flow::FrozenPhysicalPhaseIdentity
phase_id(
    std::string key) {
    return {
        "fixture/phase-continuation",
        std::move(key)};
}

flow::PhaseIdentityContinuationSnapshot
continuation() {
    flow::PhaseSetTransitionCandidate
        candidate;
    candidate.source_phase_count = 1U;
    candidate.target_phase_count = 2U;
    candidate.status =
        flow::
            PhaseSetTransitionCandidateStatus::
                target_resolved;
    candidate.trigger =
        flow::
            PhaseSetTransitionTrigger::
                stability_witness;
    candidate.evidence_profile =
        "fixture/fresh-target-resolve/v1";
    candidate.diagnostic =
        "explicit test continuation";

    return flow::
        make_phase_identity_continuation_snapshot(
            candidate,
            flow::FrozenActivePhaseIdentityMap{
                {phase_id("phase-a")}},
            flow::FrozenActivePhaseIdentityMap{
                {
                    phase_id("phase-a"),
                    phase_id("phase-b")}});
}

flow::NaturalVariableStateIdentity3P
host_identity(
    std::vector<std::string> component_ids,
    double pressure_pa,
    double temperature_k,
    std::vector<double> composition) {
    flow::NaturalVariableLayout1P
        layout{
            component_ids.size()};
    return {
        layout.descriptor(),
        std::move(component_ids),
        pressure_pa,
        temperature_k,
        {1.0, 0.0, 0.0},
        {
            std::move(composition),
            {},
            {}}};
}

flow::AbsentPhaseThermodynamicCoordinateExtension
coordinates(
    flow::FrozenPhysicalPhaseIdentity identity,
    flow::NaturalVariableStateIdentity3P host,
    std::vector<double> composition) {
    const std::size_t q =
        host.layout.unknown_count();
    const std::size_t n =
        host.layout.component_count();
    flow::AbsentPhaseThermodynamicCoordinateExtension
        result;
    result.identity =
        std::move(identity);
    result.host_state_identity =
        std::move(host);
    result.phase_pressure_pa =
        result.host_state_identity
            .reference_pressure_pa;
    result.phase_pressure_gradient.assign(
        q,
        0.0);
    result.phase_pressure_gradient[
        result.host_state_identity
            .layout.pressure_unknown_index()] =
        1.0;
    result.hypothetical_composition =
        std::move(composition);
    result.hypothetical_composition_jacobian
        .assign(
            n * q,
            0.0);

    if (n > 1U) {
        const flow::NaturalVariableLayout1P
            layout{
                n};
        for (std::size_t component = 0U;
             component < n;
             ++component) {
            const auto column =
                layout
                    .independent_composition_unknown_index(
                        component);
            if (!column) {
                continue;
            }
            result.hypothetical_composition_jacobian[
                component * q +
                *column] = 1.0;
            result.hypothetical_composition_jacobian[
                layout
                        .dependent_composition_component() *
                    q +
                *column] = -1.0;
        }
    }
    result.provenance =
        "fixture/hypothetical-composition-continuation/v1";
    result.validate();
    return result;
}

void require_extension(
    const flow::
        AbsentPhasePotentialExtensionLinearization&
            extension,
    const flow::
        AbsentPhaseThermodynamicCoordinateExtension&
            source,
    double expected_mass_density) {
    require_provider(
        extension.identity() ==
                source.identity &&
            extension.phase_pressure_pa() ==
                source.phase_pressure_pa &&
            extension.input_count() ==
                source.input_count() &&
            extension
                    .phase_pressure_gradient()
                    .size() ==
                source.input_count() &&
            extension
                    .mass_density_gradient()
                    .size() ==
                source.input_count() &&
            std::isfinite(
                extension
                    .mass_density_kg_per_m3()) &&
            extension
                    .mass_density_kg_per_m3() >
                0.0,
        "absent-phase provider result metadata mismatch");
    require_provider(
        std::abs(
            extension
                .mass_density_kg_per_m3() -
            expected_mass_density) <=
            1.0e-11 *
                std::max(
                    1.0,
                    std::abs(
                        expected_mass_density)),
        "absent-phase provider mass-density primal mismatch");
    for (double derivative :
         extension.mass_density_gradient()) {
        require_provider(
            std::isfinite(derivative),
            "absent-phase provider published non-finite mass-density derivative");
    }
    require_provider(
        extension.provenance().find(
            "fixture/fresh-target-resolve/v1") !=
            std::string_view::npos,
        "transition evidence profile was not retained in provider provenance");
}

void pr76_provider_contract() {
    const auto parameters =
        pr_parameters();
    const auto model =
        th::Pr76Phase<double>::
            from_parameters(parameters);
    constexpr double pressure = 1.0e6;
    constexpr double temperature = 450.0;
    const std::vector<double>
        x{0.4, 0.6};
    th::Pr76PhaseWorkspace<double>
        roots_workspace;
    const auto roots =
        model.roots_full(
            pressure,
            temperature,
            x,
            roots_workspace);
    require_provider(
        roots.status ==
                th::Pr76RootStatus::success &&
            roots.count > 0U,
        "PR76 provider fixture has no selected branch");
    const std::size_t selected_root =
        roots.count - 1U;

    const auto cont =
        continuation();
    auto registry =
        flow::
            make_transition_selected_phase_branch_registry(
                cont,
                std::vector<
                    flow::
                        FrozenSelectedPhaseBranchBinding<
                            th::Pr76SelectedPhase>>{
                    {
                        phase_id("phase-a"),
                        {0U, {}},
                        "fixture/source-root"},
                    {
                        phase_id("phase-b"),
                        {selected_root, {}},
                        "fixture/target-root"}});

    flow::
        Pr76AbsentPhasePotentialExtensionProvider<double>
        provider{
            model,
            std::move(registry)};
    const auto coordinate =
        coordinates(
            phase_id("phase-b"),
            host_identity(
                {"A", "B"},
                pressure,
                temperature,
                x),
            x);
    const auto extension =
        provider.evaluate(
            coordinate);

    th::Pr76PhaseWorkspace<double>
        density_workspace;
    const auto molar_density =
        th::evaluate_selected_phase_molar_density(
            model,
            pressure,
            temperature,
            std::span<const double>{x},
            th::Pr76SelectedPhase{
                selected_root, {}},
            density_workspace)
            .molar_density_mol_per_m3;
    const double mixture_molar_mass =
        0.4 * 0.020 +
        0.6 * 0.040;
    require_extension(
        extension,
        coordinate,
        molar_density *
            mixture_molar_mass);
}

void sw92_provider_contract() {
    auto prepared =
        swt::binary_input(
            swt::co2);
    prepared.catalog[0].molar_mass =
        swt::scalar(
            0.04401,
            th::Unit::kilogram_per_mole,
            swt::synthetic(
                "CO2 molar mass for provider test"));
    prepared.catalog[1].molar_mass =
        swt::scalar(
            0.018015,
            th::Unit::kilogram_per_mole,
            swt::synthetic(
                "water molar mass for provider test"));
    const auto parameters =
        th::Sw92ParameterSet::create(
            prepared.catalog,
            prepared.order,
            prepared.input,
            th::DataPolicy::
                allow_synthetic_tests);
    const auto model =
        th::Sw92Phase<double>::
            from_parameters(parameters);

    constexpr double pressure = 3.0e6;
    constexpr double temperature = 340.0;
    constexpr double molality = 0.0;
    const std::vector<double>
        x{0.7, 0.3};
    th::Sw92PhaseWorkspace<double>
        root_workspace;
    const auto roots =
        model.roots(
            pressure,
            temperature,
            x,
            molality,
            th::SwPhaseFamily::
                nonaqueous,
            root_workspace);
    require_provider(
        roots.status ==
                th::Sw92RootStatus::success &&
            roots.count == 3U,
        "SW92 provider fixture lost nonaqueous root topology");

    const auto cont =
        continuation();
    auto registry =
        flow::
            make_transition_selected_phase_branch_registry(
                cont,
                std::vector<
                    flow::
                        FrozenSelectedPhaseBranchBinding<
                            th::Sw92SelectedPhase<double>>>{
                    {
                        phase_id("phase-a"),
                        {
                            molality,
                            th::SwPhaseFamily::aqueous,
                            0U,
                            {}},
                        "fixture/aqueous-family"},
                    {
                        phase_id("phase-b"),
                        {
                            molality,
                            th::SwPhaseFamily::nonaqueous,
                            2U,
                            {}},
                        "fixture/nonaqueous-family-root2"}});

    flow::
        Sw92AbsentPhasePotentialExtensionProvider<double>
        provider{
            model,
            std::move(registry)};
    const auto coordinate =
        coordinates(
            phase_id("phase-b"),
            host_identity(
                {
                    "carbon-dioxide",
                    "water"},
                pressure,
                temperature,
                x),
            x);
    const auto extension =
        provider.evaluate(
            coordinate);

    th::Sw92PhaseWorkspace<double>
        density_workspace;
    const auto molar_density =
        th::evaluate_selected_phase_molar_density(
            model,
            pressure,
            temperature,
            std::span<const double>{x},
            th::Sw92SelectedPhase<double>{
                molality,
                th::SwPhaseFamily::nonaqueous,
                2U,
                {}},
            density_workspace)
            .molar_density_mol_per_m3;
    require_extension(
        extension,
        coordinate,
        molar_density *
            (0.7 * 0.04401 +
             0.3 * 0.018015));
}

void cpa_provider_contract() {
    const auto parameters =
        cpa_parameters();
    const auto model =
        th::CpaPtPhase::
            from_parameters(parameters);
    constexpr double pressure = 1.0e5;
    constexpr double temperature = 200.0;
    const std::vector<double> x{0.6, 0.4};
    const auto roots =
        model.roots(
            pressure,
            temperature,
            x);
    require_provider(
        roots.status ==
                th::CpaPtRootStatus::success &&
            !roots.roots.empty(),
        "CPA provider fixture has no resolved density root");
    const std::size_t selected_root =
        roots.roots.size() - 1U;

    const auto cont =
        continuation();
    auto registry =
        flow::
            make_transition_selected_phase_branch_registry(
                cont,
                std::vector<
                    flow::
                        FrozenSelectedPhaseBranchBinding<
                            th::CpaSelectedPhase>>{
                    {
                        phase_id("phase-a"),
                        {0U, {}},
                        "fixture/source-density-root"},
                    {
                        phase_id("phase-b"),
                        {selected_root, {}},
                        "fixture/target-density-root"}});

    flow::
        CpaAbsentPhasePotentialExtensionProvider
        provider{
            model,
            std::move(registry)};
    const auto coordinate =
        coordinates(
            phase_id("phase-b"),
            host_identity(
                {"A", "B"},
                pressure,
                temperature,
                x),
            x);
    const auto extension =
        provider.evaluate(
            coordinate);
    const auto molar_density =
        th::evaluate_selected_phase_molar_density(
            model,
            pressure,
            temperature,
            std::span<const double>{x},
            th::CpaSelectedPhase{
                selected_root, {}})
            .molar_density_mol_per_m3;
    require_extension(
        extension,
        coordinate,
        molar_density *
            (0.6 * 0.044 +
             0.4 * 0.030));
}

void identity_continuation_rejection() {
    flow::PhaseSetTransitionCandidate
        candidate;
    candidate.source_phase_count = 1U;
    candidate.target_phase_count = 2U;
    candidate.status =
        flow::
            PhaseSetTransitionCandidateStatus::
                target_resolved;
    candidate.evidence_profile =
        "fixture/identity-rejection";

    bool caught = false;
    try {
        (void)flow::
            make_phase_identity_continuation_snapshot(
                candidate,
                flow::
                    FrozenActivePhaseIdentityMap{
                        {phase_id("phase-a")}},
                flow::
                    FrozenActivePhaseIdentityMap{
                        {
                            phase_id("renamed-phase"),
                            phase_id("phase-b")}});
    } catch (const std::invalid_argument&) {
        caught = true;
    }
    require_provider(
        caught,
        "phase appearance silently renamed an existing physical identity");
}

} // namespace

int main() {
    try {
        identity_continuation_rejection();
        pr76_provider_contract();
        sw92_provider_contract();
        cpa_provider_contract();
        return 0;
    } catch (const std::exception& error) {
        std::cerr
            << "[FAIL] "
            << error.what()
            << '\n';
        return 1;
    }
}
