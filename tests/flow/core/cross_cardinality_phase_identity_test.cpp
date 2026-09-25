#include <mpmc/flow/cross_cardinality_phase_identity.hpp>

#include <cmath>
#include <cstddef>
#include <source_location>
#include <stdexcept>
#include <string>
#include <string_view>
#include <utility>
#include <vector>

namespace {

namespace flow = mpmc::flow;

void require_phase_identity(
    bool condition,
    std::string_view message,
    std::source_location where =
        std::source_location::current()) {
    if (!condition) {
        throw std::runtime_error(
            std::string{where.file_name()} +
            ":" +
            std::to_string(where.line()) +
            ": " +
            std::string{message});
    }
}

template <class Function>
void expect_invalid_phase_identity(
    Function&& function,
    std::string_view fragment) {
    try {
        function();
    } catch (const std::invalid_argument& error) {
        require_phase_identity(
            std::string_view{error.what()}.find(
                fragment) !=
                std::string_view::npos,
            "phase-identity invalid diagnostic changed");
        return;
    }
    throw std::runtime_error(
        "expected invalid_argument in phase-identity contract");
}

flow::FrozenPhysicalPhaseIdentity
id(
    std::string key) {
    return {
        "fixture/explicit-identity",
        std::move(key)};
}

} // namespace

void cross_cardinality_phase_identity_contract() {
    const flow::FrozenActivePhaseIdentityMap
        one{
            {id("aqueous")}};
    const flow::FrozenActivePhaseIdentityMap
        two{
            {
                id("aqueous"),
                id("hydrocarbon-0")}};
    const flow::FrozenActivePhaseIdentityMap
        three{
            {
                id("aqueous"),
                id("hydrocarbon-0"),
                id("hydrocarbon-1")}};

    const auto one_two =
        flow::
            make_cross_cardinality_face_phase_identity_plan(
                one,
                two);
    require_phase_identity(
        one_two.owner_phase_count() == 1U &&
            one_two.neighbour_phase_count() == 2U &&
            !one_two.same_cardinality() &&
            !one_two.same_active_identity_set() &&
            !one_two.slot_aligned_same_active_set() &&
            one_two.requires_absent_phase_extension() &&
            one_two.bindings().size() == 2U,
        "1P/2P phase-identity plan metadata mismatch");

    const auto* aqueous =
        one_two.find(
            id("aqueous"));
    const auto* hydrocarbon0 =
        one_two.find(
            id("hydrocarbon-0"));
    require_phase_identity(
        aqueous != nullptr &&
            aqueous->presence ==
                flow::
                    CrossCardinalityPhasePresence::
                        active_on_both_sides &&
            aqueous->owner_active_phase_index ==
                std::optional<std::size_t>{0U} &&
            aqueous->neighbour_active_phase_index ==
                std::optional<std::size_t>{0U} &&
            hydrocarbon0 != nullptr &&
            hydrocarbon0->presence ==
                flow::
                    CrossCardinalityPhasePresence::
                        neighbour_only_active &&
            !hydrocarbon0
                 ->owner_active_phase_index
                 .has_value() &&
            hydrocarbon0
                    ->neighbour_active_phase_index ==
                std::optional<std::size_t>{1U},
        "1P/2P shared/absent phase mapping mismatch");

    const auto reversed =
        flow::
            make_cross_cardinality_face_phase_identity_plan(
                two,
                one);
    require_phase_identity(
        reversed.bindings().size() ==
                one_two.bindings().size() &&
            reversed.bindings()[0].identity ==
                one_two.bindings()[0].identity &&
            reversed.bindings()[1].identity ==
                one_two.bindings()[1].identity &&
            reversed.bindings()[1].presence ==
                flow::
                    CrossCardinalityPhasePresence::
                        owner_only_active,
        "phase binding order/orientation reversal is not deterministic");

    const auto two_three =
        flow::
            make_cross_cardinality_face_phase_identity_plan(
                two,
                three);
    require_phase_identity(
        two_three.bindings().size() == 3U &&
            two_three.requires_absent_phase_extension() &&
            two_three.find(
                id("hydrocarbon-1")) !=
                nullptr &&
            two_three
                    .find(
                        id("hydrocarbon-1"))
                    ->presence ==
                flow::
                    CrossCardinalityPhasePresence::
                        neighbour_only_active,
        "2P/3P absent-phase mapping mismatch");

    const flow::FrozenActivePhaseIdentityMap
        swapped{
            {
                id("hydrocarbon-0"),
                id("aqueous")}};
    const auto same_set_swapped_slots =
        flow::
            make_cross_cardinality_face_phase_identity_plan(
                two,
                swapped);
    require_phase_identity(
        same_set_swapped_slots.same_cardinality() &&
            same_set_swapped_slots
                .same_active_identity_set() &&
            !same_set_swapped_slots
                 .requires_absent_phase_extension() &&
            !same_set_swapped_slots
                 .slot_aligned_same_active_set(),
        "same active set with permuted slots was treated as direct TPFA compatible");

    const auto same =
        flow::
            make_cross_cardinality_face_phase_identity_plan(
                three,
                three);
    require_phase_identity(
        same.same_cardinality() &&
            same.same_active_identity_set() &&
            same.slot_aligned_same_active_set() &&
            !same.requires_absent_phase_extension(),
        "slot-aligned identical active phase set was not recognized");

    require_phase_identity(
        flow::
                CrossCardinalityAbsentPhaseSemantics::
                    mobility_per_pa_s ==
            0.0 &&
            flow::
                CrossCardinalityAbsentPhaseSemantics::
                    mobility_gradient_is_structural_zero &&
            flow::
                CrossCardinalityAbsentPhaseSemantics::
                    phase_pressure_requires_explicit_extension &&
            flow::
                CrossCardinalityAbsentPhaseSemantics::
                    mass_density_requires_explicit_extension &&
            flow::
                CrossCardinalityAbsentPhaseSemantics::
                    absent_composition_is_undefined &&
            flow::
                CrossCardinalityAbsentPhaseSemantics::
                    absent_enthalpy_is_undefined &&
            flow::
                CrossCardinalityAbsentPhaseSemantics::
                    nonzero_advective_payload_requires_active_upwind_phase,
        "absent-phase semantics changed");

    const flow::
        AbsentPhasePotentialExtensionLinearization
        extension{
            id("hydrocarbon-1"),
            18.0e6,
            420.0,
            4U,
            {1.0, 0.0, 2.0, 0.0},
            {0.0, -0.1, 0.0, 0.2},
            "fixture/metastable-potential-extension/v1"};
    require_phase_identity(
        extension.identity() ==
                id("hydrocarbon-1") &&
            extension.phase_pressure_pa() ==
                18.0e6 &&
            extension.mass_density_kg_per_m3() ==
                420.0 &&
            extension.input_count() == 4U &&
            extension
                    .phase_pressure_gradient()
                    .size() ==
                4U &&
            extension
                    .mass_density_gradient()
                    .size() ==
                4U &&
            !extension.provenance().empty(),
        "explicit absent-phase potential extension carrier mismatch");

    expect_invalid_phase_identity(
        [] {
            (void)flow::
                FrozenActivePhaseIdentityMap{
                    {
                        id("aqueous"),
                        id("aqueous")}};
        },
        "duplicate");
    expect_invalid_phase_identity(
        [] {
            (void)flow::
                FrozenActivePhaseIdentityMap{
                    {
                        {
                            "",
                            "phase"}}};
        },
        "nonempty provenance");
    expect_invalid_phase_identity(
        [] {
            (void)flow::
                AbsentPhasePotentialExtensionLinearization{
                    id("hydrocarbon-1"),
                    0.0,
                    420.0,
                    4U,
                    {0.0, 0.0, 0.0, 0.0},
                    {0.0, 0.0, 0.0, 0.0},
                    "fixture"};
        },
        "invalid explicit");
    expect_invalid_phase_identity(
        [] {
            (void)flow::
                AbsentPhasePotentialExtensionLinearization{
                    id("hydrocarbon-1"),
                    18.0e6,
                    420.0,
                    4U,
                    {0.0, 0.0},
                    {0.0, 0.0, 0.0, 0.0},
                    "fixture"};
        },
        "invalid explicit");
}
