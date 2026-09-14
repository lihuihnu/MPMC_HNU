#include <mpmc/model_configuration_grpc/model_wire.hpp>
#include <google/rpc/status.pb.h>

#include <limits>
#include <utility>

namespace mpmc::model_configuration_grpc {
grpc::Status model_error_status(grpc::StatusCode code, std::string_view message,
    const model_configuration::v1::ModelServiceError& detail) {
    google::rpc::Status envelope;
    envelope.set_code(static_cast<int>(code));
    envelope.set_message(std::string(message));
    if (!envelope.add_details()->PackFrom(detail)) {
        return {code, std::string(message)};
    }
    return {code, std::string(message), envelope.SerializeAsString()};
}
namespace {
namespace mc = ::mpmc::model_configuration;
namespace wire = ::mpmc::model_configuration::v1;
namespace old = ::mpmc::runtime::v1;
namespace fl = ::mpmc::flash;

void required(bool present, const char* field) {
    if (!present) {
        throw mc::ModelConfigurationError(mc::ModelConfigurationErrorCode::missing_field,
                                           field, "explicit wire field required");
    }
}
[[noreturn]] void invalid_enum(const char* field) {
    throw mc::ModelConfigurationError(mc::ModelConfigurationErrorCode::invalid_value,
                                       field, "unknown enum value");
}
mc::ThermodynamicModelFamily family(wire::ModelFamily value) {
    switch (value) {
    case wire::MODEL_FAMILY_UNSPECIFIED: return mc::ThermodynamicModelFamily::unspecified;
    case wire::MODEL_FAMILY_PR76: return mc::ThermodynamicModelFamily::peng_robinson_1976;
    case wire::MODEL_FAMILY_SW92: return mc::ThermodynamicModelFamily::soreide_whitson_1992;
    case wire::MODEL_FAMILY_CPA: return mc::ThermodynamicModelFamily::cubic_plus_association;
    default: invalid_enum("family");
    }
}
wire::ModelFamily family(mc::ThermodynamicModelFamily value) {
    switch (value) {
    case mc::ThermodynamicModelFamily::unspecified: return wire::MODEL_FAMILY_UNSPECIFIED;
    case mc::ThermodynamicModelFamily::peng_robinson_1976: return wire::MODEL_FAMILY_PR76;
    case mc::ThermodynamicModelFamily::soreide_whitson_1992: return wire::MODEL_FAMILY_SW92;
    case mc::ThermodynamicModelFamily::cubic_plus_association: return wire::MODEL_FAMILY_CPA;
    default: invalid_enum("family");
    }
}
mc::ComponentKind component_kind(wire::ComponentKind value) {
    switch (value) {
    case wire::COMPONENT_KIND_UNSPECIFIED: return mc::ComponentKind::unspecified;
    case wire::COMPONENT_KIND_PURE: return mc::ComponentKind::pure;
    case wire::COMPONENT_KIND_PSEUDO: return mc::ComponentKind::pseudo;
    default: invalid_enum("component_kind");
    }
}
wire::ComponentKind component_kind(mc::ComponentKind value) {
    switch (value) {
    case mc::ComponentKind::unspecified: return wire::COMPONENT_KIND_UNSPECIFIED;
    case mc::ComponentKind::pure: return wire::COMPONENT_KIND_PURE;
    case mc::ComponentKind::pseudo: return wire::COMPONENT_KIND_PSEUDO;
    default: invalid_enum("component_kind");
    }
}
mc::SourceKind source_kind(wire::SourceKind value) {
    switch (value) {
    case wire::SOURCE_KIND_UNSPECIFIED: return mc::SourceKind::unspecified;
    case wire::SOURCE_KIND_LITERATURE: return mc::SourceKind::literature;
    case wire::SOURCE_KIND_DATABASE: return mc::SourceKind::database;
    case wire::SOURCE_KIND_USER_SUPPLIED: return mc::SourceKind::user_supplied;
    case wire::SOURCE_KIND_ASSUMPTION: return mc::SourceKind::assumption;
    case wire::SOURCE_KIND_SYNTHETIC_TEST: return mc::SourceKind::synthetic_test;
    default: invalid_enum("source_kind");
    }
}
wire::SourceKind source_kind(mc::SourceKind value) {
    switch (value) {
    case mc::SourceKind::unspecified: return wire::SOURCE_KIND_UNSPECIFIED;
    case mc::SourceKind::literature: return wire::SOURCE_KIND_LITERATURE;
    case mc::SourceKind::database: return wire::SOURCE_KIND_DATABASE;
    case mc::SourceKind::user_supplied: return wire::SOURCE_KIND_USER_SUPPLIED;
    case mc::SourceKind::assumption: return wire::SOURCE_KIND_ASSUMPTION;
    case mc::SourceKind::synthetic_test: return wire::SOURCE_KIND_SYNTHETIC_TEST;
    default: invalid_enum("source_kind");
    }
}
mc::PtSolverSettingsKind settings_kind(wire::SolverSettingsKind value) {
    switch (value) {
    case wire::SOLVER_SETTINGS_KIND_UNSPECIFIED: return mc::PtSolverSettingsKind::unspecified;
    case wire::SOLVER_SETTINGS_KIND_PRESET: return mc::PtSolverSettingsKind::preset;
    case wire::SOLVER_SETTINGS_KIND_CUSTOM: return mc::PtSolverSettingsKind::custom;
    default: invalid_enum("settings_kind");
    }
}
wire::SolverSettingsKind settings_kind(mc::PtSolverSettingsKind value) {
    switch (value) {
    case mc::PtSolverSettingsKind::unspecified: return wire::SOLVER_SETTINGS_KIND_UNSPECIFIED;
    case mc::PtSolverSettingsKind::preset: return wire::SOLVER_SETTINGS_KIND_PRESET;
    case mc::PtSolverSettingsKind::custom: return wire::SOLVER_SETTINGS_KIND_CUSTOM;
    default: invalid_enum("settings_kind");
    }
}
mc::ModelProvenance read_provenance(const wire::ModelProvenance& w) {
    required(w.has_kind(), "provenance.kind");
    return {source_kind(w.kind()), w.reference(), w.revision(), w.locator(), w.note(),
            w.acquisition(), w.usage_terms()};
}
void write_provenance(const mc::ModelProvenance& n, wire::ModelProvenance& w) {
    w.set_kind(source_kind(n.kind));
    w.set_reference(n.reference); w.set_revision(n.revision); w.set_locator(n.locator);
    w.set_note(n.note); w.set_acquisition(n.acquisition); w.set_usage_terms(n.usage_terms);
}
mc::ModelScalar read_scalar(const wire::ModelScalar& w) {
    required(w.has_value(), "scalar.value");
    required(w.has_provenance(), "scalar.provenance");
    return {w.value(), read_provenance(w.provenance()), w.original_unit(), w.conversion()};
}
void write_scalar(const mc::ModelScalar& n, wire::ModelScalar& w) {
    w.set_value(n.value); write_provenance(n.provenance, *w.mutable_provenance());
    w.set_original_unit(n.original_unit); w.set_conversion(n.conversion);
}
mc::PtEosRootSettings read_settings(const wire::PtEosRootSettings& w) {
    mc::PtEosRootSettings n;
    if (w.has_max_iterations()) { n.max_iterations = w.max_iterations(); }
    return n;
}
void write_settings(const mc::PtEosRootSettings& n, wire::PtEosRootSettings& w) {
    if (n.max_iterations) { w.set_max_iterations(*n.max_iterations); }
}
mc::PtStabilitySettings read_settings(const wire::PtStabilitySettings& w) {
    mc::PtStabilitySettings n;
    if (w.has_tpd_tolerance()) { n.tpd_tolerance = w.tpd_tolerance(); }
    if (w.has_stationarity_tolerance()) { n.stationarity_tolerance = w.stationarity_tolerance(); }
    if (w.has_line_search_armijo_coefficient()) { n.line_search_armijo_coefficient = w.line_search_armijo_coefficient(); }
    if (w.has_max_log_composition_step()) { n.max_log_composition_step = w.max_log_composition_step(); }
    if (w.has_max_iterations()) { n.max_iterations = w.max_iterations(); }
    if (w.has_max_backtracks()) { n.max_backtracks = w.max_backtracks(); }
    if (w.has_max_property_evaluations()) { n.max_property_evaluations = w.max_property_evaluations(); }
    if (w.has_automatic_multistart()) { n.automatic_multistart = w.automatic_multistart(); }
    if (w.has_max_starts()) { n.max_starts = w.max_starts(); }
    return n;
}
void write_settings(const mc::PtStabilitySettings& n, wire::PtStabilitySettings& w) {
    if (n.tpd_tolerance) { w.set_tpd_tolerance(*n.tpd_tolerance); }
    if (n.stationarity_tolerance) { w.set_stationarity_tolerance(*n.stationarity_tolerance); }
    if (n.line_search_armijo_coefficient) { w.set_line_search_armijo_coefficient(*n.line_search_armijo_coefficient); }
    if (n.max_log_composition_step) { w.set_max_log_composition_step(*n.max_log_composition_step); }
    if (n.max_iterations) { w.set_max_iterations(*n.max_iterations); }
    if (n.max_backtracks) { w.set_max_backtracks(*n.max_backtracks); }
    if (n.max_property_evaluations) { w.set_max_property_evaluations(*n.max_property_evaluations); }
    if (n.automatic_multistart) { w.set_automatic_multistart(*n.automatic_multistart); }
    if (n.max_starts) { w.set_max_starts(*n.max_starts); }
}
mc::PtTwoPhaseSettings read_settings(const wire::PtTwoPhaseSettings& w) {
    mc::PtTwoPhaseSettings n;
    if (w.has_fugacity_equilibrium_tolerance()) { n.fugacity_equilibrium_tolerance = w.fugacity_equilibrium_tolerance(); }
    if (w.has_absolute_mass_balance_tolerance()) { n.absolute_mass_balance_tolerance = w.absolute_mass_balance_tolerance(); }
    if (w.has_relative_mass_balance_tolerance()) { n.relative_mass_balance_tolerance = w.relative_mass_balance_tolerance(); }
    if (w.has_minimum_phase_fraction()) { n.minimum_phase_fraction = w.minimum_phase_fraction(); }
    if (w.has_minimum_log_composition_separation()) { n.minimum_log_composition_separation = w.minimum_log_composition_separation(); }
    if (w.has_minimum_relative_z_separation()) { n.minimum_relative_z_separation = w.minimum_relative_z_separation(); }
    if (w.has_max_log_equilibrium_ratio_step()) { n.max_log_equilibrium_ratio_step = w.max_log_equilibrium_ratio_step(); }
    if (w.has_residual_progress_coefficient()) { n.residual_progress_coefficient = w.residual_progress_coefficient(); }
    if (w.has_gibbs_progress_coefficient()) { n.gibbs_progress_coefficient = w.gibbs_progress_coefficient(); }
    if (w.has_max_iterations()) { n.max_iterations = w.max_iterations(); }
    if (w.has_max_backtracks()) { n.max_backtracks = w.max_backtracks(); }
    if (w.has_max_property_evaluations()) { n.max_property_evaluations = w.max_property_evaluations(); }
    if (w.has_rachford_rice_max_iterations()) { n.rachford_rice_max_iterations = w.rachford_rice_max_iterations(); }
    if (w.has_max_split_attempts()) { n.max_split_attempts = w.max_split_attempts(); }
    return n;
}
void write_settings(const mc::PtTwoPhaseSettings& n, wire::PtTwoPhaseSettings& w) {
    if (n.fugacity_equilibrium_tolerance) { w.set_fugacity_equilibrium_tolerance(*n.fugacity_equilibrium_tolerance); }
    if (n.absolute_mass_balance_tolerance) { w.set_absolute_mass_balance_tolerance(*n.absolute_mass_balance_tolerance); }
    if (n.relative_mass_balance_tolerance) { w.set_relative_mass_balance_tolerance(*n.relative_mass_balance_tolerance); }
    if (n.minimum_phase_fraction) { w.set_minimum_phase_fraction(*n.minimum_phase_fraction); }
    if (n.minimum_log_composition_separation) { w.set_minimum_log_composition_separation(*n.minimum_log_composition_separation); }
    if (n.minimum_relative_z_separation) { w.set_minimum_relative_z_separation(*n.minimum_relative_z_separation); }
    if (n.max_log_equilibrium_ratio_step) { w.set_max_log_equilibrium_ratio_step(*n.max_log_equilibrium_ratio_step); }
    if (n.residual_progress_coefficient) { w.set_residual_progress_coefficient(*n.residual_progress_coefficient); }
    if (n.gibbs_progress_coefficient) { w.set_gibbs_progress_coefficient(*n.gibbs_progress_coefficient); }
    if (n.max_iterations) { w.set_max_iterations(*n.max_iterations); }
    if (n.max_backtracks) { w.set_max_backtracks(*n.max_backtracks); }
    if (n.max_property_evaluations) { w.set_max_property_evaluations(*n.max_property_evaluations); }
    if (n.rachford_rice_max_iterations) { w.set_rachford_rice_max_iterations(*n.rachford_rice_max_iterations); }
    if (n.max_split_attempts) { w.set_max_split_attempts(*n.max_split_attempts); }
}
mc::PtThreePhaseSettings read_settings(const wire::PtThreePhaseSettings& w) {
    mc::PtThreePhaseSettings n;
    if (w.has_chemical_potential_tolerance()) { n.chemical_potential_tolerance = w.chemical_potential_tolerance(); }
    if (w.has_absolute_mass_balance_tolerance()) { n.absolute_mass_balance_tolerance = w.absolute_mass_balance_tolerance(); }
    if (w.has_relative_mass_balance_tolerance()) { n.relative_mass_balance_tolerance = w.relative_mass_balance_tolerance(); }
    if (w.has_generalized_rr_balance_tolerance()) { n.generalized_rr_balance_tolerance = w.generalized_rr_balance_tolerance(); }
    if (w.has_minimum_phase_fraction()) { n.minimum_phase_fraction = w.minimum_phase_fraction(); }
    if (w.has_minimum_log_composition_separation()) { n.minimum_log_composition_separation = w.minimum_log_composition_separation(); }
    if (w.has_max_log_step()) { n.max_log_step = w.max_log_step(); }
    if (w.has_residual_progress_coefficient()) { n.residual_progress_coefficient = w.residual_progress_coefficient(); }
    if (w.has_max_iterations()) { n.max_iterations = w.max_iterations(); }
    if (w.has_max_line_search_backtracks()) { n.max_line_search_backtracks = w.max_line_search_backtracks(); }
    if (w.has_max_balance_iterations()) { n.max_balance_iterations = w.max_balance_iterations(); }
    if (w.has_max_balance_backtracks()) { n.max_balance_backtracks = w.max_balance_backtracks(); }
    if (w.has_max_property_evaluations()) { n.max_property_evaluations = w.max_property_evaluations(); }
    if (w.has_max_three_phase_attempts()) { n.max_three_phase_attempts = w.max_three_phase_attempts(); }
    if (w.has_new_phase_seed_fraction()) { n.new_phase_seed_fraction = w.new_phase_seed_fraction(); }
    return n;
}
void write_settings(const mc::PtThreePhaseSettings& n, wire::PtThreePhaseSettings& w) {
    if (n.chemical_potential_tolerance) { w.set_chemical_potential_tolerance(*n.chemical_potential_tolerance); }
    if (n.absolute_mass_balance_tolerance) { w.set_absolute_mass_balance_tolerance(*n.absolute_mass_balance_tolerance); }
    if (n.relative_mass_balance_tolerance) { w.set_relative_mass_balance_tolerance(*n.relative_mass_balance_tolerance); }
    if (n.generalized_rr_balance_tolerance) { w.set_generalized_rr_balance_tolerance(*n.generalized_rr_balance_tolerance); }
    if (n.minimum_phase_fraction) { w.set_minimum_phase_fraction(*n.minimum_phase_fraction); }
    if (n.minimum_log_composition_separation) { w.set_minimum_log_composition_separation(*n.minimum_log_composition_separation); }
    if (n.max_log_step) { w.set_max_log_step(*n.max_log_step); }
    if (n.residual_progress_coefficient) { w.set_residual_progress_coefficient(*n.residual_progress_coefficient); }
    if (n.max_iterations) { w.set_max_iterations(*n.max_iterations); }
    if (n.max_line_search_backtracks) { w.set_max_line_search_backtracks(*n.max_line_search_backtracks); }
    if (n.max_balance_iterations) { w.set_max_balance_iterations(*n.max_balance_iterations); }
    if (n.max_balance_backtracks) { w.set_max_balance_backtracks(*n.max_balance_backtracks); }
    if (n.max_property_evaluations) { w.set_max_property_evaluations(*n.max_property_evaluations); }
    if (n.max_three_phase_attempts) { w.set_max_three_phase_attempts(*n.max_three_phase_attempts); }
    if (n.new_phase_seed_fraction) { w.set_new_phase_seed_fraction(*n.new_phase_seed_fraction); }
}
[[nodiscard]] std::uint32_t checked_uint32(std::size_t value) {
    if (value > std::numeric_limits<std::uint32_t>::max()) {
        throw std::overflow_error("PT service index does not fit uint32");
    }
    return static_cast<std::uint32_t>(value);
}

[[nodiscard]] old::PtPhaseTransitionSupport transition_support(
    flash::PtPhaseTransitionSupport value) {
    switch (value) {
    case flash::PtPhaseTransitionSupport::detection_only:
        return old::PT_PHASE_TRANSITION_SUPPORT_DETECTION_ONLY;
    case flash::PtPhaseTransitionSupport::fresh_target_resolve:
        return old::PT_PHASE_TRANSITION_SUPPORT_FRESH_TARGET_RESOLVE;
    }
    throw std::logic_error("unknown phase transition support");
}

[[nodiscard]] old::PtPhaseTransitionTrigger transition_trigger(
    flash::PtPhaseTransitionTrigger value) {
    switch (value) {
    case flash::PtPhaseTransitionTrigger::initial_stability_witness:
        return old::PT_PHASE_TRANSITION_TRIGGER_INITIAL_STABILITY_WITNESS;
    case flash::PtPhaseTransitionTrigger::final_phase_set_instability:
        return old::PT_PHASE_TRANSITION_TRIGGER_FINAL_PHASE_SET_INSTABILITY;
    case flash::PtPhaseTransitionTrigger::phase_disappearance:
        return old::PT_PHASE_TRANSITION_TRIGGER_PHASE_DISAPPEARANCE;
    case flash::PtPhaseTransitionTrigger::provider_topology_witness:
        return old::PT_PHASE_TRANSITION_TRIGGER_PROVIDER_TOPOLOGY_WITNESS;
    case flash::PtPhaseTransitionTrigger::provider_boundary_route:
        return old::PT_PHASE_TRANSITION_TRIGGER_PROVIDER_BOUNDARY_ROUTE;
    }
    throw std::logic_error("unknown phase transition trigger");
}

[[nodiscard]] old::PtPhaseTransitionResolution transition_resolution(
    flash::PtPhaseTransitionResolution value) {
    switch (value) {
    case flash::PtPhaseTransitionResolution::accepted_target:
        return old::PT_PHASE_TRANSITION_RESOLUTION_ACCEPTED_TARGET;
    case flash::PtPhaseTransitionResolution::target_resolve_required:
        return old::PT_PHASE_TRANSITION_RESOLUTION_TARGET_RESOLVE_REQUIRED;
    case flash::PtPhaseTransitionResolution::target_resolve_failed:
        return old::PT_PHASE_TRANSITION_RESOLUTION_TARGET_RESOLVE_FAILED;
    case flash::PtPhaseTransitionResolution::candidate_not_accepted:
        return old::PT_PHASE_TRANSITION_RESOLUTION_CANDIDATE_NOT_ACCEPTED;
    case flash::PtPhaseTransitionResolution::broader_topology_required:
        return old::PT_PHASE_TRANSITION_RESOLUTION_BROADER_TOPOLOGY_REQUIRED;
    case flash::PtPhaseTransitionResolution::indeterminate:
        return old::PT_PHASE_TRANSITION_RESOLUTION_INDETERMINATE;
    }
    throw std::logic_error("unknown phase transition resolution");
}

void map_transition_capability(
    const flash::PtPhaseTransitionCapability& source,
    old::PtPhaseTransitionCapability& target) {
    target.set_convention(flash::PtPhaseTransitionCapability::convention.data(),
                          flash::PtPhaseTransitionCapability::convention.size());
    for (const auto& edge : source.edges) {
        auto* output = target.add_edges();
        output->set_source_phase_count(checked_uint32(edge.source_phase_count));
        output->set_target_phase_count(checked_uint32(edge.target_phase_count));
        output->set_support(transition_support(edge.support));
        output->set_requires_fresh_target_solve(
            edge.requires_fresh_target_solve);
    }
}

void map_capability(const flash::PtFlashBackendCapability& source,
                    old::PtBackendCapability& target) {
    target.set_convention(flash::PtFlashBackendCapability::convention.data(),
                          flash::PtFlashBackendCapability::convention.size());
    target.set_backend_id(source.backend_id);
    target.set_model_profile(source.model_profile);
    target.set_algorithm_profile(source.algorithm_profile);
    target.set_publication_profile(source.publication_profile);
    target.set_configuration_profile(source.configuration_profile);
    target.set_dataset_id(source.dataset_id);
    target.set_revision(source.revision);
    for (const auto& component_id : source.component_ids) {
        target.add_component_ids(component_id);
    }
    for (const std::size_t count : source.supported_phase_counts) {
        target.add_supported_phase_counts(checked_uint32(count));
    }
    for (const auto& setting : source.scalar_settings) {
        auto* output = target.add_scalar_settings();
        output->set_id(setting.id);
        output->set_value(setting.value);
        output->set_unit(setting.unit);
    }
    map_transition_capability(source.transition_capability,
                              *target.mutable_transition_capability());
    target.set_performs_initial_stability_search(
        source.performs_initial_stability_search);
    target.set_performs_final_phase_set_review(
        source.performs_final_phase_set_review);
    target.set_performs_boundary_neighbor_resolve(
        source.performs_boundary_neighbor_resolve);
    target.set_global_stability_proven(source.global_stability_proven);
    target.set_phase_metadata_namespace(source.phase_metadata_namespace);
}

void map_transition_report(const flash::PtPhaseTransitionReport& source,
                           old::PtPhaseTransitionReport& target) {
    target.set_convention(flash::PtPhaseTransitionReport::convention.data(),
                          flash::PtPhaseTransitionReport::convention.size());
    for (const auto& evidence : source.evidence) {
        auto* output = target.add_evidence();
        output->set_source_phase_count(
            checked_uint32(evidence.source_phase_count));
        if (evidence.target_phase_count) {
            output->set_target_phase_count(
                checked_uint32(*evidence.target_phase_count));
        }
        output->set_trigger(transition_trigger(evidence.trigger));
        output->set_resolution(transition_resolution(evidence.resolution));
        output->set_fresh_target_solve_attempted(
            evidence.fresh_target_solve_attempted);
        output->set_target_topology_closed(evidence.target_topology_closed);
        output->set_provider_evidence_profile(
            evidence.provider_evidence_profile);
        output->set_diagnostic(evidence.diagnostic);
    }
}
} // namespace

mc::ThermodynamicModelDefinition decode_definition(const wire::ThermodynamicModelDefinition& w,
                                                       const mc::ModelConfigurationLimits& limits) {
    required(w.has_version(), "definition.version");
    required(w.has_family(), "definition.family");
    required(w.has_provenance(), "definition.provenance");
    required(w.has_applicability(), "definition.applicability");
    required(w.applicability().has_provenance(), "applicability.provenance");
    mc::ThermodynamicModelDefinition n;
    n.family = family(w.family());
    if (n.family != mc::ThermodynamicModelFamily::peng_robinson_1976) {
        throw mc::ModelConfigurationError(mc::ModelConfigurationErrorCode::unsupported_family,
                                           "family", "custom model family is unsupported");
    }
    required(w.has_pr76(), "definition.pr76");
    if (static_cast<std::size_t>(w.components_size()) > limits.max_components ||
        static_cast<std::size_t>(w.pr76().pure_size()) > limits.max_components ||
        static_cast<std::size_t>(w.pr76().binary_size()) > limits.max_pair_records) {
        throw mc::ModelConfigurationError(mc::ModelConfigurationErrorCode::resource_limit,
                                           "definition", "wire shape exceeds host limits");
    }
    n.version = w.version(); n.display_name = w.display_name();
    n.dataset_id = w.dataset_id(); n.revision = w.revision();
    n.provenance = read_provenance(w.provenance());
    for (const auto& c : w.components()) {
        required(c.has_kind(), "components.kind");
        required(c.has_provenance(), "components.provenance");
        mc::ComponentDefinition item{c.component_id(), c.display_name(), component_kind(c.kind()),
                                     read_provenance(c.provenance()), {}};
        if (c.has_molar_mass_kg_per_mol()) { item.molar_mass_kg_per_mol = read_scalar(c.molar_mass_kg_per_mol()); }
        n.components.push_back(std::move(item));
    }
    const auto& bounds = w.applicability();
    n.applicability.provenance = read_provenance(bounds.provenance());
    if (bounds.has_temperature_lower_k()) { n.applicability.temperature_lower_k = bounds.temperature_lower_k(); }
    if (bounds.has_temperature_upper_k()) { n.applicability.temperature_upper_k = bounds.temperature_upper_k(); }
    if (bounds.has_pressure_lower_pa()) { n.applicability.pressure_lower_pa = bounds.pressure_lower_pa(); }
    if (bounds.has_pressure_upper_pa()) { n.applicability.pressure_upper_pa = bounds.pressure_upper_pa(); }
    mc::Pr76ParameterDefinition pr;
    for (const auto& record : w.pr76().pure()) {
        mc::Pr76PureParameters item; item.component_id = record.component_id();
        if (record.has_critical_temperature_k()) { item.critical_temperature_k = read_scalar(record.critical_temperature_k()); }
        if (record.has_critical_pressure_pa()) { item.critical_pressure_pa = read_scalar(record.critical_pressure_pa()); }
        if (record.has_acentric_factor()) { item.acentric_factor = read_scalar(record.acentric_factor()); }
        pr.pure.push_back(std::move(item));
    }
    for (const auto& record : w.pr76().binary()) {
        mc::Pr76BinaryInteraction item{record.first_component_id(), record.second_component_id(), {}};
        if (record.has_kij()) { item.kij = read_scalar(record.kij()); }
        pr.binary.push_back(std::move(item));
    }
    n.parameters = std::move(pr);
    return n;
}
void encode_definition(const mc::ThermodynamicModelDefinition& n, wire::ThermodynamicModelDefinition& w) {
    w.Clear();
    w.set_version(n.version); w.set_family(family(n.family)); w.set_display_name(n.display_name);
    w.set_dataset_id(n.dataset_id); w.set_revision(n.revision);
    write_provenance(n.provenance, *w.mutable_provenance());
    for (const auto& c : n.components) {
        auto& item = *w.add_components();
        item.set_component_id(c.component_id); item.set_display_name(c.display_name);
        item.set_kind(component_kind(c.kind)); write_provenance(c.provenance, *item.mutable_provenance());
        if (c.molar_mass_kg_per_mol) { write_scalar(*c.molar_mass_kg_per_mol, *item.mutable_molar_mass_kg_per_mol()); }
    }
    auto& bounds = *w.mutable_applicability();
    write_provenance(n.applicability.provenance, *bounds.mutable_provenance());
    if (n.applicability.temperature_lower_k) { bounds.set_temperature_lower_k(*n.applicability.temperature_lower_k); }
    if (n.applicability.temperature_upper_k) { bounds.set_temperature_upper_k(*n.applicability.temperature_upper_k); }
    if (n.applicability.pressure_lower_pa) { bounds.set_pressure_lower_pa(*n.applicability.pressure_lower_pa); }
    if (n.applicability.pressure_upper_pa) { bounds.set_pressure_upper_pa(*n.applicability.pressure_upper_pa); }
    if (const auto* pr = std::get_if<mc::Pr76ParameterDefinition>(&n.parameters)) {
        auto& payload = *w.mutable_pr76();
        for (const auto& record : pr->pure) {
            auto& item = *payload.add_pure(); item.set_component_id(record.component_id);
            if (record.critical_temperature_k) { write_scalar(*record.critical_temperature_k, *item.mutable_critical_temperature_k()); }
            if (record.critical_pressure_pa) { write_scalar(*record.critical_pressure_pa, *item.mutable_critical_pressure_pa()); }
            if (record.acentric_factor) { write_scalar(*record.acentric_factor, *item.mutable_acentric_factor()); }
        }
        for (const auto& record : pr->binary) {
            auto& item = *payload.add_binary();
            item.set_first_component_id(record.first_component_id); item.set_second_component_id(record.second_component_id);
            if (record.kij) { write_scalar(*record.kij, *item.mutable_kij()); }
        }
    }
}
mc::PtSolverSettings decode_settings(const wire::PtSolverSettings& w) {
    required(w.has_version(), "settings.version"); required(w.has_kind(), "settings.kind");
    mc::PtSolverSettings n;
    n.version = w.version(); n.kind = settings_kind(w.kind()); n.preset_id = w.preset_id();
    required(w.has_eos_root(), "settings.eos_root"); n.eos_root = read_settings(w.eos_root());
    required(w.has_initial_stability(), "settings.initial_stability"); n.initial_stability = read_settings(w.initial_stability());
    required(w.has_two_phase(), "settings.two_phase"); n.two_phase = read_settings(w.two_phase());
    required(w.has_final_two_phase_stability(), "settings.final_two_phase_stability"); n.final_two_phase_stability = read_settings(w.final_two_phase_stability());
    required(w.has_three_phase(), "settings.three_phase"); n.three_phase = read_settings(w.three_phase());
    required(w.has_final_three_phase_stability(), "settings.final_three_phase_stability"); n.final_three_phase_stability = read_settings(w.final_three_phase_stability());
    return n;
}
void encode_settings(const mc::PtSolverSettings& n, wire::PtSolverSettings& w) {
    w.Clear(); w.set_version(n.version); w.set_kind(settings_kind(n.kind)); w.set_preset_id(n.preset_id);
    write_settings(n.eos_root, *w.mutable_eos_root());
    write_settings(n.initial_stability, *w.mutable_initial_stability());
    write_settings(n.two_phase, *w.mutable_two_phase());
    write_settings(n.final_two_phase_stability, *w.mutable_final_two_phase_stability());
    write_settings(n.three_phase, *w.mutable_three_phase());
    write_settings(n.final_three_phase_stability, *w.mutable_final_three_phase_stability());
}
void encode_snapshot(const mc::Pr76RegisteredModelSnapshot& n, wire::ModelSnapshot& w) {
    w.Clear();
    encode_definition(n.definition, *w.mutable_definition());
    encode_settings(n.settings, *w.mutable_settings());
    map_capability(n.capability, *w.mutable_capability());
    w.mutable_parameter_limits()->set_max_components(n.parameter_limits.max_components);
    w.mutable_parameter_limits()->set_max_pair_records(n.parameter_limits.max_pair_records);
    w.mutable_parameter_limits()->set_max_matrix_entries(n.parameter_limits.max_matrix_entries);
    w.mutable_parameter_limits()->set_max_identifier_bytes(n.parameter_limits.max_identifier_bytes);
    w.mutable_parameter_limits()->set_max_display_name_bytes(n.parameter_limits.max_display_name_bytes);
    w.mutable_parameter_limits()->set_max_provenance_field_bytes(n.parameter_limits.max_provenance_field_bytes);
    w.mutable_parameter_limits()->set_max_total_text_bytes(n.parameter_limits.max_total_text_bytes);
    w.mutable_solver_limits()->set_max_root_iterations(n.solver_limits.max_root_iterations);
    w.mutable_solver_limits()->set_max_iterations(n.solver_limits.max_iterations);
    w.mutable_solver_limits()->set_max_backtracks(n.solver_limits.max_backtracks);
    w.mutable_solver_limits()->set_max_property_evaluations(n.solver_limits.max_property_evaluations);
    w.mutable_solver_limits()->set_max_stability_starts(n.solver_limits.max_stability_starts);
    w.mutable_solver_limits()->set_max_split_attempts(n.solver_limits.max_split_attempts);
    w.mutable_solver_limits()->set_max_three_phase_attempts(n.solver_limits.max_three_phase_attempts);
    w.mutable_solver_limits()->set_max_components(n.solver_limits.max_components);
    w.mutable_solver_limits()->set_max_stability_start_entries(n.solver_limits.max_stability_start_entries);
    w.mutable_solver_limits()->set_max_three_phase_starts(n.solver_limits.max_three_phase_starts);
    w.mutable_solver_limits()->set_max_three_phase_start_entries(n.solver_limits.max_three_phase_start_entries);
}
void encode_result(const fl::PtFlashBackendResult& n, wire::FullPtResult& w) {
    if (!n.structurally_valid()) { throw std::logic_error("malformed native model result"); }
    w.Clear();
    w.set_backend_result_convention(std::string(fl::PtFlashBackendResult::convention));
    w.set_phase_set_convention(std::string(fl::PtPhaseSetResult::convention));
    map_capability(n.capability, *w.mutable_capability());
    switch (n.solution.status) {
    case fl::PtPhaseSetStatus::accepted: w.set_outcome(old::PT_COMPUTATION_OUTCOME_ACCEPTED); break;
    case fl::PtPhaseSetStatus::phase_set_unstable: w.set_outcome(old::PT_COMPUTATION_OUTCOME_PHASE_SET_UNSTABLE); break;
    case fl::PtPhaseSetStatus::indeterminate: w.set_outcome(old::PT_COMPUTATION_OUTCOME_INDETERMINATE); break;
    default: throw std::logic_error("unknown native model outcome");
    }
    w.set_maximum_phase_count(checked_uint32(n.solution.capability.maximum_phase_count));
    w.set_pressure_pa(n.solution.pressure_pa); w.set_temperature_k(n.solution.temperature_k);
    for (double value : n.solution.feed) { w.add_feed(value); }
    if (n.solution.candidate_phase_set) {
        auto& set = *w.mutable_candidate_phase_set();
        for (const auto& phase : n.solution.candidate_phase_set->phases) {
            auto& item = *set.add_phases();
            item.set_mole_phase_fraction(phase.mole_phase_fraction);
            for (double value : phase.composition) { item.add_composition(value); }
            for (double value : phase.activity.ln_phi) { item.add_ln_fugacity_coefficient(value); }
            item.set_provider_branch(phase.activity.branch); item.set_provider_branch_smooth(phase.activity.smooth);
            if (phase.compressibility_factor) { item.set_compressibility_factor(*phase.compressibility_factor); }
        }
    }
    w.set_global_stability_proven(n.solution.global_stability_proven);
    w.set_diagnostic(n.solution.diagnostic);
    map_transition_report(n.transition_report, *w.mutable_transition_report());
    w.set_provider_result_convention(n.provider_result_convention);
    for (const auto& metadata : n.phase_metadata) {
        auto& item = *w.add_phase_metadata(); item.set_role_id(metadata.role_id); item.set_family_id(metadata.family_id);
    }
    w.set_morphology_resolved(n.morphology_resolved);
}
} // namespace mpmc::model_configuration_grpc
