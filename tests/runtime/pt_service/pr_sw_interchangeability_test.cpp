#include <mpmc/flash/pr76_pt_flash_backend.hpp>
#include <mpmc/flash/sw92_profile_c_pt_flash_backend.hpp>
#include <mpmc/runtime/pt_service.hpp>

#include "../../thermodynamics/sw92/test_support.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <limits>
#include <optional>
#include <stdexcept>
#include <string>
#include <string_view>
#include <utility>
#include <vector>

namespace {
namespace fl = mpmc::flash;
namespace rt = mpmc::runtime;
namespace th = mpmc::thermodynamics;

void require(bool condition, std::string_view message) {
    if (!condition) { throw std::runtime_error(std::string(message)); }
}

void require_near(double first, double second, double relative_tolerance,
                  std::string_view message) {
    const double scale = std::max({1.0, std::abs(first), std::abs(second)});
    require(std::isfinite(first) && std::isfinite(second) &&
                std::abs(first - second) <= relative_tolerance * scale,
            message);
}

th::Provenance synthetic_pr_binary_source() {
    return {th::SourceKind::synthetic_test,
            "test://runtime/pt-service/pr-sw-interchangeability",
            "v1",
            "Classic PR binary interaction fixture",
            "Zero binary interaction is used only for interface conformance and "
            "trend comparison; no physical-accuracy claim is made",
            "tests/runtime/pt_service/pr_sw_interchangeability_test.cpp",
            "Synthetic test datum; not experimental or fitted data"};
}

th::Pr76Phase<double> classic_pr_model_with_sw_inventory() {
    const auto pure_source = sw92_test::paper("Table 3 physical properties");
    const auto binary_source = synthetic_pr_binary_source();
    const std::array<const sw92_test::Spec*, 2> specs{
        &sw92_test::co2, &sw92_test::water};

    std::vector<th::Component> catalog;
    std::vector<std::string> order;
    th::PrParameterInput input;
    input.model_id = std::string(th::pr76_profile);
    input.dataset_id = "pt-service-pr-sw-interface-fixture";
    input.revision = "v1";
    input.applicability = {std::nullopt, std::nullopt, binary_source};

    for (const auto* spec : specs) {
        catalog.push_back({spec->id, spec->display, th::ComponentKind::pure,
                           pure_source, std::nullopt});
        order.emplace_back(spec->id);
        input.pure.push_back({
            spec->id,
            sw92_test::scalar(spec->tc, th::Unit::kelvin, pure_source,
                              "K", "identity"),
            sw92_test::scalar(spec->pc_bar * 100000.0, th::Unit::pascal,
                              pure_source, "bar", "bar * 100000 -> Pa"),
            sw92_test::scalar(spec->omega, th::Unit::dimensionless,
                              pure_source)});
    }
    input.binary.push_back({
        sw92_test::co2.id, sw92_test::water.id,
        sw92_test::scalar(0.0, th::Unit::dimensionless, binary_source,
                          "dimensionless", "identity")});

    const auto parameters = th::PrParameterSet::create(
        catalog, order, input, th::DataPolicy::allow_synthetic_tests);
    return th::Pr76Phase<double>::from_parameters(parameters);
}

void require_role_neutral_response(
    const rt::PtServiceResponse& response,
    std::string_view configured_backend_id) {
    require(response.structurally_valid(),
            "service returned a structurally invalid response");
    require(response.outcome == rt::PtServiceOutcome::accepted,
            "shared PR/SW conformance state must be accepted");
    require(response.accepted_phase_count() == 1U,
            "shared pure-CO2 conformance state must remain one phase");
    require(response.result.has_value(),
            "accepted response lost its result envelope");

    const auto& result = *response.result;
    require(result.provenance.backend.configured_backend_id ==
                configured_backend_id,
            "service did not preserve the selected configured backend ID");
    require(result.provenance.backend.capability.phase_metadata_namespace.empty(),
            "PR/SW external contract must not require provider phase identities");
    require(!result.morphology_resolved,
            "PR/SW external contract must not publish phase morphology");
    require(result.pressure_pa == 3.0e6 && result.temperature_k == 340.0,
            "service changed the common PT request");
    require(result.feed.size() == 2U &&
                result.feed[0].component_id == sw92_test::co2.id &&
                result.feed[0].mole_fraction == 1.0 &&
                result.feed[1].component_id == sw92_test::water.id &&
                result.feed[1].mole_fraction == 0.0,
            "service did not map the common component-ID feed to backend order");

    for (const auto& phase : result.phases) {
        require(!phase.provider_metadata.has_value(),
                "PR/SW public phase unexpectedly carries EOS-specific identity metadata");
        require(phase.components.size() == 2U &&
                    phase.components[0].component_id == sw92_test::co2.id &&
                    phase.components[1].component_id == sw92_test::water.id,
                "PR/SW public phase changed the shared component envelope");
    }
}

void require_same_inventory(const rt::PtRuntimeComponentInventory& first,
                            const rt::PtRuntimeComponentInventory& second) {
    require(first.components.size() == second.components.size(),
            "switching EOS backend changed component inventory size");
    for (std::size_t i = 0; i < first.components.size(); ++i) {
        require(first.components[i].component_id == second.components[i].component_id &&
                    first.components[i].feed_index == second.components[i].feed_index,
                "switching EOS backend changed the public component inventory");
    }
}

void require_pure_co2_common_limit(
    const th::Pr76Phase<double>& pr_model,
    const th::Sw92Phase<double>& sw_model) {
    const std::array<double, 2> composition{1.0, 0.0};
    const std::array<std::pair<double, double>, 4> states{{
        {1.0e6, 300.0}, {5.0e6, 340.0}, {10.0e6, 400.0}, {20.0e6, 450.0}}};

    for (const auto& [pressure_pa, temperature_k] : states) {
        th::Pr76PhaseWorkspace<double> pr_workspace;
        th::Sw92PhaseWorkspace<double> sw_workspace;
        const auto pr_roots = pr_model.roots_full(
            pressure_pa, temperature_k, composition, pr_workspace);
        const auto sw_roots = sw_model.roots(
            pressure_pa, temperature_k, composition, 0.0,
            th::SwPhaseFamily::nonaqueous, sw_workspace);
        require(pr_roots.status == sw_roots.status &&
                    pr_roots.count == sw_roots.count,
                "PR/SW pure-CO2 root topology diverged");

        for (std::size_t root = 0; root < pr_roots.count; ++root) {
            const auto pr = pr_model.evaluate_full(
                pressure_pa, temperature_k, composition, root, pr_workspace);
            const auto sw = sw_model.evaluate(
                pressure_pa, temperature_k, composition, 0.0,
                th::SwPhaseFamily::nonaqueous, root, sw_workspace);
            require_near(pr.z, sw.z, 2.0e-13,
                         "PR/SW pure-CO2 Z common limit diverged");
            require_near(pr.ln_phi.at(0), sw.ln_phi.at(0), 2.0e-13,
                         "PR/SW pure-CO2 fugacity common limit diverged");
        }
    }
}

std::string_view outcome_name(rt::PtServiceOutcome outcome) {
    switch (outcome) {
    case rt::PtServiceOutcome::accepted: return "accepted";
    case rt::PtServiceOutcome::phase_set_unstable: return "phase_set_unstable";
    case rt::PtServiceOutcome::indeterminate: return "indeterminate";
    case rt::PtServiceOutcome::error: return "error";
    }
    return "unknown";
}

struct TrendSummary {
    std::size_t phase_count{};
    double min_water_x{std::numeric_limits<double>::quiet_NaN()};
    double max_water_x{std::numeric_limits<double>::quiet_NaN()};
    double beta_at_min_water{std::numeric_limits<double>::quiet_NaN()};
    double beta_at_max_water{std::numeric_limits<double>::quiet_NaN()};
    double min_z{std::numeric_limits<double>::quiet_NaN()};
    double max_z{std::numeric_limits<double>::quiet_NaN()};
    double material_balance_residual{std::numeric_limits<double>::quiet_NaN()};
    double fugacity_log_residual{std::numeric_limits<double>::quiet_NaN()};
};

TrendSummary summarize_accepted(const rt::PtServiceComputationResult& result) {
    TrendSummary summary;
    summary.phase_count = result.phases.size();
    if (result.phases.empty()) { return summary; }

    std::size_t water_index = result.feed.size();
    for (std::size_t i = 0; i < result.feed.size(); ++i) {
        if (result.feed[i].component_id == sw92_test::water.id) {
            water_index = i;
            break;
        }
    }
    require(water_index < result.feed.size(), "trend audit lost water component");

    summary.min_water_x = std::numeric_limits<double>::infinity();
    summary.max_water_x = -std::numeric_limits<double>::infinity();
    summary.min_z = std::numeric_limits<double>::infinity();
    summary.max_z = -std::numeric_limits<double>::infinity();
    bool all_z = true;
    for (const auto& phase : result.phases) {
        const double water_fraction =
            phase.components.at(water_index).mole_fraction;
        if (water_fraction < summary.min_water_x) {
            summary.min_water_x = water_fraction;
            summary.beta_at_min_water = phase.mole_phase_fraction;
        }
        if (water_fraction > summary.max_water_x) {
            summary.max_water_x = water_fraction;
            summary.beta_at_max_water = phase.mole_phase_fraction;
        }
        if (phase.compressibility_factor.has_value()) {
            summary.min_z = std::min(summary.min_z,
                                     *phase.compressibility_factor);
            summary.max_z = std::max(summary.max_z,
                                     *phase.compressibility_factor);
        } else {
            all_z = false;
        }
    }
    if (!all_z) {
        summary.min_z = std::numeric_limits<double>::quiet_NaN();
        summary.max_z = std::numeric_limits<double>::quiet_NaN();
    }

    summary.material_balance_residual = 0.0;
    for (std::size_t component = 0; component < result.feed.size(); ++component) {
        double reconstructed = 0.0;
        for (const auto& phase : result.phases) {
            reconstructed += phase.mole_phase_fraction *
                             phase.components.at(component).mole_fraction;
        }
        summary.material_balance_residual = std::max(
            summary.material_balance_residual,
            std::abs(reconstructed - result.feed[component].mole_fraction));
    }

    summary.fugacity_log_residual = 0.0;
    if (result.phases.size() > 1U) {
        for (std::size_t component = 0; component < result.feed.size(); ++component) {
            double minimum = std::numeric_limits<double>::infinity();
            double maximum = -std::numeric_limits<double>::infinity();
            bool comparable = true;
            for (const auto& phase : result.phases) {
                const double mole_fraction =
                    phase.components.at(component).mole_fraction;
                if (!(mole_fraction > 1.0e-14)) {
                    comparable = false;
                    break;
                }
                const double log_f_over_p =
                    std::log(mole_fraction) +
                    phase.components.at(component).ln_fugacity_coefficient;
                minimum = std::min(minimum, log_f_over_p);
                maximum = std::max(maximum, log_f_over_p);
            }
            if (comparable) {
                summary.fugacity_log_residual = std::max(
                    summary.fugacity_log_residual, maximum - minimum);
            }
        }
    }
    return summary;
}

void emit_trend_row(rt::PtService& service, std::string_view axis,
                    std::string_view backend_id, double pressure_pa,
                    double temperature_k, double z_co2, double z_water) {
    rt::PtServiceRequest request;
    request.configured_backend_id = std::string(backend_id);
    request.pressure_pa = pressure_pa;
    request.temperature_k = temperature_k;
    request.feed = {{sw92_test::water.id, z_water},
                    {sw92_test::co2.id, z_co2}};
    const auto response = service.solve(request);
    require(response.structurally_valid(),
            "trend audit received structurally invalid service response");

    std::cout << std::setprecision(17)
              << "PR_SW_TREND axis=" << axis
              << " backend=" << backend_id
              << " p_pa=" << pressure_pa
              << " t_k=" << temperature_k
              << " z_co2=" << z_co2
              << " z_water=" << z_water
              << " outcome=" << outcome_name(response.outcome);
    if (response.result.has_value() &&
        response.outcome == rt::PtServiceOutcome::accepted) {
        const auto summary = summarize_accepted(*response.result);
        require(summary.material_balance_residual <= 1.0e-12,
                "trend audit material balance regression");
        require(summary.fugacity_log_residual <= 1.0e-9,
                "trend audit fugacity equilibrium regression");
        std::cout << " phases=" << summary.phase_count
                  << " min_x_water=" << summary.min_water_x
                  << " max_x_water=" << summary.max_water_x
                  << " beta_min_x_water=" << summary.beta_at_min_water
                  << " beta_max_x_water=" << summary.beta_at_max_water
                  << " min_z=" << summary.min_z
                  << " max_z=" << summary.max_z
                  << " mb_residual=" << summary.material_balance_residual
                  << " fug_log_residual=" << summary.fugacity_log_residual;
    }
    std::cout << '\n';
}

void emit_trend_audit(rt::PtService& service) {
    const std::array<std::pair<double, double>, 4> pure_co2_states{{
        {1.0e6, 300.0}, {5.0e6, 340.0}, {10.0e6, 400.0}, {20.0e6, 450.0}}};
    for (const auto& [pressure_pa, temperature_k] : pure_co2_states) {
        emit_trend_row(service, "pure-co2", "classic-pr", pressure_pa,
                       temperature_k, 1.0, 0.0);
        emit_trend_row(service, "pure-co2", "sw92", pressure_pa,
                       temperature_k, 1.0, 0.0);
    }

    // Classic PR deliberately uses kij=0 here. It is a structural/trend
    // comparator, not an aqueous-solubility reference model.
    for (const double pressure_mpa :
         std::array{1.0, 2.0, 5.0, 10.0, 15.0, 20.0}) {
        emit_trend_row(service, "pressure", "classic-pr",
                       pressure_mpa * 1.0e6, 350.0, 0.5, 0.5);
        emit_trend_row(service, "pressure", "sw92",
                       pressure_mpa * 1.0e6, 350.0, 0.5, 0.5);
    }
    for (const double temperature_k :
         std::array{300.0, 325.0, 350.0, 375.0, 400.0, 425.0, 450.0}) {
        emit_trend_row(service, "temperature", "classic-pr", 10.0e6,
                       temperature_k, 0.5, 0.5);
        emit_trend_row(service, "temperature", "sw92", 10.0e6,
                       temperature_k, 0.5, 0.5);
    }
}

} // namespace

int main() {
    try {
        const auto pr_model = classic_pr_model_with_sw_inventory();
        const auto sw_model = th::Sw92Phase<double>::from_parameters(
            sw92_test::binary_parameters(sw92_test::co2));

        // Non-water common limit: SW92 deliberately retains the PR76 pure kernel,
        // cubic and fugacity expression. Zero water must therefore preserve the
        // CO2 root topology, Z and CO2 fugacity coefficient at the same state.
        require_pure_co2_common_limit(pr_model, sw_model);

        fl::Pr76VleEvaluator pr_evaluator(pr_model);
        fl::Pr76PtFlashBackend pr_backend(pr_evaluator);
        fl::Sw92ProfileCPtFlashBackend sw_backend(sw_model);

        require(pr_backend.capability().component_ids ==
                    sw_backend.capability().component_ids,
                "fixture no longer gives PR and SW the same component inventory");
        require(pr_backend.capability().phase_metadata_namespace.empty() &&
                    sw_backend.capability().phase_metadata_namespace.empty(),
                "PR and SW no longer share the role-neutral external phase contract");

        const std::array<rt::PtServiceBackendRegistration, 2> backends{{
            {"classic-pr", &pr_backend},
            {"sw92", &sw_backend}}};
        rt::PtService service(backends);

        rt::PtServiceRequest request;
        request.pressure_pa = 3.0e6;
        request.temperature_k = 340.0;
        request.feed = {{sw92_test::water.id, 0.0},
                        {sw92_test::co2.id, 1.0}};

        request.configured_backend_id = "classic-pr";
        const auto pr_response = service.solve(request);
        require_role_neutral_response(pr_response, "classic-pr");

        request.configured_backend_id = "sw92";
        const auto sw_response = service.solve(request);
        require_role_neutral_response(sw_response, "sw92");

        require_same_inventory(
            pr_response.result->provenance.backend.component_inventory,
            sw_response.result->provenance.backend.component_inventory);

        emit_trend_audit(service);
        return 0;
    } catch (const std::exception& error) {
        std::cerr << "PR_SW_TREND_FAIL " << error.what() << '\n';
        return 1;
    }
}
