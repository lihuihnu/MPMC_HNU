#include <mpmc/physics/pr76_thermodynamic_closure.hpp>

#include "closure_references.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <exception>
#include <iostream>
#include <limits>
#include <optional>
#include <source_location>
#include <span>
#include <stdexcept>
#include <string>
#include <string_view>
#include <utility>
#include <vector>

namespace {
namespace fl = mpmc::flash;
namespace ph = mpmc::physics;
namespace th = mpmc::thermodynamics;
using Vec = std::vector<double>;

void require(bool value, std::string_view message,
             std::source_location where = std::source_location::current()) {
    if (!value) {
        throw std::runtime_error(
            std::string(where.file_name()) + ":" + std::to_string(where.line()) +
            ": " + std::string(message));
    }
}

void near_primal(double actual, long double expected) {
    const long double tolerance = 5.0e-11L * (1.0L + std::abs(expected));
    require(std::isfinite(actual) &&
                std::abs(static_cast<long double>(actual) - expected) <= tolerance,
            "closure primal differs from independent Decimal reference");
}

void near_derivative(double actual, long double expected) {
    const long double scale = std::max(1.0e-7L, std::abs(expected));
    const long double tolerance = 3.0e-8L * scale;
    require(std::isfinite(actual) &&
                std::abs(static_cast<long double>(actual) - expected) <= tolerance,
            "closure derivative differs from independent Decimal reference");
}

th::SourcedScalar sourced_scalar(double value, th::Unit unit,
                                  const th::Provenance& source,
                                  std::string original_unit,
                                  std::string conversion) {
    return {value, unit, source, std::move(original_unit), std::move(conversion)};
}

th::Pr76Phase<double> binary_model(bool carbon) {
    const th::Provenance source{
        th::SourceKind::literature,
        "https://academicweb.nd.edu/~markst/zm97a.pdf",
        "revised October 1997",
        carbon ? "Sections 4.1/4.4 and Table 5" : "Section 4.3 and Table 4",
        "PR numerical regression, not experimental validation",
        "Read author-hosted PDF pages",
        "Limited attributed factual parameters; no article text or code reproduced"
    };
    const std::array<std::string, 2> ids = carbon
        ? std::array<std::string, 2>{"carbon-dioxide", "methane"}
        : std::array<std::string, 2>{"nitrogen", "ethane"};
    const std::array<std::array<double, 3>, 2> specs = carbon
        ? std::array<std::array<double, 3>, 2>{{
              {304.2, 7.38e6, 0.225}, {190.6, 4.60e6, 0.008}}}
        : std::array<std::array<double, 3>, 2>{{
              {126.2, 3.39e6, 0.04}, {305.4, 4.88e6, 0.098}}};

    std::vector<th::Component> catalog;
    th::PrParameterInput input;
    input.model_id = std::string(th::pr76_profile);
    input.dataset_id = carbon ? "Hua1997-CO2-methane" : "Hua1997-nitrogen-ethane";
    input.revision = "split-regression-v1";
    input.applicability = {std::nullopt, std::nullopt, source};
    for (std::size_t i = 0; i < ids.size(); ++i) {
        catalog.push_back({ids[i], ids[i], th::ComponentKind::pure, source, {}});
        input.pure.push_back({
            ids[i],
            sourced_scalar(specs[i][0], th::Unit::kelvin, source, "K", "identity"),
            sourced_scalar(specs[i][1], th::Unit::pascal, source, "bar",
                           "bar * 100000 -> Pa"),
            sourced_scalar(specs[i][2], th::Unit::dimensionless, source,
                           "dimensionless", "identity")
        });
    }
    input.binary.push_back({
        ids[0], ids[1],
        sourced_scalar(carbon ? 0.095 : 0.08, th::Unit::dimensionless,
                       source, "dimensionless", "identity")
    });
    const std::vector<std::string> order{ids.begin(), ids.end()};
    return th::Pr76Phase<double>::from_parameters(
        th::PrParameterSet::create(catalog, order, input));
}

th::Pr76Phase<double> ternary_model() {
    const th::Provenance pure{
        th::SourceKind::literature,
        "https://doi.org/10.1002/aic.16730",
        "AIChE Journal 65(11), e16730; first published 26 July 2019",
        "Table 1, methane/ethane/propane Peng-Robinson parameters",
        "Tc/Pc/omega only; fitted methane-propane binary kij is not used here",
        "Read open Wiley full text and Table 1",
        "Limited attributed factual parameters; no article text reproduced"
    };
    const th::Provenance zero_kij{
        th::SourceKind::literature,
        "https://doi.org/10.1021/i160057a011",
        "Ind. Eng. Chem. Fundam. 15(1), 1976",
        "journal p.62, methane-ethane-propane ternary discussion and Figure 4",
        "Original ternary example used no interaction coefficients; explicit zero kij",
        "Read project-supplied Peng-Robinson 1976 PDF",
        "Limited attributed factual model statement; no article text reproduced"
    };
    constexpr std::array<std::string_view, 3> ids{
        "methane", "ethane", "propane"};
    constexpr std::array<std::array<double, 3>, 3> specs{{
        {190.555, 4.595e6, 0.0},
        {305.4, 4.88e6, 0.099},
        {369.825, 4.248e6, 0.15308}
    }};

    std::vector<th::Component> catalog;
    th::PrParameterInput input;
    input.model_id = std::string(th::pr76_profile);
    input.dataset_id = "DeitersBell-aic16730-PengRobinson1976-ternary";
    input.revision = "ternary-regression-v1";
    input.applicability = {std::nullopt, std::nullopt, zero_kij};
    for (std::size_t i = 0; i < ids.size(); ++i) {
        catalog.push_back({std::string(ids[i]), std::string(ids[i]),
                           th::ComponentKind::pure, pure, {}});
        input.pure.push_back({
            std::string(ids[i]),
            sourced_scalar(specs[i][0], th::Unit::kelvin, pure, "K", "identity"),
            sourced_scalar(specs[i][1], th::Unit::pascal, pure, "MPa",
                           "MPa * 1e6 -> Pa"),
            sourced_scalar(specs[i][2], th::Unit::dimensionless, pure,
                           "dimensionless", "identity")
        });
    }
    const auto zero = sourced_scalar(
        0.0, th::Unit::dimensionless, zero_kij, "dimensionless", "identity");
    input.binary.push_back({"methane", "ethane", zero});
    input.binary.push_back({"methane", "propane", zero});
    input.binary.push_back({"ethane", "propane", zero});
    const std::vector<std::string> order{"methane", "ethane", "propane"};
    return th::Pr76Phase<double>::from_parameters(
        th::PrParameterSet::create(catalog, order, input));
}

template <std::size_t Q>
void check_reference(
    const ph::ThermodynamicClosureSnapshot& closure,
    const std::array<long double, 4>& primal_reference,
    const std::array<long double, Q>& liquid_z_reference,
    const std::array<long double, Q>& vapor_z_reference,
    const std::array<long double, Q>& liquid_density_reference,
    const std::array<long double, Q>& vapor_density_reference) {
    require(closure.primal_status == ph::ThermodynamicClosurePrimalStatus::valid,
            closure.diagnostic);
    require(closure.linearization_status ==
                ph::ThermodynamicClosureLinearizationStatus::available,
            closure.diagnostic);
    require(closure.linearization_reason ==
                ph::ThermodynamicClosureLinearizationReason::none,
            "successful closure has a failure reason");
    require(closure.residual_available() && closure.can_seed_newton(),
            "successful closure availability helpers");
    require(closure.primal.has_value() && closure.linearization.has_value(),
            "successful closure missing atomic payload");

    const auto& primal = *closure.primal;
    const auto& linearization = *closure.linearization;
    require(linearization.input_count == Q &&
                linearization.component_count + 1 == Q,
            "closure linearization shape");
    near_primal(primal.liquid.compressibility_factor, primal_reference[0]);
    near_primal(primal.vapor.compressibility_factor, primal_reference[1]);
    near_primal(primal.liquid.molar_density_mol_per_m3, primal_reference[2]);
    near_primal(primal.vapor.molar_density_mol_per_m3, primal_reference[3]);

    const double r = th::Pr76Pure<double>::gas_constant();
    near_primal(
        primal.liquid.molar_density_mol_per_m3,
        static_cast<long double>(closure.pressure_pa /
            (primal.liquid.compressibility_factor * r * closure.temperature_k)));
    near_primal(
        primal.vapor.molar_density_mol_per_m3,
        static_cast<long double>(closure.pressure_pa /
            (primal.vapor.compressibility_factor * r * closure.temperature_k)));

    for (std::size_t column = 0; column < Q; ++column) {
        near_derivative(linearization.d_liquid_compressibility(column),
                        liquid_z_reference[column]);
        near_derivative(linearization.d_vapor_compressibility(column),
                        vapor_z_reference[column]);
        near_derivative(linearization.d_liquid_molar_density(column),
                        liquid_density_reference[column]);
        near_derivative(linearization.d_vapor_molar_density(column),
                        vapor_density_reference[column]);
    }
}

void binary_reference() {
    const auto model = binary_model(false);
    fl::Pr76VleEvaluator evaluator(model);
    const Vec feed{0.30, 0.70};
    const auto split = fl::solve_pr76_pt_vle(7.6e6, 270.0, feed, evaluator);
    const auto closure = ph::build_pr76_pt_vle_closure(split, evaluator);
    check_reference(
        closure, closure_reference::binary_primal,
        closure_reference::binary_liquid_z_gradient,
        closure_reference::binary_vapor_z_gradient,
        closure_reference::binary_liquid_density_gradient,
        closure_reference::binary_vapor_density_gradient);
    require(std::abs(closure.primal->vapor.mole_phase_fraction -
                     split.solution.candidate()->fractions.vapor_fraction) < 1e-14,
            "vapor mole phase fraction changed");
}

void ternary_reference() {
    const auto model = ternary_model();
    fl::Pr76VleEvaluator evaluator(model);
    const Vec feed{
        0.69717974137952967145,
        0.21360224237319418100,
        0.08921801624727614754};
    const auto split = fl::solve_pr76_pt_vle(4.0e6, 220.0, feed, evaluator);
    const auto closure = ph::build_pr76_pt_vle_closure(split, evaluator);
    check_reference(
        closure, closure_reference::ternary_primal,
        closure_reference::ternary_liquid_z_gradient,
        closure_reference::ternary_vapor_z_gradient,
        closure_reference::ternary_liquid_density_gradient,
        closure_reference::ternary_vapor_density_gradient);
}

void unavailable_policy() {
    const auto model = binary_model(false);
    fl::Pr76VleEvaluator evaluator(model);
    const Vec feed{0.30, 0.70};
    const auto split = fl::solve_pr76_pt_vle(7.6e6, 270.0, feed, evaluator);

    auto snapshot = ph::build_pr76_pt_vle_closure(split, evaluator);
    require(snapshot.can_seed_newton(), "fixture is not linearizable");

    ph::Pr76PtVleClosureOptions boundary_options;
    boundary_options.sensitivity.minimum_derivative_phase_fraction = 0.45;
    snapshot = ph::build_pr76_pt_vle_closure(
        split, evaluator, boundary_options);
    require(snapshot.primal_status == ph::ThermodynamicClosurePrimalStatus::valid &&
                snapshot.residual_available(),
            "phase-boundary derivative failure discarded a valid primal");
    require(snapshot.linearization_status ==
                ph::ThermodynamicClosureLinearizationStatus::unavailable &&
                snapshot.linearization_reason ==
                    ph::ThermodynamicClosureLinearizationReason::phase_boundary &&
                !snapshot.linearization.has_value() && !snapshot.can_seed_newton(),
            "phase-boundary derivative was exposed as Newton-usable");

    ph::Pr76PtVleClosureOptions conditioning_options;
    conditioning_options.sensitivity.minimum_reciprocal_condition = 0.90;
    snapshot = ph::build_pr76_pt_vle_closure(
        split, evaluator, conditioning_options);
    require(snapshot.primal_status == ph::ThermodynamicClosurePrimalStatus::valid &&
                snapshot.residual_available(),
            "ill-conditioned derivative failure discarded a valid primal");
    require(snapshot.linearization_status ==
                ph::ThermodynamicClosureLinearizationStatus::unavailable &&
                snapshot.linearization_reason ==
                    ph::ThermodynamicClosureLinearizationReason::ill_conditioned_equilibrium &&
                !snapshot.linearization.has_value() && !snapshot.can_seed_newton(),
            "ill-conditioned derivative was exposed as Newton-usable");
}

void indeterminate_policy() {
    const auto model = binary_model(true);
    fl::Pr76VleEvaluator evaluator(model);
    const Vec feed{0.30, 0.70};
    const auto split = fl::solve_pr76_pt_vle(
        2339623.8796839858, 220.0, feed, evaluator);
    require(split.solution.status == fl::PtSplitStatus::indeterminate,
            "near-dew fixture no longer exercises indeterminate flash semantics");
    const auto closure = ph::build_pr76_pt_vle_closure(split, evaluator);
    require(closure.primal_status == ph::ThermodynamicClosurePrimalStatus::indeterminate,
            "indeterminate flash became a closure primal");
    require(!closure.primal.has_value() && !closure.linearization.has_value() &&
                !closure.residual_available() && !closure.can_seed_newton(),
            "indeterminate closure exposed stale/unsafe payload");
}

using Test = std::pair<std::string_view, void (*)()>;
constexpr Test tests[]{
    {"binary_reference", binary_reference},
    {"ternary_reference", ternary_reference},
    {"unavailable_policy", unavailable_policy},
    {"indeterminate_policy", indeterminate_policy}
};

} // namespace

int main(int argc, char** argv) {
    try {
        if (argc != 2) {
            throw std::invalid_argument("one test name required");
        }
        for (const auto& [name, run] : tests) {
            if (name == argv[1]) {
                run();
                std::cout << "[PASS] " << name << '\n';
                return 0;
            }
        }
        throw std::invalid_argument("unknown test");
    } catch (const std::exception& error) {
        std::cerr << "[FAIL] " << error.what() << '\n';
        return 1;
    }
}
