#include <mpmc/flow/detail/validation.hpp>

#include <mpmc/flow/phase_transport.hpp>

#include <array>
#include <cmath>
#include <limits>
#include <stdexcept>
#include <string>
#include <vector>

namespace {

namespace validation = mpmc::flow::validation_detail;

void require(bool condition, const char* message) {
    if (!condition) {
        throw std::runtime_error(message);
    }
}

template <std::size_t Multiplier>
void roundoff_boundary() {
    const auto near = validation::near_roundoff<Multiplier>;
    const double threshold = static_cast<double>(Multiplier) *
        std::numeric_limits<double>::epsilon();
    const double infinity = std::numeric_limits<double>::infinity();
    const double nan = std::numeric_limits<double>::quiet_NaN();
    require(near(0.0, threshold, 0.0), "exact roundoff boundary must be accepted");
    require(!near(0.0, std::nextafter(threshold, infinity), 0.0),
            "one ULP beyond the roundoff boundary must be rejected");
    require(near(0.0, -threshold, 0.0), "roundoff allowance must be sign symmetric");
    require(near(0.0, 8.0 * threshold, 8.0), "explicit closure scale must be retained");
    require(!near(0.0, 8.0 * threshold, 0.0), "closure scale must not be implicit");
    require(!near(infinity, infinity, 0.0) && !near(nan, 0.0, 0.0) &&
                !near(0.0, 0.0, infinity) && !near(0.0, 0.0, nan) &&
                !near(0.0, 0.0, -1.0),
            "non-finite inputs and invalid closure scales must be rejected");
    require(near(-0.0, 0.0, -0.0), "signed zero retains its existing meaning");
    require(!near(std::numeric_limits<double>::max(),
                  -std::numeric_limits<double>::max(), 0.0),
            "overflowing finite difference must not compare equal");
}

mpmc::flow::NaturalVariableStateIdentity make_identity(std::size_t phases) {
    mpmc::flow::NaturalVariableStateIdentity identity;
    identity.layout = mpmc::flow::NaturalVariableLayoutDescriptor{
        3U, phases, std::vector<std::size_t>(phases, 2U)};
    identity.component_ids = {"water", "carbon-dioxide", "hydrocarbon"};
    identity.reference_pressure_pa = 1.0e7;
    identity.temperature_k = 350.0;
    for (std::size_t phase = 0U; phase < phases; ++phase) {
        identity.saturation[phase] = 1.0 / static_cast<double>(phases);
        identity.phase_composition[phase] = {0.2, 0.3, 0.5};
    }
    return identity;
}

} // namespace

void shared_validation_contract() {
    roundoff_boundary<4096U>();
    roundoff_boundary<8192U>();
    roundoff_boundary<16384U>();

    const std::vector<std::string> names{"water", "carbon-dioxide"};
    require(validation::same_component_ids(names, names), "identity equality");
    require(!validation::same_component_ids(names, std::vector<std::string>{"water"}),
            "identity length is significant");
    require(!validation::same_component_ids(
                names, std::vector<std::string>{"carbon-dioxide", "water"}),
            "identity order is significant");

    const auto near = [](double first, double second) {
        return validation::near_roundoff<4096U>(first, second);
    };
    for (std::size_t phases = 1U; phases <= 3U; ++phases) {
        const auto original = make_identity(phases);
        const auto same = [&](const auto& candidate) {
            return validation::same_state_identity(original, candidate, near);
        };
        require(same(original), "valid 1P/2P/3P identities must match");
        auto changed = original;
        changed.layout = mpmc::flow::NaturalVariableLayoutDescriptor{
            3U, phases, std::vector<std::size_t>(phases, 0U)};
        require(!same(changed), "equal dimensions do not imply the same pivot");
        changed = original;
        changed.component_ids[0] = "different";
        require(!same(changed), "component identity mismatch");
        changed = original;
        changed.temperature_k += 1.0;
        require(!same(changed), "temperature identity mismatch");
        changed = original;
        changed.reference_pressure_pa = std::numeric_limits<double>::quiet_NaN();
        require(!same(changed), "non-finite pressure mismatch");
        changed = original;
        changed.phase_composition[0].pop_back();
        require(!same(changed), "composition shape mismatch");
        changed = original;
        changed.phase_composition[0][0] += 0.01;
        require(!same(changed), "composition value mismatch");
        changed = original;
        changed.saturation[0] += 0.01;
        require(!same(changed), "active saturation mismatch");
        if (phases < 3U) {
            changed = original;
            changed.saturation[phases] = std::numeric_limits<double>::denorm_min();
            require(!same(changed), "inactive saturation must be exactly zero");
            changed = original;
            changed.phase_composition[phases] = {0.0};
            require(!same(changed), "inactive composition must be empty");
        }
    }
    require(!validation::same_layout(make_identity(1U).layout, make_identity(2U).layout),
            "phase count is part of the chart identity");
}
