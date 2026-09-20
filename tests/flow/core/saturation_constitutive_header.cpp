#include <mpmc/flow/saturation_constitutive.hpp>

#include <type_traits>

struct HeaderRelativePermeability3P {
    template <typename Number>
    mpmc::flow::RelativePermeabilityEvaluation3P<Number>
    operator()(
        const mpmc::flow::
            ThreePhaseSaturationState3P<Number>& state) const {
        return {
            std::array<Number, 3>{
                state.saturation[0],
                state.saturation[1],
                state.saturation[2]}};
    }
};

static_assert(
    mpmc::flow::RelativePermeabilityEvaluator3P<
        HeaderRelativePermeability3P,
        double>);

static_assert(
    mpmc::flow::CapillaryPressureEvaluator3P<
        mpmc::flow::NoCapillaryPressure3P,
        double>);

static_assert(
    mpmc::flow::
        ThreePhaseSaturationConstitutiveEvaluation3P<double>::
            convention ==
    mpmc::flow::
        three_phase_saturation_constitutive_convention);

bool saturation_constitutive_header() {
    return true;
}
