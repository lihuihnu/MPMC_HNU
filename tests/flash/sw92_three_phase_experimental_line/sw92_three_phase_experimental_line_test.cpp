#include <mpmc/flash/sw92_stability.hpp>

#include "test_support.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <iostream>
#include <source_location>
#include <stdexcept>
#include <string>
#include <string_view>
#include <vector>

namespace {
namespace fl = mpmc::flash;
namespace th = mpmc::thermodynamics;
using Vec = std::vector<double>;

void require(bool condition, std::string_view message,
             std::source_location where = std::source_location::current()) {
    if (!condition) {
        throw std::runtime_error(std::string(where.file_name()) + ":" +
                                 std::to_string(where.line()) + ": " +
                                 std::string(message));
    }
}

void near(double actual, long double expected, long double relative = 5e-9L,
          long double absolute = 5e-12L,
          std::source_location where = std::source_location::current()) {
    if (!std::isfinite(actual) ||
        std::abs(static_cast<long double>(actual) - expected) >
            absolute + relative * std::abs(expected)) {
        std::cerr << "actual=" << actual << " expected=" << expected << '\n';
        require(false, "Decimal reference mismatch", where);
    }
}

struct TripleLineGolden {
    long double temperature_k;
    long double pressure_pa;
    long double x_butane_w;
    long double x_butane_h0;
    long double x_butane_h1;
    long double z_w;
    long double z_h0;
    long double z_h1;
    long double experimental_pressure_bar;
    long double experimental_x_butane_h0;
    long double experimental_x_butane_h1;
};

// Model values are independently solved by reference_decimal.py. Experimental
// pressure and hydrocarbon-phase compositions are transcribed from Reamer et al.
// 1944, DOI 10.1021/ie50412a024. H0/H1 are numerical names only; this test does
// not publish physical liquid/vapor morphology.
constexpr TripleLineGolden golden[] = {
    {310.93L, 366115.2242123761867209414980249077L, 0.00006145142935965007151728057515805494L, 0.9992890086701407690271617340432109L, 0.9834709095671742044581485208424232L, 0.003029822413207653255262496622640287L, 0.01410358072573536399266895145222988L, 0.9067845134705541850350445113973057L, 3.62L, 0.9995L, 0.9833L},
    {327.59L, 577023.6622234274946889777804063806L, 0.00006854484141433789690960193974856582L, 0.9984864940104285452360422463032948L, 0.975183324678663947357170654454763L, 0.004587630466254574190109124380022269L, 0.0220863698343264901055381831941121L, 0.8705607490043180018793099872437074L, 5.71L, 0.999L, 0.9755L},
    {344.26L, 872813.4285485332777315917743723425L, 0.00008004535308235784989663406568337875L, 0.997001583118795661765790542003923L, 0.9651775013716443866812248131038978L, 0.006689279839705071335055579858541577L, 0.03357429748940142502794147453241662L, 0.8261556480296558127857395838283497L, 8.65L, 0.9979L, 0.9662L},
    {360.93L, 1276364.659140507478194752351523234L, 0.00009729058410623347720233512035998895L, 0.9943772366573021195204550950587515L, 0.9539569286162522859267063435825542L, 0.009460439156507461024082704111757394L, 0.05006848709604471398038915302187447L, 0.7722973773109762120565071897943959L, 12.6L, 0.9957L, 0.9561L},
    {377.59L, 1815178.719588864623776827865166275L, 0.0001223249697185207665316108079145566L, 0.9898556893884632822852138430321794L, 0.9424063897008130986044357866486701L, 0.0130531087487992309146359106319712L, 0.07413307501498326169433182829371512L, 0.7066833529302514058172399788815064L, 17.88L, 0.9915L, 0.9459L},
    {394.26L, 2524082.720795113322311437705399813L, 0.000158210976400834422553218523540456L, 0.9820203008341390905011354755058195L, 0.9320370831509181788802833903574147L, 0.01766421527080147757362141215197109L, 0.1110768934715609303167863562277189L, 0.6241595383321420573558388798658149L, 24.82L, 0.9843L, 0.9361L},
    {410.93L, 3445400.12621569395196516229129189L, 0.0002092675935300077162838444659057195L, 0.9673734570274493306642727872678245L, 0.926362641981655197909001404161218L, 0.02353769440816449033268510418784217L, 0.1766336987346601237250636879099011L, 0.5091486568099046872237656024517462L, 33.85L, 0.9732L, 0.9292L},
    {416.48L, 3808203.984644166090035644283266387L, 0.0002306125006318315221347679932426901L, 0.9589614056289235296150618240339575L, 0.9275083696974733284529249158180607L, 0.02582599553065140033602006864035551L, 0.2158066022854759719719101611939148L, 0.4526599588887746195970223573454892L, 37.4L, 0.9683L, 0.9302L},
};

th::Sw92Phase<double> model(bool reverse = false) {
    return th::Sw92Phase<double>::from_parameters(
        sw92_test::binary_parameters(sw92_test::n_butane, reverse));
}

Vec composition(long double x_butane, bool reverse = false) {
    const double hydrocarbon = static_cast<double>(x_butane);
    const double water = 1.0 - hydrocarbon;
    return reverse ? Vec{water, hydrocarbon} : Vec{hydrocarbon, water};
}

double common_mu_norm(const fl::Sw92FamilySelectedPhase& first,
                      std::span<const double> first_x,
                      const fl::Sw92FamilySelectedPhase& second,
                      std::span<const double> second_x) {
    require(first_x.size() == second_x.size() &&
                first.activity.ln_phi.size() == first_x.size() &&
                second.activity.ln_phi.size() == second_x.size(),
            "phase/activity dimension mismatch");
    double norm = 0.0;
    for (std::size_t i = 0; i < first_x.size(); ++i) {
        const double first_mu = std::log(first_x[i]) + first.activity.ln_phi[i];
        const double second_mu = std::log(second_x[i]) + second.activity.ln_phi[i];
        norm = std::max(norm, std::abs(first_mu - second_mu));
    }
    return norm;
}

void model_matches_decimal_and_experiment() {
    const auto m = model(false);
    fl::Sw92FamilyStabilityEvaluator aq(
        m, 0.0, th::SwPhaseFamily::aqueous);
    fl::Sw92FamilyStabilityEvaluator na(
        m, 0.0, th::SwPhaseFamily::nonaqueous);

    double max_pressure_relative = 0.0;
    double max_h0_absolute = 0.0;
    double max_h1_absolute = 0.0;
    for (const auto& point : golden) {
        const double pressure = static_cast<double>(point.pressure_pa);
        const double temperature = static_cast<double>(point.temperature_k);
        const Vec w = composition(point.x_butane_w);
        const Vec h0 = composition(point.x_butane_h0);
        const Vec h1 = composition(point.x_butane_h1);
        const auto phase_w = aq.evaluate_selected(pressure, temperature, w);
        const auto phase_h0 = na.evaluate_selected(pressure, temperature, h0);
        const auto phase_h1 = na.evaluate_selected(pressure, temperature, h1);

        near(phase_w.compressibility_factor, point.z_w);
        near(phase_h0.compressibility_factor, point.z_h0);
        near(phase_h1.compressibility_factor, point.z_h1);
        require(phase_w.activity.smooth && phase_h0.activity.smooth &&
                    phase_h1.activity.smooth,
                "experimental triple-line point lost smooth minimum-Gibbs roots");
        require(common_mu_norm(phase_w, w, phase_h0, h0) <= 2.0e-8 &&
                    common_mu_norm(phase_w, w, phase_h1, h1) <= 2.0e-8,
                "experimental triple-line oracle no longer lies on one common tangent");

        const double predicted_bar = pressure / 1.0e5;
        max_pressure_relative = std::max(
            max_pressure_relative,
            std::abs(predicted_bar -
                     static_cast<double>(point.experimental_pressure_bar)) /
                static_cast<double>(point.experimental_pressure_bar));
        max_h0_absolute = std::max(
            max_h0_absolute,
            std::abs(static_cast<double>(point.x_butane_h0 -
                                         point.experimental_x_butane_h0)));
        max_h1_absolute = std::max(
            max_h1_absolute,
            std::abs(static_cast<double>(point.x_butane_h1 -
                                         point.experimental_x_butane_h1)));
        require(point.x_butane_w < 0.001L && point.x_butane_h0 > point.x_butane_h1,
                "triple-line relative phase-composition ordering changed");
    }

    require(max_pressure_relative <= 0.02,
            "SW92 n-butane triple-line pressure left declared 2% experimental envelope");
    require(max_h0_absolute <= 0.01,
            "SW92 dense-NA composition left declared 0.01 experimental envelope");
    require(max_h1_absolute <= 0.005,
            "SW92 second-NA composition left declared 0.005 experimental envelope");
}

void component_permutation() {
    constexpr std::size_t index = 4U;
    const auto& point = golden[index];
    const auto normal_model = model(false);
    const auto reverse_model = model(true);
    fl::Sw92FamilyStabilityEvaluator normal_aq(
        normal_model, 0.0, th::SwPhaseFamily::aqueous);
    fl::Sw92FamilyStabilityEvaluator normal_na(
        normal_model, 0.0, th::SwPhaseFamily::nonaqueous);
    fl::Sw92FamilyStabilityEvaluator reverse_aq(
        reverse_model, 0.0, th::SwPhaseFamily::aqueous);
    fl::Sw92FamilyStabilityEvaluator reverse_na(
        reverse_model, 0.0, th::SwPhaseFamily::nonaqueous);

    const double pressure = static_cast<double>(point.pressure_pa);
    const double temperature = static_cast<double>(point.temperature_k);
    const auto check = [&](long double x_butane,
                           fl::Sw92FamilyStabilityEvaluator& forward,
                           fl::Sw92FamilyStabilityEvaluator& reverse) {
        const auto a = forward.evaluate_selected(
            pressure, temperature, composition(x_butane, false));
        const auto b = reverse.evaluate_selected(
            pressure, temperature, composition(x_butane, true));
        near(a.compressibility_factor, b.compressibility_factor, 2e-12L, 2e-13L);
        require(a.activity.ln_phi.size() == 2U && b.activity.ln_phi.size() == 2U,
                "permuted phase activity dimension changed");
        near(a.activity.ln_phi[0], b.activity.ln_phi[1], 2e-11L, 2e-12L);
        near(a.activity.ln_phi[1], b.activity.ln_phi[0], 2e-11L, 2e-12L);
    };
    check(point.x_butane_w, normal_aq, reverse_aq);
    check(point.x_butane_h0, normal_na, reverse_na);
    check(point.x_butane_h1, normal_na, reverse_na);
}

} // namespace

int main(int argc, char** argv) {
    if (argc != 2) {
        std::cerr << "expected one test-case name\n";
        return 2;
    }
    const std::string_view name{argv[1]};
    try {
        if (name == "experimental_line") model_matches_decimal_and_experiment();
        else if (name == "permutation") component_permutation();
        else throw std::invalid_argument("unknown test-case name");
    } catch (const std::exception& error) {
        std::cerr << error.what() << '\n';
        return 1;
    }
    return 0;
}
