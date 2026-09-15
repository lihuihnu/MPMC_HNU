#define main mpmc_pr76_c10_kappa_original_audit_main
#include "pr76_c10_kappa_mismatch_audit_test.cpp"
#undef main

int main() {
    try {
        const auto formula_source = ufc_source(
            "Section 4.4.1, Eqs.20-21: piecewise Peng-Robinson form factor m(omega)",
            "Eq.20 is the original quadratic for omega < 0.49; Eq.21 is the cubic branch for omega >= 0.49.");
        require(formula_source.kind == th::SourceKind::literature &&
                    formula_source.reference.find("riufc/80725") != std::string::npos,
                "high-omega kappa audit lost its UFC literature provenance");
        require(decane_omega >= 0.49,
                "n-decane no longer falls on the dissertation high-omega branch");

        const double kappa_strict = strict_kappa(decane_omega);
        const double kappa_high = soria_kappa(decane_omega);
        const double equivalent_omega = model_decane_omega(KappaMode::soria_high_omega);
        const double alpha_strict_323 = alpha(323.15, decane_tc_k, kappa_strict);
        const double alpha_high_323 = alpha(323.15, decane_tc_k, kappa_high);

        const auto strict_fit = fit_kij(KappaMode::strict_pr76);
        const auto high_same_kij = evaluate_kij(
            strict_kij_from_pr101, KappaMode::soria_high_omega);
        const auto high_fit = fit_kij(KappaMode::soria_high_omega);
        const auto high_published_kij = evaluate_kij(
            soria_table12_kij, KappaMode::soria_high_omega);

        std::cout << std::setprecision(12)
                  << "C10 kappa audit: omega_physical=" << decane_omega
                  << " kappa_strict=" << kappa_strict
                  << " kappa_Soria_Eq21=" << kappa_high
                  << " relative_kappa_shift=" << (kappa_high / kappa_strict - 1.0)
                  << " equivalent_test_omega=" << equivalent_omega << '\n'
                  << "alpha_323.15_strict=" << alpha_strict_323
                  << " alpha_323.15_high_omega=" << alpha_high_323
                  << " relative_alpha_shift="
                  << (alpha_high_323 / alpha_strict_323 - 1.0) << '\n';
        print_fit("strict-PR76 refit", strict_fit);
        if (high_same_kij) {
            print_fit("high-omega with PR101 kij", *high_same_kij);
        } else {
            std::cout << "high-omega with PR101 kij=0.05226578047: "
                         "not representable on one or more of the same four incipient-VLE branches\n";
        }
        print_fit("high-omega refit", high_fit);
        if (high_published_kij) {
            print_fit("high-omega with Soria Table12 kij", *high_published_kij);
        } else {
            std::cout << "high-omega with Soria Table12 kij=0.09670: "
                         "not directly comparable on all four pressure-only VLE branches\n";
        }
        std::cout << "AARD delta high-omega-refit minus strict-refit="
                  << high_fit.aard - strict_fit.aard << '\n';

        require(std::abs(strict_fit.kij - strict_kij_from_pr101) <= 5.0e-7,
                "audit no longer reproduces the merged PR101 strict C10 kij baseline");
        require(std::abs(strict_fit.aard - 0.1324992938) <= 5.0e-6,
                "audit no longer reproduces the merged PR101 strict C10 AARD baseline");
        require(std::isfinite(high_fit.aard) && std::isfinite(high_fit.relative_sse),
                "high-omega refit produced nonfinite diagnostics");
        if (high_same_kij) {
            require(high_fit.relative_sse <= high_same_kij->relative_sse,
                    "high-omega refit is worse than keeping the independently fitted strict kij");
        }
        require(decane_323_data.size() == 4U,
                "C10 mismatch audit dataset shape changed unexpectedly");

        std::cout << "[PASS] pr76_c10_high_omega_kappa_mismatch_audit\n";
        return 0;
    } catch (const std::exception& error) {
        std::cerr << "[FAIL] " << error.what() << '\n';
        return 1;
    }
}
