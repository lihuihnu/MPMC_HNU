"""Independent Decimal(80) anchors for SW92 Xu Gate 3B.3 final acceptance.

This script imports no production code. It reuses the previously independent
Gate-3B.1 Decimal implementation of corrected SW92/PR equations and fixed-pair
Newton solve, then independently constructs the selected AQ+AQ pair common
reduced tangent, compares pair Gibbs against the lower-envelope feed Gibbs, and
evaluates prescribed AQ/NA TPD values against that exact tangent.

The values are model/numerical anchors. The finite probes below are not an
interval/global phase-stability proof.
"""

from decimal import Decimal as D, localcontext
from pathlib import Path
import re
import runpy

FIXED_PAIR_REFERENCE = (
    Path(__file__).resolve().parents[1]
    / "sw92_asymmetric_fixed_pair"
    / "reference_decimal.py"
)
reference = runpy.run_path(
    str(FIXED_PAIR_REFERENCE), run_name="sw92_fixed_pair_reference"
)

minimum_phase = reference["minimum_phase"]
pair_anchor = reference["pair_anchor"]
FEED = reference["FEED"]


def reduced_gibbs(composition_gas: D, phase) -> D:
    x = [composition_gas, D(1) - composition_gas]
    return sum(x[i] * (x[i].ln() + phase[1][i]) for i in range(2))


def pair_common_tangent():
    values = pair_anchor("AQ", "AQ", "0.006", "0.987")
    x, y, beta = values[0], values[1], values[2]
    phase0 = minimum_phase(x, "AQ")
    phase1 = minimum_phase(y, "AQ")
    mu0 = [x.ln() + phase0[1][0],
           (D(1) - x).ln() + phase0[1][1]]
    mu1 = [y.ln() + phase1[1][0],
           (D(1) - y).ln() + phase1[1][1]]
    residual = [mu0[i] - mu1[i] for i in range(2)]
    common = [(mu0[i] + mu1[i]) / D(2) for i in range(2)]
    pair_gibbs = ((D(1) - beta) *
                  (x * mu0[0] + (D(1) - x) * mu0[1]) +
                  beta * (y * mu1[0] + (D(1) - y) * mu1[1]))
    return x, y, beta, pair_gibbs, common, residual


def imposed_tpd(gas_fraction: D, family: str, common) -> D:
    phase = minimum_phase(gas_fraction, family)
    w = [gas_fraction, D(1) - gas_fraction]
    return sum(
        w[i] * (w[i].ln() + phase[1][i] - common[i])
        for i in range(2)
    )


def final_anchors():
    x, y, _, pair_gibbs, common, residual = pair_common_tangent()
    feed_phase = minimum_phase(FEED, "AQ")
    feed_gibbs = reduced_gibbs(FEED, feed_phase)
    return (
        feed_gibbs,
        pair_gibbs,
        pair_gibbs - feed_gibbs,
        common[0],
        common[1],
        imposed_tpd(x, "NA", common),
        imposed_tpd(y, "NA", common),
        imposed_tpd(D("0.5"), "AQ", common),
        imposed_tpd(D("0.5"), "NA", common),
        x,
        y,
        residual,
        common,
    )


if __name__ == "__main__":
    with localcontext() as context:
        context.prec = 80
        values = final_anchors()
        anchors = values[:9]
        x, y, residual, common = values[9], values[10], values[11], values[12]

        if max(abs(value) for value in residual) > D("1e-50"):
            raise AssertionError("independent AQ+AQ common chemical potentials drifted")
        if abs(imposed_tpd(x, "AQ", common)) > D("1e-50") or \
           abs(imposed_tpd(y, "AQ", common)) > D("1e-50"):
            raise AssertionError("independent selected AQ phases are not on common tangent")
        if anchors[2] >= D(0):
            raise AssertionError("selected reference pair is no longer below feed Gibbs")

        # Coarse deterministic probe only catches gross tangent/sign regressions.
        # It is deliberately bounded and is not a global stability certificate.
        for family in ("AQ", "NA"):
            lowest = None
            for k in range(1, 100):
                value = imposed_tpd(D(k) / D(100), family, common)
                lowest = value if lowest is None else min(lowest, value)
            if lowest < D("-1e-28"):
                raise AssertionError(
                    f"independent coarse {family} probe found unexpected negative TPD: {lowest}"
                )

        cpp = Path(__file__).with_name("sw92_asymmetric_max2_test.cpp").read_text(
            encoding="utf-8"
        )
        block = cpp.split("constexpr FinalGolden golden{", 1)[1].split("};", 1)[0]
        literals = [D(value) for value in re.findall(
            r"([-+]?[0-9]+\.[0-9]*(?:[eE][-+]?[0-9]+)?)L", block
        )]
        if len(literals) != len(anchors):
            raise AssertionError(
                f"Gate-3B.3 golden count changed: {len(literals)} != {len(anchors)}"
            )
        for literal, value in zip(literals, anchors):
            if abs(literal - value) > D("1e-31") * max(D(1), abs(value)):
                raise AssertionError(
                    f"Gate-3B.3 anchor mismatch: {literal} vs {value}"
                )

        print("SW92 Gate 3B.3 Decimal(80) reference passed")
        print("feed G/RT =", format(anchors[0], ".36g"))
        print("pair G/RT =", format(anchors[1], ".36g"))
        print("pair-feed =", format(anchors[2], ".36g"))
        print("common d =", format(anchors[3], ".36g"),
              format(anchors[4], ".36g"))
