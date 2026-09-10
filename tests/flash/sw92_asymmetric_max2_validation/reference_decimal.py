"""Independent Decimal(80) validation expansion for SW92 Xu max2.

No production code is imported. The script reuses only the already-independent
stdlib Decimal SW92 equation implementation from the one-family VLE reference,
then adds minimum-Gibbs root selection, arbitrary AQ/NA fixed-family pair solves,
Xu lower-envelope family comparisons, pair/feed Gibbs and common-tangent values.

Sources for constants/correlations remain Soreide-Whitson 1992 Table 3, Table 5,
corrected Eq. (12), and Eq. (9). These are model-regression anchors, not direct
experimental measurements and not a global-stability certificate.
"""

from decimal import Decimal as D, localcontext
from pathlib import Path
import re
import runpy

BASE = Path(__file__).resolve().parents[1] / "sw92_family_vle" / "reference_decimal.py"
ref = runpy.run_path(str(BASE), run_name="sw92_family_vle_reference")
R = ref["R"]
SPEC = ref["SPEC"]
NA_KIJ = ref["NA_KIJ"]
pure_ab = ref["pure_ab"]
aqueous_kij = ref["aqueous_kij"]
cubic_roots = ref["cubic_roots"]


def minimum_phase(gas_fraction: D, gas: str, family: str,
                  pressure: D, temperature: D, molality: D):
    composition = [gas_fraction, D(1) - gas_fraction]
    pure = [pure_ab(gas, temperature, molality), pure_ab("H2O", temperature, molality)]
    ai = [item[0] for item in pure]
    bi = [item[1] for item in pure]
    kij = NA_KIJ[gas] if family == "NA" else aqueous_kij(gas, temperature, molality)

    aij = [[D(0), D(0)], [D(0), D(0)]]
    for i in range(2):
        for j in range(2):
            pair_kij = D(0) if i == j else kij
            aij[i][j] = (ai[i] * ai[j]).sqrt() * (D(1) - pair_kij)
    rows = [sum(composition[j] * aij[i][j] for j in range(2)) for i in range(2)]
    attraction = sum(composition[i] * rows[i] for i in range(2))
    covolume = sum(composition[i] * bi[i] for i in range(2))
    aa = attraction * pressure / (R * temperature) ** 2
    bb = covolume * pressure / (R * temperature)
    c2 = bb - D(1)
    c1 = aa - D(3) * bb * bb - D(2) * bb
    sqrt2 = D(2).sqrt()

    candidates = []
    for z in cubic_roots(aa, bb):
        derivative = (D(3) * z + D(2) * c2) * z + c1
        if derivative <= 0:
            continue
        log_ratio = ((z + (D(1) + sqrt2) * bb) /
                     (z + (D(1) - sqrt2) * bb)).ln()
        ln_phi = []
        for i in range(2):
            ratio = bi[i] / covolume
            expression = D(2) * rows[i] / attraction - ratio
            ln_phi.append(ratio * (z - D(1)) - (z - bb).ln() -
                          aa / (D(2) * sqrt2 * bb) * expression * log_ratio)
        root_gibbs = sum(composition[i] * ln_phi[i] for i in range(2))
        candidates.append((root_gibbs, z, ln_phi))
    if not candidates:
        raise ArithmeticError("no mechanically admissible minimum-Gibbs reference root")
    return min(candidates, key=lambda item: item[0])


def pair_residuals(x: D, y: D, gas: str, family0: str, family1: str,
                   pressure: D, temperature: D, molality: D):
    phase0 = minimum_phase(x, gas, family0, pressure, temperature, molality)
    phase1 = minimum_phase(y, gas, family1, pressure, temperature, molality)
    return (x.ln() + phase0[2][0] - y.ln() - phase1[2][0],
            (D(1) - x).ln() + phase0[2][1] -
            (D(1) - y).ln() - phase1[2][1])


def solve_pair(gas: str, family0: str, family1: str,
               pressure: D, temperature: D, molality: D,
               x_seed: str, y_seed: str):
    x = D(x_seed)
    y = D(y_seed)
    h0 = D("1e-24")
    for _ in range(70):
        r0, r1 = pair_residuals(x, y, gas, family0, family1,
                                pressure, temperature, molality)
        if max(abs(r0), abs(r1)) < D("1e-55"):
            break
        hx = min(h0, x / D(10), (D(1) - x) / D(10))
        hy = min(h0, y / D(10), (D(1) - y) / D(10))
        xp = pair_residuals(x + hx, y, gas, family0, family1,
                            pressure, temperature, molality)
        xm = pair_residuals(x - hx, y, gas, family0, family1,
                            pressure, temperature, molality)
        yp = pair_residuals(x, y + hy, gas, family0, family1,
                            pressure, temperature, molality)
        ym = pair_residuals(x, y - hy, gas, family0, family1,
                            pressure, temperature, molality)
        j00 = (xp[0] - xm[0]) / (D(2) * hx)
        j10 = (xp[1] - xm[1]) / (D(2) * hx)
        j01 = (yp[0] - ym[0]) / (D(2) * hy)
        j11 = (yp[1] - ym[1]) / (D(2) * hy)
        determinant = j00 * j11 - j01 * j10
        dx = (-r0 * j11 + j01 * r1) / determinant
        dy = (j10 * r0 - j00 * r1) / determinant
        damping = D(1)
        while not (D(0) < x + damping * dx < D(1) and
                   D(0) < y + damping * dy < D(1)):
            damping /= D(2)
        x += damping * dx
        y += damping * dy
    else:
        raise ArithmeticError("fixed-family-pair Decimal Newton iteration limit")
    residual = pair_residuals(x, y, gas, family0, family1,
                              pressure, temperature, molality)
    if max(abs(residual[0]), abs(residual[1])) > D("1e-50"):
        raise ArithmeticError("fixed-family-pair Decimal residual did not converge")
    return x, y


def family_difference(gas_fraction: D, gas: str,
                      pressure: D, temperature: D, molality: D) -> D:
    composition = [gas_fraction, D(1) - gas_fraction]
    aq = minimum_phase(gas_fraction, gas, "AQ", pressure, temperature, molality)
    na = minimum_phase(gas_fraction, gas, "NA", pressure, temperature, molality)
    return sum(composition[i] * (aq[2][i] - na[2][i]) for i in range(2))


def case(pressure_value: str, temperature_value: str,
         molality_value: str, feed_value: str):
    gas = "CH4"
    pressure = D(pressure_value)
    temperature = D(temperature_value)
    molality = D(molality_value)
    feed = D(feed_value)

    x, y = solve_pair(gas, "AQ", "AQ", pressure, temperature, molality,
                      "0.0005", "0.99")
    beta = (feed - x) / (y - x)
    phase0 = minimum_phase(x, gas, "AQ", pressure, temperature, molality)
    phase1 = minimum_phase(y, gas, "AQ", pressure, temperature, molality)
    mu0 = [x.ln() + phase0[2][0], (D(1) - x).ln() + phase0[2][1]]
    mu1 = [y.ln() + phase1[2][0], (D(1) - y).ln() + phase1[2][1]]
    common = [(mu0[i] + mu1[i]) / D(2) for i in range(2)]
    pair_gibbs = ((D(1) - beta) *
                  (x * mu0[0] + (D(1) - x) * mu0[1]) +
                  beta * (y * mu1[0] + (D(1) - y) * mu1[1]))

    feed_difference = family_difference(feed, gas, pressure, temperature, molality)
    if feed_difference >= 0:
        raise AssertionError("validation reference unexpectedly lost lower AQ feed family")
    feed_phase = minimum_phase(feed, gas, "AQ", pressure, temperature, molality)
    feed_gibbs = (feed * (feed.ln() + feed_phase[2][0]) +
                  (D(1) - feed) * ((D(1) - feed).ln() + feed_phase[2][1]))

    aqna_x, aqna_y = solve_pair(gas, "AQ", "NA", pressure, temperature, molality,
                                "0.0005", "0.995")

    return (x, y, beta, phase0[1], phase1[1],
            feed_gibbs, pair_gibbs, pair_gibbs - feed_gibbs,
            common[0], common[1],
            family_difference(x, gas, pressure, temperature, molality),
            family_difference(y, gas, pressure, temperature, molality),
            aqna_x, aqna_y)


def anchors():
    return [
        case("1e7", "350", "0", "0.5"),
        case("2e7", "376.15", "4", "0.5"),
    ]


if __name__ == "__main__":
    with localcontext() as context:
        context.prec = 80
        computed = anchors()

        # Finite family-envelope audit only: at these traceable binary states AQ
        # is lower at every sampled interior composition. This is not a proof for
        # all states and does not turn model-family identity into a physical label.
        state_args = [(D("1e7"), D("350"), D("0")),
                      (D("2e7"), D("376.15"), D("4"))]
        for pressure, temperature, molality in state_args:
            for k in range(1, 100):
                gap = family_difference(D(k) / D(100), "CH4",
                                        pressure, temperature, molality)
                if not gap < 0:
                    raise AssertionError(
                        f"sampled CH4 family envelope changed sign at x={k}/100: {gap}")

        cpp = Path(__file__).with_name("sw92_asymmetric_max2_validation_test.cpp").read_text(
            encoding="utf-8")
        block = cpp.split("constexpr BinaryGolden golden[] = {", 1)[1].split("};", 1)[0]
        literals = [D(value) for value in re.findall(
            r"([-+]?[0-9]+\.[0-9]*(?:[eE][-+]?[0-9]+)?)L", block)]
        expected_count = 14 * len(computed)
        if len(literals) != expected_count:
            raise AssertionError(
                f"max2 validation anchor count changed: {len(literals)} != {expected_count}")
        cursor = 0
        for index, values in enumerate(computed):
            for value in values:
                literal = literals[cursor]
                cursor += 1
                if abs(literal - value) > D("2e-31") * max(D(1), abs(value)):
                    raise AssertionError(
                        f"max2 validation anchor mismatch: {literal} vs {value}")
            print("SW92 max2 validation case", index + 1,
                  "x0=", format(values[0], ".35g"),
                  "x1=", format(values[1], ".35g"),
                  "beta1=", format(values[2], ".35g"),
                  "Gpair-Gfeed=", format(values[7], ".35g"),
                  "AQ+NA x0/x1=", format(values[12], ".35g"),
                  format(values[13], ".35g"))
        print("Independent Decimal(80) CH4/fresh+brine max2 validation passed")
