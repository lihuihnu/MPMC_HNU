"""Independent Decimal(80) anchors for SW92 Xu Gate 3B.1 fixed pairs.

Stdlib only; no production imports. The script independently rebuilds the
corrected-original CO2/H2O SW92 pure/mixing equations, PR cubic roots and
minimum-Gibbs root choice, then solves the two common-reduced-chemical-potential
equations for three explicit family assignments: AQ+AQ, AQ+NA and NA+NA.

The returned values are equation/model anchors. A solved fixed-family pair is
not an accepted global phase set: the script separately evaluates AQ-minus-NA
Gibbs ordering at each converged phase so the C++ lower-envelope dominance
checks can reject a solved but dominated assignment.
"""

from decimal import Decimal as D, localcontext
import math
from pathlib import Path
import re

R = D("8.31446261815324")
TC_CO2 = D("304.2")
PC_CO2 = D("73.8e5")
OMEGA_CO2 = D("0.2273")
TC_WATER = D("647.3")
PC_WATER = D("221.2e5")
OMEGA_WATER = D("0.3434")
NA_KIJ = D("0.1896")
PRESSURE = D("3e6")
TEMPERATURE = D("340")
MOLALITY = D(0)
FEED = D("0.7")


def powd(value: D, exponent: D) -> D:
    if value == 0:
        if exponent > 0:
            return D(0)
        raise ValueError("zero base outside reference domain")
    return (exponent * value.ln()).exp()


def water_alpha(temperature: D, molality: D) -> D:
    tr = temperature / TC_WATER
    root = (D(1) + D("0.4530") *
            (D(1) - tr * (D(1) - D("0.0103") * powd(molality, D("1.1")))) +
            D("0.0034") * (tr ** D(-3) - D(1)))
    return root * root


def co2_alpha(temperature: D) -> D:
    kappa = (D("0.37464") + D("1.54226") * OMEGA_CO2 -
             D("0.26992") * OMEGA_CO2 * OMEGA_CO2)
    return (D(1) + kappa * (D(1) - (temperature / TC_CO2).sqrt())) ** 2


def pure_ab(temperature: D, molality: D):
    co2_a = (D("0.45724") * R * R * TC_CO2 * TC_CO2 /
             PC_CO2 * co2_alpha(temperature))
    co2_b = D("0.07780") * R * TC_CO2 / PC_CO2
    water_a = (D("0.45724") * R * R * TC_WATER * TC_WATER /
               PC_WATER * water_alpha(temperature, molality))
    water_b = D("0.07780") * R * TC_WATER / PC_WATER
    return [(co2_a, co2_b), (water_a, water_b)]


def aqueous_kij(temperature: D, molality: D) -> D:
    tr = temperature / TC_CO2
    return (-D("0.31092") *
            (D(1) + D("0.15587") * powd(molality, D("0.7505"))) +
            D("0.23580") *
            (D(1) + D("0.17837") * powd(molality, D("0.979"))) * tr -
            D("21.2566") * (-D("6.7222") * tr - molality).exp())


def cubic_roots(a: D, b: D):
    c2 = b - D(1)
    c1 = a - D(3) * b * b - D(2) * b
    c0 = b ** 3 + b * b - a * b
    q = (D(3) * c1 - c2 * c2) / D(9)
    r = (D(9) * c2 * c1 - D(27) * c0 - D(2) * c2 ** 3) / D(54)
    discriminant = q ** 3 + r * r
    if discriminant >= 0:
        sd = math.sqrt(float(discriminant))

        def cbrt(value: float) -> float:
            return math.copysign(abs(value) ** (1.0 / 3.0), value)

        seeds = [cbrt(float(r) + sd) + cbrt(float(r) - sd) - float(c2) / 3.0]
    else:
        theta = math.acos(float(r / (-q ** 3).sqrt()))
        seeds = [2.0 * math.sqrt(float(-q)) *
                 math.cos((theta + 2.0 * k * math.pi) / 3.0) - float(c2) / 3.0
                 for k in range(3)]

    roots = []
    for seed in seeds:
        z = D(str(seed))
        for _ in range(180):
            residual = ((z + c2) * z + c1) * z + c0
            derivative = (D(3) * z + D(2) * c2) * z + c1
            next_z = z - residual / derivative
            if abs(next_z - z) < D("1e-70"):
                z = next_z
                break
            z = next_z
        residual = ((z + c2) * z + c1) * z + c0
        scale = D(1) + abs(z) ** 3 + abs(c2 * z * z) + abs(c1 * z) + abs(c0)
        if abs(residual) > D("1e-67") * scale:
            raise ArithmeticError("reference cubic root did not converge")
        if z > b and all(abs(z - old) > D("1e-48") for old in roots):
            roots.append(z)
    return sorted(roots)


def family_candidates(gas_fraction: D, family: str, pressure=PRESSURE,
                      temperature=TEMPERATURE, molality=MOLALITY,
                      na_kij=NA_KIJ):
    composition = [gas_fraction, D(1) - gas_fraction]
    pure = pure_ab(temperature, molality)
    ai = [entry[0] for entry in pure]
    bi = [entry[1] for entry in pure]
    kij = aqueous_kij(temperature, molality) if family == "AQ" else na_kij

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
        candidates.append((z, ln_phi, root_gibbs))
    if not candidates:
        raise ArithmeticError("no mechanically admissible reference root")
    return candidates


def minimum_phase(gas_fraction: D, family: str, na_kij=NA_KIJ):
    return min(family_candidates(gas_fraction, family, na_kij=na_kij),
               key=lambda item: item[2])


def pair_residuals(x: D, y: D, family0: str, family1: str):
    phase0 = minimum_phase(x, family0)
    phase1 = minimum_phase(y, family1)
    return (x.ln() + phase0[1][0] - y.ln() - phase1[1][0],
            (D(1) - x).ln() + phase0[1][1] -
            (D(1) - y).ln() - phase1[1][1])


def solve_pair(family0: str, family1: str, x_seed: str, y_seed: str):
    x = D(x_seed)
    y = D(y_seed)
    finite_difference_step = D("1e-24")
    for _ in range(60):
        r0, r1 = pair_residuals(x, y, family0, family1)
        if max(abs(r0), abs(r1)) < D("1e-55"):
            break
        hx = min(finite_difference_step, x / D(10), (D(1) - x) / D(10))
        hy = min(finite_difference_step, y / D(10), (D(1) - y) / D(10))
        x_plus = pair_residuals(x + hx, y, family0, family1)
        x_minus = pair_residuals(x - hx, y, family0, family1)
        y_plus = pair_residuals(x, y + hy, family0, family1)
        y_minus = pair_residuals(x, y - hy, family0, family1)
        j00 = (x_plus[0] - x_minus[0]) / (D(2) * hx)
        j10 = (x_plus[1] - x_minus[1]) / (D(2) * hx)
        j01 = (y_plus[0] - y_minus[0]) / (D(2) * hy)
        j11 = (y_plus[1] - y_minus[1]) / (D(2) * hy)
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
        raise ArithmeticError("reference fixed-pair Newton iteration limit")

    residual = pair_residuals(x, y, family0, family1)
    if max(abs(residual[0]), abs(residual[1])) > D("1e-50"):
        raise ArithmeticError("reference fixed-pair residual did not converge")
    return x, y


def family_difference(gas_fraction: D) -> D:
    composition = [gas_fraction, D(1) - gas_fraction]
    aq = minimum_phase(gas_fraction, "AQ")
    na = minimum_phase(gas_fraction, "NA")
    return sum(composition[i] * (aq[1][i] - na[1][i]) for i in range(2))


def pair_anchor(family0: str, family1: str, x_seed: str, y_seed: str):
    x, y = solve_pair(family0, family1, x_seed, y_seed)
    beta = (FEED - x) / (y - x)
    if not (D(0) < beta < D(1)):
        raise AssertionError("reference pair does not bracket the feed")
    phase0 = minimum_phase(x, family0)
    phase1 = minimum_phase(y, family1)
    mu0 = [x.ln() + phase0[1][0],
           (D(1) - x).ln() + phase0[1][1]]
    mu1 = [y.ln() + phase1[1][0],
           (D(1) - y).ln() + phase1[1][1]]
    gibbs = ((D(1) - beta) *
             (x * mu0[0] + (D(1) - x) * mu0[1]) +
             beta * (y * mu1[0] + (D(1) - y) * mu1[1]))
    return (x, y, beta, phase0[0], phase1[0], gibbs,
            family_difference(x), family_difference(y))


def anchors():
    return [
        pair_anchor("AQ", "AQ", "0.006", "0.987"),
        pair_anchor("AQ", "NA", "0.006", "0.989"),
        pair_anchor("NA", "NA", "0.00028", "0.989"),
    ]


if __name__ == "__main__":
    with localcontext() as context:
        context.prec = 80
        computed = anchors()
        if not (computed[0][6] < 0 and computed[0][7] < 0):
            raise AssertionError("AQ+AQ reference lost lower-envelope AQ assignment")
        if not (computed[1][6] < 0 and computed[1][7] < 0):
            raise AssertionError("AQ+NA reference lost expected NA dominance evidence")
        if not (computed[2][6] < 0 and computed[2][7] < 0):
            raise AssertionError("NA+NA reference lost expected AQ dominance evidence")

        text = Path(__file__).with_name("sw92_asymmetric_fixed_pair_test.cpp").read_text(
            encoding="utf-8")
        block = text.split("constexpr PairGolden pair_golden[] = {", 1)[1].split("};", 1)[0]
        literals = [D(value) for value in re.findall(
            r"([-+]?[0-9]+\.[0-9]*(?:[eE][-+]?[0-9]+)?)L", block)]
        expected_count = 8 * len(computed)
        if len(literals) != expected_count:
            raise AssertionError(
                f"fixed-pair golden count changed: {len(literals)} != {expected_count}")
        cursor = 0
        for index, values in enumerate(computed):
            for value in values:
                literal = literals[cursor]
                cursor += 1
                if abs(literal - value) > D("1e-31") * max(D(1), abs(value)):
                    raise AssertionError(
                        f"fixed-pair anchor mismatch: {literal} vs {value}")
            print("SW92 fixed-pair case", index + 1,
                  "x0=", format(values[0], ".35g"),
                  "x1=", format(values[1], ".35g"),
                  "beta1=", format(values[2], ".35g"),
                  "G/RT=", format(values[5], ".35g"))

        kij_match = re.search(
            r"aq_co2_kij_340k_fresh\s*=\s*\n?\s*([-+]?[0-9.]+);", text)
        if kij_match is None:
            raise AssertionError("synthetic family-tie kij literal missing")
        kij_literal = D(kij_match.group(1))
        kij_reference = aqueous_kij(TEMPERATURE, MOLALITY)
        if abs(kij_literal - kij_reference) > D("1e-16"):
            raise AssertionError(
                f"synthetic family-tie kij drift: {kij_literal} vs {kij_reference}")

        print(f"Independent Decimal(80) SW92 fixed-pair anchors: {len(computed)}/{len(computed)} checked")
