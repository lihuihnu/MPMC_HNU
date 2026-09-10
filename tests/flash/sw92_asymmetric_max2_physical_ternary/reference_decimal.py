"""Independent Decimal(80) CH4/CO2/H2O SW92/Xu max2 reference.

Stdlib only; no production imports. SW92 Table 3 supplies CH4/CO2/H2O pure
properties; corrected SW92 water-pair rules/Table 5 supply CH4-H2O and
CO2-H2O interactions. Fateen, Khalil & Elnabawy (2013), Table 1 row 34,
DOI 10.1016/j.jare.2012.03.004, supplies the external Peng-Robinson +
van-der-Waals constant CH4/CO2 binary interaction parameter kij=0.0919.
The same sourced gas/gas scalar is used in both AQ and NA family matrices.

The numerical state p=10 MPa, T=350 K, z=(0.35,0.15,0.50), fresh water is a
model-regression point, not a direct ternary experimental datum and not a
global-stability certificate.
"""

from decimal import Decimal as D, localcontext
import math
from pathlib import Path
import re

R = D("8.31446261815324")
PRESSURE = D("1e7")
TEMPERATURE = D("350")
MOLALITY = D("0")
FEED = [D("0.35"), D("0.15"), D("0.50")]
GAS_GAS_KIJ = D("0.0919")

SPEC = {
    "CH4": (D("190.6"), D("46.0e5"), D("0.0108")),
    "CO2": (D("304.2"), D("73.8e5"), D("0.2273")),
    "H2O": (D("647.3"), D("221.2e5"), D("0.3434")),
}
NA_WATER_KIJ = {"CH4": D("0.4850"), "CO2": D("0.1896")}
NAMES = ("CH4", "CO2", "H2O")


def powd(value: D, exponent: D) -> D:
    if value == 0:
        if exponent > 0:
            return D(0)
        raise ValueError("zero base outside reference power domain")
    return (exponent * value.ln()).exp()


def water_alpha(temperature: D, molality: D) -> D:
    tr = temperature / SPEC["H2O"][0]
    q = (D(1) + D("0.4530") *
         (D(1) - tr * (D(1) - D("0.0103") * powd(molality, D("1.1")))) +
         D("0.0034") * (tr ** D(-3) - D(1)))
    return q * q


def standard_alpha(name: str, temperature: D) -> D:
    tc, _, omega = SPEC[name]
    kappa = D("0.37464") + D("1.54226") * omega - D("0.26992") * omega * omega
    return (D(1) + kappa * (D(1) - (temperature / tc).sqrt())) ** 2


def pure_ab(name: str, temperature: D, molality: D):
    tc, pc, _ = SPEC[name]
    alpha = water_alpha(temperature, molality) if name == "H2O" else standard_alpha(name, temperature)
    return (D("0.45724") * R * R * tc * tc / pc * alpha,
            D("0.07780") * R * tc / pc)


def aqueous_water_kij(name: str, temperature: D, molality: D) -> D:
    tc, _, omega = SPEC[name]
    tr = temperature / tc
    if name == "CH4":
        a0 = D("1.1120") - D("1.7369") * powd(omega, D("-0.1"))
        a1 = D("1.1001") + D("0.8360") * omega
        a2 = -D("0.15742") - D("1.0988") * omega
        return (a0 * (D(1) + D("0.017407") * molality) +
                a1 * tr * (D(1) + D("0.033516") * molality) +
                a2 * tr * tr * (D(1) + D("0.011478") * molality))
    if name == "CO2":
        c1 = powd(molality, D("0.7505"))
        c2 = powd(molality, D("0.979"))
        return (-D("0.31092") * (D(1) + D("0.15587") * c1) +
                D("0.23580") * (D(1) + D("0.17837") * c2) * tr -
                D("21.2566") * (-D("6.7222") * tr - molality).exp())
    raise ValueError(name)


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
        for _ in range(200):
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
        if z > b and all(abs(z - old) > D("1e-45") for old in roots):
            roots.append(z)
    return sorted(roots)


def minimum_phase(composition, family: str):
    pure = [pure_ab(name, TEMPERATURE, MOLALITY) for name in NAMES]
    ai = [entry[0] for entry in pure]
    bi = [entry[1] for entry in pure]
    kij = [[D(0) for _ in range(3)] for _ in range(3)]
    kij[0][1] = kij[1][0] = GAS_GAS_KIJ
    if family == "AQ":
        kij[0][2] = kij[2][0] = aqueous_water_kij("CH4", TEMPERATURE, MOLALITY)
        kij[1][2] = kij[2][1] = aqueous_water_kij("CO2", TEMPERATURE, MOLALITY)
    elif family == "NA":
        kij[0][2] = kij[2][0] = NA_WATER_KIJ["CH4"]
        kij[1][2] = kij[2][1] = NA_WATER_KIJ["CO2"]
    else:
        raise ValueError(family)

    aij = [[(ai[i] * ai[j]).sqrt() * (D(1) - kij[i][j])
            for j in range(3)] for i in range(3)]
    rows = [sum(composition[j] * aij[i][j] for j in range(3)) for i in range(3)]
    attraction = sum(composition[i] * rows[i] for i in range(3))
    covolume = sum(composition[i] * bi[i] for i in range(3))
    aa = attraction * PRESSURE / (R * TEMPERATURE) ** 2
    bb = covolume * PRESSURE / (R * TEMPERATURE)
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
        for i in range(3):
            ratio = bi[i] / covolume
            expression = D(2) * rows[i] / attraction - ratio
            ln_phi.append(ratio * (z - D(1)) - (z - bb).ln() -
                          aa / (D(2) * sqrt2 * bb) * expression * log_ratio)
        root_gibbs = sum(composition[i] * ln_phi[i] for i in range(3))
        candidates.append((root_gibbs, z, ln_phi))
    if not candidates:
        raise ArithmeticError("no mechanically admissible minimum-Gibbs root")
    return min(candidates, key=lambda item: item[0])


def chemical_residual(unknown):
    x0 = [unknown[0], unknown[1], D(1) - unknown[0] - unknown[1]]
    x1 = [unknown[2], unknown[3], D(1) - unknown[2] - unknown[3]]
    beta = unknown[4]
    phase0 = minimum_phase(x0, "AQ")
    phase1 = minimum_phase(x1, "AQ")
    mu0 = [x0[i].ln() + phase0[2][i] for i in range(3)]
    mu1 = [x1[i].ln() + phase1[2][i] for i in range(3)]
    residual = [mu0[i] - mu1[i] for i in range(3)]
    residual.extend((D(1) - beta) * x0[i] + beta * x1[i] - FEED[i]
                    for i in range(2))
    return residual


def valid_unknown(unknown) -> bool:
    x0 = [unknown[0], unknown[1], D(1) - unknown[0] - unknown[1]]
    x1 = [unknown[2], unknown[3], D(1) - unknown[2] - unknown[3]]
    return (all(D(0) < value < D(1) for value in x0 + x1) and
            D(0) < unknown[4] < D(1))


def solve_linear(matrix, rhs):
    n = len(rhs)
    work = [list(matrix[i]) + [rhs[i]] for i in range(n)]
    for column in range(n):
        pivot = max(range(column, n), key=lambda row: abs(work[row][column]))
        if work[pivot][column] == 0:
            raise ArithmeticError("singular Decimal Newton matrix")
        work[column], work[pivot] = work[pivot], work[column]
        scale = work[column][column]
        for j in range(column, n + 1):
            work[column][j] /= scale
        for row in range(n):
            if row == column:
                continue
            factor = work[row][column]
            if factor == 0:
                continue
            for j in range(column, n + 1):
                work[row][j] -= factor * work[column][j]
    return [work[i][n] for i in range(n)]


def solve_pair():
    unknown = [D("0.001"), D("0.004"), D("0.695"), D("0.294"), D("0.5")]
    base_step = D("1e-20")
    for _ in range(40):
        residual = chemical_residual(unknown)
        norm = max(abs(value) for value in residual)
        if norm < D("1e-50"):
            break
        jacobian = [[D(0) for _ in range(5)] for _ in range(5)]
        for column in range(5):
            step = min(base_step, max(D("1e-28"), abs(unknown[column]) * D("1e-18")))
            plus = list(unknown)
            minus = list(unknown)
            plus[column] += step
            minus[column] -= step
            rp = chemical_residual(plus)
            rm = chemical_residual(minus)
            for row in range(5):
                jacobian[row][column] = (rp[row] - rm[row]) / (D(2) * step)
        delta = solve_linear(jacobian, [-value for value in residual])
        damping = D(1)
        accepted = False
        while damping > D("1e-12"):
            candidate = [unknown[i] + damping * delta[i] for i in range(5)]
            if valid_unknown(candidate):
                next_norm = max(abs(value) for value in chemical_residual(candidate))
                if next_norm < norm:
                    unknown = candidate
                    accepted = True
                    break
            damping /= D(2)
        if not accepted:
            raise ArithmeticError("Decimal Newton line search failed")
    else:
        raise ArithmeticError("Decimal Newton iteration limit")
    if max(abs(value) for value in chemical_residual(unknown)) > D("1e-48"):
        raise ArithmeticError("Decimal ternary pair residual did not converge")
    return unknown


def reduced_gibbs(composition, phase) -> D:
    return sum(composition[i] * (composition[i].ln() + phase[2][i]) for i in range(3))


def family_difference(composition) -> D:
    aq = minimum_phase(composition, "AQ")
    na = minimum_phase(composition, "NA")
    return sum(composition[i] * (aq[2][i] - na[2][i]) for i in range(3))


def anchors():
    unknown = solve_pair()
    x0 = [unknown[0], unknown[1], D(1) - unknown[0] - unknown[1]]
    x1 = [unknown[2], unknown[3], D(1) - unknown[2] - unknown[3]]
    beta = unknown[4]
    phase0 = minimum_phase(x0, "AQ")
    phase1 = minimum_phase(x1, "AQ")
    mu0 = [x0[i].ln() + phase0[2][i] for i in range(3)]
    mu1 = [x1[i].ln() + phase1[2][i] for i in range(3)]
    common = [(mu0[i] + mu1[i]) / D(2) for i in range(3)]
    pair_gibbs = ((D(1) - beta) * reduced_gibbs(x0, phase0) +
                  beta * reduced_gibbs(x1, phase1))
    feed_aq = minimum_phase(FEED, "AQ")
    feed_na = minimum_phase(FEED, "NA")
    feed_aq_gibbs = reduced_gibbs(FEED, feed_aq)
    feed_na_gibbs = reduced_gibbs(FEED, feed_na)

    for i in range(3):
        recovered = (D(1) - beta) * x0[i] + beta * x1[i]
        if abs(recovered - FEED[i]) > D("1e-50"):
            raise AssertionError("independent ternary material balance drifted")
    if max(abs(mu0[i] - mu1[i]) for i in range(3)) > D("1e-48"):
        raise AssertionError("independent ternary chemical potentials drifted")
    gap0 = family_difference(x0)
    gap1 = family_difference(x1)
    if not (feed_aq_gibbs < feed_na_gibbs and gap0 < 0 and gap1 < 0):
        raise AssertionError("independent ternary lower-envelope family assignment changed")
    if not pair_gibbs < feed_aq_gibbs:
        raise AssertionError("independent ternary pair no longer lowers feed Gibbs")

    return (
        x0[0], x0[1], x0[2],
        x1[0], x1[1], x1[2],
        beta, phase0[1], phase1[1],
        feed_aq_gibbs, feed_na_gibbs, feed_aq_gibbs - feed_na_gibbs,
        pair_gibbs, pair_gibbs - feed_aq_gibbs,
        common[0], common[1], common[2],
        gap0, gap1,
    )


if __name__ == "__main__":
    with localcontext() as context:
        context.prec = 80
        computed = anchors()
        cpp = Path(__file__).with_name(
            "sw92_asymmetric_max2_physical_ternary_test.cpp").read_text(encoding="utf-8")
        block = cpp.split("constexpr TernaryGolden golden{", 1)[1].split("};", 1)[0]
        literals = [D(value) for value in re.findall(
            r"([-+]?[0-9]+\.[0-9]*(?:[eE][-+]?[0-9]+)?)L", block)]
        if len(literals) != len(computed):
            raise AssertionError(
                f"physical ternary golden count changed: {len(literals)} != {len(computed)}")
        for literal, value in zip(literals, computed):
            if abs(literal - value) > D("2e-30") * max(D(1), abs(value)):
                raise AssertionError(f"physical ternary anchor mismatch: {literal} vs {value}")

        print("Independent Decimal(80) SW92/Xu CH4-CO2-H2O ternary reference passed")
        print("water-rich =", *(format(value, ".35g") for value in computed[0:3]))
        print("gas-rich   =", *(format(value, ".35g") for value in computed[3:6]))
        print("gas-rich beta =", format(computed[6], ".35g"))
        print("Gpair-Gfeed =", format(computed[13], ".35g"))
