"""Independent Decimal(80) SW92 fixed-family TPD reference.

No production imports. Raw CO2/water values are transcribed from the same
Soreide-Whitson 1992 corrected-original source used by the SW92 thermodynamics
regression. The script independently evaluates water alpha, CO2/water BIPs,
classical PR mixing, the original-Z cubic, same-family mechanical/Gibbs root
selection, and one explicit tangent-plane distance.
"""
from decimal import Decimal as D, localcontext
import math
from pathlib import Path
import re

R = D("8.31446261815324")
SQRT2 = D(2).sqrt()
TC_CO2 = D("304.2")
PC_CO2 = D("73.8e5")
OMEGA_CO2 = D("0.2273")
TC_WATER = D("647.3")
PC_WATER = D("221.2e5")
OMEGA_WATER = D("0.3434")
KNA_CO2_WATER = D("0.1896")


def powd(value: D, exponent: D) -> D:
    if value == 0:
        if exponent > 0:
            return D(0)
        raise ValueError("zero base outside reference domain")
    return (exponent * value.ln()).exp()


def water_alpha(t: D, molality: D) -> D:
    tr = t / TC_WATER
    root = (D(1) + D("0.4530") *
            (D(1) - tr * (D(1) - D("0.0103") * powd(molality, D("1.1")))) +
            D("0.0034") * (tr ** D(-3) - D(1)))
    return root * root


def pure_ab(name: str, t: D, molality: D):
    if name == "co2":
        tc, pc, omega = TC_CO2, PC_CO2, OMEGA_CO2
        kappa = D("0.37464") + D("1.54226") * omega - D("0.26992") * omega * omega
        alpha = (D(1) + kappa * (D(1) - (t / tc).sqrt())) ** 2
    elif name == "water":
        tc, pc, omega = TC_WATER, PC_WATER, OMEGA_WATER
        alpha = water_alpha(t, molality)
    else:
        raise ValueError(name)
    ac = D("0.45724") * R * R * tc * tc / pc
    b = D("0.07780") * R * tc / pc
    return ac * alpha, b


def all_roots(a: D, b: D):
    c2 = b - D(1)
    c1 = a - D(3) * b * b - D(2) * b
    c0 = b ** 3 + b * b - a * b
    q = (D(3) * c1 - c2 * c2) / D(9)
    r = (D(9) * c2 * c1 - D(27) * c0 - D(2) * c2 ** 3) / D(54)
    disc = q ** 3 + r * r
    if disc >= 0:
        sd = math.sqrt(float(disc))

        def cbrt(value: float) -> float:
            return math.copysign(abs(value) ** (1.0 / 3.0), value)

        seeds = [cbrt(float(r) + sd) + cbrt(float(r) - sd) - float(c2) / 3.0]
    else:
        theta = math.acos(float(r / (-q ** 3).sqrt()))
        seeds = [2.0 * math.sqrt(float(-q)) *
                 math.cos((theta + 2.0 * k * math.pi) / 3.0) - float(c2) / 3.0
                 for k in range(3)]

    values = []
    for seed in seeds:
        z = D(str(seed))
        for _ in range(80):
            f = ((z + c2) * z + c1) * z + c0
            derivative = (D(3) * z + D(2) * c2) * z + c1
            next_z = z - f / derivative
            if abs(next_z - z) < D("1e-70"):
                z = next_z
                break
            z = next_z
        residual = ((z + c2) * z + c1) * z + c0
        scale = D(1) + abs(z) ** 3 + abs(c2 * z * z) + abs(c1 * z) + abs(c0)
        if abs(residual) > D("1e-68") * scale:
            raise ArithmeticError("reference root did not converge")
        if z > b and all(abs(z - old) > D("1e-50") for old in values):
            values.append(z)
    return sorted(values)


def phase(composition, t: D, molality: D, p: D):
    names = ("co2", "water")
    x = [D(str(value)) for value in composition]
    pure = [pure_ab(name, t, molality) for name in names]
    ai = [item[0] for item in pure]
    bi = [item[1] for item in pure]
    rows = [D(0), D(0)]
    attraction = D(0)
    for i in range(2):
        for j in range(2):
            kij = D(0) if i == j else KNA_CO2_WATER
            aij = (ai[i] * ai[j]).sqrt() * (D(1) - kij)
            rows[i] += x[j] * aij
            attraction += x[i] * x[j] * aij
    covolume = x[0] * bi[0] + x[1] * bi[1]
    aa = attraction * p / (R * t) ** 2
    bb = covolume * p / (R * t)
    c2 = bb - D(1)
    c1 = aa - D(3) * bb * bb - D(2) * bb

    candidates = []
    for root_index, z in enumerate(all_roots(aa, bb)):
        derivative = (D(3) * z + D(2) * c2) * z + c1
        if derivative <= 0:
            continue
        log_ratio = ((z + (D(1) + SQRT2) * bb) /
                     (z + (D(1) - SQRT2) * bb)).ln()
        ln_phi = []
        for i in range(2):
            ratio = bi[i] / covolume
            expression = D(2) * rows[i] / attraction - ratio
            ln_phi.append(ratio * (z - D(1)) - (z - bb).ln() -
                          aa / (D(2) * SQRT2 * bb) * expression * log_ratio)
        gibbs = x[0] * ln_phi[0] + x[1] * ln_phi[1]
        candidates.append((gibbs, root_index, z, ln_phi))
    if not candidates:
        raise ArithmeticError("no mechanically admissible reference root")
    return min(candidates, key=lambda item: item[0])


def tpd(feed, trial, t: D, molality: D, p: D):
    z = [D(str(value)) for value in feed]
    w = [D(str(value)) for value in trial]
    reference = phase(feed, t, molality, p)
    candidate = phase(trial, t, molality, p)
    value = sum((w[i] * ((w[i] / z[i]).ln() + candidate[3][i] - reference[3][i])
                 for i in range(2)), D(0))
    return value, reference, candidate


if __name__ == "__main__":
    with localcontext() as context:
        context.prec = 80
        value, reference, candidate = tpd(
            ["0.7", "0.3"], ["0.5", "0.5"], D("340"), D(0), D("3e6"))
        if reference[1] != 2 or candidate[1] != 0:
            raise AssertionError(
                f"same-family Gibbs branch selection changed: {reference[1]} -> {candidate[1]}")

        text = Path(__file__).with_name("sw92_stability_test.cpp").read_text(encoding="utf-8")
        block = text.split("constexpr long double sw92_tpd_golden[] = {", 1)[1].split("};", 1)[0]
        literals = re.findall(r"([-+]?[0-9]+\.[0-9]*(?:[eE][-+]?[0-9]+)?)L", block)
        if len(literals) != 1:
            raise AssertionError("SW92 stability anchor count changed")
        literal = D(literals[0])
        if abs(literal - value) > D("1e-35") * max(D(1), abs(value)):
            raise AssertionError(f"SW92 TPD anchor mismatch: {literal} vs {value}")

        print("SW92 NA feed minimum-Gibbs root index:", reference[1])
        print("SW92 NA trial minimum-Gibbs root index:", candidate[1])
        print("SW92 NA fixed-family TPD:", format(value, ".40g"))
        print("Independent Decimal(80) SW92 stability anchors: 1/1 checked")
