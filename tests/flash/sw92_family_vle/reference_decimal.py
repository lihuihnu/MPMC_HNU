"""Independent Decimal(80) SW92 one-family VLE anchors.

Stdlib only; no production imports. Raw CO2/CH4/water properties and non-aqueous
BIPs are transcribed from Soreide-Whitson 1992 Tables 3/5 and the corrected
original AQ correlations. Each case independently solves fugacity equality for
a smallest-Z mechanically admissible liquid candidate and a largest-Z
mechanically admissible vapor candidate under ONE fixed family, then obtains the
phase fraction from material balance. These are model-regression anchors, not
experimental measurements or a global-stability proof.
"""

from decimal import Decimal as D, localcontext
import math
from pathlib import Path
import re

R = D("8.31446261815324")
SPEC = {
    "CO2": (D("304.2"), D("73.8e5"), D("0.2273")),
    "CH4": (D("190.6"), D("46.0e5"), D("0.0108")),
    "H2O": (D("647.3"), D("221.2e5"), D("0.3434")),
}
NA_KIJ = {"CO2": D("0.1896"), "CH4": D("0.4850")}


def powd(value: D, exponent: D) -> D:
    if value == 0:
        if exponent > 0:
            return D(0)
        raise ValueError("zero base outside reference domain")
    return (exponent * value.ln()).exp()


def water_alpha(temperature: D, molality: D) -> D:
    tr = temperature / SPEC["H2O"][0]
    root = (D(1) + D("0.4530") *
            (D(1) - tr * (D(1) - D("0.0103") * powd(molality, D("1.1")))) +
            D("0.0034") * (tr ** D(-3) - D(1)))
    return root * root


def standard_alpha(gas: str, temperature: D) -> D:
    tc, _, omega = SPEC[gas]
    kappa = D("0.37464") + D("1.54226") * omega - D("0.26992") * omega * omega
    return (D(1) + kappa * (D(1) - (temperature / tc).sqrt())) ** 2


def pure_ab(name: str, temperature: D, molality: D):
    tc, pc, _ = SPEC[name]
    alpha = water_alpha(temperature, molality) if name == "H2O" else standard_alpha(name, temperature)
    return (D("0.45724") * R * R * tc * tc / pc * alpha,
            D("0.07780") * R * tc / pc)


def aqueous_kij(gas: str, temperature: D, molality: D) -> D:
    tc, _, omega = SPEC[gas]
    tr = temperature / tc
    if gas == "CO2":
        return (-D("0.31092") * (D(1) + D("0.15587") * powd(molality, D("0.7505"))) +
                D("0.23580") * (D(1) + D("0.17837") * powd(molality, D("0.979"))) * tr -
                D("21.2566") * (-D("6.7222") * tr - molality).exp())
    if gas == "CH4":
        a0 = D("1.1120") - D("1.7369") * powd(omega, D("-0.1"))
        a1 = D("1.1001") + D("0.8360") * omega
        a2 = -D("0.15742") - D("1.0988") * omega
        return (a0 * (D(1) + D("0.017407") * molality) +
                a1 * tr * (D(1) + D("0.033516") * molality) +
                a2 * tr * tr * (D(1) + D("0.011478") * molality))
    raise ValueError(gas)


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
        for _ in range(100):
            residual = ((z + c2) * z + c1) * z + c0
            derivative = (D(3) * z + D(2) * c2) * z + c1
            next_z = z - residual / derivative
            if abs(next_z - z) < D("1e-70"):
                z = next_z
                break
            z = next_z
        residual = ((z + c2) * z + c1) * z + c0
        scale = D(1) + abs(z) ** 3 + abs(c2 * z * z) + abs(c1 * z) + abs(c0)
        if abs(residual) > D("1e-68") * scale:
            raise ArithmeticError("reference cubic root did not converge")
        if z > b and all(abs(z - old) > D("1e-50") for old in roots):
            roots.append(z)
    return sorted(roots)


def phase(gas_fraction: D, gas: str, family: str,
          pressure: D, temperature: D, molality: D, role: str):
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
    sqrt2 = D(2).sqrt()  # evaluated in the caller's active Decimal context

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
        candidates.append((z, ln_phi))
    if not candidates:
        raise ArithmeticError("no mechanically admissible reference root")
    return candidates[0] if role == "liquid" else candidates[-1]


def residuals(x: D, y: D, gas: str, family: str,
              pressure: D, temperature: D, molality: D):
    liquid = phase(x, gas, family, pressure, temperature, molality, "liquid")
    vapor = phase(y, gas, family, pressure, temperature, molality, "vapor")
    return (x.ln() + liquid[1][0] - y.ln() - vapor[1][0],
            (D(1) - x).ln() + liquid[1][1] -
            (D(1) - y).ln() - vapor[1][1])


def solve_case(gas: str, family: str, pressure_value: str, temperature_value: str,
               molality_value: str, feed_value: str, x_seed: str, y_seed: str):
    pressure = D(pressure_value)
    temperature = D(temperature_value)
    molality = D(molality_value)
    feed = D(feed_value)
    x = D(x_seed)
    y = D(y_seed)
    finite_difference_step = D("1e-24")

    for _ in range(50):
        r0, r1 = residuals(x, y, gas, family, pressure, temperature, molality)
        if max(abs(r0), abs(r1)) < D("1e-55"):
            break
        hx = min(finite_difference_step, x / D(10), (D(1) - x) / D(10))
        hy = min(finite_difference_step, y / D(10), (D(1) - y) / D(10))
        x_plus = residuals(x + hx, y, gas, family, pressure, temperature, molality)
        x_minus = residuals(x - hx, y, gas, family, pressure, temperature, molality)
        y_plus = residuals(x, y + hy, gas, family, pressure, temperature, molality)
        y_minus = residuals(x, y - hy, gas, family, pressure, temperature, molality)
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
        raise ArithmeticError("reference VLE Newton iteration limit")

    r0, r1 = residuals(x, y, gas, family, pressure, temperature, molality)
    if max(abs(r0), abs(r1)) > D("1e-50"):
        raise ArithmeticError("reference VLE residual did not converge")
    liquid = phase(x, gas, family, pressure, temperature, molality, "liquid")
    vapor = phase(y, gas, family, pressure, temperature, molality, "vapor")
    beta = (feed - x) / (y - x)
    return x, y, beta, liquid[0], vapor[0]


def anchors():
    return [
        solve_case("CO2", "NA", "3e6", "340", "0", "0.7", "0.00028", "0.9887"),
        solve_case("CH4", "AQ", "1e7", "350", "1", "0.5", "0.0009", "0.9904"),
    ]


if __name__ == "__main__":
    with localcontext() as context:
        context.prec = 80
        computed = anchors()
        text = Path(__file__).with_name("sw92_family_vle_test.cpp").read_text(encoding="utf-8")
        block = text.split("constexpr VleGolden golden[] = {", 1)[1].split("};", 1)[0]
        literals = [D(value) for value in re.findall(
            r"([-+]?[0-9]+\.[0-9]*(?:[eE][-+]?[0-9]+)?)L", block)]
        expected_count = 5 * len(computed)
        if len(literals) != expected_count:
            raise AssertionError(f"SW92 VLE anchor count changed: {len(literals)} != {expected_count}")
        cursor = 0
        for index, values in enumerate(computed):
            for value in values:
                literal = literals[cursor]
                cursor += 1
                if abs(literal - value) > D("1e-32") * max(D(1), abs(value)):
                    raise AssertionError(f"SW92 VLE anchor mismatch: {literal} vs {value}")
            print("SW92 one-family case", index + 1,
                  "x=", format(values[0], ".35g"),
                  "y=", format(values[1], ".35g"),
                  "beta=", format(values[2], ".35g"))
        print(f"Independent Decimal(80) SW92 one-family VLE anchors: {len(computed)}/{len(computed)} checked")
