"""Independent Decimal(80) references for SW92 Profile-C Gate C1.

Stdlib only; no production imports and no Profile-A/Profile-B result reuse.
The script transcribes the already-audited corrected-original Soreide-Whitson
binary model equations, selects the mechanically admissible minimum-Gibbs root
inside the *assigned* family, and directly solves the cross-family equations

    ln(x_i^AQ) + ln(phi_i^AQ) = ln(x_i^NA) + ln(phi_i^NA)

for CO2/H2O and CH4/H2O. The feed phase fraction is then obtained from the
single common material balance. These are model numerical references, not
experimental data, global-stability proofs, or autonomous phase-number results.
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


def minimum_phase(gas_fraction: D, gas: str, family: str,
                  pressure: D, temperature: D, molality: D):
    composition = [gas_fraction, D(1) - gas_fraction]
    pure = [pure_ab(gas, temperature, molality), pure_ab("H2O", temperature, molality)]
    ai = [entry[0] for entry in pure]
    bi = [entry[1] for entry in pure]
    kij = aqueous_kij(gas, temperature, molality) if family == "AQ" else NA_KIJ[gas]

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
        raise ArithmeticError("no mechanically admissible reference root")
    return min(candidates, key=lambda item: item[0])


def residuals(x: D, y: D, gas: str, pressure: D,
              temperature: D, molality: D):
    aqueous = minimum_phase(x, gas, "AQ", pressure, temperature, molality)
    nonaqueous = minimum_phase(y, gas, "NA", pressure, temperature, molality)
    return (x.ln() + aqueous[2][0] - y.ln() - nonaqueous[2][0],
            (D(1) - x).ln() + aqueous[2][1] -
            (D(1) - y).ln() - nonaqueous[2][1])


def solve_assigned_pair(gas: str, pressure: D, temperature: D, molality: D,
                        x_seed: str, y_seed: str):
    x = D(x_seed)
    y = D(y_seed)
    h0 = D("1e-24")
    for _ in range(80):
        r0, r1 = residuals(x, y, gas, pressure, temperature, molality)
        if max(abs(r0), abs(r1)) < D("1e-55"):
            break
        hx = min(h0, x / D(10), (D(1) - x) / D(10))
        hy = min(h0, y / D(10), (D(1) - y) / D(10))
        xp = residuals(x + hx, y, gas, pressure, temperature, molality)
        xm = residuals(x - hx, y, gas, pressure, temperature, molality)
        yp = residuals(x, y + hy, gas, pressure, temperature, molality)
        ym = residuals(x, y - hy, gas, pressure, temperature, molality)
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
        raise ArithmeticError("reference AQ/NA Newton iteration limit")

    final = residuals(x, y, gas, pressure, temperature, molality)
    if max(abs(final[0]), abs(final[1])) > D("1e-50"):
        raise ArithmeticError("reference AQ/NA residual did not converge")
    return x, y


def case(gas: str, pressure_value: str, temperature_value: str,
         molality_value: str, feed_value: str, x_seed: str, y_seed: str):
    pressure = D(pressure_value)
    temperature = D(temperature_value)
    molality = D(molality_value)
    feed = D(feed_value)
    x, y = solve_assigned_pair(gas, pressure, temperature, molality, x_seed, y_seed)
    beta = (feed - x) / (y - x)
    if not (D(0) < beta < D(1)):
        raise AssertionError("assigned AQ/NA coexistence does not bracket feed")
    aqueous = minimum_phase(x, gas, "AQ", pressure, temperature, molality)
    nonaqueous = minimum_phase(y, gas, "NA", pressure, temperature, molality)
    mu_aq = [x.ln() + aqueous[2][0],
              (D(1) - x).ln() + aqueous[2][1]]
    mu_na = [y.ln() + nonaqueous[2][0],
              (D(1) - y).ln() + nonaqueous[2][1]]
    common = [(mu_aq[i] + mu_na[i]) / D(2) for i in range(2)]
    pair_gibbs = ((D(1) - beta) *
                  (x * mu_aq[0] + (D(1) - x) * mu_aq[1]) +
                  beta * (y * mu_na[0] + (D(1) - y) * mu_na[1]))
    water_difference = (D(1) - x) - (D(1) - y)
    return (x, y, beta, aqueous[1], nonaqueous[1], pair_gibbs,
            common[0], common[1], water_difference)


def anchors():
    return [
        case("CO2", "3e6", "340", "0", "0.7", "0.006", "0.989"),
        case("CH4", "1e7", "350", "1", "0.5", "0.0009", "0.995"),
    ]


if __name__ == "__main__":
    with localcontext() as context:
        context.prec = 80
        computed = anchors()
        for values in computed:
            if not values[8] > 0:
                raise AssertionError("assigned aqueous phase is not water-richer")

        cpp = Path(__file__).with_name("sw92_phase_assigned_joint_test.cpp").read_text(
            encoding="utf-8")
        block = cpp.split("constexpr JointGolden golden[] = {", 1)[1].split("};", 1)[0]
        literals = [D(value) for value in re.findall(
            r"([-+]?[0-9]+\.[0-9]*(?:[eE][-+]?[0-9]+)?)L", block)]
        expected_count = 9 * len(computed)
        if len(literals) != expected_count:
            raise AssertionError(
                f"phase-assigned golden count changed: {len(literals)} != {expected_count}")
        cursor = 0
        for index, values in enumerate(computed):
            for value in values:
                literal = literals[cursor]
                cursor += 1
                if abs(literal - value) > D("2e-31") * max(D(1), abs(value)):
                    raise AssertionError(
                        f"phase-assigned anchor mismatch: {literal} vs {value}")
            print("SW92 phase-assigned AQ/NA case", index + 1,
                  "x_AQ=", format(values[0], ".35g"),
                  "y_NA=", format(values[1], ".35g"),
                  "beta_NA=", format(values[2], ".35g"),
                  "G/RT=", format(values[5], ".35g"))

        print("Independent Decimal(80) SW92 Profile-C Gate C1 references passed")
