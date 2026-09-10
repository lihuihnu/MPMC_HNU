"""Independent Decimal(80) anchors for SW92 Xu-style Gate 3A stability.

Stdlib only; no production imports. The script independently rebuilds the
corrected-original CO2/H2O SW92 equations, PR mixing/cubic roots, same-family
minimum-Gibbs root selection, cross-family feed Gibbs comparison, one common
feed tangent and prescribed AQ/NA TPD values.

Pure-vertex verification is deliberately for the ACTIVE pure component only.
The zero-composition component has an infinite-dilution fugacity coefficient
that may depend on the AQ/NA cross BIP and is not required to be equal.

These are model numerical anchors, not experimental validation and not a global
phase-stability proof.
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
NA_KIJ_REAL = D("0.1896")


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
    co2_a = D("0.45724") * R * R * TC_CO2 * TC_CO2 / PC_CO2 * co2_alpha(temperature)
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
        for _ in range(160):
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


def family_candidates(gas_fraction: D, family: str, pressure: D,
                      temperature: D, molality: D, na_kij=NA_KIJ_REAL):
    x = [gas_fraction, D(1) - gas_fraction]
    pure = pure_ab(temperature, molality)
    ai = [item[0] for item in pure]
    bi = [item[1] for item in pure]
    kij = aqueous_kij(temperature, molality) if family == "AQ" else na_kij

    aij = [[D(0), D(0)], [D(0), D(0)]]
    for i in range(2):
        for j in range(2):
            pair_kij = D(0) if i == j else kij
            aij[i][j] = (ai[i] * ai[j]).sqrt() * (D(1) - pair_kij)
    rows = [sum(x[j] * aij[i][j] for j in range(2)) for i in range(2)]
    attraction = sum(x[i] * rows[i] for i in range(2))
    covolume = sum(x[i] * bi[i] for i in range(2))
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
        root_gibbs = sum(x[i] * ln_phi[i] for i in range(2))
        candidates.append((z, ln_phi, root_gibbs))
    if not candidates:
        raise ArithmeticError("no mechanically admissible reference root")
    return candidates


def minimum_phase(gas_fraction: D, family: str, pressure: D,
                  temperature: D, molality: D, na_kij=NA_KIJ_REAL):
    return min(family_candidates(gas_fraction, family, pressure, temperature,
                                 molality, na_kij), key=lambda item: item[2])


def reduced_feed_gibbs(feed: D, phase) -> D:
    z = [feed, D(1) - feed]
    return sum(z[i] * (z[i].ln() + phase[1][i]) for i in range(2) if z[i] != 0)


def tpd(trial: D, feed: D, trial_phase, common_d) -> D:
    w = [trial, D(1) - trial]
    return sum(w[i] * (w[i].ln() + trial_phase[1][i] - common_d[i])
               for i in range(2) if w[i] != 0)


def feed_anchors():
    pressure = D("3e6")
    temperature = D("340")
    molality = D("0")
    feed = D("0.7")
    aq = minimum_phase(feed, "AQ", pressure, temperature, molality)
    na = minimum_phase(feed, "NA", pressure, temperature, molality)
    aq_g = reduced_feed_gibbs(feed, aq)
    na_g = reduced_feed_gibbs(feed, na)
    delta = aq_g - na_g
    z = [feed, D(1) - feed]
    common_d = [z[i].ln() + aq[1][i] for i in range(2)]
    uniform = D("0.5")
    aq_uniform = minimum_phase(uniform, "AQ", pressure, temperature, molality)
    na_uniform = minimum_phase(uniform, "NA", pressure, temperature, molality)
    return (aq_g, na_g, delta,
            aq[1][0], aq[1][1], na[1][0], na[1][1],
            common_d[0], common_d[1],
            tpd(uniform, feed, aq_uniform, common_d),
            tpd(uniform, feed, na_uniform, common_d))


def pure_active_anchors():
    pressure = D("3e6")
    temperature = D("340")
    molality = D(0)
    co2_aq = minimum_phase(D(1), "AQ", pressure, temperature, molality)
    co2_na = minimum_phase(D(1), "NA", pressure, temperature, molality)
    water_aq = minimum_phase(D(0), "AQ", pressure, temperature, molality)
    water_na = minimum_phase(D(0), "NA", pressure, temperature, molality)
    if co2_aq[0] != co2_na[0] or co2_aq[1][0] != co2_na[1][0]:
        raise AssertionError("active pure CO2 gauge mismatch")
    if water_aq[0] != water_na[0] or water_aq[1][1] != water_na[1][1]:
        raise AssertionError("active pure water gauge mismatch")
    if co2_aq[1][1] == co2_na[1][1] or water_aq[1][0] == water_na[1][0]:
        raise AssertionError("zero-component infinite-dilution check lost sensitivity")
    return co2_aq[1][0], water_aq[1][1]


def saturation_pressure() -> D:
    temperature = D("280")
    molality = D(0)
    left = D("4.1e6")
    right = D("4.2e6")

    def gap(pressure: D) -> D:
        candidates = family_candidates(D(1), "AQ", pressure, temperature, molality)
        if len(candidates) < 2:
            raise ArithmeticError("pure-CO2 reference lost two mechanical roots")
        return candidates[0][1][0] - candidates[-1][1][0]

    f_left = gap(left)
    f_right = gap(right)
    if not (f_left > 0 and f_right < 0):
        raise AssertionError("saturation bracket changed")
    for _ in range(220):
        middle = (left + right) / D(2)
        f_middle = gap(middle)
        if abs(f_middle) < D("1e-62"):
            return middle
        if f_middle > 0:
            left = middle
        else:
            right = middle
    return (left + right) / D(2)


def structural_tie_checks(k_literal: D):
    pressure = D("3e6")
    temperature = D("340")
    molality = D(0)
    feed = D("0.7")
    exact_k = aqueous_kij(temperature, molality)
    aq = minimum_phase(feed, "AQ", pressure, temperature, molality)
    na_exact = minimum_phase(feed, "NA", pressure, temperature, molality, exact_k)
    exact_delta = reduced_feed_gibbs(feed, aq) - reduced_feed_gibbs(feed, na_exact)
    if abs(exact_delta) > D("1e-60"):
        raise AssertionError("exact synthetic family tie is not exact")

    one_ulp_float = math.nextafter(float(k_literal), math.inf)
    na_near = minimum_phase(feed, "NA", pressure, temperature, molality,
                            D(str(one_ulp_float)))
    near_delta = reduced_feed_gibbs(feed, aq) - reduced_feed_gibbs(feed, na_near)
    if abs(near_delta) > D("1e-13"):
        raise AssertionError("one-ULP synthetic family tie is not numerically near")
    return exact_k, near_delta


if __name__ == "__main__":
    with localcontext() as context:
        context.prec = 80
        cpp = Path(__file__).with_name("sw92_asymmetric_stability_test.cpp").read_text(
            encoding="utf-8")

        block = cpp.split("constexpr AsymmetricGolden golden{", 1)[1].split("};", 1)[0]
        literals = [D(value) for value in re.findall(
            r"([-+]?[0-9]+\.[0-9]*(?:[eE][-+]?[0-9]+)?)L", block)]
        computed = feed_anchors()
        if len(literals) != len(computed):
            raise AssertionError(f"asymmetric golden count changed: {len(literals)} != {len(computed)}")
        for literal, value in zip(literals, computed):
            if abs(literal - value) > D("1e-31") * max(D(1), abs(value)):
                raise AssertionError(f"asymmetric anchor mismatch: {literal} vs {value}")

        pure_co2_literal = D(re.search(
            r"pure_co2_active_lnphi\s*=\s*\n?\s*([-+]?[0-9.]+)L", cpp).group(1))
        pure_water_literal = D(re.search(
            r"pure_water_active_lnphi\s*=\s*\n?\s*([-+]?[0-9.]+)L", cpp).group(1))
        pure_co2, pure_water = pure_active_anchors()
        if abs(pure_co2_literal - pure_co2) > D("1e-31"):
            raise AssertionError("pure CO2 active-component anchor mismatch")
        if abs(pure_water_literal - pure_water) > D("1e-31"):
            raise AssertionError("pure water active-component anchor mismatch")

        saturation_literal = D(re.search(
            r"co2_pr_saturation_pressure_pa\s*=\s*([0-9.]+);", cpp).group(1))
        saturation = saturation_pressure()
        if abs(saturation_literal - saturation) > D("1e-8"):
            raise AssertionError(f"pure CO2 saturation anchor mismatch: {saturation_literal} vs {saturation}")
        literal_gap_candidates = family_candidates(
            D(1), "AQ", saturation_literal, D("280"), D(0))
        literal_gap = literal_gap_candidates[0][1][0] - literal_gap_candidates[-1][1][0]
        if abs(literal_gap) > D("1e-13"):
            raise AssertionError("double saturation anchor is not a same-family Gibbs near-tie")

        kij_literal = D(re.search(
            r"aq_co2_kij_340k_fresh\s*=\s*([-+]?[0-9.]+);", cpp).group(1))
        exact_k, near_delta = structural_tie_checks(kij_literal)
        if abs(kij_literal - exact_k) > D("1e-15"):
            raise AssertionError(f"synthetic AQ/NA equality BIP literal drifted: {kij_literal} vs {exact_k}")

        print("SW92 asymmetric feed delta=", format(computed[2], ".36g"))
        print("common-reference uniform TPD AQ=", format(computed[9], ".36g"),
              "NA=", format(computed[10], ".36g"))
        print("pure active lnphi CO2=", format(pure_co2, ".36g"),
              "H2O=", format(pure_water, ".36g"))
        print("pure CO2 PR saturation pressure=", format(saturation, ".36g"),
              "literal gap=", format(literal_gap, ".12g"))
        print("synthetic one-ULP family delta=", format(near_delta, ".12g"))
        print("Independent Decimal(80) SW92 asymmetric stability anchors: PASS")
