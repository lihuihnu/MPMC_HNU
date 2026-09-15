"""Independent Decimal PR76 three-phase reference for Li/Firoozabadi acid-gas benchmark.

Scientific source:
- Li & Firoozabadi, SPE Journal 17(4), 1096-1107 (2012),
  DOI 10.2118/129844-PA.
- Table 3 supplies Tc, Pc, omega and nonzero kij for acid gas.
- Table 9 selects T=178.8 K, P=20 bar, overall CO2 mole fraction=0.5 and
  explicitly performs a three-phase split.
- The full six-component feed used here is a transparent derived P-Z state:
  z_CO2=0.5 and the other five components retain their Table-3 relative ratios.

This file independently implements the repository PR76 convention:
printed 1976 coefficients, R=8.31446261815324 SI, original quadratic kappa,
classical symmetric van-der-Waals mixing, and exact sqrt(2).
It imports no production C++ and does not fit production output.
"""

from __future__ import annotations

import argparse
from decimal import Decimal as D, localcontext
from pathlib import Path

COMPONENTS = ("carbon-dioxide", "nitrogen", "hydrogen-sulfide", "methane", "ethane", "propane")
TC = tuple(map(D, ("304.211", "126.2", "373.2", "190.564", "305.322", "369.825")))
PC_PA = tuple(D(v) * D(100000) for v in ("73.819", "33.9", "89.4", "45.992", "48.718", "42.462"))
OMEGA = tuple(map(D, (".225", ".039", ".081", ".01141", ".10574", ".15813")))
TABLE3_INITIAL = tuple(map(D, (".70592", ".07026", ".01966", ".06860", ".10559", ".02967")))
PRESSURE_PA = D("2000000")
TEMPERATURE_K = D("178.8")
R = D("8.31446261815324")

# Table 3 lists nonzero kij. Every unlisted pair is explicitly zero here.
NONZERO_KIJ = {
    (0, 1): D("-0.02"),
    (0, 2): D("0.12"),
    (1, 2): D("0.2"),
    (0, 3): D("0.125"),
    (1, 3): D("0.031"),
    (2, 3): D("0.1"),
    (0, 4): D("0.135"),
    (1, 4): D("0.042"),
    (2, 4): D("0.08"),
    (0, 5): D("0.150"),
    (1, 5): D("0.091"),
    (2, 5): D("0.08"),
}

# Rounded, deliberately low-precision numerical initialization only. These values
# are never emitted as reference values and do not bypass the nonlinear solve.
PHASE_GUESSES = (
    ("0.045", "0.650", "0.002", "0.279", "0.024", "0.0008"),
    ("0.383", "0.057", "0.040", "0.132", "0.297", "0.091"),
    ("0.782", "0.035", "0.035", "0.049", "0.084", "0.015"),
)
BETA_GUESS = ("0.12", "0.484", "0.396")
Z_GUESS = ("0.868", "0.058", "0.047")


def feed():
    remainder = sum(TABLE3_INITIAL[1:])
    values = [D("0.5")]
    values.extend(D("0.5") * value / remainder for value in TABLE3_INITIAL[1:])
    if sum(values) != D(1):
        raise ArithmeticError("derived feed does not sum exactly to one")
    return values


def normalized(values):
    total = sum(values)
    return [v / total for v in values]


def softmax5(log_ratios):
    values = [v.exp() for v in log_ratios] + [D(1)]
    return normalized(values)


def phase_fraction_map(log_ratios):
    values = [D(1), log_ratios[0].exp(), log_ratios[1].exp()]
    return normalized(values)


def log_coordinates(values):
    return [(values[i] / values[-1]).ln() for i in range(len(values) - 1)]


def beta_coordinates(values):
    return [(values[1] / values[0]).ln(), (values[2] / values[0]).ln()]


def parameter_state():
    sqrt2 = D(2).sqrt()
    n = len(COMPONENTS)
    kij = [[D(0) for _ in range(n)] for _ in range(n)]
    for (i, j), value in NONZERO_KIJ.items():
        kij[i][j] = value
        kij[j][i] = value

    pure_a = []
    pure_b = []
    for tc, pc, omega in zip(TC, PC_PA, OMEGA):
        kappa = D("0.37464") + omega * (D("1.54226") - D("0.26992") * omega)
        ac = D("0.45724") * (R * tc) * (R * tc / pc)
        b = D("0.07780") * R * tc / pc
        alpha = (D(1) + kappa * (D(1) - (TEMPERATURE_K / tc).sqrt())) ** 2
        pure_a.append(ac * alpha)
        pure_b.append(b)

    cross = [
        [(pure_a[i] * pure_a[j]).sqrt() * (D(1) - kij[i][j]) for j in range(n)]
        for i in range(n)
    ]
    return sqrt2, pure_a, pure_b, cross


def mixture(composition, pure_b, cross):
    n = len(composition)
    attraction = D(0)
    for i in range(n):
        for j in range(n):
            attraction += composition[i] * composition[j] * cross[i][j]
    covolume = sum(composition[i] * pure_b[i] for i in range(n))
    partial = [
        sum(composition[j] * cross[i][j] for j in range(n))
        for i in range(n)
    ]
    A = attraction * PRESSURE_PA / (R * TEMPERATURE_K) ** 2
    B = covolume * PRESSURE_PA / (R * TEMPERATURE_K)
    return attraction, covolume, partial, A, B


def cubic_residual(zfactor, A, B):
    return (
        zfactor ** 3
        - (D(1) - B) * zfactor ** 2
        + (A - D(3) * B * B - D(2) * B) * zfactor
        - (A * B - B * B - B ** 3)
    )


def phase_properties(composition, zfactor, params):
    sqrt2, _, pure_b, cross = params
    attraction, covolume, partial, A, B = mixture(composition, pure_b, cross)
    if zfactor <= B:
        raise ArithmeticError("PR76 candidate has Z <= B")
    logarithm = (
        (zfactor + (D(1) + sqrt2) * B)
        / (zfactor + (D(1) - sqrt2) * B)
    ).ln()
    ln_phi = []
    for i in range(len(composition)):
        b_ratio = pure_b[i] / covolume
        mixing_term = D(2) * partial[i] / attraction - b_ratio
        value = (
            b_ratio * (zfactor - D(1))
            - (zfactor - B).ln()
            - A / (D(2) * sqrt2 * B) * mixing_term * logarithm
        )
        ln_phi.append(value)
    return ln_phi, cubic_residual(zfactor, A, B)


def unpack(unknowns):
    phases = [
        softmax5(unknowns[0:5]),
        softmax5(unknowns[5:10]),
        softmax5(unknowns[10:15]),
    ]
    fractions = phase_fraction_map(unknowns[15:17])
    zfactors = [unknowns[17 + i].exp() for i in range(3)]
    return phases, fractions, zfactors


def residual(unknowns, params, overall):
    phases, fractions, zfactors = unpack(unknowns)
    ln_phi = []
    cubics = []
    for phase, zfactor in zip(phases, zfactors):
        values, cubic = phase_properties(phase, zfactor, params)
        ln_phi.append(values)
        cubics.append(cubic)

    chemical = [
        [phases[p][i].ln() + ln_phi[p][i] for i in range(len(COMPONENTS))]
        for p in range(3)
    ]

    values = []
    for p in (1, 2):
        values.extend(chemical[p][i] - chemical[0][i] for i in range(len(COMPONENTS)))
    for i in range(len(COMPONENTS) - 1):
        values.append(
            sum(fractions[p] * phases[p][i] for p in range(3)) - overall[i]
        )
    values.extend(cubics)
    return values


def solve_linear(matrix, rhs):
    augmented = [row[:] + [value] for row, value in zip(matrix, rhs)]
    n = len(augmented)
    for column in range(n):
        pivot = max(range(column, n), key=lambda row: abs(augmented[row][column]))
        if augmented[pivot][column] == 0:
            raise ArithmeticError("reference Newton Jacobian is singular")
        if pivot != column:
            augmented[column], augmented[pivot] = augmented[pivot], augmented[column]
        denominator = augmented[column][column]
        for row in range(column + 1, n):
            factor = augmented[row][column] / denominator
            if factor == 0:
                continue
            for entry in range(column, n + 1):
                augmented[row][entry] -= factor * augmented[column][entry]

    answer = [D(0)] * n
    for row in range(n - 1, -1, -1):
        value = augmented[row][n] - sum(
            augmented[row][column] * answer[column]
            for column in range(row + 1, n)
        )
        answer[row] = value / augmented[row][row]
    return answer


def initial_unknowns():
    phases = [normalized([D(value) for value in phase]) for phase in PHASE_GUESSES]
    fractions = normalized([D(value) for value in BETA_GUESS])
    return (
        log_coordinates(phases[0])
        + log_coordinates(phases[1])
        + log_coordinates(phases[2])
        + beta_coordinates(fractions)
        + [D(value).ln() for value in Z_GUESS]
    )


def solve_reference(precision):
    with localcontext() as context:
        context.prec = precision
        params = parameter_state()
        overall = feed()
        unknowns = initial_unknowns()
        step = D(10) ** (-(precision // 3))
        target = D(10) ** (-(precision // 2))
        for _ in range(20):
            current = residual(unknowns, params, overall)
            norm = max(abs(value) for value in current)
            if norm < target:
                break
            jacobian = [[D(0)] * len(unknowns) for _ in current]
            for column in range(len(unknowns)):
                upper = unknowns[:]
                lower = unknowns[:]
                upper[column] += step
                lower[column] -= step
                r_upper = residual(upper, params, overall)
                r_lower = residual(lower, params, overall)
                for row in range(len(current)):
                    jacobian[row][column] = (
                        r_upper[row] - r_lower[row]
                    ) / (D(2) * step)
            delta = solve_linear(jacobian, [-value for value in current])
            alpha = D(1)
            accepted = False
            for _ in range(24):
                trial = [
                    unknowns[index] + alpha * delta[index]
                    for index in range(len(unknowns))
                ]
                try:
                    trial_norm = max(
                        abs(value) for value in residual(trial, params, overall)
                    )
                except (ArithmeticError, ValueError):
                    trial_norm = D("Infinity")
                if trial_norm < norm:
                    unknowns = trial
                    accepted = True
                    break
                alpha /= D(2)
            if not accepted:
                raise ArithmeticError("reference Newton line search failed")
        else:
            raise ArithmeticError("reference Newton iteration limit")

        final = residual(unknowns, params, overall)
        norm = max(abs(value) for value in final)
        if norm >= target:
            raise ArithmeticError(f"reference residual {norm} exceeds target {target}")
        phases, fractions, zfactors = unpack(unknowns)
        if min(fractions) <= 0:
            raise ArithmeticError("reference contains a nonpositive phase fraction")
        for phase in phases:
            if min(phase) <= 0 or abs(sum(phase) - D(1)) >= target:
                raise ArithmeticError("reference phase is not a positive normalized composition")
        for i in range(len(COMPONENTS)):
            recovered = sum(fractions[p] * phases[p][i] for p in range(3))
            if abs(recovered - overall[i]) >= target:
                raise ArithmeticError("reference material balance failed")
        return overall, phases, fractions, zfactors, norm


def literal(value, digits=34):
    # Fixed significant digits are intentionally fewer than both 80- and 96-digit solves.
    return format(value, f".{digits}g")


def render(overall, phases, fractions, zfactors):
    text = """// Generated only by reference_sour_gas_decimal.py; no production output.
// Li & Firoozabadi (2012) acid-gas P-Z point, independently reconstructed
// with the repository PR76 printed-coefficient / SI-R / exact-sqrt2 convention.
#ifndef MPMC_TEST_PR76_SOUR_GAS_REFERENCES_HPP
#define MPMC_TEST_PR76_SOUR_GAS_REFERENCES_HPP

#include <array>

namespace pr76_sour_gas_reference {
inline constexpr double pressure_pa = 2000000.0;
inline constexpr double temperature_k = 178.8;
inline constexpr std::array<const char*, 6> component_ids{
    "carbon-dioxide", "nitrogen", "hydrogen-sulfide",
    "methane", "ethane", "propane"};
"""
    text += "inline constexpr std::array<double, 6> feed{\n    "
    text += ",\n    ".join(literal(value) for value in overall)
    text += "};\n"
    text += "inline constexpr std::array<double, 3> phase_fractions{\n    "
    text += ",\n    ".join(literal(value) for value in fractions)
    text += "};\n"
    text += "inline constexpr std::array<double, 3> compressibility_factors{\n    "
    text += ",\n    ".join(literal(value) for value in zfactors)
    text += "};\n"
    text += "inline constexpr std::array<std::array<double, 6>, 3> phases{{\n"
    for phase in phases:
        text += "    {{\n        " + ",\n        ".join(literal(value) for value in phase) + "}},\n"
    text += "}};\n} // namespace pr76_sour_gas_reference\n\n#endif\n"
    return text


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--precision", type=int, choices=(80, 96), default=80)
    parser.add_argument("--write", action="store_true")
    args = parser.parse_args()

    overall, phases, fractions, zfactors, norm = solve_reference(args.precision)
    output = render(overall, phases, fractions, zfactors)
    path = Path(__file__).with_name("sour_gas_references.hpp")
    if args.write:
        path.write_text(output, encoding="utf-8")
    elif output != path.read_text(encoding="utf-8"):
        raise AssertionError("sour-gas reference mismatch; do not fit production output")

    print(f"Decimal({args.precision}) PR76 sour-gas three-phase reference checked")
    print("max coupled equilibrium residual:", norm)
    print("phase fractions:", *(literal(value, 20) for value in fractions))
    print("Z:", *(literal(value, 20) for value in zfactors))


if __name__ == "__main__":
    main()
