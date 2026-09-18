"""Independent Decimal(80) SW92/corrected-original reference regeneration.

No production imports. Implements the paper's Eq. (9), corrected Eqs. (12)-(15),
Eq. (17), classical PR mixing Eq. (10), the original-Z PR cubic, and direct
Peng-Robinson fugacity-coefficient expression. Constants are transcribed from
the user-supplied SW92 paper and its appended authors' errata.
"""
from decimal import Decimal as D, localcontext
from pathlib import Path
import re

R = D("8.31446261815324")
S2 = D(2).sqrt()

SPECS = {
    "methane": (D("190.6"), D("46.0e5"), D("0.0108"), D("0.4850")),
    "nitrogen": (D("126.1"), D("34.0e5"), D("0.0403"), D("0.4778")),
    "co2": (D("304.2"), D("73.8e5"), D("0.2273"), D("0.1896")),
    "h2s": (D("373.2"), D("89.4e5"), D("0.1081"), None),
    "water": (D("647.3"), D("221.2e5"), D("0.3434"), None),
}


def powd(x: D, exponent: D) -> D:
    if x == 0:
        if exponent > 0:
            return D(0)
        raise ValueError("zero base outside reference domain")
    return (exponent * x.ln()).exp()


def water_alpha(t: D, c: D) -> D:
    tr = t / SPECS["water"][0]
    root = (D(1) + D("0.4530") *
            (D(1) - tr * (D(1) - D("0.0103") * powd(c, D("1.1")))) +
            D("0.0034") * (tr ** D(-3) - D(1)))
    return root * root


def aq_kij(name: str, t: D, c: D) -> D:
    tc, _, omega, _ = SPECS[name]
    tr = t / tc
    if name == "methane":
        a0 = D("1.1120") - D("1.7369") * powd(omega, D("-0.1"))
        a1 = D("1.1001") + D("0.8360") * omega
        a2 = -D("0.15742") - D("1.0988") * omega
        return (a0 * (D(1) + D("0.017407") * c) +
                a1 * tr * (D(1) + D("0.033516") * c) +
                a2 * tr * tr * (D(1) + D("0.011478") * c))
    if name == "nitrogen":
        c75 = powd(c, D("0.75"))
        return (-D("1.70235") * (D(1) + D("0.025587") * c75) +
                D("0.44338") * (D(1) + D("0.08126") * c75) * tr)
    if name == "co2":
        return (-D("0.31092") * (D(1) + D("0.15587") * powd(c, D("0.7505"))) +
                D("0.23580") * (D(1) + D("0.17837") * powd(c, D("0.979"))) * tr -
                D("21.2566") * (-D("6.7222") * tr - c).exp())
    if name == "h2s":
        return -D("0.20441") + D("0.23426") * tr
    raise ValueError(name)


def na_kij(name: str, t: D) -> D:
    tc, _, _, constant = SPECS[name]
    if name == "h2s":
        return D("0.19031") - D("0.05965") * (t / tc)
    assert constant is not None
    return constant


def pure_ab(name: str, t: D, c: D) -> tuple[D, D]:
    tc, pc, omega, _ = SPECS[name]
    ac = D("0.45724") * R * R * tc * tc / pc
    b = D("0.07780") * R * tc / pc
    if name == "water":
        alpha = water_alpha(t, c)
    else:
        kappa = D("0.37464") + D("1.54226") * omega - D("0.26992") * omega * omega
        q = D(1) + kappa * (D(1) - (t / tc).sqrt())
        alpha = q * q
    return ac * alpha, b


def mixture(names: list[str], x: list[D], t: D, c: D, family: str,
            background: dict[tuple[int, int], D] | None = None):
    ai, bi = zip(*(pure_ab(name, t, c) for name in names))
    rows = [D(0) for _ in names]
    a = D(0)
    for i in range(len(names)):
        for j in range(len(names)):
            if i == j:
                kij = D(0)
            elif names[i] == "water" or names[j] == "water":
                other = names[j] if names[i] == "water" else names[i]
                kij = aq_kij(other, t, c) if family == "aq" else na_kij(other, t)
            else:
                if background is None:
                    raise ValueError("explicit non-water pair required")
                key = tuple(sorted((i, j)))
                kij = background[key]
            aij = ai[i].sqrt() * ai[j].sqrt() * (D(1) - kij)
            rows[i] += x[j] * aij
            a += x[i] * x[j] * aij
    b = sum((x[i] * bi[i] for i in range(len(names))), D(0))
    return a, b, list(ai), list(bi), rows


def roots_and_phi(names: list[str], x: list[D], t: D, c: D, p: D, family: str):
    a, b, _, bi, rows = mixture(names, x, t, c, family)
    A = a * p / (R * t) ** 2
    B = b * p / (R * t)
    c2 = B - D(1)
    c1 = A - D(3) * B * B - D(2) * B
    c0 = B ** 3 + B * B - A * B
    found = []
    for guess in (D(".03"), D(".08"), D(".9")):
        z = guess
        for _ in range(200):
            f = z**3 + c2*z*z + c1*z + c0
            fp = D(3)*z*z + D(2)*c2*z + c1
            nz = z - f / fp
            if abs(nz-z) < D("1e-70"):
                z = nz
                break
            z = nz
        else:
            raise RuntimeError("reference Newton did not converge")
        if z > B and all(abs(z-r[0]) > D("1e-50") for r in found):
            ratio_log = ((z + (D(1)+S2)*B) / (z + (D(1)-S2)*B)).ln()
            phi = []
            for i in range(len(names)):
                ratio = bi[i] / b
                e = D(2)*rows[i]/a - ratio
                phi.append(ratio*(z-D(1)) - (z-B).ln() -
                           A/(D(2)*S2*B)*e*ratio_log)
            found.append((z, *phi))
    found.sort(key=lambda row: row[0])
    return a, b, A, B, found


def anchors():
    with localcontext() as ctx:
        ctx.prec = 80
        alpha = [water_alpha(D("377.15"), c) for c in (D(0), D(1), D(5))]
        cases = [("methane", D("377.15")), ("nitrogen", D("376.15")),
                 ("co2", D("423.15")), ("h2s", D("377.15"))]
        kij = [[aq_kij(name, t, D(1)), na_kij(name, t)] for name, t in cases]
        aqmix = mixture(["methane", "water"], [D(".003"), D(".997")],
                        D("377.15"), D(1), "aq")[:2]
        namix = mixture(["methane", "water"], [D(".003"), D(".997")],
                        D("377.15"), D(1), "na")[:2]
        na3 = roots_and_phi(["co2", "water"], [D(".7"), D(".3")],
                            D(340), D(0), D("3e6"), "na")
        aq1 = roots_and_phi(["co2", "water"], [D(".02"), D(".98")],
                            D("423.15"), D(1), D("20e6"), "aq")
        if len(na3[4]) != 3 or len(aq1[4]) != 1:
            raise RuntimeError("unexpected independent reference root topology")
        return alpha, kij, [aqmix, namix], na3[4], aq1[4][0]


def literals_after(text: str, marker: str) -> list[D]:
    block = text.split(marker, 1)[1].split("};", 1)[0]
    return [D(x) for x in re.findall(
        r"([-+]?(?:[0-9]+\.[0-9]*|\.[0-9]+)(?:[eE][-+]?[0-9]+)?)L", block)]


if __name__ == "__main__":
    alpha, kij, mix, na3, aq1 = anchors()
    expected = {
        "inline constexpr long double alpha_golden[]={": alpha,
        "inline constexpr long double kij_golden[4][2]={": [v for row in kij for v in row],
        "inline constexpr long double methane_mix_golden[2][2]={": [v for row in mix for v in row],
        "inline constexpr long double co2_na_phase_golden[3][3]={": [v for row in na3 for v in row],
        "inline constexpr long double co2_aq_phase_golden[3]={": list(aq1),
    }
    text = (Path(__file__).parents[2] / "support" / "sw92" / "test_support.hpp").read_text(encoding="utf-8")
    count = 0
    for marker, values in expected.items():
        literals = literals_after(text, marker)
        if len(literals) != len(values):
            raise RuntimeError(f"anchor count mismatch for {marker}: {len(literals)} != {len(values)}")
        for literal, value in zip(literals, values):
            scale = max(D(1), abs(value))
            if abs(literal - value) > D("2e-47") * scale:
                raise RuntimeError(f"Decimal reference mismatch: {literal} vs {value}")
        count += len(values)
    print("SW92 water-alpha:", *(format(v, ".30g") for v in alpha))
    for (name, _), row in zip((("methane",0),("nitrogen",0),("co2",0),("h2s",0)), kij):
        print(name, "kAQ/kNA:", *(format(v, ".30g") for v in row))
    print("CO2/water NA roots:", *(format(row[0], ".30g") for row in na3))
    print(f"Independent Decimal(80) SW92 anchors: {count}/{count} checked")
