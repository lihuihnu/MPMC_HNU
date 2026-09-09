"""Regenerate printed synthetic PT anchors; no production imports or source edits.

Python standard-library Decimal, 80 digits. Independent original-Z Newton
iteration and direct Eq.(19), not the production shifted polynomial/bisection.
The input constants (apart from mathematical decimals such as R) are exactly
binary-representable, so no binary64-to-decimal fixture ambiguity is present.
"""
from decimal import Decimal as D, localcontext
from pathlib import Path
import re


def anchors() -> list[list[D]]:
    with localcontext() as ctx:
        ctx.prec = 80
        pressure, temperature, r = D(1000000), D(450), D("8.31446261815324")
        specs = [(D(400), D(4000000), D("-.125")),
                 (D(500), D(3000000), D(".25")),
                 (D(600), D(5000000), D(".5"))]
        w = [D(".25"), D(".5"), D(".25")]
        kij = [[D(0), D(".125"), D("-.0625")],
               [D(".125"), D(0), D(".0625")],
               [D("-.0625"), D(".0625"), D(0)]]
        ai, bi = [], []
        for tc, pc, omega in specs:
            k = D(".37464") + D("1.54226") * omega - D(".26992") * omega**2
            ai.append(D(".45724") * r**2 * tc**2 / pc *
                      (1 + k * (1 - (temperature / tc).sqrt()))**2)
            bi.append(D(".07780") * r * tc / pc)
        pair = [[(1 - kij[i][j]) * ai[i].sqrt() * ai[j].sqrt()
                 for j in range(3)] for i in range(3)]
        rows = [sum(w[j] * pair[i][j] for j in range(3)) for i in range(3)]
        a = sum(w[i] * rows[i] for i in range(3))
        b = sum(w[i] * bi[i] for i in range(3))
        A, B = a * pressure / (r * temperature)**2, b * pressure / (r * temperature)
        c2, c1, c0 = B - 1, A - 3 * B**2 - 2 * B, B**3 + B**2 - A * B
        answer = []
        for guess in (".04", ".09", ".85"):
            z = D(guess)
            for _ in range(100):
                f = z**3 + c2 * z**2 + c1 * z + c0
                derivative = 3 * z**2 + 2 * c2 * z + c1
                next_z = z - f / derivative
                if abs(next_z - z) < D("1e-70"):
                    z = next_z
                    break
                z = next_z
            else:
                raise RuntimeError("Independent reference Newton iteration did not converge")
            if not z > B or abs(z**3 + c2 * z**2 + c1 * z + c0) > D("1e-68"):
                raise RuntimeError("Independent reference residual/domain failed")
            s2 = D(2).sqrt()
            log_ratio = ((z + (1 + s2) * B) / (z + (1 - s2) * B)).ln()
            phi = [bi[i] / b * (z - 1) - (z - B).ln() -
                   A / (2 * s2 * B) * (2 * rows[i] / a - bi[i] / b) * log_ratio
                   for i in range(3)]
            answer.append([z, *phi])
        if not answer[0][0] < answer[1][0] < answer[2][0]:
            raise RuntimeError("Reference roots are not distinct/increasing")
        return answer


if __name__ == "__main__":
    values = anchors()
    text = Path(__file__).with_name("pt_test.cpp").read_text(encoding="utf-8")
    block = text.split("constexpr long double golden[3][4] = {", 1)[1].split("};", 1)[0]
    literals = re.findall(r"([-+]?(?:[0-9]+\.[0-9]*|\.[0-9]+))L", block)
    expected = [number for row in values for number in row]
    if len(literals) != len(expected):
        raise RuntimeError("Golden anchor count changed; audit its definition")
    for literal, computed in zip(literals, expected):
        if abs(D(literal) - computed) > D("1e-32") * abs(computed):
            raise RuntimeError(f"Independent Decimal reference mismatch: {literal} vs {computed}")
    for row in values:
        print(", ".join(format(number, ".36g") for number in row))
    print("Independent Decimal(80) anchors: 12/12 checked")
