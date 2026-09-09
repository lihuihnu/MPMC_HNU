"""Regenerate printed synthetic PT anchors; no production imports or source edits.

Python standard-library Decimal, 80 digits. Independent original-Z Newton
iteration and direct Eq.(19), not the production shifted polynomial/bisection.
The input constants (apart from mathematical decimals such as R) are exactly
binary-representable, so no binary64-to-decimal fixture ambiguity is present.
"""
from decimal import Decimal as D, localcontext
from pathlib import Path
import re


def anchors() -> tuple[list[list[D]], list[list[list[D]]]]:
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
        ai, bi, slopes = [], [], []
        for tc, pc, omega in specs:
            k = D(".37464") + D("1.54226") * omega - D(".26992") * omega**2
            ai.append(D(".45724") * r**2 * tc**2 / pc *
                      (1 + k * (1 - (temperature / tc).sqrt()))**2)
            bi.append(D(".07780") * r * tc / pc)
            ac = D(".45724") * r**2 * tc**2 / pc
            q = 1 + k * (1 - (temperature / tc).sqrt())
            slopes.append(-ac * k * q / (temperature * tc).sqrt())
        pair = [[(1 - kij[i][j]) * ai[i].sqrt() * ai[j].sqrt()
                 for j in range(3)] for i in range(3)]
        rows = [sum(w[j] * pair[i][j] for j in range(3)) for i in range(3)]
        a = sum(w[i] * rows[i] for i in range(3))
        b = sum(w[i] * bi[i] for i in range(3))
        A, B = a * pressure / (r * temperature)**2, b * pressure / (r * temperature)
        c2, c1, c0 = B - 1, A - 3 * B**2 - 2 * B, B**3 + B**2 - A * B
        row_slopes = [sum(w[j] * pair[i][j] *
                          (slopes[i] / ai[i] + slopes[j] / ai[j]) / 2
                          for j in range(3)) for i in range(3)]
        da_dt = sum(w[i] * row_slopes[i] for i in range(3))
        answer, reduced = [], []
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
            derivatives = [[D(0)] * 4 for _ in range(4)]
            g = A / (2 * s2 * B)
            plus, minus = z + (1 + s2) * B, z + (1 - s2) * B
            fz = 3 * z**2 + 2 * (B - 1) * z + A - 3 * B**2 - 2 * B
            fa, fb = z - B, z**2 - (6 * B + 2) * z + 3 * B**2 + 2 * B - A
            for column in range(4):
                # Direct tangent coordinates, not subtraction of rounded full gradients.
                da = da_dt if column == 1 else (2 * (rows[column - 2] - rows[2])
                                                if column >= 2 else D(0))
                db = bi[column - 2] - bi[2] if column >= 2 else D(0)
                dA = (A / pressure if column == 0 else
                      A * (da_dt / a - 2 / temperature) if column == 1 else A * da / a)
                dB = (B / pressure if column == 0 else
                      -B / temperature if column == 1 else B * db / b)
                dz = -(fa * dA + fb * dB) / fz
                derivatives[0][column] = dz
                dg = (dA * B - A * dB) / (2 * s2 * B**2)
                dlog = (dz + (1 + s2) * dB) / plus - (dz + (1 - s2) * dB) / minus
                for i in range(3):
                    ds = (row_slopes[i] if column == 1 else
                          pair[i][column - 2] - pair[i][2] if column >= 2 else D(0))
                    ratio, dratio = bi[i] / b, -bi[i] * db / b**2
                    e = 2 * rows[i] / a - ratio
                    de = 2 * (ds * a - rows[i] * da) / a**2 - dratio
                    derivatives[i + 1][column] = (dratio * (z - 1) + ratio * dz -
                        (dz - dB) / (z - B) - (dg * e + g * de) * log_ratio - g * e * dlog)
            reduced.append(derivatives)
        if not answer[0][0] < answer[1][0] < answer[2][0]:
            raise RuntimeError("Reference roots are not distinct/increasing")
        return answer, reduced


if __name__ == "__main__":
    values, derivatives = anchors()
    text = Path(__file__).with_name("pt_test.cpp").read_text(encoding="utf-8")
    groups = [("constexpr long double golden[3][4] = {", [x for row in values for x in row]),
              ("constexpr long double reduced_golden[3][4][4] = {",
               [x for root in derivatives for row in root for x in row])]
    count = 0
    for marker, expected in groups:
        block = text.split(marker, 1)[1].split("};", 1)[0]
        literals = re.findall(r"([-+]?(?:[0-9]+\.[0-9]*|\.[0-9]+)(?:[eE][-+]?[0-9]+)?)L", block)
        if len(literals) != len(expected):
            raise RuntimeError("Golden anchor count changed; audit its definition")
        for literal, computed in zip(literals, expected):
            if abs(D(literal) - computed) > D("1e-32") * abs(computed):
                raise RuntimeError(f"Independent Decimal reference mismatch: {literal} vs {computed}")
        count += len(literals)
    for row in values:
        print(", ".join(format(number, ".36g") for number in row))
    print("Reduced derivative diagnostic (root 1, component 1, w0):", derivatives[1][2][2])
    print(f"Independent Decimal(80) anchors: {count}/{count} checked (12 values + 48 derivatives)")
