"""Independent Decimal(80) structural reference for Profile-C Gate C2b.2.

Stdlib only; no production imports. This builds on the already-independent C2b.1
SW92 transcription and re-evaluates the nested C1 and W+H0+H1 states. The
CH4/CO2 non-water pair remains the explicitly tagged synthetic kij=0 fixture,
so these checks are structural numerical evidence only, not experimental or
physical validation.
"""

from decimal import Decimal as D, localcontext
from pathlib import Path
import runpy

BASE = (Path(__file__).parents[1] /
        "sw92_phase_assigned_three_phase" / "reference_decimal.py")
ref = runpy.run_path(str(BASE), run_name="sw92_c2b2_base_reference")
anchors = ref["anchors"]
mu = ref["mu"]
tpd = ref["tpd"]


if __name__ == "__main__":
    with localcontext() as context:
        context.prec = 80
        values = anchors()

        c1_mu_w, _ = mu(values["c1_w"], "AQ")
        c1_mu_h, _ = mu(values["c1_h"], "NA")
        common_c1 = [(c1_mu_w[i] + c1_mu_h[i]) / D(2) for i in range(3)]
        common_three = values["common"]

        g_c1 = sum(values["feed"][i] * common_c1[i] for i in range(3))
        g_three = sum(values["feed"][i] * common_three[i] for i in range(3))
        if not g_three < g_c1:
            raise AssertionError("C2b.1 nested three-phase Gibbs is not below C1 source")

        h0_tpd = tpd(values["h0"], common_three, "NA")
        h1_tpd = tpd(values["h1"], common_three, "NA")
        old_witness_tpd = tpd(values["h1"], common_three, "NA")
        for label, value in (("H0", h0_tpd), ("H1", h1_tpd),
                             ("old C2a1 witness", old_witness_tpd)):
            if abs(value) > D("3e-47"):
                raise AssertionError(f"{label} is not on the final C2b.1 tangent: {value}")

        na_at_w = tpd(values["w"], common_three, "NA")
        if values["w"][2] <= max(values["h0"][2], values["h1"][2]):
            raise AssertionError("structural W phase lost the relative-water topology guard")

        print("C1 reduced G/RT=", format(g_c1, ".35g"))
        print("C2b.1 reduced G/RT=", format(g_three, ".35g"))
        print("nested delta G/RT=", format(g_three - g_c1, ".35g"))
        print("final TPD(H0)=", format(h0_tpd, ".35g"))
        print("final TPD(H1/source witness)=", format(h1_tpd, ".35g"))
        print("diagnostic final NA TPD(W composition)=", format(na_at_w, ".35g"))
        print("Independent Decimal(80) Profile-C C2b.2 structural reference passed")
