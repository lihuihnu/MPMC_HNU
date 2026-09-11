# SW92 experimental three-phase coexistence validation

## Purpose

The physical Sample-6 regression validates a complete multicomponent
`W(AQ)+H0(NA)+H1(NA)` feed flash with unique positive phase fractions. This
second physical regression asks a different question: does the exact project
SW92/PR76 profile reproduce an **experimental three-phase coexistence line**
without fitting any parameter to that line?

The validation uses the binary n-butane/water system because it is fully covered
by the existing corrected-original SW92 parameter contract and does not require
an external hydrocarbon/hydrocarbon BIP.

## Sources

Thermodynamic model data are from Soreide & Whitson 1992, Fluid Phase
Equilibria 77, 217-240, DOI `10.1016/0378-3812(92)85105-H`:

- Table 3 n-C4: `Tc = 425.2 K`, `Pc = 38.0 bar`, `omega = 0.1931`;
- Table 5 n-C4/water non-aqueous BIP: `kij = 0.5091`;
- the aqueous pair uses the corrected-original SW92 hydrocarbon/water
  correlation already implemented and independently tested by the project.

Experimental three-phase pressures and the two hydrocarbon-phase n-butane mole
fractions are transcribed from:

H. H. Reamer, R. H. Olds, B. H. Sage and W. N. Lacey,
*Phase Equilibria in Hydrocarbon Systems. n-Butane-Water System in Three-Phase
Region*, Industrial & Engineering Chemistry 36 (1944) 381-383,
DOI `10.1021/ie50412a024`.

No experimental value is used as a fitted SW92 parameter.

## Independent oracle

`tests/flash/sw92_three_phase_experimental_line/reference_decimal.py` is a
stdlib-only Decimal(80) calculation with no production imports. For each fixed
experimental temperature it solves four unknowns:

```text
P,
x_nC4^W,
x_nC4^H0,
x_nC4^H1
```

from the four component chemical-potential equalities

```text
mu_nC4^W(AQ)  = mu_nC4^H0(NA)
mu_H2O^W(AQ)  = mu_H2O^H0(NA)
mu_nC4^W(AQ)  = mu_nC4^H1(NA)
mu_H2O^W(AQ)  = mu_H2O^H1(NA).
```

Each phase independently chooses the mechanically admissible minimum-Gibbs
cubic root inside its fixed AQ/NA family. H0/H1 remain unordered numerical phase
instances; the regression does not use their density/Z ordering to publish a
liquid/vapor morphology label.

Binary three-phase **phase fractions are intentionally not part of this
regression**. At a fixed point on a binary three-phase line the coexistence
phase properties are fixed, but a general overall composition does not uniquely
determine all three positive phase fractions. The multicomponent Sample-6
regression remains the phase-fraction validation.

## Experimental model-regression envelope

Eight points from 310.93 K through 416.48 K are used. The fixed acceptance
limits are:

```text
max |P_model - P_exp| / P_exp <= 0.02
max |x_nC4,H0(model) - x_nC4,H0(exp)| <= 0.01
max |x_nC4,H1(model) - x_nC4,H1(exp)| <= 0.005
```

These are **model-validation regression limits**, not flash convergence
 tolerances. They do not modify EOS equations, BIPs, root selection, TPD
thresholds, mass-balance tolerances or chemical-potential tolerances.

The independently solved exact-profile line has a maximum pressure relative
error of about 1.82%, maximum H0 composition error about 0.00934 and maximum H1
composition error about 0.00406 over those eight points.

Near the upper experimental endpoint the two non-aqueous phases coalesce. The
current exact project profile reaches its numerical coalescence somewhat before
the final experimental critical endpoint. Points above 416.48 K are therefore
reserved for a later explicit near-critical/coalescence diagnostic rather than
forcing the cubic solver to retain two distinguishable NA phases.

## Production regression

The C++ regression evaluates the existing production
`Sw92FamilyStabilityEvaluator` at every independently solved coexistence state
and requires:

- the AQ and both NA minimum-Gibbs roots to remain smooth;
- compressibility factors to match the Decimal oracle tightly;
- common reduced chemical potentials across all three phase instances;
- the experimental error envelope above;
- component-order permutation invariance.

This validation changes no production thermodynamic or flash code.

## Role in authoritative three-phase flash

The two physical validation pillars are now complementary:

1. Mortezazadeh-Rasaei Sample-6: multicomponent physical three-phase topology,
   common chemical potentials, unique material-balanced phase fractions and
   top-level routing;
2. Reamer et al. n-butane/water: direct experimental three-phase coexistence
   pressure and phase-composition validation of the exact SW92/PR76 profile.

After disappearance-neighbor re-solve is merged, these physical regressions are
sufficient to proceed to the explicit Profile-C `PtPhaseSetResult` publication
adapter. `global_stability_proven=false` must remain explicit because the
production stability review is finite, and H morphology remains unresolved
until a separately validated morphology classifier exists.

Frontend work remains frozen until that authoritative publication layer is
implemented and validated.
