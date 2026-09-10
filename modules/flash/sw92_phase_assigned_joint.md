# SW92 Profile-C fixed AQ/NA joint two-phase primitive

## Scope

This increment implements Gate C1 from `sw92_equilibrium_profile_separation_audit.md` under the
reserved equilibrium profile

```text
SW92-equilibrium/phase-assigned-aq-na-joint/v1
```

and numerical primitive identity

```text
SW92-equilibrium/phase-assigned-aq-na-joint/fixed-two-phase-logK-SSI-RR/v1
```

It solves exactly one constrained two-phase topology:

- physical aqueous phase -> `SwPhaseFamily::aqueous`;
- physical nonaqueous phase -> `SwPhaseFamily::nonaqueous`.

It does **not** perform autonomous phase-number selection, single-phase classification, final
constrained two-family stability, three-phase flash, salt-inventory conservation, derivatives, or
physics closure. `global_stability_proven` and `final_stability_checked` are therefore both false.

## Separation from Profile A and Profile B

This primitive is not the Whitson dual-model observable route. Both phases participate in one
material balance and one set of cross-phase reduced chemical-potential equalities; no AQ-pass and
NA-pass results are stitched together.

It is also not the Xu lower-envelope route. The physical role fixes which family evaluates each
phase. The solver does not evaluate the opposite family at a converged phase in order to reassign
that phase and does not use Profile-B `family_assignment_dominated` as an acceptance criterion.

The same mathematical AQ+NA equation solution may therefore have different publication semantics:

- Profile B can retain the equation point but reject it when the NA assignment is above the AQ
  lower-envelope surface;
- Profile C can retain it as a C1 candidate when the fixed physical roles, material balance,
  cross-phase equilibrium, root smoothness and physical-role topology checks pass.

That difference is intentional and is covered by a focused regression.

## Continuous problem

For active EOS components, the physical aqueous composition is `x`, the physical nonaqueous
composition is `y`, and `beta` is the nonaqueous mole fraction. Define

```text
logK_i = ln(y_i / x_i).
```

For fixed `logK`, ordinary two-phase material balance is solved with the existing Rachford-Rice
kernel:

```text
x_i = z_i / (1 - beta + beta K_i)
y_i = K_i x_i
sum_i z_i (K_i - 1) / (1 - beta + beta K_i) = 0.
```

The legacy Rachford-Rice field names are private implementation details. Public Profile-C output
immediately maps them to `aqueous_phase`, `nonaqueous_phase`, and `nonaqueous_phase_fraction`.

At each state the solver evaluates the mechanically admissible minimum-Gibbs root **within the
assigned family** and forms

```text
m_i^AQ = ln(x_i) + ln(phi_i^AQ(x))
m_i^NA = ln(y_i) + ln(phi_i^NA(y))
r_i    = m_i^AQ - m_i^NA.
```

The pair equations are converged only if `max_i |r_i|`, absolute material balance and relative
material balance all pass their declared thresholds.

The SSI/logK direction and Gibbs/residual line search follow the same audited fixed-family-pair
thermodynamic derivation used in Gate 3B.1, but this implementation has an independent Profile-C
result/status contract and contains no lower-envelope family-selection stage.

## Physical role and family are separate metadata

Every retained phase records both

```text
physical_role
thermodynamic_family
```

C1 fixes their mapping as

```text
aqueous    -> SwPhaseFamily::aqueous
nonaqueous -> SwPhaseFamily::nonaqueous.
```

Compressibility factor, root index, vector slot, and density ordering are diagnostics, not physical
role selectors.

For the C1 two-phase topology, water must be robustly richer in the aqueous-role phase:

```text
x_water^AQ > x_water^NA.
```

No absolute water-fraction cutoff is used. The comparison uses a roundoff-scaled guard. A resolved
reversal returns `phase_role_reversed`; an unresolved difference returns
`phase_role_indeterminate`. Neither may publish a C1 candidate.

This relative ordering is a constrained two-phase topology check only. It does not solve the
separate single-phase role problem identified as blocked by the profile-separation audit.

## Candidate versus accepted equilibrium

`Sw92PhaseAssignedJointResult::candidate()` is deliberately narrower than an accepted phase-set
API. It is non-null only when C1 equations, balance, phase distinction, same-family root smoothness,
and physical-role ordering pass.

C1 does not call this an accepted equilibrium phase set because Gate C2 final constrained
common-tangent stability is not yet implemented. A caller must not infer global stability from a C1
candidate.

## Independent Decimal(80) references

`tests/flash/sw92_phase_assigned_joint/reference_decimal.py` is stdlib-only and imports no
production code or Profile-A/Profile-B result. It independently rebuilds the corrected-original
binary SW92/PR equations and directly solves the cross-family AQ+NA equations.

### CO2/H2O

```text
p = 3 MPa
T = 340 K
m_NaCl = 0 mol/kg H2O
z_CO2 = 0.7

x_CO2^AQ = 0.00594506516618610144073756823151865...
y_CO2^NA = 0.988787734309055310658555988544970...
beta_NA  = 0.706170943350571824136142980382328...
Z_AQ      = 0.0232691618225347655031695609311701...
Z_NA      = 0.886575325995761419715861484043613...
G_pair/RT = -1.49504682831804812553073260552711...
common gas   = -0.121620504033995572408892284829747...
common water = -4.69970825164750408281502668715430...
```

### CH4/H2O

```text
p = 10 MPa
T = 350 K
m_NaCl = 1 mol/kg H2O
z_CH4 = 0.5

x_CH4^AQ = 0.000910799547271048457583638344068579...
y_CH4^NA = 0.994678544467031754705810261260436...
beta_NA  = 0.502219158353772754133143176267965...
Z_AQ      = 0.0755672461713566851072962391061499...
Z_NA      = 0.910625639821077058120594969177001...
G_pair/RT = -2.78260864163829488520641350877109...
common gas   = -0.113632542098874397480096323695159...
common water = -5.45158474117771537293273069384702...
```

These are direct model-equation references. They are not obtained by combining the two independent
Profile-A flashes and are not Profile-B AQ+AQ anchors.

## Failure and resource semantics

The primitive retains explicit statuses for no interior Rachford-Rice state, phase disappearance,
indistinguishable phases, family-root nonsmoothness, physical-role reversal/indeterminacy, property
failure, iteration/evaluation limits, and line-search failure. A phase-disappearance state is
classified only after the chemical-potential and material-balance thresholds have passed.

Property budgets count assigned-family property evaluations only. No opposite-family Gibbs check is
performed after convergence, so Profile-B lower-envelope diagnostics cannot consume or bias the C1
resource budget.

## Next gate

Gate C2 should take a C1 candidate, rebuild its common tangent, include both converged phase
compositions in both assigned-family finite-search start sets, and perform constrained AQ and NA
final stability against that common tangent. Only after C2 may a Profile-C API publish an accepted
two-phase phase set; `global_stability_proven=false` must remain explicit for finite multistart
searches.
