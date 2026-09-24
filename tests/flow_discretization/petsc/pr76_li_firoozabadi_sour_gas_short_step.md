# Li-Firoozabadi sour-gas three-phase non-isothermal serial short-step

## Scope

This regression closes the first source-complete three-phase multicomponent flow
benchmark for PR #118 without adding a new physical model or control mode.

The production solve is deliberately stationary:

- one Cartesian cell;
- six components: CO2, N2, H2S, CH4, C2H6, C3H8;
- three frozen PR76 equilibrium phases;
- P = 2.0 MPa, T = 178.8 K;
- dt = 1 s;
- no internal faces;
- fixed BHP exactly equal to the common phase pressure.

Therefore the independent Peaceman well-rate oracle is exactly zero and the
backward-Euler component/energy inventories must remain unchanged.  The test is
still non-isothermal: phase h/u and stationary-rock thermal storage are evaluated
and differentiated with respect to T.

## Thermodynamic state

The equilibrium anchor is the existing repository benchmark:

Li, Z. and Firoozabadi, A. (2012), *SPE Journal* 17(4), 1096-1107,
DOI `10.2118/129844-PA`.

`tests/flash/pr76_three_phase/reference_sour_gas_decimal.py` independently
solves the three-phase PR76 equilibrium using Python standard-library Decimal and
does not import production C++.  The short-step oracle imports that independent
solve rather than copying production phase compositions.

## Transport and caloric model boundary

No new transport/caloric model is introduced.  The six-component provider uses
the same model chain already used by the repository methane/ethane/propane
provider:

- selected PR76 density/fugacity branch;
- Stiel-Thodos dilute viscosity;
- Herning-Zipperer dilute-mixture blending;
- Kay pseudo-critical mixing;
- Lohrenz-Bray-Clark dense-fluid correction;
- sourced ideal-gas Cp integration;
- PR76 residual-enthalpy departure;
- `u = h - p/rho_mass`.

LBC reference: Lohrenz, Bray and Clark (1964), *JPT* 16(10), 1171-1176,
DOI `10.2118/915-PA`.

Low-temperature ideal-gas Cp data are restricted to the explicit common
100-298.15 K interval:

- CO2: NIST-JANAF table C-095;
- N2: NIST Chemistry WebBook SRD 69 Shomate equation, 100-500 K;
- H2S: NIST-JANAF table H-080;
- CH4/C2H6/C3H8: NIST Chemistry WebBook SRD 69 recommended gas Cp tables.

Molar masses and the critical volumes/densities required by LBC are likewise
taken from NIST Chemistry WebBook SRD 69.  The provider rejects temperatures
outside the sourced 100-298.15 K caloric interval.

The nonreactive enthalpy reference is
`h_i^ig(298.15 K) = 0` for every component.  This is a flow energy reference,
not a heat-of-formation model.

## Rock thermal storage

The stationary rock is the 99.93% quartz sample in NIST Structural Ceramics
Database SRD 30, citation Z00788:

Anderson, C. T. (1936), *J. Am. Chem. Soc.* 58, 568-570,
DOI `10.1021/ja01295a008`.

The source reports density 2.6378 g/cm3 and Cp = 469.3 J/(kg K) at 169.1 K and
510.6 J/(kg K) at 184.8 K.  The regression linearly interpolates Cp only inside
that measured interval.  Rock internal energy is referenced to zero at 178.8 K;
the production Jacobian retains the sourced `rho_rock * Cp(T)` derivative.

## Cell and well geometry

The engineering geometry is taken from the OPM copy of SPE1/Odeh input,
`OPM/opm-data@00472476647ccd99cb70bb585300b5d282ed368b`,
`spe1/SPE1CASE1.DATA`:

- dx = dy = 1000 ft;
- first-layer dz = 20 ft;
- porosity = 0.3;
- first-layer Kx = Ky = 500 mD;
- the deck explicitly records a 0.5 ft well-bore internal diameter, hence
  `r_w = 0.25 ft`.

For this single-cell benchmark the same 500 mD value is supplied in the axial
slot because the Peaceman z-well index consumes only transverse Kx/Ky; no claim
about Odeh's unavailable Kz measurement is made.

Peaceman geometry continues to use the repository's existing 1983 contract,
DOI `10.2118/10528-PA`.

## Relative-permeability boundary

There is no literature V-L1-L2 relative-permeability dataset for this equilibrium
benchmark.  The PETSc production evaluator requires a three-phase saturation
constitutive callback, so the test supplies a local linear `kr_alpha=S_alpha`
carrier solely to satisfy that existing software interface.

It is **not** part of the external physical oracle:

- there are no internal faces;
- capillary pressure is disabled;
- `p_bhp = p_alpha` for every phase, so every well phase rate is exactly zero;
- all validated state/property/inventory quantities are independent of the
  structural kr value.

The test must never be cited as relative-permeability validation.

## Independent oracle and acceptance

`tests/flow/core/reference_li_firoozabadi_sour_gas_flow.py` recomputes at
Decimal(80) and Decimal(96):

- the existing independent three-phase equilibrium;
- phase volume fractions converted to saturations;
- PR76 molar/mass density;
- LBC viscosity;
- ideal + residual enthalpy and internal energy;
- cell component inventories and total internal energy;
- quartz volumetric heat capacity.

The generated header is frozen only after both precisions render the same
fixed-digit result.

Production acceptance requires:

1. the six-component provider matches the independent oracle;
2. the serial production SNES uses the repository NewtonLS + GMRES + ASM path;
3. the converged state remains on the independent three-phase stationary state;
4. Peaceman phase/component/energy rates remain zero;
5. all six component inventories and total energy remain at the independent
   reference values;
6. a fresh final residual evaluation closes without changing production
   tolerances.

This benchmark does not add rate control, minimum-BHP switching, schedules,
phase transitions, multiple cells, MPI decomposition, or any new EOS,
viscosity, caloric, capillary or relative-permeability model.


## Serial ↔ 2-rank decomposition-invariance extension

The same source-complete stationary state is also used in a 2×1 two-cell
decomposition regression.  This is not a second physical benchmark: both cells
reuse the exact 20-bar / 178.8-K Li-Firoozabadi equilibrium, the same sourced
transport/caloric closure, the same quartz rock storage and the same SPE1/Odeh
cell dimensions/porosity/permeability.

The two cells are connected by one explicitly materialized x-normal TPFA face.
For identical Cartesian cells,

`T_f = K_x * (dy * dz) / dx`

is derived directly from the already-cited SPE1 geometry/permeability.  Gravity
is zero, both sides have identical pressure/temperature/phase state, and thermal
face conductance is explicitly disabled rather than assigned an unsourced
conductivity.  The independent stationary face oracle is therefore

- phase volumetric flux = 0;
- every component molar flux = 0;
- total internal-face energy rate = 0.

Cell 10 retains the same fixed-BHP Peaceman completion with
`p_bhp = p_cell = 2 MPa`, so its independent well-rate oracle remains zero.
Cell 20 has no external source.

The regression solves the *same* 38-scalar two-cell nonlinear problem twice:

1. serial: both cells owned on `PETSC_COMM_SELF`;
2. distributed: rank 0 owns cell 10 and the authoritative internal face, rank 1
   owns cell 20, and each rank retains the other cell as the local ghost copy.

Before the nonlinear solve, the test assembles the physical Jacobian and requires
a nonzero cross-cell block.  This prevents the zero-flux stationary state from
passing without exercising the internal-face coupling.  In the serial solve the
coupling resides in the local diagonal MPIAIJ block; in the 2-rank solve the same
coupling crosses the rank boundary.

After NewtonLS + GMRES + ASM convergence, the regression compares:

- both stable-cell natural-variable states `p/T/S/x`;
- the frozen three physical phase identities;
- phase/component internal-face rates and total face energy rate;
- fixed-BHP phase/component/energy well rates;
- all six global component inventories;
- global total internal energy;
- fresh final residual L2 norm;
- the aggregate absolute cross-cell Jacobian coupling.

All physical zero-rate and inventory references remain those of the existing
source-complete oracle.  No new EOS, transport/caloric, relative-permeability,
capillary, well-control, phase-transition or solver model is introduced by this
decomposition test.
