# PT component inventory and local Jacobian

## Scope

`mpmc/physics/component_inventory.hpp` consumes an already owned thermodynamic-closure snapshot and constructs a **model-neutral, local fluid-volume component inventory**. It does not create pore volume, saturation, mesh, flux, source/sink terms, time stepping or a global residual/Jacobian assembly.

Two existing closure payloads are supported without changing either upstream contract:

- fixed liquid/vapor `ThermodynamicClosureSnapshot` (current PR76 path);
- variable-cardinality `PtPhaseSetThermodynamicClosureSnapshot` (current SW92 Profile-C path).

The public convention is:

```text
PT/component-inventory/phase-resolved-reduced-feed/v1
```

## Primal definition

For accepted phases `alpha=0..P-1`, the thermodynamic closure already owns mole phase fractions `beta_alpha`, ordered phase compositions `x_alpha,i`, and phase molar densities `c_alpha [mol/m^3]`.

Define the mixture molar volume per total mole by

```text
v_bar = sum_alpha beta_alpha / c_alpha
```

and the total fluid molar density by

```text
c_mix = 1 / v_bar.
```

The phase-resolved component mole fraction is

```text
s_i = sum_alpha beta_alpha x_alpha,i.
```

The component inventory per **total fluid volume** is

```text
a_i = c_mix s_i  [mol/m^3 fluid].
```

Because the accepted flash state must satisfy material balance, the implementation independently requires

```text
s_i ~= z_i
```

and checks the equivalent feed form

```text
a_i ~= c_mix z_i.
```

This is not pore-volume accumulation. No porosity or saturation is implied.

## Reduced-feed coordinates and Jacobian

The local derivative coordinates are inherited unchanged from the source closure:

```text
q = (p_Pa, T_K, z_0, ..., z_{N-2})
z_{N-1} = 1 - sum(z_0, ..., z_{N-2}).
```

For each column `q_j`,

```text
dv_bar/dq_j = sum_alpha [
    (d beta_alpha/dq_j) / c_alpha
    - beta_alpha (d c_alpha/dq_j) / c_alpha^2
]

dc_mix/dq_j = -c_mix^2 dv_bar/dq_j

ds_i/dq_j = sum_alpha [
    (d beta_alpha/dq_j) x_alpha,i
    + beta_alpha (d x_alpha,i/dq_j)
]

da_i/dq_j = (dc_mix/dq_j) s_i + c_mix (ds_i/dq_j).
```

The implementation also constructs the independent feed-form derivative

```text
da_i/dq_j = (dc_mix/dq_j) z_i + c_mix (dz_i/dq_j)
```

and requires the two forms to agree. Therefore phase-fraction, composition and phase-density Jacobians all participate in the production validation; the inventory layer does not silently bypass them with the overall feed.

## Availability and failure policy

Primal inventory availability and local-Jacobian availability remain separate decisions.

- If the source thermodynamic primal is valid, a component inventory may be published even when the source linearization is unavailable.
- If the source linearization is unavailable because of a phase-boundary guard, ill conditioning, unsupported feed support, property failure or derivative arithmetic failure, the inventory remains usable but the inventory Jacobian remains unavailable with the propagated reason.
- If the source payload is malformed, phase/feed material balance is inconsistent, or differentiated phase/feed identities fail, no zero/stale/fallback Jacobian is published.
- Invalid inventory integrity options are caller errors and throw `std::invalid_argument`.

The layer exposes `inventory_available()` and `linearization_available()`. It deliberately does **not** claim that a local inventory Jacobian can seed a future global Newton solve; that decision belongs to the future residual/discretization/solver layers.

## 0D fixed-volume integration test

The focused regression uses a test-only isothermal, fixed-fluid-volume component residual

```text
R_i = V a_i - N_i_target
```

with unknown coordinates

```text
u = (p, z_0, ..., z_{N-2}).
```

The Jacobian is assembled from the inventory linearization by selecting the pressure column and reduced-feed columns while holding temperature fixed:

```text
dR_i/du_j = V da_i/du_j.
```

This 0D residual is an integration test only. It is not a production mesh/discretization API and does not introduce porosity, Darcy flux, transmissibility, relative permeability, time accumulation or a global nonlinear system.

## Validation

The focused suite covers:

1. synthetic two-phase algebra and exact phase-slot permutation invariance;
2. current PR76 accepted two-phase fixed-VLE closure;
3. SW92 Profile-C one-phase and physical Sample-6 three-phase closures;
4. runtime component permutation for Sample-6 inventory and p/T Jacobian rows;
5. derivative-unavailable propagation while retaining a valid inventory primal;
6. phase/feed tamper rejection and invalid-option rejection;
7. PR76 and SW92 binary 0D fixed-volume Jacobians cross-checked against fresh symmetric re-solves in pressure and reduced feed;
8. public-header self containment.

The local Jacobian remains valid only inside the source closure's fixed differentiable state. It does not cross phase appearance/disappearance, selected-root/family switching, critical/near-multiple-root states or Profile-C topology switching.
