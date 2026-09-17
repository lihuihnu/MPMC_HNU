# CPA no-cold-reference association observable certificate audit

## Status

This document derives a **test-only a-posteriori certificate** for the production
safety contract

```text
MPMC_HNU/CPA/association-continuation-safety/v1
```

The purpose is to decide whether a warm association fixed-point candidate can be
accepted **without first running the current cold association solve at the same
density**.

Nothing in this slice enables continuation in production. The production
association equations, `1e-12` stopping tolerance, density scan, root solver,
stability solver, PT flash solver and failure semantics remain unchanged.

## 1. Association fixed-point map

For the configured site classes, index sites by `a,b`. For site `b`, define

```text
w_b = x_component(b) * multiplicity_b
```

and, at fixed `(T,rho,x)`,

```text
A_ab = rho * w_b * Delta_ab >= 0.
```

`Delta_ab` is the same explicit CPA association strength used by production.
The production mass-action equation is then

```text
F_a(X) = 1 / (1 + sum_b A_ab X_b)
X*     = F(X*).
```

The existing damping changes the iteration path only; it does not change this
fixed point.

For a converged warm candidate `Xhat`, define the fresh, post-convergence
residual

```text
h(Xhat) = Xhat - F(Xhat).
```

The certificate recomputes this residual once. It does not use a cold solution.

## 2. Local root enclosure

Define

```text
H(X) = X - F(X).
```

Its Jacobian is

```text
J_ab(X) = delta_ab + A_ab / (1 + sum_c A_ac X_c)^2.
```

For a positive box

```text
B = [Xhat-r, Xhat+r],
```

each Jacobian entry has a direct enclosure because every `A_ab` is
nonnegative. If `l=Xhat-r` and `u=Xhat+r`, then

```text
Dmin_a = 1 + sum_c A_ac l_c
Dmax_a = 1 + sum_c A_ac u_c

Jlow_ab  = delta_ab + A_ab / Dmax_a^2
Jhigh_ab = delta_ab + A_ab / Dmin_a^2.
```

The test implementation uses the inverse of `J(Xhat)` as a preconditioner `B0`
and forms a Krawczyk-style enclosure

```text
K = Xhat - B0 H(Xhat) + (I - B0 J(B)) (B-Xhat).
```

With midpoint/radius Jacobian bounds, a conservative componentwise radius for
the second term is

```text
M = abs(I - B0 Jmid) + abs(B0) Jrad.
```

The box is accepted only when

```text
abs(-B0 H(Xhat)) + M r < r
```

component by component and every lower site-fraction bound remains positive.
The implementation enlarges `r` monotonically until this inclusion succeeds or
certification fails.

This produces a candidate-independent enclosure

```text
|X*_a - Xhat_a| <= r_a.
```

The current implementation evaluates the certificate in `long double` and uses
explicit numerical inflation. It is an engineering a-posteriori certificate,
not a formal directed-rounding interval proof. Before production enablement the
same algebra must either retain demonstrably conservative floating-point
inflation on all supported platforms or use outward-rounded interval arithmetic.

## 3. Pressure-error certificate

At fixed `(T,rho,x)`, the physical SRK pressure does not depend on `X`. The
association pressure is

```text
P_assoc = -C_P * sum_a w_a (1-X_a)

C_P = 0.5 * R*T*rho * (1 + rho*dln(g)/drho).
```

Therefore the site enclosure gives the direct no-cold bound

```text
|Delta P| <= C_P * sum_a w_a r_a.
```

A warm candidate can satisfy the v1 observable pressure scale only if this bound
is no larger than the frozen

```text
5e-6 Pa
```

contract value.

## 4. `ln(phi)` certificate

At fixed target pressure, temperature, density and composition, the cubic
chemical-potential term and `-ln(Z)` do not change with the association site
fractions. Production uses

```text
mu_assoc_i = sum_(a in i) m_a ln(X_a)
             - K_i * sum_a w_a (1-X_a)

K_i = (1.9/8) * rho * b_i * g.
```

For the certified box, define

```text
l_a = Xhat_a-r_a
u_a = Xhat_a+r_a.
```

Because `l_a>0`, the logarithmic contribution has the exact box bound

```text
L_a = max( ln(Xhat_a/l_a), ln(u_a/Xhat_a) ).
```

Thus, for component `i`,

```text
|Delta ln(phi_i)| <=
    sum_(a in i) m_a L_a
    + K_i * sum_a w_a r_a.
```

The no-cold candidate satisfies the v1 fugacity scale only when the maximum of
these component bounds is no larger than

```text
1e-10.
```

## 5. What the certificate does and does not prove

If the local enclosure, pressure bound and `ln(phi)` bound all pass, the warm
candidate has a certified nearby association fixed point whose affected phase
observables lie inside the predeclared v1 numerical scales **without comparing to
a cold solve**.

This certificate is deliberately local. It does not by itself prove that the
association equations have no other fixed point outside the certified box.
Production continuation therefore still relies on the v1 ownership/continuity
chain: every root-search cache starts empty, the first density uses the existing
cold path, and only accepted states may seed later densities in the same root
search.

The certificate also does not by itself authorize a change in root topology.
When production integration is considered, pressure-sign decisions near the
root-search tolerance boundary must either be certified from the pressure
interval or fall back to the current cold evaluation.

## 6. Validation protocol

The existing 333.15 K methanol-water five-state continuity replay is reused.
For every seeded warm candidate the audit will:

1. compute the certificate from the warm state only;
2. record enclosure success/failure;
3. record pressure-bound and `ln(phi)`-bound pass/fail counts;
4. separately compute the existing cold result **for validation only**;
5. verify that observed warm/cold site, pressure and root-`ln(phi)` differences
   are consistent with the independently computed enclosure/bounds whenever the
   corresponding candidate is certified;
6. report how many of the 305,431 seeded density evaluations could be accepted
   without a cold reference under the frozen v1 scales.

Cold values are not inputs to the certificate decision. They are retained only
as an independent validation oracle for this audit.

No performance claim will be made from certificate-audit wall time. The first
question is whether the mathematical acceptance test is sufficiently safe and
sufficiently selective to justify a later production design.
