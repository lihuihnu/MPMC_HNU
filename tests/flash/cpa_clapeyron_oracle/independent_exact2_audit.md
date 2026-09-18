# Pinned compressed 2x2 X_exact2 independent audit

## Scope

This audit isolates only the pinned compressed two-variable association path for
`state=4 phase=liquid`.

Pinned external implementation:

`ClapeyronThermo/Clapeyron.jl@229b09452f36c2f812486150df0bb43b197bb4e5`

Frozen profile:

`CPA/SRK-physical/simplified-rdf-1.9eta/explicit-site-pairs/v1`

No pinned Clapeyron source, MPMC_HNU production code, CPA parameter, or frozen Gate-E threshold is modified.

Evidence workflow run:

`35300162361`

Evidence code head:

`5e7e167448566b0eca742893f42609ec1e3768ce`

Independent exact2 artifact digest:

`sha256:ea267112289aa2e92ec4437ccec2516483c2d27f27974852ba59a13450e62ed5`

## Frozen independent solver

At each evaluated volume:

1. build the same original four-site association matrix using the pinned
   `delta_assoc` and `assoc_site_matrix`;
2. apply the pinned unmodified matrix compression;
3. require the same reduced dimension 2 and mapping `[1,1,2,2]`;
4. solve the two reduced equations independently, without calling
   `X_exact2!`, `assoc_matrix_x0!`, `assoc_matrix_solve`, or another Clapeyron nonlinear solver.

The equations are

`R_i(x)=x_i*(1+sum_j K_ij*x_j)-1=0`.

The independent analytic Jacobian is

`J_ij=delta_ij*(1+(Kx)_i)+x_i*K_ij`.

Frozen Newton settings:

- precisions: 256 and 512 BigFloat bits;
- initial point: `(0.5,0.5)`;
- residual target: `max(abs(R)) <= 1e-30`;
- maximum iterations: 100;
- line search: `lambda=1,1/2,...,2^-40`;
- accepted steps must stay in `0 < x_i <= 1` and strictly reduce the residual.

The derivative probe is also frozen:

- centered 7-point O(h^6);
- `h/V=1e-6`.

## Independent Newton convergence

Both 256- and 512-bit runs converge in 8 Newton iterations.

At 512 bits the reduced solution is

`x1 = 0.0721735494520661912926687691021456023169856014113249370263879500885608082682036558...`

`x2 = 0.1830794007502389035605951321528209992301732176136659846492497120149274240695468955...`

and the expanded four-site solution is

`[x1,x1,x2,x2]`.

The reduced and expanded maximum mass-action residual are both

`3.5986384120953434785537598891531440039e-41`,

far below the frozen `1e-30` target.

The pinned compressed-`X_exact2!` solution differs by only

`max |delta X| = 1.4872759452715762818352456175940970e-14`.

Thus the pinned site fractions look numerically extremely close, but they are not at the high-precision stationary solution.

## Independent Q / Helmholtz reconstruction

Using the original four-site matrix and the expanded independent X, the audit evaluates the same stationary association functional independently:

`Q2 = sum_i w_i*(ln(X_i)+1-X_i)`

`Q1 = -0.5*sum_i w_i*X_i*(KX)_i`

`Q = Q1+Q2`.

For the high-precision solution:

`Q = -4.5078513457715525951250774798653055438585150859253612443813373322069802362382149...`

The direct equilibrium association expression

`sum_i w_i*(ln(X_i)-X_i/2+1/2)`

agrees with Q to

`3.1016310870987952238568140563351255e-41`.

This independently confirms stationarity.

For the pinned `X_exact2!` site fractions, the same independent formula has a stationary/direct mismatch of

`9.88150354437284680389900288076314e-14`.

That is consistent with the previously observed `~1.147e-13` expanded mass-action residual.

### Center scalar difference

Replacing pinned X with the independent stationary X changes the center association residual Helmholtz energy by only

`2.7371401388441839350367789199255574e-10 J`.

The full residual Helmholtz energy changes by exactly the same amount because the cubic term is unchanged.

This is why the scalar `F_res` comparison can remain excellent even when the derivative is affected.

## V-derivative result

At every shifted volume `V +/- {1,2,3}h`, the compressed association equations are independently re-solved to approximately `3.6e-41` residual before evaluating the scalar.

At 512 bits:

Independent association derivative:

`dA_assoc/dV = 1.26111591029416424063346531235030199692061240434590495801757070317555762517749e8 Pa`

Pinned-`X_exact2!` association derivative:

`dA_assoc/dV = 1.26111591029367485713055978332145607858889485587363004143934487704623483519134e8 Pa`

Difference:

`4.8938350290552902884591833171754847e-5 Pa`.

The cubic scalar path is unchanged, so the full residual derivative changes by the same amount:

`abs(delta dA_res/dV) = 4.8938350290552902884591833171754847e-5 Pa`.

## Pressure consequence

With the independent stationary solution:

`P_independent = 72831.99999998462551138330133107220996838668980034932837... Pa`.

Frozen target pressure:

`P_target = 72832 Pa`.

Therefore

`abs(P_independent-P_target) = 1.537448861669867e-8 Pa`.

This is more than two orders of magnitude inside the frozen Gate-E pressure threshold

`5e-6 Pa`.

By contrast, the same 7-point scalar derivative using the pinned `X_exact2!` path gives

`P_pinned_fd = 72832.0000489229758019362042156640431401415... Pa`

and

`abs(P_pinned_fd-P_target) = 4.8922975801936206e-5 Pa`.

The independent correction relative to the pinned scalar derivative is

`4.8938350290552902884591833171754847e-5 Pa`.

This accounts for essentially the entire state-4 pressure discrepancy.

## Precision reproducibility

The 256- and 512-bit runs produce the same result at every scale relevant to Gate E:

- 8 Newton iterations;
- reduced/expanded residual approximately `3.5986e-41`;
- `max |delta X| ~= 1.4873e-14`;
- center association energy delta `~=2.73714e-10 J`;
- derivative correction `~=4.893835e-5 Pa`;
- corrected pressure error vs target `~=1.53745e-8 Pa`.

The result is therefore not a precision-selection artifact.

## Answer to the audit question

Yes: the pinned compressed `X_exact2!` numerical path is sufficient to explain essentially all of the remaining state-4 pressure offset.

The causal chain supported by the audit is:

1. the pinned 4-site problem is compressed to 2x2;
2. `X_exact2!` returns a site solution with about `1.5e-14` error in X and about `1e-13` stationarity residual;
3. that error changes the scalar association energy by only about `2.7e-10 J`;
4. because the error varies with volume, it changes `dA_assoc/dV` by about `4.893835e-5 Pa`;
5. independently solving the same compressed equations to `~3.6e-41` residual removes that derivative offset and brings pressure to within `1.54e-8 Pa` of the frozen target.

This directly explains the earlier pattern in which scalar Helmholtz and residual chemical potential evidence were green while pressure and the pressure-dependent Z contribution were outside the extremely tight frozen external-oracle envelope.

## Gate consequence

This audit does **not** modify the pinned Clapeyron oracle or the already-frozen Gate-E acceptance contract.

Under the literal frozen contract, the unmodified pinned implementation still produces the old pressure / ln(phi) values, so Gate E remains **BLOCKED** and Gate F remains **NOT READY**.

However, the blocker is now attributable to a specific numerical implementation path in the pinned external oracle rather than to:

- MPMC_HNU's scalar residual Helmholtz formulation;
- Float64 cancellation;
- ordinary finite-difference error;
- configurable association `rtol/atol/max_iters`;
- or a mismatch in the frozen CPA formulation/parameters.

No threshold widening and no production source-of-truth switch are authorized by this audit.
