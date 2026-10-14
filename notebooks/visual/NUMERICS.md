# Numerical curves

Chlorination uses the [generated event model's finite forward equations](../../src/mfobs/numerics/context_nac.py) for the mean, variance and blocking probability. Total chlorine is four. Tail growth is unbounded but has zero jump on these coordinates, so the observable dynamics need no tail truncation.

Rémy integrates the closed coefficients produced by linear discovery with [the shared linear solver](../../src/mfobs/numerics/linear.py). DOP853 uses relative tolerance `1e-10` and absolute tolerance `1e-12`. Rémy starts with one leaf and its exterior stem.

## ABA: expectation and variance

The growing-path notebook uses the **exported event model** for both SSA and deterministic numerical integration. Changing the motif reruns support discovery and supplies the new events to the same solver. No equation or suffix list is supplied to the numerical solver.

[Notebook adapter: `integrate_events`](../paper_view.py) → [solver: `solve_moments`](../../src/mfobs/numerics/finite_state.py).

## What is integrated

Let `q` be the selected motif count and `z` all remaining coordinates. The solver checks that no rate depends on `q` and every jump of `q` is nonnegative. It then sums out `q` exactly. For each remaining state `z`, it evolves its probability `p`, its probability-weighted first query moment `m`, and its probability-weighted second query moment `s`:

```
p' = A p
m' = A m + B1 p
s' = A s + 2 B1 m + B2 p
```

`A` is the forward matrix constructed from the computed events. A transition `j → i` contributes `rate(j) * Δq` to `B1[i,j]` and `rate(j) * (Δq)^2` to `B2[i,j]`. Self-transitions in `z` that change `q` also contribute. The displayed mean is `sum(m)` and the variance is `sum(s) - sum(m)^2`. The initial path fixes `p`, `m` and `s`.

These are forward equations over states, not a claim that finitely many unconditional moments close. Sparse matrix exponentiation solves the truncated equations. The EMGF PDE is not discretised.

## Truncation and its bounds

Boolean coordinates and their reachable combinations are obtained from the event model. For every other state coordinate, the solver verifies that only constant-rate events can increase it. Their event counts give Poisson upper bounds. The solver also searches the event jumps for a query bound of the form `q ≤ q(0) + N(t)` or `q ≤ q(0) + z_i(0) + N(t)` (with a step multiplier when needed). The returned `curves['numerical']` record contains the verified envelope.

If `E` is the event of ever leaving the chosen finite region before the final time `T`, a union bound gives `P(E) ≤ ε`. If `B(T)` is the Poisson envelope of the query, Cauchy–Schwarz gives

```
missing mean          ≤ sqrt(E[B(T)^2] * ε) = e1
missing second moment ≤ sqrt(E[B(T)^4] * ε) = e2
absolute variance error ≤ e2 + 2 E[B(T)] e1.
```

This does not assume independence of the query and the exit event. Cutoffs are increased until both mean and variance bounds meet the requested tolerance, uniformly on `[0,T]`. Outgoing probability is killed, never reflected or renormalised. A small probability deficit alone would not bound unbounded count moments.

The figure reports **state-truncation bounds**. Floating-point matrix-exponential error is separate and is not included in these bounds. Positivity and probability-mass checks detect numerical failures. Tests compare against analytical Binomial/Poisson laws and independently enumerated finite word chains; the notebook's default run is also checked at a tighter truncation tolerance.

For a long time horizon, the a-priori Poisson box can be unnecessarily large. The solver then starts with a smaller box and integrates **first-exit flux** explicitly, rather than subtracting its mass from one. On an exited path, the query is bounded by its value plus the witnessed resource immediately after exit, followed by future Poisson immigration. Three additional equations accumulate the exited probability and the first two moments of this envelope. Their first/second moments bound the omitted query moments directly, including paths that later re-enter the box. The solver increases the boundary with the greatest outflow until both error bounds meet the requested tolerance. These bounds depend on numerical integration; its floating-point error remains separate.

The default tolerance is `1e-6` for mean and variance, with a budget of 150,000 retained states. Models outside the checked envelope conditions, or an exhausted state budget, produce an explicit error instead of an uncertified curve. This numerical consumer is generic within that supported class; it is not a solver for every possible polynomial jump model.

[Finite State Projection: Munsky and Khammash (2006)](https://doi.org/10.1063/1.2145882).

## Reading the plots

The red curves are computed from the forward equations. The blue mean and sample variance come from the same set of SSA paths. The shaded band around the SSA mean is an approximate pointwise Monte Carlo interval (`±1.96` standard errors); it is **not the process variance** or a simultaneous confidence band. Variance is shown separately in the third panel.

Set `T` (end time) and `RUNS` (independent SSA runs) in the notebook's settings cell. The rate sliders inherit both settings, as well as the random seed and numerical tolerance, from the current plot. They do not reset the time range or sample count. Larger `RUNS` reduces Monte Carlo uncertainty; it does not change the computed FSP curves. Only the first eight individual trajectories are drawn. Long horizons such as `T=100` take longer to compute and may need more states.

[Back to the notebook](../examples/02_aba.ipynb) · [Code map](CORE.md)
