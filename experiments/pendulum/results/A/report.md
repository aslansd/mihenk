# mihenk report: pendulum A: informed library, large swings

*mihenk 0.1.1*

4 warning(s): the equation may be real, but parts of it are not yet supported

## Equation

```
theta' = 0.134 - 0.04273 theta + 0.9999 omega
omega' = 0.4616 - 0.1941 theta + 63.82 sin(theta) - 0.2435 cos(theta) + 0.4428 sin(omega)
```

## Findings

### WARNING `DERIVATIVE_ARTIFACT` (tier 0)
*where:* omega'  
Coefficients that move by more than 30% when estimated derivative-free: the constant term +0.4616 -> +0.1457, theta -0.1941 -> -0.08566, sin(omega) +0.4428 -> +0.2656.

*fix:* Refit with smoothed or weak-form derivatives (e.g. PySINDy's SmoothedFiniteDifference) and report the derivative-free coefficients next to the discovered ones.

### WARNING `FORKING_PATHS` (tier 1)
*where:* analysis  
The equation's structure comes back under only 0 of 5 defensible analysis variants (2 distinct structures in total). Which law you report depends on analysis choices.

*fix:* Report the variants, or pick the analysis on held-out data before looking at the equation.

### WARNING `UNSTABLE_STRUCTURE` (tier 1)
*where:* theta'  
Refitting on 30 resamples of your data changes this equation's terms. Terms in the equation that often disappear: the constant term (43%), theta (43%). Terms absent from the equation that often appear: sin(theta) (80%) (share of refits). The structure depends on which part of the data the discovery saw.

*fix:* Treat these terms as tentative; check them on independent trajectories, or with more or less noisy data.

### WARNING `WEAK_TERM` (tier 0)
*where:* omega'  
Dropping the constant term (F = 0.0842), theta (F = 0.291) barely changes the fit (F below 10): the data cannot tell the equation with these terms from the equation without them.

*fix:* Do not interpret these terms without more or less noisy data; a sparser equation fits as well.

## Checks

| check | tier | status | note |
|---|---|---|---|
| EQUIFINALITY | 0 | clean |  |
| DERIVATIVE_ARTIFACT | 0 | flagged |  |
| WEAK_TERM | 0 | flagged |  |
| UNDERFIT | 0 | clean |  |
| ROLLOUT_TRAP | 0 | clean |  |
| UNSTABLE_STRUCTURE | 1 | flagged |  |
| FORKING_PATHS | 1 | flagged |  |
| RECOVERY | 2 | not_assessed | no true equations (synthetic data only) |

## Metrics

| metric | value |
|---|---|
| variables | 2 |
| trajectories | 6 |
| samples | 55001 |
| active_terms | 8 |
| derivative_r2 | 0.985 |
| integral_r2 | 0.9997 |
| library_relations | 0 |
| rollout_r2_short | 0.9998 |
| structure_reproduced | 0.1 |
| structure_under_variants | 0 |
| distinct_structures | 2 |

## Notes

- held-out data: held-out trajectories (val)
- derivative-free reference: integral form over windows of 50 samples
- candidate library for alternatives: 7 terms
