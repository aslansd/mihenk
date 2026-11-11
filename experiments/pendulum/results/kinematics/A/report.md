# mihenk report: pendulum A with theta' = omega declared

*mihenk 0.1.1*

1 error(s), 4 warning(s): do not interpret this equation as the law yet

## Equation

```
theta' = 0.134 - 0.04273 theta + 0.9999 omega
omega' = 0.4616 - 0.1941 theta + 63.82 sin(theta) - 0.2435 cos(theta) + 0.4428 sin(omega)
```

## Findings

### ERROR `KINEMATIC_MISMATCH` (tier 0)
*where:* theta' = omega  
In the data, theta does not change by the integral of omega: omega runs late by 0.00085 s (0.85 samples). These explain 0.62% of the change in theta (F = 2.77e+03). The discovered equation is not theta' = omega: extra terms 1 (+0.134), theta (-0.04273) (they change theta' by 0.58% of omega). A omega late by 0.00085 s predicts extra terms equal to 0.00085 times the discovered omega' equation; that accounts for 98% of them. Correcting omega for this changes the law: omega' 1 +0.1457 -> +0.2283 (derivative-free refits of the equation's own terms).

*fix:* Check how omega is measured or filtered: estimate it from theta with a matched, zero-phase filter, or shift it by the reported lag before discovery. Alternatively, fix theta' = omega instead of discovering it, and interpret the other equations with the reported lag in mind.

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
| KINEMATIC_MISMATCH | 0 | flagged |  |
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
| kin_theta_lag | 0.00085 |
| kin_theta_gain | 1 |
| structure_reproduced | 0.1 |
| structure_under_variants | 0 |
| distinct_structures | 2 |

## Notes

- held-out data: held-out trajectories (val)
- derivative-free reference: integral form over windows of 50 samples
- candidate library for alternatives: 7 terms
