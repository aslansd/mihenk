# mihenk report: pendulum C with theta' = omega declared

*mihenk 0.1.1*

3 error(s), 3 warning(s): do not interpret this equation as the law yet

## Equation

```
theta' = 0.9916 omega + 0.008537 sin(omega)
omega' = 23.15 + 3.933 theta + 15.36 omega + 68.86 sin(theta) + 35.36 cos(theta) - 16.31 sin(omega) - 0.1821 cos(omega)
```

## Findings

### ERROR `EQUIFINALITY` (tier 0)
*where:* theta'  
A sparser equation fits your data exactly as well (1 term instead of 2): theta' = 0.9998 omega. On the observed states some candidate terms are combinations of others, so the data cannot choose between these laws.

*fix:* Collect data that leaves the current set of states (other initial conditions, transients, perturbations) so the alternatives make different predictions.

### ERROR `EQUIFINALITY` (tier 0)
*where:* omega'  
A sparser equation fits your data exactly as well (3 terms instead of 7): omega' = 16.06 omega + 64.78 sin(theta) - 17.02 sin(omega). On the observed states some candidate terms are combinations of others, so the data cannot choose between these laws.

*fix:* Collect data that leaves the current set of states (other initial conditions, transients, perturbations) so the alternatives make different predictions.

### ERROR `UNSTABLE_STRUCTURE` (tier 1)
*where:* theta'  
Refitting on 30 resamples of your data changes this equation's terms. Terms in the equation that often disappear: sin(omega) (23%) (share of refits). The structure depends on which part of the data the discovery saw.

*fix:* Treat these terms as tentative; check them on independent trajectories, or with more or less noisy data.

### WARNING `DERIVATIVE_ARTIFACT` (tier 0)
*where:* omega'  
Coefficients that move by more than 30% when estimated derivative-free: sin(theta) +68.86 -> +1337.

*fix:* Refit with smoothed or weak-form derivatives (e.g. PySINDy's SmoothedFiniteDifference) and report the derivative-free coefficients next to the discovered ones.

### WARNING `FORKING_PATHS` (tier 1)
*where:* analysis  
The equation's structure comes back under only 0 of 5 defensible analysis variants (1 distinct structures in total). Which law you report depends on analysis choices.

*fix:* Report the variants, or pick the analysis on held-out data before looking at the equation.

### WARNING `WEAK_TERM` (tier 0)
*where:* omega'  
Dropping sin(theta) (F = 1.26) barely changes the fit (F below 10): the data cannot tell the equation with these terms from the equation without them.

*fix:* Do not interpret these terms without more or less noisy data; a sparser equation fits as well.

## Checks

| check | tier | status | note |
|---|---|---|---|
| EQUIFINALITY | 0 | flagged |  |
| DERIVATIVE_ARTIFACT | 0 | flagged |  |
| WEAK_TERM | 0 | flagged |  |
| UNDERFIT | 0 | clean |  |
| ROLLOUT_TRAP | 0 | clean |  |
| KINEMATIC_MISMATCH | 0 | clean |  |
| UNSTABLE_STRUCTURE | 1 | flagged |  |
| FORKING_PATHS | 1 | flagged |  |
| RECOVERY | 2 | not_assessed | no true equations (synthetic data only) |

## Metrics

| metric | value |
|---|---|
| variables | 2 |
| trajectories | 1 |
| samples | 23001 |
| active_terms | 9 |
| derivative_r2 | 0.1583 |
| integral_r2 | 0.9889 |
| library_relations | 3 |
| rollout_r2_short | 1 |
| kin_theta_lag | 5e-05 |
| kin_theta_gain | 0.9998 |
| structure_reproduced | 0.2333 |
| structure_under_variants | 0 |
| distinct_structures | 1 |

## Notes

- held-out data: last 20% of each trajectory
- derivative-free reference: integral form over windows of 50 samples
- candidate library for alternatives: 7 terms
