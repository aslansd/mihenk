# mihenk report: pysindy

*mihenk 0.1.1*

5 warning(s): the equation may be real, but parts of it are not yet supported

## Equation

```
x' = 0.1173 - 9.815 x + 9.851 y
y' = 27.12 x - 0.7343 y - 0.9799 x*z
z' = -2.653 z + 0.9943 x*y
```

## Findings

### WARNING `DERIVATIVE_ARTIFACT` (tier 0)
*where:* x'  
Coefficients that move by more than 30% when estimated derivative-free: the constant term +0.1173 -> +0.03256.

*fix:* Refit with smoothed or weak-form derivatives (e.g. PySINDy's SmoothedFiniteDifference) and report the derivative-free coefficients next to the discovered ones.

### WARNING `FORKING_PATHS` (tier 1)
*where:* analysis  
The equation's structure comes back under only 1 of 5 defensible analysis variants (3 distinct structures in total). Which law you report depends on analysis choices.

*fix:* Report the variants, or pick the analysis on held-out data before looking at the equation.

### WARNING `UNSTABLE_STRUCTURE` (tier 1)
*where:* y'  
Refitting on 30 resamples of your data changes this equation's terms. Terms absent from the equation that often appear: the constant term (97%) (share of refits). The structure depends on which part of the data the discovery saw.

*fix:* Treat these terms as tentative; check them on independent trajectories, or with more or less noisy data.

### WARNING `UNSTABLE_STRUCTURE` (tier 1)
*where:* z'  
Refitting on 30 resamples of your data changes this equation's terms. Terms absent from the equation that often appear: the constant term (90%), x (50%), y (53%) (share of refits). The structure depends on which part of the data the discovery saw.

*fix:* Treat these terms as tentative; check them on independent trajectories, or with more or less noisy data.

### WARNING `WEAK_TERM` (tier 0)
*where:* x'  
Dropping the constant term (F = 0.00631) barely changes the fit (F below 10): the data cannot tell the equation with these terms from the equation without them.

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
| variables | 3 |
| trajectories | 1 |
| samples | 5000 |
| active_terms | 8 |
| derivative_r2 | 0.2745 |
| integral_r2 | 0.9903 |
| library_relations | 0 |
| rollout_r2_short | 0.982 |
| structure_reproduced | 0 |
| structure_under_variants | 0.2 |
| distinct_structures | 3 |

## Notes

- held-out data: last 20% of each trajectory
- derivative-free reference: integral form over windows of 50 samples
- candidate library for alternatives: 10 terms
