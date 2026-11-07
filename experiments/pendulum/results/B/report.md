# mihenk report: pendulum B: polynomial library, large swings

*mihenk 0.1.1*

1 error(s), 2 warning(s): do not interpret this equation as the law yet

## Equation

```
theta' = 0.134 - 0.04273 theta + 0.9999 omega
omega' = -92.43 + 214.9 theta - 88.51 theta^2 + 9.38 theta^3
```

## Findings

### ERROR `EQUIFINALITY` (tier 0)
*where:* theta'  
A different equation of the same size fits your data exactly as well: theta' = 0.07188 - 0.007016 theta^2 + 0.9999 omega. On the observed states some candidate terms are combinations of others, so the data cannot choose between these laws.

*fix:* Collect data that leaves the current set of states (other initial conditions, transients, perturbations) so the alternatives make different predictions.

### WARNING `FORKING_PATHS` (tier 1)
*where:* analysis  
The equation's structure comes back under only 2 of 5 defensible analysis variants (2 distinct structures in total). Which law you report depends on analysis choices.

*fix:* Report the variants, or pick the analysis on held-out data before looking at the equation.

### WARNING `UNSTABLE_STRUCTURE` (tier 1)
*where:* omega'  
Refitting on 30 resamples of your data changes this equation's terms. Terms absent from the equation that often appear: omega (63%), omega*theta (63%), omega*theta^2 (63%) (share of refits). The structure depends on which part of the data the discovery saw.

*fix:* Treat these terms as tentative; check them on independent trajectories, or with more or less noisy data.

## Checks

| check | tier | status | note |
|---|---|---|---|
| EQUIFINALITY | 0 | flagged |  |
| DERIVATIVE_ARTIFACT | 0 | clean |  |
| WEAK_TERM | 0 | clean |  |
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
| active_terms | 7 |
| derivative_r2 | 0.985 |
| integral_r2 | 0.9997 |
| library_relations | 1 |
| rollout_r2_short | 0.9999 |
| structure_reproduced | 0.3667 |
| structure_under_variants | 0.4 |
| distinct_structures | 2 |

## Notes

- held-out data: held-out trajectories (val)
- derivative-free reference: integral form over windows of 50 samples
- candidate library for alternatives: 10 terms
