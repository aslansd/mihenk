# mihenk 0.3 checks: calibration and pre-registered test

mihenk 0.3 adds or changes four things, each answering a miss in the
experimental pendulum ([../pendulum/RESULTS.md](../pendulum/RESULTS.md)):

| change | the miss it answers |
|---|---|
| `PHYSICS_VIOLATION`: declared `dissipates=H` or `conserves=H` | an equation that injected energy passed |
| `ROLLOUT_TRAP`: twin rollouts decide whether full-length R² is fair | a long rollout with R² 0.73 passed |
| `UNDERFIT`: searches two polynomial degrees beyond the library (a warning from 0.1% of the signal, an error from 1% as before) | sin θ's missing degree-5 term was invisible |
| centred relation search | "1 ≈ θ" was reported instead of "sin θ ≈ π − θ" |

Because the pendulum motivated all four, **it is a development case, not a
test.** The test is this synthetic sweep. Thresholds were chosen on **seeds
0–9** (`results/calibration_seeds_0-9.txt`). The predictions below were
committed before **seeds 10–19** were run.

## Setup

`scripts/sweep.py` runs 15 cases at noise levels 0, 0.5% and 1% of each
variable's standard deviation. Each row counts the outcome of the case's
*target* check: error, warning or clean.

| case | system | equation scanned | target |
|---|---|---|---|
| P1 | pendulum | the true law | `PHYSICS_VIOLATION`, dissipates H = ½ω² − 64 cos θ |
| P2 | pendulum | PySINDy (poly 1 + Fourier 1) | same |
| P3 | pendulum | the true law with damping replaced by +0.03 ω (injects energy) | same |
| P4 | pendulum, ω recorded 0.7 ms late | PySINDy | same |
| P5 | pendulum | the true law with three times the damping | same |
| C1 | Lotka–Volterra | the true law | `PHYSICS_VIOLATION`, conserves H = ½x − ln x + ½y − ln y |
| C2 | Lotka–Volterra | PySINDy (poly 2) | same |
| C3 | Lotka–Volterra | the true law with 5% damping on x | same |
| L1 | pendulum | the true law without damping | `ROLLOUT_TRAP` |
| L2–L4 | pendulum, Van der Pol, Duffing | PySINDy | `ROLLOUT_TRAP` |
| L5 | Lorenz | PySINDy | `ROLLOUT_TRAP`, plus how the twins classify it |
| U1 | pendulum (swings of 1–2 rad) | PySINDy with a cubic polynomial library | `UNDERFIT` |
| R1 | pendulum, swings of 0.05–0.1 rad recorded about θ = π | PySINDy | relations table |

## Thresholds, fixed on seeds 0–9

These are the defaults in `mihenk/scan.py`:

| threshold | value | why |
|---|---|---|
| `phys_inject` | 0.1 | An equation gaining H at 10% of the data's loss rate is an error. The injecting law sat at −0.43 of the data's rate (a gain), and correct laws at +0.91 to +1.09. |
| `phys_ratio_low` / `phys_ratio_high` | 0.5 / 2.0 | A warning when the equation loses H at under half, or over twice, the data's rate. Correct laws stayed within 0.91–1.09, and the over-damped law sat at 3.0. |
| `phys_rollout_gain` | 0.1 | A simulation raising H by 10% of H's scale is an error. |
| `phys_conserve_rel` / `phys_conserve_drift` | 0.01 / 0.05 | Conservation tolerances. Correct laws drifted at most 0.8% of H's scale, and the 5%-damped law 3.7–11%. |
| `rollout_long_r2` | 0.8 | Full-length R² below this, for a non-chaotic equation, is drift. Clean non-chaotic fits had R² ≥ 0.91, and the undamped law 0.06 down to −0.88. |
| `twin_growth` | 1000 | Twins whose gap grows more than 1000-fold mark a chaotic equation. Clean non-chaotic fits grew it at most 18.7-fold, a wrong noisy law 108-fold, and Lorenz typically 4×10⁴–10⁶-fold. The exception is one trajectory that stays on a single lobe (5.5-fold, R² 1.0), which is genuinely non-chaotic on this time scale. |
| `underfit_degrees` | 2 | Search two degrees beyond the library. |
| `underfit_share` / `underfit_share_error` | 0.001 / 0.01 | A warning from 0.1% of the signal, an error from 1% (the 0.1 value). No clean validation case had a left-out term passing the gain and F gates at any size. |

**A known gap, recorded before testing.** The 0.1 long-run spread-and-mean
comparison still warns on some single Lorenz trajectories, which can spend
more of 10 time units on one lobe. That is unchanged behaviour, not part of
this test.

## Predictions for seeds 10–19

The code and thresholds are frozen at the commit that adds this file.

- **Q1. Correct laws pass the energy check.** No `PHYSICS_VIOLATION` in any of
  the 60 runs of P1 and P2.
- **Q2. An energy-injecting law is an error** in 30 of 30 runs (P3).
- **Q3. Over-damping is a warning** in 30 of 30 runs (P5), with a median
  rate ratio between 2.8 and 3.2.
- **Q4. The energy check does not see a sensor lag** (P4). There is no flag
  in at least 27 of 30 runs. The lagged data and SINDy's fit to them agree on
  the energy budget; the lag biases the damping coefficient, which is
  `KINEMATIC_MISMATCH`'s job.
- **Q5. Conservation.**
  - No flag in any of the 60 runs of C1 and C2.
  - Flagged in 30 of 30 runs of C3.
- **Q6. Long rollouts.**
  - `ROLLOUT_TRAP` fires in 30 of 30 runs of L1.
  - It does not fire in at least 88 of the 90 runs of L2–L4.
- **Q7. Lorenz** (L5).
  - Classified chaotic in at least 8 of 10 runs at each noise level.
  - The new drift warning ("drifts away") fires in none of the 30 runs.
- **Q8. Degree-5 term** (U1).
  - `UNDERFIT` fires in at least 7 of 10 runs at noise 0, and θ⁵ is the best
    left-out term in at least 8 of 10.
  - No claim with noise. Calibration found no power at 0.5% and 1%.
- **Q9. Relations** (R1). A relation containing both θ and sin θ is reported
  in 30 of 30 runs.
- **Q10. No regressions.** In `mihenk validate --seeds 20`, every case passes
  (360 of 360), and no run where SINDy found the true law gains a flag it did
  not have in mihenk 0.2.0.

Results go in `RESULTS.md`, with every prediction scored as written here.
