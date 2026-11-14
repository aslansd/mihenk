# mihenk 0.3 checks: test results (seeds 10–19)

The predictions and thresholds were committed in [PROTOCOL.md](PROTOCOL.md)
(`61ca977`) before these seeds were run. The raw outputs are in
`results/calibration_seeds_0-9.txt` and `results/test_seeds_10-19.txt`, and
the script is `scripts/sweep.py`.

## Score

| | prediction | outcome |
|---|---|---|
| Q1 | correct pendulum laws pass the energy check, 0 of 60 flagged | **held** (0 of 60) |
| Q2 | an energy-injecting law is an error, 30 of 30 | **held** (30 of 30) |
| Q3 | over-damping is a warning, 30 of 30, with median ratio 2.8–3.2 | **held** (30 of 30; medians 2.99, 2.98, 2.96) |
| Q4 | the energy check does not see a velocity lag, no flag in ≥ 27 of 30 | **held** (0 of 30 flagged) |
| Q5 | conservation: correct laws 0 of 60, the damped law 30 of 30 | **held** (0 of 60; 30 of 30) |
| Q6 | long rollouts: the undamped law 30 of 30; clean fits ≥ 88 of 90 unflagged | **held** (30 of 30; 90 of 90) |
| Q7 | Lorenz classified chaotic in ≥ 8 of 10 per noise level; no drift warning in 30 | **failed**. Chaotic in 6, 7 and 6 of 10, and 2 drift warnings. |
| Q8 | `UNDERFIT` in ≥ 7 of 10 runs at noise 0; θ⁵ best in ≥ 8 of 10 | **held** (8 of 10; 10 of 10) |
| Q9 | relation with θ and sin θ reported, 30 of 30 | **held** (30 of 30) |
| Q10 | no regressions: 360/360 pass, and no correct-law run gains a flag | **failed, then fixed after the test.** 359/360 passed; no correct-law run gained a flag |

## Reading the results

**The energy check behaves as designed.** Correct laws, from the truth or from
SINDy, never tripped it, with rate ratios between 0.84 and 1.04. A law with
damping replaced by anti-damping was an error in every run, at about −0.43
of the data's rate (a gain). Triple damping was a warning in every run, at
about 3.0.

**It does not see a sensor lag, as predicted (Q4).** With ω recorded 0.7 ms
late, SINDy halves the damping coefficient. But its law and the lagged data
agree on the energy budget (ratio 0.84–1.02): the lag moves dissipation
between terms rather than removing it. That is `KINEMATIC_MISMATCH`'s case,
and the two checks complement each other.

**Conservation is tight.** The true Lotka–Volterra law drifts by about 10⁻⁷
of H's scale, and SINDy's clean fits by at most 0.7%. 5% damping drifts by
3–15%, a warning in every run.

**Long rollouts.** The undamped pendulum's full-length R² ranged from 0.30
down to −0.91, flagged in all 30 runs. No clean non-chaotic fit was flagged.

**Where Q7 failed.** Four of the ten test trajectories (seeds 11, 12, 15 and
18) spend all 10 time units spiralling around one of Lorenz's two fixed
points. On that time scale their twins really do stay together (6- to
229-fold growth), so "not chaotic" is the right call, and the prediction
assumed a more typical draw of seeds. Calibration had one such trajectory in
ten.

Two of those runs got the new drift warning:

| seed | noise | R² | twin growth | discovered law |
|---|---|---|---|---|
| 15 | 0 | 0.42 | 39× | correct |
| 11 | 1% | −3.26 | 229× | correct structure |

A spiral around an unstable fixed point eventually jumps to the other lobe,
and *when* it jumps is extremely sensitive to tiny coefficient errors, not
just to the starting state. Twin rollouts perturb only the start, so they
miss this kind of sensitivity. **Seed 15 is a false warning on a correct
law.** We did not change any threshold after seeing it. The fix belongs in
the check itself (for example, twins that also perturb the coefficients) and
is on the roadmap.

Two other Lorenz flags in the test runs were already there in 0.2.0, and are
unchanged behaviour:

- an escape error (seed 11, 0.5% noise);
- a spread/mean warning (seed 12, 0.5% noise).

**`UNDERFIT` two degrees out works without noise and has no power with it.**
At noise 0, θ⁵ was the best left-out term in 10 of 10 runs, and it was
flagged (as a warning) in 8. At 0.5% and 1% noise it was never flagged. A
term worth 0.1% of the signal sits below the noise, exactly as calibration
showed.

## Q10: the regression the test found, and its fix

We compared the 20-seed validation (18 cases, 360 runs) run by run against
mihenk 0.2.0.

**As committed, the new code failed one run.** On `lotka_volterra_one_orbit_noisy`,
seed 1, 0.2.0 had found the near-ellipse that a single noisy orbit lies on.
From it, 0.2.0 offered the *true* law (x' = 0.9995 x − 0.4999 x·y) as a
sparser, equally good alternative, an `EQUIFINALITY` error. The centred
relation search judges that ellipse against the variables' *variation*
rather than their raw size, which sits well above zero for x and y around 2.
So it missed it. The wrong law was still flagged, but only with warnings,
and the case failed. Two other runs also lost `EQUIFINALITY`, while picking
up other flags. No run where SINDy found the true law gained a flag.

**The fix, made after the test:** search both ways. Centred relations come
first, and an uncentred relation (0.2's method) is added unless its terms are
a subset of a centred relation's terms. That keeps "1 ≈ θ" out of the
small-swing pendulum's report and keeps the ellipse in.

**Re-checked after the fix** against 0.2.0, run by run:

| seeds | 0.2.0 | 0.3.0 | runs gaining a flag | gains on correct laws | flags lost | wrong laws unflagged |
|---|---|---|---|---|---|---|
| 0–19 (used above, no longer blind) | 360/360 | 360/360 | 32 | 0 | 0 | 0 of 171 |
| 20–29 (fresh) | 180/180 | 180/180 | 14 | 0 | 0 | 0 of 87 |

The physics sweep's test seeds give identical output before and after the
fix. All the gained flags are on wrong laws, mostly the new drift warning
(`duffing_small_library`, `van_der_pol_noise10`) and `UNDERFIT` on the
lagged pendulum.

**Score as committed: 8 held, 2 failed (Q7 and Q10).** Q10's failure was
fixed after the test, and the fix was checked on fresh seeds. That does not
turn it into a pass.
