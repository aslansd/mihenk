# Experimental pendulum: results

The protocol and predictions are in [PROTOCOL.md](PROTOCOL.md), committed
(`59fbb4b`) before any fit or scan. Everything below was run afterwards with
the analyses and thresholds fixed there.

- `results/logs/run_mihenk-0.1.0.txt`: the pre-registered run, with mihenk 0.1.0.
- `results/logs/run_mihenk-0.1.1.txt`: the same analyses after the one bug fix this
  test exposed (see C). The fix changes no flag in any analysis.
- `results/0.1-analyses/A`, `B`, `C`, `E1`: full reports (HTML to read, JSON with every number).

Reproduce:

```bash
git clone --depth 1 https://github.com/dynamicslab/MultiArm-Pendulum
python experiments/pendulum/scripts/prepare.py MultiArm-Pendulum   # optional: data/ holds the bundles
python experiments/pendulum/scripts/run.py
python experiments/pendulum/scripts/explore_lag.py
```

## Score

| | prediction | outcome |
|---|---|---|
| A1 | SINDy finds θ' = ω; ω' has sin θ within 10% of 64.22 and a negative ω term | **partly held**. sin θ = 63.82 (0.6% off). But θ' has two extra terms, and there is no damping term. |
| A2 | no `EQUIFINALITY` or `ROLLOUT_TRAP` error | **held** |
| A3 | extra ω' terms are flagged, or are stable across resamples (a candidate real effect) | **held** |
| A4 | the sin θ term is not flagged | **held** |
| B1 | `UNDERFIT` or `ROLLOUT_TRAP` on ω' (clean counts as a miss) | **failed**. B was flagged, but not for this. |
| C1 | `EQUIFINALITY` on ω' | **held** (error) |
| C2 | ω' not reported clean | **held** |

5 held, 1 partly held, 1 failed. Two things mihenk *missed* that no prediction
covered are listed under "What mihenk missed".

## A: informed library, large swings

```
SINDy:   theta' = 0.134 - 0.04273 theta + 0.9999 omega
         omega' = 0.4616 - 0.1941 theta + 63.82 sin(theta) - 0.2435 cos(theta) + 0.4428 sin(omega)
mihenk:  4 warnings (DERIVATIVE_ARTIFACT, WEAK_TERM, UNSTABLE_STRUCTURE, FORKING_PATHS), no errors
```

**The gravity term is right.** sin θ has coefficient 63.82, against the lab's
64.22. It is the same under every derivative method (63.92–63.96), its F
statistic for dropping it is above 10⁵, and it appears in 30 of 30 refits.
mihenk flags nothing about it.

**The extra θ' terms were a real flaw in the data, and mihenk pointed at them.**
θ' should be exactly ω. The extra `0.134 − 0.0427 θ` is the same under every
derivative method (0.141, −0.045), so it is not differentiation noise.
`UNSTABLE_STRUCTURE` and `FORKING_PATHS` showed that it trades places with
sin θ across resamples and analysis variants. Those are the signs of a term
that fits a real but unmodelled effect.

That effect is a **lagging velocity signal**. If the measured ω trails the
true velocity by D, then θ' = ω_meas − 64·D·(θ − π). This matches both
coefficients:

- −64·D = −0.0427, so D ≈ 0.67 ms;
- 0.0427·π = 0.134, the constant.

Measured directly (`scripts/explore_lag.py`), the lab's filtered ω lags dθ/dt by:

| segment | 1 | 2 | 3 | 4 | 5 | 6 |
|---|---|---|---|---|---|---|
| lag (ms) | 0.73 | 0.79 | 0.87 | 0.97 | 1.10 | 1.38 |

The lag grows as the swings shrink. That fits the lab's non-linear denoising
of ω (median filter plus wavelet thresholding) better than a fixed delay. We
do not claim the exact mechanism.

**The lag also explains the missing damping.** A lag D adds about +64·D·ω
(≈ +0.043·ω) of apparent *anti*-damping to ω'. That cancels most of the true
−0.067·ω, so the net damping is too small to survive the sparsity step.
SINDy put a +0.44·sin ω term in its place.

**Exploratory E1 (not pre-registered).** We refit analysis A with ω
differentiated from the angle itself, so it has no filter lag:

```
SINDy:   theta' = omega
         omega' = 0.1025 - 0.07709 theta - 0.05721 omega + 63.95 sin(theta) - 0.2234 cos(theta) - 0.05026 sin(omega)
```

θ' = ω exactly, and damping comes back with the right sign. The ω and sin ω
terms together give −0.107 at small speeds, against the lab's linear −0.067.
The lab fitted linear damping only, so the difference may be real non-linear
friction. mihenk now reports the leftover constant as a `DERIVATIVE_ARTIFACT`
**error**: it has the opposite sign without differentiation.

**cos θ** (−0.24 in A, −0.22 in E1) is in 30 of 30 refits and is stable
under every derivative method. It is what an angle offset of about 3.8 mrad
looks like: sin(θ + δ) ≈ sin θ + δ cos θ. Per A3, we record it as a candidate
real effect (encoder zero slightly off "straight down"), not a miss.

## B: polynomial library, large swings

```
SINDy:   omega' = -92.43 + 214.9 theta - 88.51 theta^2 + 9.38 theta^3
mihenk:  EQUIFINALITY error on theta' (the lag terms again); UNSTABLE_STRUCTURE on omega'
         (omega, omega*theta, omega*theta^2 appear in 63% of refits); FORKING_PATHS
```

**B1 failed.** We predicted that mihenk would notice that the library lacks
sin θ. It did not, for two reasons:

1. **Within the observed range, the cubic is a good law.** Around θ = π,
   sin θ = −(θ − π) + (θ − π)³/6 − …, so a cubic captures it to within the
   5th-order term. The short held-out rollouts have R² 0.9999. No check can
   prove the library wrong without data where the cubic and sin θ differ,
   such as full rotations.
2. **`UNDERFIT` looks only one degree beyond the library.** Here the first
   missing term is degree **5**: sin θ is odd about π, so a degree-4 term
   does not help. Searching two degrees, or by parity, would have found it.
   This is on the roadmap.

mihenk did not report B clean. Its findings concern the lag in θ', and the
damping trying to enter ω' (the 63% of refits that add ω terms).

## C: informed library, small swings

```
SINDy:   omega' = 23.15 + 3.933 theta + 15.36 omega + 68.86 sin(theta) + 35.36 cos(theta)
                  - 16.31 sin(omega) - 0.1821 cos(omega)
mihenk:  EQUIFINALITY errors on theta' and omega'; UNSTABLE_STRUCTURE error; DERIVATIVE_ARTIFACT,
         WEAK_TERM, FORKING_PATHS warnings
```

**C1 and C2 held.** With ±0.07 rad swings, the seven-term ω' is meaningless.
mihenk says so: it finds a sparser law that fits as well,
`omega' = 16.06 omega + 64.78 sin(theta) - 17.02 sin(omega)`. In that law,
sin θ is within 1% of the lab's value, and ω and sin ω are near-copies of
each other at these small speeds.

**This test exposed a mihenk bug (fixed in 0.1.1).** The held-out last 20%
of this recording is the pendulum **at rest**: θ is constant to the encoder's
resolution. In 0.1.0, the short-rollout R² was computed against that
segment's own spread, which was zero, and reported −2×10²⁰. Rollouts are now
scored against the variance of the whole data set. Every flag in A, B and C,
and the whole 20-seed synthetic validation, is unchanged by the fix.

That the pendulum comes to a complete stop is itself physical information.
Linear damping never stops a pendulum in finite time. Dry (Coulomb) friction
does, and the library has no `sign(ω)` term for it.

## What mihenk missed

1. **An energy-injecting law passed without a word.** A's +0.44·sin ω term
   adds energy at small speeds, which is impossible for a free pendulum, and
   the equation has no damping at all. mihenk flagged the term's coefficient
   (`DERIVATIVE_ARTIFACT`: it moves to +0.27 without differentiation), but
   no check asks whether the law dissipates energy. That needs a declared
   physical property (for example "this system loses energy"), which mihenk
   cannot infer.
2. **Long rollouts that drift were not flagged.** A's full-length rollout of
   segment 2 has R² 0.73. Under 0.1.0's per-segment scaling, segments 2 and 6
   scored 0.83 and 0.84. `ROLLOUT_TRAP` checks only short horizons
   and long-run spread and mean, a deliberately coarse test so that chaotic
   systems are not punished. For a system that is not chaotic, full-length R²
   is a fair test. Telling the two apart needs a chaos check, such as twin
   rollouts from nearby starts.

## What this changes

| finding | where it goes |
|---|---|
| R² on a resting segment | **fixed in 0.1.1** |
| velocity lag shows up as θ' ≠ ω | **done in 0.2**: `KINEMATIC_MISMATCH` (see below) |
| an energy-injecting law passes | roadmap: declared physical properties (dissipative, conserved quantities) |
| `UNDERFIT` misses odd functions | roadmap: search two degrees beyond the library, or by parity |
| drifting long rollouts | roadmap: a chaos-aware long-rollout check |
| relations reported as "1 ≈ θ" when θ barely moves around a large offset (C) | roadmap: centre columns before looking for relations |

No threshold was changed.

## mihenk 0.2: `KINEMATIC_MISMATCH` (development case, not a test)

This check was **designed after** these data showed the lag, so applying it
here only shows that it does what it was built for. Its test is the
pre-registered synthetic sweep in
[../kinematic_synthetic](../kinematic_synthetic/RESULTS.md), where it scored
5 held and 1 partly held. Reproduce with
`python experiments/pendulum/scripts/run_kinematics.py`. The reports are in
`results/0.2-kinematics/`.

| analysis | ω comes from | result with `theta = omega` declared |
|---|---|---|
| A | the lab's filtered velocity | **error**: ω runs late by 0.85 ms (F = 2.8×10³); the lag accounts for 98% of the extra θ' terms |
| B | the lab's filtered velocity | warning: the same lag and the same 98% |
| C | Savitzky–Golay derivative of θ | clean |
| E1 | Savitzky–Golay derivative of θ | clean |

C and E1 are the controls. Their ω is computed from θ itself, so the pair
holds by construction, and the check stays quiet.

The per-segment lags are fitted derivative-free, on window increments of θ.
They agree with `explore_lag.py`, which matched ω against dθ/dt instead:

| segment | 1 | 2 | 3 | 4 | 5 | 6 |
|---|---|---|---|---|---|---|
| `KINEMATIC_MISMATCH` (ms) | 0.75 | 0.80 | 0.85 | 0.95 | 1.10 | 1.35 |
| `explore_lag.py` (ms) | 0.73 | 0.79 | 0.87 | 0.97 | 1.10 | 1.38 |

**What the check cannot do here.** Correcting ω changes A's ω' constant by
36%. But A's ω' has no ω term, so the correction cannot show the damping
coming back. The check only re-weighs terms the equation already has. We also
rediscovered with SINDy on ω shifted by the fitted lag, both one global lag
and one lag per segment. A damping term appears (−0.054 and −0.052 ω), but θ'
still keeps extra terms. A constant shift cannot undo a lag that grows as the
swings shrink. Computing ω from θ (E1) is the better repair, and it is the
first one the finding suggests.

**A synthetic confirmation of the explanation.** `synth.simulate("pendulum",
sensor={"omega": {"delay": 7e-4}})` puts a 0.7 ms lag on a pendulum with known
law. SINDy then returns θ' = 0.9999 ω − 0.0448 sin θ (−64 × 0.0007 = −0.0448),
and the damping falls from −0.070 to −0.032. That is the experiment's pattern,
now with ground truth. There, correcting ω restores the damping to −0.070.

## Data license

The data come from <https://github.com/dynamicslab/MultiArm-Pendulum>
(MIT license, © 2022 dynamicslab). The license text is in
[data/LICENSE.txt](data/LICENSE.txt). The bundles in `data/` are derived from
`SinglePendulumDataForParameterEstimation.mat` and
`SingleDataFreeSwing_1_Dt_0_0001.mat` by `scripts/prepare.py`.
