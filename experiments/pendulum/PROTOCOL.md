# Experimental pendulum: protocol and predictions

**This file was committed before any SINDy fit or mihenk scan was run on these
data.** Results go in `RESULTS.md`, written afterwards, and every prediction
below is scored there as it stands, whether it held or not.

## Data

The experimental multi-arm pendulum from the Brunton lab, the same group that
develops PySINDy:

- Kaheman, Fasel, Bramburger, Strom, Kutz & Brunton (2022), *The Experimental
  Multi-Arm Pendulum on a Cart: A Benchmark System for Chaos, Learning, and
  Control*.
- Repository: <https://github.com/dynamicslab/MultiArm-Pendulum> (MIT license).

Only the **single pendulum** (cart fixed) is used here. Its law is known from
mechanics, and it is written out in the repository's
`ParameterEstimation/SinglePendulum/SinglePendulumODE.m`:

```
theta' = omega
omega' = (a1 g m1 / J) sin(theta) - (k1 / J) omega,      J = m1 a1^2 + I1
```

θ = π is hanging straight down. The lab's own parameter estimate is in
`EstimatedParameters_SinglePendulum.mat`
(`a1, m1, I1, k1, g = 0.14775, 0.14758, 1.0912e-4, 2.2394e-4, 9.8100`). It gives
the **reference law**:

```
theta' = omega
omega' = 64.22 sin(theta) - 0.0672 omega
```

The coefficients are the lab's estimates, fitted by simulation to the same
recordings, not exact truths. The *structure* (`sin θ` and linear damping) is the
physics. Real friction in bearings is not purely linear, so a faithful
discovery may legitimately find more than these two terms.

Two data sets, both built by `prepare.py`:

| bundle | source file | states | amplitude | split |
|---|---|---|---|---|
| `pendulum_large.npz` | `SinglePendulumDataForParameterEstimation.mat` | θ and ω, both measured by the encoder (ω denoised by the lab) | ±1.6 rad decaying to ±0.4 | the lab's own split: 4 fitting segments, 2 held-out (`val`) |
| `pendulum_small.npz` | `Datas/SinglePendulum/SingleDataFreeSwing_1_Dt_0_0001.mat`, 60–83 s | θ measured; ω computed from θ (Savitzky–Golay, 101 samples at 10 kHz, cubic) and then both downsampled to 1 kHz | ±0.07 rad around π | mihenk default: last 20% held out |

## Analyses (fixed in advance)

All fits use PySINDy 2.1.0 with STLSQ at **threshold 0.02**. That is below the
reference damping of 0.067, so a correct fit is allowed to keep the damping
term. The differentiation is PySINDy's default (`FiniteDifference`).

| analysis | data | library |
|---|---|---|
| **A** informed | large | `PolynomialLibrary(degree=1) + FourierLibrary(n_frequencies=1)`: `1, θ, ω, sin θ, cos θ, sin ω, cos ω` |
| **B** naive | large | `PolynomialLibrary(degree=3)`: no trigonometric terms |
| **C** small swings | small | same library as A |

mihenk 0.1.0 runs with all default thresholds and seed 0, through
`mihenk.scan_pysindy` (tiers 0 and 1). **No threshold will be changed after
seeing these results.** If a threshold looks wrong, that goes in `RESULTS.md`
and `ROADMAP.md`, not into the code.

## Predictions

**A (informed library, large swings).**

- A1. SINDy finds θ' = ω, with no other term above the threshold. For ω', it
  finds `sin θ` with a coefficient within 10% of 64.22, and an `ω` term with a
  negative coefficient.
- A2. mihenk reports no `EQUIFINALITY` **error** and no `ROLLOUT_TRAP`
  **error**. Large swings separate `sin θ` from `θ`, and the law reproduces the
  held-out segments.
- A3. If ω' contains terms beyond `sin θ` and `ω`, each of them is flagged by at
  least one of `WEAK_TERM`, `DERIVATIVE_ARTIFACT` or `UNSTABLE_STRUCTURE`, or
  is consistent across resamples. In that last case it is reported as a
  candidate real effect (for example non-linear friction), not as a miss.
- A4. The `sin θ` term itself is not flagged by any check.

**B (polynomial library, large swings).** The true law is not in the library.

- B1. mihenk flags the ω' equation with `UNDERFIT`, or `ROLLOUT_TRAP`, or both.
  A clean report for B counts as a **miss**.

**C (informed library, small swings).** Near the bottom, sin θ ≈ π − θ, so
`sin θ`, `θ` and the constant are nearly interchangeable.

- C1. mihenk reports `EQUIFINALITY` (error or warning) on ω'. No
  `EQUIFINALITY` counts as a **miss**.
- C2. Whatever SINDy returns for ω', mihenk does not report it clean.

## Scoring

Each prediction is marked held, partly held, or failed in `RESULTS.md`. A
failure means one of two things: either mihenk missed something (a detector
problem, for the roadmap), or the prediction was wrong about the physics or
about SINDy. `RESULTS.md` says which, with the evidence.
