# Roadmap

mihenk's question stays the same: *is the equation you discovered the law of
your system, or an artifact?* Each release should make the answer more
trustworthy, not just longer.

## Done

### 0.1 — the touchstone

- **Tier 0, data + equation.** These checks never call the discovery tool.
  They refit only the coefficients of term sets the equation already has,
  against a derivative-free integral form:
  - `EQUIFINALITY`, with the sparsest equally good law and one-term swaps;
  - `DERIVATIVE_ARTIFACT`;
  - `WEAK_TERM`;
  - `UNDERFIT`, searching the library and one polynomial degree beyond it;
  - `ROLLOUT_TRAP`.
- **Tier 1, refits through PySINDy.**
  - `UNSTABLE_STRUCTURE`, using a block bootstrap;
  - `FORKING_PATHS`, over the derivative method (finite vs smoothed
    differences) and the sparsity threshold (×0.5, ×1, ×1.5).
- **Tier 2, ground truth** for synthetic data.
- **Inputs:**
  - a PySINDy 2.x adapter;
  - equations as text (including `model.print()` output), SymPy or JSON;
  - `.npz` and `.csv` data.
- **Outputs and tools:**
  - CLI (`scan`, `fit`, `synth`, `validate`, `info`) that exits with status 1
    on errors;
  - HTML, JSON, Markdown and text reports;
  - optional daftar provenance.
- **Validation.** 7 synthetic systems and 14 validation cases. Over 20 seeds,
  all 91 wrong laws were flagged, no correct law got an error, and 3 of 155
  correct-law runs got a warning.

### What 0.1 development taught

Three detector designs failed on clean data and were replaced before release:

- **Trapezoid integration** left an error in the derivative-free reference.
  Stray cubic terms could absorb it, so `UNDERFIT` fired on clean Lorenz. The
  fix: Simpson integration, and a left-out term must also explain at least 1%
  of the signal, not just of an already tiny residual.
- **Random 60% windows** as the bootstrap asked "would a shorter recording
  give the same law?". A Lorenz trajectory that stays on one lobe fails that
  question legitimately. The fix: a block bootstrap, which keeps the amount
  and spread of data the same.
- **STLSQ leftovers.** PySINDy's STLSQ can stop at `max_iter` with
  coefficients far below its own threshold still non-zero. In refits these
  looked like new terms. They are now treated as zero in refits, never in
  your own equation.

A near relation among library terms (one Lotka–Volterra orbit lies within
0.14% of an ellipse) now warns only when it actually leaves the equation's
coefficients loosely determined.

### 0.1.1 — first contact with experimental data

- **A pre-registered test on a real pendulum**
  ([experiments/pendulum](../experiments/pendulum/RESULTS.md)). It uses the
  Brunton lab's experimental single pendulum. Three analyses were fixed, with
  seven predictions, and committed before any fit:
  - **Results:** 5 held, 1 partly held, 1 failed.
  - **The gravity term** came out within 0.6% of the lab's estimate, and
    mihenk did not flag it.
  - **Small swings:** mihenk reported `EQUIFINALITY` errors, as predicted.
  - **A real sensor problem.** mihenk's warnings on θ' led to a 0.7–1.4 ms
    lag in the lab's filtered velocity signal. That lag hid the damping term.
- **Fixed a bug the test exposed.** Rollout R² was scored against a held-out
  segment's own spread, which is meaningless when the system is at rest.
- **What it missed:**
  - a law that injects energy;
  - a long rollout that drifts;
  - a missing degree-5 term (sin θ around π).

  These lead the 0.2 plan below.

### 0.2 — `KINEMATIC_MISMATCH`

- **A new tier-0 check for declared derivative pairs** (`theta' = omega`),
  given with `kinematics=`, `--kinematic theta=omega`, or the bundle's meta.
  It does four things:
  - fits `p' = gain * v(t + lag) + offset` derivative-free, overall and per
    trajectory;
  - compares the discovered `p'` with `v`;
  - tests whether the lag predicts the extra terms (lag times the discovered
    `v'` equation);
  - refits the other equations with `v` corrected, so it can show that a
    sensor flaw changes the law.

  Undeclared pairs are suggested.
- **Calibrated on seeds 0–9, tested on seeds 10–19,** with predictions
  committed in between
  ([experiments/kinematic_synthetic](../experiments/kinematic_synthetic/RESULTS.md)):
  5 held, 1 partly held. There were no false alarms in 80 runs. A 0.7 ms lag
  is an error up to 0.5% noise. At 1% position noise only half of the lags
  are caught.
- **On the experimental pendulum** (a development case, since the check was
  designed after seeing it): ω runs 0.85 ms late, and the lag accounts for
  98% of SINDy's extra θ' terms. The per-segment lags match last release's
  independent estimate to within 0.03 ms.
- **`synth`** gains a `pendulum` system and a `sensor=` option (lag, gain
  and offset for any variable). `mihenk fit` gains `--fourier N`.

### 0.3 — physics, honest long rollouts, a wider `UNDERFIT`, readable relations

- **`PHYSICS_VIOLATION`** (tier 0). Declare `dissipates=H` or `conserves=H`,
  with `properties=` or `--dissipates/--conserves`. mihenk compares H's
  derivative-free rate of change in the data with the equation's
  `grad H . f` along the data, and follows H along the equation's own
  simulations. Gaining energy where the data lose it is an error; losing it
  at under half or over twice the data's rate is a warning.
- **`ROLLOUT_TRAP` uses twin rollouts.** Each full simulation gets a twin
  started 10⁻⁶ of the range away. If the gap grows less than 1000-fold, the
  equation is not chaotic on this time scale, and a full-length R² below
  0.8 is reported as drift.
- **`UNDERFIT` searches two polynomial degrees beyond the library.** It warns
  from 0.1% of the signal, and an error still needs 1%.
- **Relation search, centred and uncentred.** The relations table reports
  θ + sin θ ≈ π instead of "1 ≈ θ".
- **`synth`** gains `pendulum_small` and `pendulum_dry`.
- **Tested** with predictions committed before the test seeds
  ([experiments/physics_synthetic](../experiments/physics_synthetic/RESULTS.md)):
  8 of 10 held. The two that failed:
  - **Lorenz:** trajectories that stay on one lobe are non-chaotic on this
    time scale, and one of them got a false drift warning.
  - **A regression from the centred relation search:** one wrong law lost
    its `EQUIFINALITY` error. It was fixed after the test by combining the
    two searches, and the fix was checked against 0.2.0 on fresh seeds.
- **On the experimental pendulum** (a development case): the drifting
  rollouts are now flagged, and the relations read correctly. The energy
  check shows that A's law is not a net energy injector, which corrects the
  0.1.1 write-up. The degree-5 term stays below the noise, as the synthetic
  test predicted.

## Next: 0.4

Ordered by what the 0.3 tests showed, then by what was already planned.

- **Twins that also perturb the coefficients.** A twin started a hair away
  cannot see sensitivity to the *law*. A Lorenz spiral that will jump lobes
  at a time set by tiny coefficient errors looks non-chaotic and then
  "drifts". Perturbing the coefficients within their uncertainty would tell
  the two apart. This is the cause of 0.3's one false warning.
- **Local energy injection.** `PHYSICS_VIOLATION` judges dissipation on
  average. A law can gain energy in one region of state space and lose more
  elsewhere. A local test has to separate a real gain from H's conservative
  part being slightly off. Binning by speed alone does not (pendulum A).
- **Held-out segments without motion.** Warn when held-out data barely move
  (the pendulum at rest), since they test nothing.
- **A second experimental system** with a different structure, such as a
  chemical oscillator or coupled oscillators, pre-registered the same way.
- **A calibration/test split for the 0.1 checks.** 0.2 and 0.3 did this for
  their own checks. The 0.1 checks still have thresholds adjusted while
  looking at seeds 0–9 of the same suite they are scored on.
- **More systems in `synth`:** a Rössler attractor, a forced Duffing
  oscillator, an angle offset for the `pendulum`, and a higher-dimensional
  system.
- **`mihenk design`.** Before an experiment, simulate a known law at your
  planned sampling rate, duration, noise level and number of trajectories,
  and ask whether the discovery can recover it at all. This is the idea of
  `miras identify`.
- **`mihenk hybrid`.** Put a known law into your own data's noise and
  sampling, then score how well the flags match the truth. This is the idea
  of `ayar hybrid`.

## Later

- **More of PySINDy.**
  - **Control inputs (`u`).** Rejected with a clear error in 0.1.
  - **Other optimizers.** SR3, `EnsembleOptimizer` and MIOSR need tests.
    Their forking-path variants need their own parameters, not STLSQ's
    `threshold`.
  - **Weak-form and PDE libraries.** `WeakPDELibrary` and `PDELibrary` are
    rejected in 0.1. PDEs need space as well as time in the bundle, and an
    integral form over both.
  - **Non-polynomial libraries.** Fourier libraries work: the pendulum test
    used `PolynomialLibrary + FourierLibrary`. Custom libraries parse, but are
    not validated.
- **`DIMENSIONALLY_INCONSISTENT`.** Optional units per variable, checked
  term by term.
- **`ILL_CONDITIONED`.** Structural stability of the formula under small
  perturbations, after the Effective Information Criterion (ICML 2026).
- **A threshold-path view for `FORKING_PATHS`.** Instead of three fixed
  factors, sweep the sparsity threshold and report how wide the plateau is on
  which the reported structure holds.
- **A PySR reader.** Read PySR's saved hall-of-fame CSV without importing PySR
  or Julia, so mihenk checks symbolic-regression output too. Installing Julia
  can be hard on restricted networks, and the checker should not need it.
- **Exact F tests.** The integral form's windows overlap by half, so 0.1's F
  statistics are screening numbers. Non-overlapping windows or a corrected
  variance would make them proper tests.

## Known limits in 0.3

- One mild-noise Lorenz trajectory in 20 lies close to a quadric surface. A
  constant term then appears in some refits, and mihenk warns. The law is
  right. The warning is correct about that recording, but it reads like a
  false alarm.
- `ROLLOUT_TRAP` compares long-run spread and mean for chaotic systems. That
  is a coarse check, and a model with the wrong attractor shape but the right
  moments passes it.
- `ROLLOUT_TRAP`'s twins perturb only the start, not the law (see 0.4). A
  single Lorenz trajectory can also spend most of 10 time units on one lobe,
  and then the long-run spread-and-mean comparison warns although the law is
  right.
- `PHYSICS_VIOLATION` checks only what you declare, and judges dissipation on
  average over the data, not locally.
- `UNDERFIT`'s two-degree search has no power once a missing term is below
  the noise. A degree-5 term worth 0.1% of the signal was not found at 0.5%
  noise.
- A recording that ends at rest is used as held-out data, although it tests
  nothing. mihenk does not yet warn about held-out segments without motion.
- `KINEMATIC_MISMATCH` corrects `v` with one lag, gain and offset per pair.
  A lag that changes within a recording (the pendulum's grows from 0.75 to
  1.35 ms as swings shrink) is reported per trajectory, but it is not
  corrected per trajectory. The coefficient-change estimate also cannot add
  terms the equation lacks, such as the pendulum's missing damping.
- `KINEMATIC_MISMATCH` searches lags up to 10 samples either way
  (`kin_max_lag`). A longer lag is reported as at the search limit.
