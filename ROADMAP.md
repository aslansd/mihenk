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
  all 91 wrong laws were flagged, no correct law got an error, and 3 of 175
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

## Next: 0.2 — first contact with experimental data

Every 0.1 threshold was set on toy systems. The next release should change
that before it adds many features.

- **Run mihenk on public experimental data with a known law.** Good
  candidates have well-understood physics and real noise, such as a damped
  pendulum, coupled oscillators or a chemical oscillator. Report every miss,
  as ayna's VALIDATION.md does.
- **A calibration/test split for validation.** Fix the thresholds on one set
  of seeds and systems. Record that split in a commit before running the
  other set, and quote only the other set. 0.1 thresholds were adjusted while
  looking at seeds 0–9.
- **More systems in `synth`.** A Rössler attractor, a forced Duffing
  oscillator, a pendulum (with `sin` terms, outside polynomial libraries) and
  a higher-dimensional system, so the detectors are not tuned to
  two-variable polynomials.
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
  - **Non-polynomial libraries.** Fourier and custom libraries parse today,
    but are not validated.
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

## Known limits in 0.1

- One mild-noise Lorenz trajectory in 20 lies close to a quadric surface. A
  constant term then appears in some refits, and mihenk warns. The law is
  right. The warning is correct about that recording, but it reads like a
  false alarm.
- `ROLLOUT_TRAP` compares long-run spread and mean for chaotic systems. That
  is a coarse check, and a model with the wrong attractor shape but the right
  moments passes it.
