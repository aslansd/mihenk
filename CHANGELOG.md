# Changelog

## 0.2.0 — 2026-11-11

### Added
- **`KINEMATIC_MISMATCH`** (tier 0). Declare derivative pairs such as
  `theta' = omega`, with `kinematics={"theta": "omega"}` in `scan` and
  `scan_pysindy`, `--kinematic theta=omega` on the command line, or
  `"kinematics"` in the bundle's meta. mihenk then:
  - fits the measured `omega` as `gain * omega(t + lag) + offset`,
    derivative-free, overall and per trajectory;
  - checks the discovered `theta'` against `omega`;
  - estimates how much of the extra `theta'` terms the lag predicts;
  - refits the other equations with `omega` corrected, so it can report when
    the sensor flaw changes the law.

  Without a declaration the check is *not assessed*, with suggested pairs.
  Nine new thresholds start with `kin_`. All of them are documented in
  `Thresholds` and can be overridden.
- **`synth.simulate(..., sensor=...)`** records a variable with a lag, gain
  and/or offset, so sensor artifacts can be studied with a known law.
- **A `pendulum` system in `synth`**, with coefficients like the Brunton-lab
  experiment and the pair `theta' = omega` declared.
- **`mihenk fit --fourier N`** and `fit_bundle(..., fourier=N)` add
  `FourierLibrary(n_frequencies=N)` to the polynomial library.
- **Four validation cases** (now 18): a clean pendulum (the kinematic check
  must stay quiet), and velocity lag, gain and offset.
- **`experiments/kinematic_synthetic/`.** Thresholds were calibrated on seeds
  0–9, and predictions for seeds 10–19 committed before running them:
  5 held, 1 partly held.
- **`experiments/pendulum/run_kinematics.py`**: the development case on the
  real pendulum. ω runs 0.85 ms late, and the lag explains 98% of the extra
  θ' terms.
- 10 new tests (35 in total).

### Changed
- The validation `Case` gains `sensor` and `forbid` (flags that must not
  fire).



First contact with experimental data, and the bug it exposed.

### Fixed
- **Rollout R² on segments where the system barely moves.** Short and full
  rollouts were scored against the held-out segment's own spread. On a
  recording that ends with the pendulum at rest, that gave R² = −2×10²⁰.
  Rollouts are now scored against the variance of the whole data set. No flag
  changed in any synthetic case: the 20-seed validation is identical to
  0.1.0. No flag changed in the pendulum analyses either. A regression test
  was added.

### Added
- **`experiments/pendulum/`: a pre-registered test on the Brunton lab's
  experimental single pendulum** (MIT-licensed data, attributed). It
  contains:
  - `PROTOCOL.md`, with three analyses and seven predictions, committed
    before any fit;
  - `prepare.py`, `run.py` and an exploratory `explore_lag.py`;
  - the derived bundles;
  - `RESULTS.md`, scoring every prediction: 5 held, 1 partly held, 1 failed.

  The tests found a 0.7–1.4 ms lag in the lab's filtered velocity signal, and
  three things mihenk misses. Those now lead the roadmap.
- `ROADMAP.md` reordered for 0.2 around what the experiment showed.

### Earlier, on GitHub only (since 0.1.0)

#### Added
- `INSTALL.md`: extras, versions, installing on an unreliable connection
  (including offline installs from downloaded wheels), and editable installs.
- `TESTING.md`: the test suite, the synthetic validation and its pass rules,
  and what to do when a validation run fails.
- `PUBLISHING.md`: release steps, GitHub-only changes, and safe re-uploads.
- `ROADMAP.md`: what 0.1 did and learned, and the plan for 0.2 (experimental
  data, a calibration/test split, `mihenk design`, `mihenk hybrid`).
- A line in the README linking these documents.
- Project URLs in `pyproject.toml`. They reach PyPI with the next release.

#### Fixed
- **Tests no longer skip the whole module when PySINDy is missing.** Only the
  6 tests that fit with PySINDy are skipped (`18 passed, 6 skipped`), so the
  core checks are tested on numpy, scipy and sympy alone. Before this, a fresh
  environment without PySINDy reported `1 skipped` and tested nothing.

## 0.1.0 — 2026-11-04

First release.

- Tier-0 checks from data and equation alone: `EQUIFINALITY`, `DERIVATIVE_ARTIFACT`, `WEAK_TERM`, `UNDERFIT`, `ROLLOUT_TRAP`.
- Tier-1 checks through PySINDy refits: `UNSTABLE_STRUCTURE` (bootstrap) and `FORKING_PATHS` (derivative method × threshold).
- Tier-2 ground-truth comparison for synthetic data.
- PySINDy 2.x adapter; equations from PySINDy models, text (including `model.print()` output), SymPy or JSON.
- CLI: `scan`, `fit`, `synth`, `validate` (with `--seeds`), `info`. HTML, JSON, Markdown and text reports.
- Optional daftar provenance (`--daftar`).
- 7 synthetic systems and 14 self-validation cases.
