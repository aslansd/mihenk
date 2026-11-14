# Changelog

## 0.3.0 — 2026-11-14

### Added
- **`PHYSICS_VIOLATION`** (tier 0). Declare `{"dissipates": "H"}` or
  `{"conserves": "H"}`, with `properties=` in `scan` and `scan_pysindy`,
  `--dissipates` / `--conserves` on the command line, or `"properties"` in
  the bundle's meta. mihenk compares H's derivative-free rate of change in
  the data with the equation's `grad H . f` along the data, and follows H
  along full simulations.
  - Dissipation: gaining H is an error; losing it at under half or over
    twice the data's rate is a warning.
  - Conservation: a rate or drift above tolerance is flagged.

  New thresholds: `phys_inject`, `phys_ratio_low`, `phys_ratio_high`,
  `phys_data_min`, `phys_rollout_gain`, `phys_conserve_rel` and
  `phys_conserve_drift`.
- `synth`: `pendulum_small` (swings of 0.05–0.1 rad) and `pendulum_dry`
  (smoothed Coulomb friction).
- `experiments/physics_synthetic/`: thresholds calibrated on seeds 0–9,
  predictions for seeds 10–19 committed before running them. 8 held,
  2 failed (one, a regression, fixed afterwards).
- `experiments/pendulum/scripts/run_physics.py`: the 0.3 development case on
  the real pendulum.
- 9 new tests (45 in total).

### Changed
- **`ROLLOUT_TRAP` uses twin rollouts.** Each full simulation gets a twin
  started `twin_eps` (10⁻⁶) of the range away. If their gap grows less than
  `twin_growth` (1000) times, the equation is not chaotic on this time
  scale, and a full-length R² below `rollout_long_r2` (0.8) is a drift
  warning. The rollouts table gains "regime" and "twin gap growth" columns.
- **`UNDERFIT` searches `underfit_degrees` (2) polynomial degrees beyond the
  library, not 1.** A left-out term now needs 0.1% of the signal for a
  warning (`underfit_share`, previously 1%). An error still needs 1%
  (`underfit_share_error`).
- **Relations are searched twice: centred, then uncentred.** The centred
  search reports a variable sitting on an offset correctly (θ + sin θ ≈ π
  instead of "1 ≈ θ"). The uncentred search keeps 0.2's sensitivity.
  - The pre-registered test found that centred-only lost one wrong law's
    `EQUIFINALITY` error.
  - The combined search fixes that, and was checked on seeds 0–19 and
    fresh seeds 20–29.

### Corrected
- **A counting error since 0.1.0.** The README, ROADMAP and TESTING said the
  synthetic validation has "9" non-pendulum clean cases with "175"
  correct-law runs. It has 8 cases and 155 runs. The outcomes were as
  reported: no errors, 3 warnings.
- `experiments/pendulum/RESULTS.md` said analysis A's law "injects energy".
  Its +0.44·sin ω term does, taken alone, at |ω| < π. But the whole law
  loses energy along the data, at 74% of the data's rate. A dated correction
  was added; the original text is kept.

## Repository tidy-up (2026-11-13, shipped with 0.3.0)

The package code changed in one place (`Report.save`).

### Changed: a tidier repository
- `INSTALL.md`, `TESTING.md`, `PUBLISHING.md` and `ROADMAP.md` moved to
  `docs/`. `README.md`, `CHANGELOG.md` and `LICENSE` stay at the root.
- Each experiment has the same shape: `PROTOCOL.md` and `RESULTS.md`, plus
  `data/`, `scripts/` and `results/`. The pendulum results are split by
  release (`0.1-analyses/`, `0.2-kinematics/`, `logs/`). The protocols were
  moved, not edited.
- `experiments/README.md` indexes the experiments and their layout.
- Experiment reports keep `report.html` and `report.json`. The `.md` and
  `.txt` copies were duplicates and are gone.
- The 0.2 kinematic reports were regenerated: they were labelled 0.1.1
  because they ran just before the version bump. Their findings are
  unchanged.
- Scratch files at the root (`lorenz.npz`, `model.json`, `report.html`,
  `report/`, `mihenk_report/`) were removed, and `.gitignore` now keeps them
  out.
- `dist/` is no longer tracked. Built files live on PyPI, and
  `docs/PUBLISHING.md` says so.

### Added
- `Report.save(out_dir, formats=...)` writes only the chosen formats. The
  default is still all four.

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
