# Testing mihenk

There are two layers. Each catches a different kind of mistake.

| layer | what it is | needs | catches |
|---|---|---|---|
| unit and integration tests | `tests/test_mihenk.py` | mihenk + pytest (+ PySINDy for 6 tests) | broken code, broken contracts, broken CLI |
| synthetic validation | `mihenk validate` | mihenk + PySINDy | detectors that stopped finding the failures they name, or started flagging correct laws |
| experimental data | `experiments/pendulum/` | mihenk + PySINDy | what toy systems cannot show: sensor artifacts, unmodelled friction, data that ends at rest |

The thresholds are calibrated on toy systems. The experimental test is scored
against predictions committed in advance, and no threshold was changed after
it.

## 1. The test suite

```bash
pip install -e ".[test]"
python -m pytest -q
```

35 tests, under a minute on a laptop CPU.

- With PySINDy installed, all 35 must pass.
- Without PySINDy, 29 pass and 6 are skipped with the reason `needs PySINDy`.
  That is expected. The core checks need only numpy, scipy and sympy, and they
  are tested on their own.
- If you see `1 skipped` and nothing else, you have a test file older than
  this one, which skipped the whole module when PySINDy was missing. Update
  `tests/test_mihenk.py`.

What the suite pins down:

- **Model contract.**
  - Text equations in every accepted form: `x' =`, `dx/dt =`, `d(x)/dt =`,
    `x_dot =`, `(x)' =`, and PySINDy's own `model.print()` output.
  - Implicit products such as `x z`, and powers written `x^2`.
  - Multi-letter names are not split into letters (`theta`, `omega`), and a
    name starting with `d` (`delta`) is not read as a derivative.
  - Unknown symbols are rejected.
  - JSON round-trip.
- **Bundle contract.** `.npz` round-trip with several trajectories and the
  true law in `meta`. Non-increasing times are rejected.
- **PySINDy adapter.** A model read through mihenk predicts exactly what
  PySINDy predicts.
- **Detectors.**
  - Clean Lorenz is clean, and ground truth agrees.
  - Limit-cycle data gives an `EQUIFINALITY` error that names `x' = -y`.
  - A library without `x^3` gives `UNDERFIT` on Duffing's `y'`.
  - A hand-edited sign gives `ROLLOUT_TRAP`.
  - A spurious constant from noisy Lorenz is flagged.
- **Tiers.** Without a refit callable, the tier-1 checks are *not assessed*,
  never clean.
- **Thresholds.** Overrides work, and an unknown threshold name is an error.
- **Reports.** `report.json` carries the version, and `report.html` is a full
  page with the findings.
- **Kinematic pairs.**
  - Declarations parse, and bad ones are rejected.
  - A synthetic lag, gain and offset are each recovered to within 2%.
  - An exact sensor is clean, and no declaration means *not assessed*, with
    a suggested pair.
  - A lag is flagged as an error. The lag explains the extra θ' term, and
    correcting ω restores the damping to −0.07.
  - The `kinematics=` argument overrides the bundle's meta.
  - The CLI `--kinematic` flag works.
- **Rollouts on resting data.** A held-out segment where the system is at
  rest gives a finite R² (the 0.1.0 bug found by the pendulum test).
- **CLI.**
  - `synth` → `fit` → `scan` end to end.
  - Exit status 1 on errors, 0 when clean, 2 on a usage error (an equation
    whose variables do not match the data).

## 2. Synthetic validation

```bash
mihenk validate                 # 18 cases, about 30 seconds
mihenk validate --quick         # 17 cases, skips the slow heavy-noise case
mihenk validate --seeds 20      # every case over seeds 0-19, about 10 minutes
```

Each case simulates a system with a known law and fits it with PySINDy (or
uses a hand-written equation). It then scans the result and compares mihenk's
tier-0/1 flags with what the case was built to show. The ground-truth
comparison (tier 2) is printed next to it. The pass rules are:

| case kind | passes when |
|---|---|
| built to be clean, and the fit found the true law | no **errors** (warnings are counted separately as false warnings) |
| built to be clean, but the fit found a wrong law | at least one flag. A wrong law that nothing flags is the worst outcome |
| built to show a failure | the expected flag fires (`a \| b` in the output means any one of them) |
| any case with `forbid` | none of the forbidden flags fires (the clean pendulum must not raise `KINEMATIC_MISMATCH`) |

The multi-seed summary has three columns:

- `passed`: runs that passed out of all runs.
- `warned (clean)`: runs where a correct law received a warning.
- `wrong law, flagged`: runs where SINDy returned a wrong law and mihenk
  flagged it.

The README quotes the 20-seed table for 0.1.0:

- 280/280 runs pass.
- All 91 wrong laws are flagged.
- No correct law gets an error.
- 3 of 175 correct-law runs get a warning.

**Caveat.** Thresholds were adjusted while looking at seeds 0–9, so seeds
10–19 are the only part of that table that played no role in calibration.
A cleaner calibration/test split is on the [roadmap](ROADMAP.md).

### When a validation run fails

1. Rerun the single case with its seed and read the full report:
   ```python
   from mihenk import validate
   C = {c.name: c for c in validate.CASES}
   r = validate.run_case(C["lorenz_noise1"], seed=6)
   print(r.report.render_text())
   ```
2. Decide whether mihenk is wrong, or whether the case expectation is wrong.
   In 0.1 development, three apparent failures were mihenk being right:
   - SINDy found a wrong law on "clean" data, and mihenk flagged it.
   - A Lorenz trajectory lay near a quadric surface.
   - A noisy fit was flagged by different detectors than the ones listed.
3. Change a threshold only with a reason that holds beyond the failing seed,
   and rerun `--seeds 20` afterwards.

## 3. Pre-registered synthetic tests

```bash
python experiments/kinematic_synthetic/sweep.py 10 20   # the KINEMATIC_MISMATCH test seeds, about 1.5 minutes
```

Compare the output with `experiments/kinematic_synthetic/test_seeds_10-19.txt`.
Thresholds come from seeds 0–9 only (`sweep.py 0 10`). If you change a
`kin_` threshold, the 10–19 numbers stop being a test. Calibrate on new
seeds, and pre-register again before running fresh test seeds.

## 4. Experimental data

```bash
python experiments/pendulum/run.py          # analyses A, B, C of PROTOCOL.md, about 10 seconds
python experiments/pendulum/explore_lag.py  # the exploratory follow-up, E1
python experiments/pendulum/run_kinematics.py   # 0.2: KINEMATIC_MISMATCH on A, B, C, E1
```

The bundles are committed, so these run offline. `prepare.py` rebuilds them
from a clone of <https://github.com/dynamicslab/MultiArm-Pendulum>. The
expected flags are recorded in `experiments/pendulum/results/summary.json`:

| analysis | flags | errors |
|---|---|---|
| A | `DERIVATIVE_ARTIFACT`, `FORKING_PATHS`, `UNSTABLE_STRUCTURE`, `WEAK_TERM` | none |
| B | `EQUIFINALITY`, `FORKING_PATHS`, `UNSTABLE_STRUCTURE` | `EQUIFINALITY` |
| C | `DERIVATIVE_ARTIFACT`, `EQUIFINALITY`, `FORKING_PATHS`, `UNSTABLE_STRUCTURE`, `WEAK_TERM` | `EQUIFINALITY`, `UNSTABLE_STRUCTURE` |

If a code change alters these flags, say why in `RESULTS.md` and keep the
old run's log next to the new one. **Never edit the predictions in
`PROTOCOL.md`.** A new experiment gets its own protocol, committed before
its first fit.

## Before a release

1. `python -m pytest -q` is green, with PySINDy installed (35 passed).
2. `mihenk validate --seeds 20` is green. If a detector or threshold changed,
   paste the new table into the README.
3. `python experiments/pendulum/run.py` gives the flags in the table above,
   or the change is explained in `RESULTS.md`.
4. Build the wheel, install it in a fresh virtualenv **without** PySINDy, and
   run a tier-0 scan there. Then check the wheel again with PySINDy. See
   [PUBLISHING.md](PUBLISHING.md).
