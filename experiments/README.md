# Experiments

Each experiment tests mihenk against something it cannot see in its own
synthetic validation. Each one is **pre-registered**: `PROTOCOL.md` fixes the
analyses, thresholds and predictions, and is committed before the first fit.
`RESULTS.md` scores every prediction afterwards, whether it held or not.
**A protocol is never edited after its commit.** A new question gets a new
protocol.

| experiment | question | protocol committed | outcome |
|---|---|---|---|
| [pendulum](pendulum/RESULTS.md) | Does mihenk behave on a real instrument with a known law (the Brunton-lab single pendulum)? | `59fbb4b` (mihenk 0.1.0) | 5 held, 1 partly held, 1 failed. It found a 0.7–1.4 ms lag in the velocity signal. |
| [kinematic_synthetic](kinematic_synthetic/RESULTS.md) | Does `KINEMATIC_MISMATCH` find lags, gains and offsets without false alarms (calibration on seeds 0–9, test on seeds 10–19)? | `6a9e672` (mihenk 0.2.0) | 5 held, 1 partly held |
| [physics_synthetic](physics_synthetic/RESULTS.md) | Do the 0.3 checks (declared physics, twin long rollouts, two-degree `UNDERFIT`, centred relations) work without false alarms, and do they cause regressions? | `61ca977` (mihenk 0.3.0) | 8 held, 2 failed; one failure was a regression, fixed afterwards |

## Layout

Every experiment folder has the same shape:

```
PROTOCOL.md      analyses and predictions, committed before any fit
RESULTS.md       the scored outcome
data/            input data, with its license (if the experiment uses outside data)
scripts/         everything needed to reproduce, run from the repository root
results/         outputs: reports (report.html to read, report.json with every number),
                 summary.json, and logs/ with raw run output
```

The pendulum's `results/` is split by mihenk release:

- `0.1-analyses/`: the pre-registered analyses A, B and C, and the exploratory
  E1.
- `0.2-kinematics/`: the same analyses with θ' = ω declared, a development
  case for `KINEMATIC_MISMATCH`.
- `0.3-physics/`: the same analyses with θ' = ω and dissipation declared, a
  development case for the 0.3 checks.

Its `logs/` folder holds the console output of each run, one file per mihenk
version.

## Running

From the repository root, with `pip install -e ".[test]"`:

```bash
python experiments/pendulum/scripts/run.py               # analyses A, B, C (about 10 s)
python experiments/pendulum/scripts/explore_lag.py       # exploratory E1
python experiments/pendulum/scripts/run_kinematics.py    # 0.2 development case
python experiments/pendulum/scripts/run_physics.py       # 0.3 development case
python experiments/kinematic_synthetic/scripts/sweep.py 10 20   # 0.2 test seeds (about 1.5 min)
python experiments/physics_synthetic/scripts/sweep.py 10 20     # 0.3 test seeds (about 2 min)
```

These overwrite the files in `results/`. Run them on a branch, and compare
with `git diff` before committing anything.
