# KINEMATIC_MISMATCH: test results (seeds 10–19)

The predictions and thresholds were committed in [PROTOCOL.md](PROTOCOL.md)
(`6a9e672`) before these seeds were run. The raw outputs are
`results/calibration_seeds_0-9.txt` and `results/test_seeds_10-19.txt`, and
the script is `scripts/sweep.py`.

```
noise 0.0    none        error  0  warning  0  clean 10 | lag est +0.000 samp, share 0.0000, F median 6.13
noise 0.0    lag0.7ms    error 10  warning  0  clean  0 | lag est +0.350 samp, share 0.0049, F median 7.51e+16
noise 0.0    lag0.3ms    error  0  warning 10  clean  0 | lag est +0.150 samp, share 0.0021, F median 1.32e+16
noise 0.0    gain1.02    error 10  warning  0  clean  0 | lag est +0.000 samp, share 0.0200, F median 1.11e+18
noise 0.0    offset0.05  error  0  warning 10  clean  0 | lag est +0.000 samp, share 0.0075, F median 1.72e+17
noise 0.002  none        error  0  warning  0  clean 10 | lag est +0.000 samp, share 0.0001, F median 0.0518
noise 0.002  lag0.7ms    error 10  warning  0  clean  0 | lag est +0.350 samp, share 0.0050, F median 98.1
noise 0.002  lag0.3ms    error  2  warning  8  clean  0 | lag est +0.150 samp, share 0.0022, F median 18
noise 0.002  gain1.02    error 10  warning  0  clean  0 | lag est +0.000 samp, share 0.0201, F median 1.58e+03
noise 0.002  offset0.05  error  0  warning 10  clean  0 | lag est +0.000 samp, share 0.0076, F median 221
noise 0.005  none        error  0  warning  0  clean 10 | lag est +0.000 samp, share 0.0004, F median 0.1
noise 0.005  lag0.7ms    error  9  warning  0  clean  1 | lag est +0.350 samp, share 0.0051, F median 15.8
noise 0.005  lag0.3ms    error  0  warning  0  clean 10 | lag est +0.150 samp, share 0.0022, F median 3.12
noise 0.005  gain1.02    error 10  warning  0  clean  0 | lag est +0.000 samp, share 0.0202, F median 256
noise 0.005  offset0.05  error  0  warning 10  clean  0 | lag est +0.000 samp, share 0.0076, F median 36
noise 0.01   none        error  0  warning  0  clean 10 | lag est +0.000 samp, share 0.0008, F median 0.122
noise 0.01   lag0.7ms    error  0  warning  5  clean  5 | lag est +0.350 samp, share 0.0051, F median 4.17
noise 0.01   lag0.3ms    error  0  warning  0  clean 10 | lag est +0.150 samp, share 0.0023, F median 0.861
noise 0.01   gain1.02    error 10  warning  0  clean  0 | lag est +0.000 samp, share 0.0204, F median 65.2
noise 0.01   offset0.05  error  0  warning 10  clean  0 | lag est +0.000 samp, share 0.0077, F median 9.29
```

## Score

| | prediction | outcome |
|---|---|---|
| K1 | no flag without an imperfection, 0 of 40 | **held** (0 of 40) |
| K2 | 0.7 ms lag: error in ≥ 9/10 at noise 0–0.5%; flagged in ≥ 7/10 at 1% | **partly held**. Errors in 10, 10 and 9 of 10. At 1% noise, flagged in only **5** of 10. |
| K3 | 0.3 ms lag: 10/10 at noise 0; ≥ 8/10 at 0.2% | **held** (10 and 10) |
| K4 | gain 1.02: error 10/10 at every noise level | **held** (40 of 40) |
| K5 | offset 0.05: flagged 10/10 at every noise level | **held** (40 of 40) |
| K6 | median lag within ±0.05 samples of the truth up to 0.5% noise | **held** (0.350 and 0.150 at every level) |

5 held, 1 partly held.

## Reading the results

**No false alarms.** Over calibration and test together, 80 runs had no
imperfection and none raised `KINEMATIC_MISMATCH`. The largest F seen without
an imperfection was 0.5, against a threshold of 10.

**The lag estimate holds up even when it is not significant.** Up to 0.5%
noise, the median estimate equals the truth in both sets. At 1% noise it is
within 0.025 samples (0.375 and 0.175 in calibration; 0.350 and 0.150 in the
test), even though F falls below the threshold there. What noise removes is the *confidence*, not the *estimate*. With
1% noise on θ (about 0.01 rad), each window's change in θ carries more noise
than a 0.35-sample timing error produces.

**Where K2 fell short.** At 1% noise, calibration flagged the 0.7 ms lag in 9
of 10 runs and the test in 5 of 10, so the prediction of at least 7 was too
optimistic. We checked where those flags came from. In both sets every
one came from the *equation side* (the discovered θ' differing from ω by more
than 0.5%), and none from the data fit: 9 of 9 in calibration, 5 of 5 in the
test. How
far SINDy's θ' strays at that noise level varies a lot from seed to seed. We
did not change any threshold after seeing this.

**What the power means in practice.** The experimental pendulum's angle comes
from an optical encoder, whose noise is far below 0.2% of the swing, and there
the check is decisive (F = 2.8×10³). With position noise of 1% or more, a
sub-sample lag cannot be shown from the data alone. mihenk then says nothing
rather than guessing.
