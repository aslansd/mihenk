# KINEMATIC_MISMATCH: calibration and pre-registered test

`KINEMATIC_MISMATCH` was designed after the experimental pendulum
([../pendulum/RESULTS.md](../pendulum/RESULTS.md)) showed a velocity lag.
**That makes the pendulum a development case for this check, not a test.** The
test is this synthetic sweep. Its thresholds were chosen on **seeds 0–9**
only, and the predictions below were committed before **seeds 10–19** were
run.

## Setup

- **System:** `synth.simulate("pendulum")`, i.e. θ' = ω and
  ω' = −64 sin θ − 0.07 ω. It starts from |θ₀| between 1 and 2 rad, is sampled
  every 2 ms for 10 s, and is declared `{"theta": "omega"}`.
- **Discovery:** PySINDy, with `PolynomialLibrary(degree=1) +
  FourierLibrary(n_frequencies=1)`, STLSQ threshold 0.02, finite differences.
- **Imperfections**, applied to the recorded ω only:
  - none;
  - a lag of 0.7 ms (0.35 samples, the size found in the experiment) or
    0.3 ms;
  - a gain of 1.02;
  - an offset of 0.05 rad/s.
- **Noise:** 0, 0.2%, 0.5% or 1% of each variable's standard deviation, on
  both θ and ω.
- **Script:** `sweep.py`.

## Thresholds, fixed on seeds 0–9

These are the defaults in `mihenk/scan.py`:

| threshold | value | why |
|---|---|---|
| `kin_f` | 10 | With no imperfection, F never exceeded 0.5 in 30 noisy runs. With a 0.7 ms lag, it was at least 11.4 up to 0.5% noise. It started at 100, which left power unused. |
| `kin_min_share` | 0.001 | Keeps clean data's numerical noise out (share 0.0000 at no noise). |
| `kin_eq_share` | 0.005 | The discovered θ' may differ from ω by 0.5% of ω before it counts on its own. |
| `kin_error_share` | 0.01 | The share at which a mismatch becomes an error. |
| `kin_coef_change` | 0.25 | Correcting ω moves a coefficient of another equation by 25% or more, so the error is about the law, not only the kinematic equation. |
| `kin_coef_min_share` | 0.001 | Ignore coefficients of terms that carry almost nothing. |

Calibration result for seeds 0–9 (`sweep.py 0 10`):

```
noise 0.0    none        error  0  warning  0  clean 10 | lag est +0.000 samp, share 0.0000, F median 5.17
noise 0.0    lag0.7ms    error 10  warning  0  clean  0 | lag est +0.350 samp, share 0.0051, F median 8.75e+16
noise 0.0    lag0.3ms    error  0  warning 10  clean  0 | lag est +0.150 samp, share 0.0022, F median 1.61e+16
noise 0.0    gain1.02    error 10  warning  0  clean  0 | lag est +0.000 samp, share 0.0200, F median 1.12e+18
noise 0.0    offset0.05  error  0  warning 10  clean  0 | lag est +0.000 samp, share 0.0083, F median 1.98e+17
noise 0.002  none        error  0  warning  0  clean 10 | lag est +0.000 samp, share 0.0001, F median 0.0723
noise 0.002  lag0.7ms    error 10  warning  0  clean  0 | lag est +0.350 samp, share 0.0051, F median 102
noise 0.002  lag0.3ms    error  3  warning  7  clean  0 | lag est +0.150 samp, share 0.0022, F median 18.3
noise 0.002  gain1.02    error 10  warning  0  clean  0 | lag est +0.000 samp, share 0.0199, F median 1.6e+03
noise 0.002  offset0.05  error  0  warning 10  clean  0 | lag est +0.000 samp, share 0.0083, F median 283
noise 0.005  none        error  0  warning  0  clean 10 | lag est +0.000 samp, share 0.0005, F median 0.162
noise 0.005  lag0.7ms    error 10  warning  0  clean  0 | lag est +0.350 samp, share 0.0051, F median 16.3
noise 0.005  lag0.3ms    error  0  warning  0  clean 10 | lag est +0.150 samp, share 0.0023, F median 3.5
noise 0.005  gain1.02    error 10  warning  0  clean  0 | lag est +0.000 samp, share 0.0197, F median 258
noise 0.005  offset0.05  error  1  warning  9  clean  0 | lag est +0.000 samp, share 0.0083, F median 45.3
noise 0.01   none        error  0  warning  0  clean 10 | lag est +0.025 samp, share 0.0011, F median 0.209
noise 0.01   lag0.7ms    error  0  warning  9  clean  1 | lag est +0.375 samp, share 0.0053, F median 4.65
noise 0.01   lag0.3ms    error  0  warning  1  clean  9 | lag est +0.175 samp, share 0.0026, F median 1.02
noise 0.01   gain1.02    error 10  warning  0  clean  0 | lag est +0.025 samp, share 0.0195, F median 65.5
noise 0.01   offset0.05  error  1  warning  9  clean  0 | lag est +0.025 samp, share 0.0083, F median 11.4
```

## Predictions for seeds 10–19

The thresholds and code are frozen at the commit that adds this file.

- **K1. No false alarms.** Without an imperfection, `KINEMATIC_MISMATCH` is
  not flagged in any of the 40 runs (4 noise levels × 10 seeds).
- **K2. The experiment-sized lag (0.7 ms).**
  - Flagged as an **error** in at least 9 of 10 runs at noise 0, 0.2% and
    0.5%.
  - Flagged at all (error or warning) in at least 7 of 10 runs at noise 1%.
- **K3. A smaller lag (0.3 ms).**
  - Flagged in 10 of 10 runs at noise 0.
  - Flagged in at least 8 of 10 runs at noise 0.2%.
  - No claim at 0.5% and 1% noise, where calibration found little power
    (0 and 1 of 10).
- **K4. Gain 1.02.** An error in 10 of 10 runs at every noise level.
- **K5. Offset 0.05.** Flagged in 10 of 10 runs at every noise level.
- **K6. The lag estimate is unbiased.** For both lags, the median estimate
  is within ±0.05 samples of the true value (0.35 and 0.15) at noise up to
  0.5%.

Results go in `RESULTS.md`, with every prediction scored as written here.
