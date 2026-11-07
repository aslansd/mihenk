# mihenk

**A touchstone for discovered equations.** Before you interpret the law that [PySINDy](https://github.com/dynamicslab/pysindy) found in your data, mihenk checks whether it could be an artifact: of noise amplified by numerical differentiation, of the particular states your data happened to visit, of a library that was too small, or of analysis choices.

*Mihenk* (Turkish *mihenk taşı*, Azerbaijani *mehək daşı*, Persian *mahak* محک) is the touchstone used to test whether gold is real.

[Install](INSTALL.md) · [Testing](TESTING.md) · [Roadmap](ROADMAP.md) ·
[Publishing](PUBLISHING.md) · [Changelog](CHANGELOG.md)

```bash
pip install "mihenk[pysindy]"
```

The core needs only numpy, scipy and sympy, runs on a CPU, and works offline. It needs no account, API or GPU.

## Why

A SINDy fit always returns *an* equation. Whether that equation is the law of your system is a separate question, and the fit's R² does not answer it. Here is SINDy on clean data from the Hopf normal form, observed on its limit cycle:

```
true law:        x' = x - y - x^3 - x*y^2
SINDy returns:   x' = -0.6667 y - 0.3333 x^2*y - 0.3333 y^3      (derivative R² = 1.000)
```

On the circle `x² + y² = 1`, both equations (and `x' = -y`) fit perfectly. mihenk finds this from the data and the equation alone:

```
ERROR     EQUIFINALITY  [tier 0: data + equation]
          where: x'
          A sparser equation fits your data exactly as well (1 term instead of 3): x' = -y. On the
          observed states some candidate terms are combinations of others, so the data cannot choose
          between these laws.
          fix: Collect data that leaves the current set of states (other initial conditions,
          transients, perturbations) so the alternatives make different predictions.
```

Following the fix works: with transients from six starting points, SINDy recovers the true law and mihenk reports no findings.

## Quick start

```python
import pysindy as ps
import mihenk
from mihenk import synth

# example data with a known law; use your own x (T, n) and t (T,) instead
b = synth.simulate("lorenz", noise=0.05, seed=0)
x, t = b.xs[0], b.ts[0]

model = ps.SINDy(optimizer=ps.STLSQ(threshold=0.1))
model.fit(x, t=t, feature_names=["x", "y", "z"])

report = mihenk.scan_pysindy(model, x, t)
print(report.render_text())
report.save("mihenk_report")       # then open mihenk_report/report.html
```

From the command line:

```bash
mihenk synth lorenz --noise 0.05 -o lorenz.npz          # or bring your own .npz / .csv
mihenk fit lorenz.npz -o model.json --threshold 0.1      # PySINDy, polynomial library
mihenk scan lorenz.npz --model model.json --out report   # exit status 1 on errors
```

`--model` also takes a plain text file, including PySINDy's own `model.print()` output:

```
(x)' = -9.999 x + 9.999 y
(y)' = 27.992 x + -0.999 y + -1.000 x z
(z)' = -2.666 z + 1.000 x y
```

## What mihenk checks

| Finding | Tier | The failure it names | How it is tested |
|---|---|---|---|
| `EQUIFINALITY` | 0 | Another law fits the data exactly as well | Relations among candidate terms on the observed states (SVD); sparsest equally good subset; one-term swaps |
| `DERIVATIVE_ARTIFACT` | 0 | A term or coefficient produced by differentiated noise | Refit the equation's own terms with a **derivative-free integral form** and compare signs and values |
| `WEAK_TERM` | 0 | The data can't tell the equation with a term from without it | F statistic for dropping each term, in the integral form |
| `UNDERFIT` | 0 | A left-out term explains much of what the equation misses | Add each library term, and each term one polynomial degree beyond the library |
| `ROLLOUT_TRAP` | 0 | Fits derivatives but does not reproduce the system when simulated | Short rollouts from held-out states; full rollouts; escapes; long-run spread and mean |
| `UNSTABLE_STRUCTURE` | 1 | Terms come and go when the data are resampled | Bootstrap (trajectories, or blocks of one trajectory), refit through PySINDy |
| `FORKING_PATHS` | 1 | The law depends on defensible analysis choices | Refit under other derivative methods and sparsity thresholds |
| `MISSED_TERM`, `SPURIOUS_TERM`, `COEF_ERROR` | 2 | (synthetic data only) The discovered law against the true one | Term-by-term comparison |

**Tiers** describe what evidence a check needs:

- Tier 0: only your data and the equation. mihenk never calls the discovery tool; it refits only the coefficients of term sets the equation already has.
- Tier 1: a refit callable, which you get automatically with a PySINDy model.
- Tier 2: ground truth, which only synthetic data has.

Every report lists which checks ran, which were clean, and which were **not assessed** and why. A check that did not run is not a clean result.

Every threshold is named and can be overridden: `mihenk scan ... --set weak_term_f=5`, or `mihenk.scan(..., thresholds=mihenk.Thresholds(weak_term_f=5))`.

## Validation

`mihenk validate` runs 14 synthetic cases built to show each failure, plus cases that must come out clean. `mihenk validate --seeds 20` repeats them over 20 random seeds (different initial conditions and noise):

```
case                               passed  warned (clean)   wrong law, flagged
lorenz_clean                       20/20             0/20                  1/1
lorenz_noise1                      20/20             3/20                  3/3
van_der_pol_clean                  20/20             0/20                    -
duffing_clean                      20/20             0/20                  1/1
linear_clean                       20/20             0/20                    -
lotka_volterra_many                20/20             0/20                    -
stuart_landau_many                 20/20             0/20                    -
stuart_landau_on_cycle             20/20                -                20/20
lotka_volterra_one_orbit           20/20             0/20                    -
lotka_volterra_one_orbit_noisy     20/20                -                  6/6
duffing_small_library              20/20                -                20/20
lorenz_noise5_low_threshold        20/20                -                20/20
van_der_pol_noise10                20/20                -                20/20
linear_wrong_sign                  20/20                -                    -

280/280 case runs pass over seeds 0-19.
```

Over these runs:
- Every one of the 91 times SINDy returned a wrong law, mihenk flagged it, including the 5 times this happened in cases meant to be clean.
- In the 9 cases built to come out clean, none of the 175 correct-law runs received an error.
- 3 of those 175 runs received a warning. All three are single mild-noise Lorenz trajectories lying close to a quadric surface, where a constant term appears in some refits.

## Experimental data

mihenk was tested on the Brunton lab's experimental single pendulum, with the
predictions committed before any fit
([protocol](experiments/pendulum/PROTOCOL.md),
[results](experiments/pendulum/RESULTS.md)). The reference law is
`omega' = 64.22 sin(theta) - 0.0672 omega`. Of seven predictions, 5 held, 1
partly held and 1 failed. In brief:

- **Large swings, informed library.** SINDy found sin θ = 63.82, and mihenk
  left that term alone. It warned about extra terms in θ' (which should be
  exactly ω). Following that up showed the lab's filtered velocity signal
  lags the angle by 0.7–1.4 ms. That lag also cancelled most of the damping,
  which is why SINDy missed it.
- **Small swings.** sin θ ≈ π − θ, so mihenk reported `EQUIFINALITY` errors
  and a sparser equally good law, as predicted.
- **What mihenk missed:**
  - an equation that injects energy;
  - a long rollout that drifts (R² 0.73);
  - that a cubic library lacks sin θ. The first missing term is degree 5,
    beyond `UNDERFIT`'s one-degree search, and within the observed swings
    the cubic is a fine law anyway.

  These lead the [roadmap](ROADMAP.md).

## Limits (read these)

What comes next is in [ROADMAP.md](ROADMAP.md).


- **Research prototype (alpha).** Thresholds are calibrated on mihenk's own toy systems and were not changed after the first experimental test. Treat them as starting points.
- **mihenk rules out named failures. It cannot prove that an equation is the law.** No findings means that no tested failure explains the equation.
- ODE models only, without control inputs. PDEs and weak-form libraries (`WeakPDELibrary`) are not supported yet.
- The polynomial library is assumed when you pass a text equation. With a PySINDy model, mihenk uses its real library (polynomial and Fourier libraries are tested).
- mihenk checks only physics you can see in the data. It does not know that a free pendulum must lose energy unless a check for that exists (see the roadmap).
- The integral form uses overlapping windows, so F statistics are screening statistics, not exact tests.
- For chaotic systems a low full-length rollout R² is expected. mihenk compares long-run spread and mean instead.

## Provenance

With `pip install "mihenk[provenance]"`, `mihenk scan ... --daftar` records each scan (inputs, thresholds, environment, findings) with [daftar](https://pypi.org/project/daftar/), so `daftar diff` can show what changed between two scans.

## Related work

- Brunton, Proctor & Kutz (2016), *PNAS*: SINDy.
- de Silva et al. (2020) and Kaptanoglu et al. (2022), *JOSS*: PySINDy.
- Reinbold, Gurevich & Grigoriev (2020/2021): integral (weak) formulations, and ensembles that flag missing or spurious terms.
- Hauger, Scholl & Kutyniok (2024): robust identifiability for symbolic recovery of differential equations.

## License

Apache-2.0. Author: Aslan Satary Dizaji.
