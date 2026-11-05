# Installing mihenk

## The package

```bash
pip install mihenk                 # core: tier-0 checks from data + equation
pip install "mihenk[pysindy]"      # adds PySINDy: `mihenk fit` and the tier-1 refits
```

Python 3.10 or newer. The core needs three packages: numpy, scipy and sympy.
It runs on a CPU, offline, with no account, API or GPU. It never imports
PySINDy unless you use the PySINDy adapter.

Optional extras:

| extra | adds | for |
|---|---|---|
| `mihenk[pysindy]` | PySINDy ≥ 2.0 (and scikit-learn) | `mihenk fit`, `mihenk.scan_pysindy`, tier-1 checks (`UNSTABLE_STRUCTURE`, `FORKING_PATHS`) |
| `mihenk[provenance]` | daftar | `mihenk scan --daftar` |
| `mihenk[test]` | pytest, PySINDy | running the test suite |

Check what is available with:

```bash
mihenk info
```

Without PySINDy, `mihenk scan` still runs every tier-0 check and lists the
tier-1 checks as *not assessed*. It never reports them as clean.

### From source

```bash
git clone https://github.com/aslansd/mihenk.git
cd mihenk
pip install -e ".[test]"
python -m pytest -q
```

## Versions it was built and tested with

| package | version |
|---|---|
| Python | 3.13 (Linux). It also installs on 3.11 (macOS, Apple Silicon). |
| PySINDy | 2.1.0 |
| numpy / scipy / sympy | 2.4 / 1.17–1.18 / 1.14 |

The PySINDy adapter is written against PySINDy **2.x**. In 2.x, feature names
go to `fit()` rather than the `SINDy(...)` constructor. mihenk does not support
PySINDy 1.x.

The derivative-free integral form uses `scipy.integrate.cumulative_simpson`
(SciPy ≥ 1.12). On older SciPy it falls back to the trapezoid rule, which is
less accurate on clean data. `UNDERFIT` then needs its share-of-signal gate to
stay quiet.

## Installing on an unreliable connection

pip downloads can break partway through, for example over an unstable VPN,
with `Connection interrupted while downloading` or `Read timed out`. Try
these in order.

**1. Give pip more patience.** pip resumes partial downloads, and with more
retries a slow connection often gets through:

```bash
pip install --retries 20 --timeout 120 --resume-retries 50 "mihenk[pysindy]"
```

(`--resume-retries` needs a recent pip; drop it if pip does not recognise it.)

**2. Install from files.** On any machine with a working connection, download
the wheels for *your* platform. Then copy them over and install without the
network. For example, for an Apple Silicon Mac with Python 3.11:

```bash
# on the machine with a good connection
pip download "mihenk[pysindy]" -d mihenk-wheels \
    --only-binary=:all: --platform macosx_12_0_arm64 \
    --python-version 3.11 --implementation cp --abi cp311

# on your machine, after copying the folder
pip install --no-index --find-links mihenk-wheels "mihenk[pysindy]"
```

For an Intel Mac use `--platform macosx_10_9_x86_64`, and for Linux use
`--platform manylinux2014_x86_64`. Match `--python-version` and `--abi` to
`python --version`. Only scikit-learn, numpy and scipy are platform-specific.
Everything else mihenk needs is pure Python.

To see which files pip still needs, run the install once without
`--no-index`. Every `Using cached ...` line is already on your machine, and
every `Downloading ...` line is a file to fetch.

**3. Conda.** If your Python environment comes from conda, conda-forge's
servers are sometimes reachable when PyPI's are not:

```bash
conda install -c conda-forge numpy scipy sympy scikit-learn pysindy
pip install --no-deps mihenk
```

## Editable installs and stale copies

If `pip show mihenk` reports a location inside `site-packages` rather than your
clone, the tests and the CLI run against that installed copy, and your edits
will not show. Switch to an editable install once:

```bash
pip install -e .
# or, with no network at all:
pip install --no-index --no-deps --no-build-isolation -e .
```

The offline form needs `setuptools` already present in the environment.
