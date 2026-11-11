"""PySINDy adapter (``pip install "mihenk[pysindy]"``), written against PySINDy 2.x.

Turns a fitted ``pysindy.SINDy`` into a :class:`~mihenk.model.Model`, and builds
the tier-1 callables: ``refit`` (the same configuration on resampled data) and
``variants`` (the same data under other defensible derivative methods and
sparsity thresholds).
"""

from __future__ import annotations

import numpy as np

from ..bundle import Bundle
from ..model import Model, parse_term


def _require():
    try:
        import pysindy  # noqa: F401
    except ImportError as exc:  # pragma: no cover
        raise ImportError('this needs PySINDy: pip install "mihenk[pysindy]"') from exc


def from_pysindy(sindy, variables=None, source="pysindy") -> Model:
    """Read variables, candidate library and coefficients from a fitted PySINDy model."""
    _require()
    if getattr(sindy, "n_control_features_", 0):
        raise NotImplementedError("models with control inputs (u) are not supported in mihenk 0.1")
    lib_name = type(sindy.feature_library).__name__
    if "Weak" in lib_name or "PDE" in lib_name:
        raise NotImplementedError(f"{lib_name} is not supported in mihenk 0.1 (ODE libraries only)")
    variables = list(variables or sindy.feature_names)
    names = sindy.get_feature_names()
    try:
        terms = [parse_term(n, variables) for n in names]
    except Exception as exc:
        raise ValueError(f"cannot read PySINDy feature names as expressions ({names[:5]}...): {exc}") from exc
    coefs = np.asarray(sindy.coefficients(), dtype=float)
    m = Model(variables, terms, coefs, library=list(terms), source=source)
    m.meta = {"library": lib_name, "optimizer": type(sindy.optimizer).__name__,
              "differentiation": type(sindy.differentiation_method).__name__}
    return m


def _fit_clone(sindy, pairs, variables, **overrides):
    from sklearn.base import clone

    m = clone(sindy)
    if overrides:
        m.set_params(**overrides)
    xs = [x for x, _ in pairs]
    ts = [t for _, t in pairs]
    import warnings

    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        if len(xs) == 1:
            m.fit(xs[0], t=ts[0], feature_names=list(variables))
        else:
            m.fit(xs, t=ts, feature_names=list(variables))
    model = from_pysindy(m, variables)
    # STLSQ can stop at max_iter with coefficients far below its own threshold still non-zero;
    # in a refit those are optimizer leftovers, not structure, so they are not counted as terms.
    thr = m.optimizer.get_params().get("threshold")
    if thr:
        small = (np.abs(model.coefficients) < thr) & (model.coefficients != 0)
        model.meta["pruned_below_threshold"] = int(small.sum())
        model.coefficients[small] = 0.0
    return model


def make_refit(sindy, variables=None):
    variables = list(variables or sindy.feature_names)
    return lambda pairs: _fit_clone(sindy, pairs, variables)


def make_variants(sindy, variables=None, threshold_factors=(0.5, 1.0, 1.5)) -> dict:
    """Other defensible analyses of the same data: derivative method x sparsity threshold."""
    import pysindy as ps

    variables = list(variables or sindy.feature_names)
    diffs = {"FiniteDifference": ps.FiniteDifference, "SmoothedFiniteDifference": ps.SmoothedFiniteDifference}
    current_diff = type(sindy.differentiation_method).__name__
    params = sindy.optimizer.get_params()
    th0 = params.get("threshold")
    out = {}
    for dname, dcls in diffs.items():
        for f in threshold_factors if th0 is not None else (1.0,):
            if dname == current_diff and f == 1.0:
                continue  # that is the original analysis
            over = {"differentiation_method": dcls()}
            label = dname
            if th0 is not None:
                over["optimizer__threshold"] = th0 * f
                label += f", threshold {th0 * f:.3g}"
            out[label] = (lambda pairs, over=over: _fit_clone(sindy, pairs, variables, **over))
    return out


def scan_pysindy(sindy, x, t, *, feature_names=None, name="pysindy", val=None, thresholds=None,
                 variants=True, seed=0, meta=None, kinematics=None):
    from ..scan import scan

    variables = list(feature_names or sindy.feature_names)
    xs = x if isinstance(x, (list, tuple)) else [x]
    if isinstance(t, (int, float)):
        ts = [np.arange(len(xi)) * float(t) for xi in xs]
    else:
        ts = t if isinstance(t, (list, tuple)) else [t]
    b = Bundle(list(xs), list(ts), variables, val=val, meta={"name": name, **(meta or {})})
    model = from_pysindy(sindy, variables)
    return scan(b, model, refit=make_refit(sindy, variables),
                variants=make_variants(sindy, variables) if variants else None,
                thresholds=thresholds, seed=seed, kinematics=kinematics)


# --------------------------------------------------------------- CLI helpers
DIFFS = {"fd": "FiniteDifference", "smoothed": "SmoothedFiniteDifference"}


def build_sindy(degree: int = 2, threshold: float = 0.1, diff: str = "fd", fourier: int = 0):
    """A PySINDy model from a small, serializable configuration (used by ``mihenk fit``).

    ``fourier > 0`` adds ``FourierLibrary(n_frequencies=fourier)`` to the polynomial library.
    """
    _require()
    import pysindy as ps

    if diff not in DIFFS:
        raise ValueError(f"diff must be one of {', '.join(DIFFS)}")
    library = ps.PolynomialLibrary(degree=degree)
    if fourier:
        library = library + ps.FourierLibrary(n_frequencies=int(fourier))
    return ps.SINDy(optimizer=ps.STLSQ(threshold=threshold), feature_library=library,
                    differentiation_method=getattr(ps, DIFFS[diff])())


def fit_bundle(bundle: Bundle, degree: int = 2, threshold: float = 0.1, diff: str = "fd", fourier: int = 0):
    """Fit PySINDy to a bundle's trajectories; returns ``(sindy, model)``."""
    import warnings

    sindy = build_sindy(degree, threshold, diff, fourier)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        if bundle.n_traj == 1:
            sindy.fit(bundle.xs[0], t=bundle.ts[0], feature_names=bundle.variables)
        else:
            sindy.fit(bundle.xs, t=bundle.ts, feature_names=bundle.variables)
    model = from_pysindy(sindy, bundle.variables, source="pysindy (mihenk fit)")
    model.meta["pysindy_config"] = {"degree": degree, "threshold": threshold, "diff": diff}
    if fourier:
        model.meta["pysindy_config"]["fourier"] = int(fourier)
    return sindy, model
