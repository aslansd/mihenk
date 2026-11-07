"""Derivatives, derivative-free (integral) regression targets, refits and rollouts.

mihenk deliberately refits only with *fixed* term sets here (ordinary least
squares on the terms an equation already uses). Changing the structure is the
discovery tool's job and happens only in the tier-1 checks, through the tool.
"""

from __future__ import annotations

import numpy as np
from scipy.integrate import cumulative_trapezoid, solve_ivp

try:  # scipy >= 1.12: fourth-order accurate, so the integral form's own error stays below real effects
    from scipy.integrate import cumulative_simpson as _cumint
except ImportError:  # pragma: no cover
    _cumint = cumulative_trapezoid
from scipy.signal import savgol_filter


# ------------------------------------------------------------- derivatives
def finite_difference(x: np.ndarray, t: np.ndarray) -> np.ndarray:
    return np.gradient(x, t, axis=0, edge_order=2)


def savgol(x: np.ndarray, t: np.ndarray, window: int, order: int = 3) -> np.ndarray:
    window = min(window, len(t) - (1 - len(t) % 2))
    if window % 2 == 0:
        window -= 1
    order = min(order, window - 1)
    dt = np.median(np.diff(t))
    return savgol_filter(x, window, order, deriv=1, delta=dt, axis=0, mode="interp")


def derivative_methods(n_samples: int) -> dict:
    """The derivative estimators mihenk compares, sized to the trajectory length."""
    w_small = 7 if n_samples >= 70 else 5
    w_large = int(np.clip(n_samples // 50, w_small + 4, 41)) | 1
    return {
        "finite_difference": lambda x, t: finite_difference(x, t),
        f"savgol_{w_small}": lambda x, t, w=w_small: savgol(x, t, w),
        f"savgol_{w_large}": lambda x, t, w=w_large: savgol(x, t, w),
    }


def stack_derivatives(pairs, method, trim: int = 0):
    """Concatenate states and derivative estimates over trajectories, trimming edges."""
    X, D = [], []
    for x, t in pairs:
        d = method(x, t)
        sl = slice(trim, len(t) - trim if trim else None)
        X.append(x[sl])
        D.append(d[sl])
    return np.vstack(X), np.vstack(D)


# ------------------------------------------------------------ integral form
def integral_targets(pairs, theta_fn, window: int):
    """Derivative-free regression: ``x(t+W) - x(t) = sum_j c_j * integral(theta_j)``.

    Noise in ``x`` enters the target once instead of being amplified by
    differentiation, which is what lets mihenk tell a derivative artifact from
    a real term.
    """
    A, B = [], []
    for x, t in pairs:
        th = theta_fn(x)
        cum = _cumint(th, x=t, axis=0, initial=0)
        w = min(window, len(t) - 1)
        starts = np.arange(0, len(t) - w, max(1, w // 2))
        A.append(cum[starts + w] - cum[starts])
        B.append(x[starts + w] - x[starts])
    return np.vstack(A), np.vstack(B)


def integral_window(n_samples: int) -> int:
    return int(np.clip(n_samples // 100, 5, 50))


# ------------------------------------------------------------- regression
def lstsq(A: np.ndarray, b: np.ndarray):
    """Least squares with column scaling; returns (coef, residual norm)."""
    if A.shape[1] == 0:
        return np.zeros(0), float(np.linalg.norm(b))
    scale = np.linalg.norm(A, axis=0)
    scale[scale == 0] = 1.0
    coef, *_ = np.linalg.lstsq(A / scale, b, rcond=None)
    coef = coef / scale
    return coef, float(np.linalg.norm(A @ coef - b))


def r2(y: np.ndarray, yhat: np.ndarray, var=None) -> float:
    """Coefficient of determination.

    ``var`` (per-variable variance) sets the scale instead of ``y``'s own
    spread. Rollouts pass the variance of the whole data set, so a short
    segment where the system barely moves (a pendulum at rest) is judged on
    the scale of the data, not on its own near-zero spread.
    """
    y = np.asarray(y)
    if var is None:
        ss = np.sum((y - y.mean(axis=0)) ** 2)
    else:
        ss = float(np.sum(np.asarray(var, float) * y.shape[0]))
    if ss == 0:
        return float("nan")
    return float(1 - np.sum((y - yhat) ** 2) / ss)


def normalized_columns(A: np.ndarray) -> np.ndarray:
    n = np.linalg.norm(A, axis=0)
    n[n == 0] = 1.0
    return A / n


# ---------------------------------------------------------------- rollout
def rollout(model, x0, t, bounds=None, max_seconds_hint=None):
    """Integrate a model from ``x0`` on times ``t``.

    Returns ``(trajectory, status)`` where status is ``"ok"``, ``"escaped"``
    (left ``bounds``), ``"nonfinite"`` or ``"failed"``. The trajectory is
    padded with NaN after an escape so it always has ``len(t)`` rows.
    """
    f = model.rhs_function()
    events = []
    if bounds is not None:
        lo, hi = bounds

        def escape(_t, x):
            return float(np.min(np.concatenate([x - lo, hi - x])))

        escape.terminal = True
        events.append(escape)
    try:
        with np.errstate(all="ignore"):
            sol = solve_ivp(f, (t[0], t[-1]), np.asarray(x0, float), t_eval=t, method="LSODA",
                            rtol=1e-8, atol=1e-10, events=events or None)
    except Exception:  # pragma: no cover - integrator internals
        return np.full((len(t), len(x0)), np.nan), "failed"
    y = sol.y.T
    out = np.full((len(t), len(x0)), np.nan)
    out[: len(y)] = y
    if not np.all(np.isfinite(y)):
        return out, "nonfinite"
    if sol.status == 1:
        return out, "escaped"
    if len(y) < len(t) or sol.status < 0:
        return out, "failed"
    return out, "ok"
