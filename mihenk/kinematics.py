"""KINEMATIC_MISMATCH: does a declared derivative pair hold in the data and in the equation?

Many systems have variables that are, by definition, derivatives of others:
position and velocity, angle and angular velocity. Declare such a pair as
``{"theta": "omega"}`` (meaning theta' = omega). mihenk then checks two things.

1. **The data.** Derivative-free, over the same windows as the integral form:
   is ``theta(t+W) - theta(t)`` equal to the integral of ``omega`` over the
   window? If not, it fits the measured omega as::

       theta' = gain * omega_meas(t + lag) + offset

   A positive ``lag`` means the recorded omega runs late. That is typical of
   a causal filter, and it was what the Brunton-lab pendulum showed. A gain
   other than 1 suggests a scale or unit error. An offset suggests a bias.

2. **The equation.** Is the discovered ``theta'`` equal to ``omega``? If it
   has extra terms and the data show a lag, those terms are predicted: a
   velocity late by ``lag`` gives ``theta' = omega_meas + lag * omega'``. So
   mihenk compares the extra terms with ``lag`` times the discovered
   ``omega'`` equation, and reports how much of them the lag explains.

This check needs knowledge mihenk cannot infer from the data: that the pair
is kinematic. Without a declaration it reports *not assessed*, and suggests
pairs whose discovered equations look kinematic.
"""

from __future__ import annotations

import numpy as np
import sympy as sp
from scipy.interpolate import CubicSpline

from . import numerics as N
from .model import Model, term_label
from .report import Check, Finding


def validate_pairs(kinematics, variables) -> dict:
    """Normalize ``{"theta": "omega"}`` / ``["theta=omega"]`` and check the names."""
    if not kinematics:
        return {}
    if isinstance(kinematics, (list, tuple)):
        pairs = {}
        for item in kinematics:
            if "=" not in item:
                raise ValueError(f"kinematic pair must look like theta=omega, got {item!r}")
            p, v = (s.strip() for s in item.split("=", 1))
            pairs[p.rstrip("'")] = v
        kinematics = pairs
    out = {}
    for p, v in dict(kinematics).items():
        p, v = str(p).strip().rstrip("'"), str(v).strip()
        for name in (p, v):
            if name not in variables:
                raise ValueError(f"kinematic pair {p}'={v}: '{name}' is not a variable "
                                 f"({', '.join(variables)})")
        if p == v:
            raise ValueError(f"kinematic pair {p}'={v} pairs a variable with itself")
        out[p] = v
    return out


def suggest_pairs(model: Model, tol: float = 0.05) -> dict:
    """Pairs whose discovered equation is dominated by another variable with coefficient ~1."""
    out = {}
    syms = {v: sp.Symbol(v) for v in model.variables}
    for i, p in enumerate(model.variables):
        for v, s in syms.items():
            if v == p or s not in model.terms:
                continue
            c = model.coefficients[i, model.terms.index(s)]
            if abs(c - 1) <= tol:
                out[p] = v
    return out


# ------------------------------------------------------------------ data fit
def _windows(t, window, margin):
    """Start/end indices of half-overlapping windows that stay ``margin`` seconds inside."""
    w = min(window, len(t) - 1)
    starts = np.arange(0, len(t) - w, max(1, w // 2))
    ends = starts + w
    keep = (t[starts] - margin >= t[0]) & (t[ends] + margin <= t[-1])
    return starts[keep], ends[keep]


def fit_pair(pairs, ip, iv, window, max_lag_samples=10, steps_per_sample=20):
    """Best ``gain``, ``lag`` and ``offset`` for ``x[:, ip]' = gain * x[:, iv](t + lag) + offset``.

    Returns a dict with the fit, the baseline (gain 1, lag 0, offset 0), and
    sums of squares, all on derivative-free window increments.
    """
    dt = float(np.median(np.concatenate([np.diff(t) for _, t in pairs])))
    L = max_lag_samples * dt
    lags = np.linspace(-L, L, 2 * max_lag_samples * steps_per_sample + 1)
    segs = []
    for x, t in pairs:
        s, e = _windows(t, window, L)
        if len(s) < 3:
            continue
        S = CubicSpline(t, x[:, iv]).antiderivative()
        segs.append((x[e, ip] - x[s, ip], t[s], t[e], S))
    if not segs:
        return None
    dp = np.concatenate([g[0] for g in segs])
    W = np.concatenate([g[2] - g[1] for g in segs])
    I0 = np.concatenate([S(te) - S(ts) for _, ts, te, S in segs])
    rss0 = float(np.sum((dp - I0) ** 2))
    best = None
    for lag in lags:
        I = np.concatenate([S(te + lag) - S(ts + lag) for _, ts, te, S in segs])
        A = np.column_stack([I, W])
        coef, *_ = np.linalg.lstsq(A, dp, rcond=None)
        rss = float(np.sum((A @ coef - dp) ** 2))
        if best is None or rss < best["rss"]:
            best = {"rss": rss, "lag": float(lag), "gain": float(coef[0]), "offset": float(coef[1])}
    # lag alone (gain 1, offset 0): how much does timing by itself explain?
    lag_only = min(
        (float(np.sum((np.concatenate([S(te + lag) - S(ts + lag) for _, ts, te, S in segs]) - dp) ** 2)), lag)
        for lag in lags)
    n = len(dp)
    best.update({
        "rss0": rss0, "rss_lag_only": lag_only[0], "lag_only": float(lag_only[1]), "n_windows": n,
        "signal": float(np.linalg.norm(dp)), "dt": dt, "max_lag": float(L),
        "at_edge": abs(abs(best["lag"]) - L) < dt / steps_per_sample,
    })
    best["share"] = float(np.sqrt(max(rss0 - best["rss"], 0.0)) / max(best["signal"], 1e-300))
    best["share_unexplained"] = float(np.sqrt(best["rss"]) / max(best["signal"], 1e-300))
    best["F"] = float(((rss0 - best["rss"]) / 3) / max(best["rss"] / max(n - 3, 1), 1e-300))
    return best


# ------------------------------------------------------------- the check
def check(ctx, kinematics, findings, checks):
    m, th = ctx.model, ctx.th
    variables = list(m.variables)
    pairs = validate_pairs(kinematics, variables)
    if not pairs:
        hint = suggest_pairs(m)
        reason = "no kinematic pairs declared"
        if hint:
            reason += "; the equations suggest " + ", ".join(f"{p}={v}" for p, v in hint.items())
        checks.append(Check("KINEMATIC_MISMATCH", 0, "not_assessed", reason))
        return
    rows, traj_rows, impact_rows = [], [], []
    flagged = False
    for p, v in pairs.items():
        ip, iv = variables.index(p), variables.index(v)
        fit = fit_pair(ctx.pairs, ip, iv, ctx.window, th.kin_max_lag)
        if fit is None:
            checks.append(Check("KINEMATIC_MISMATCH", 0, "not_assessed", f"{p}'={v}: trajectories too short"))
            continue
        for k, pair in enumerate(ctx.pairs):
            fk = fit_pair([pair], ip, iv, ctx.window, th.kin_max_lag)
            if fk is not None:
                traj_rows.append([f"{p}' = {v}", k, fk["lag"] / fk["dt"], fk["lag"], fk["gain"], fk["offset"],
                                  fk["share"]])
        rows.append([f"{p}' = {v}", fit["lag"] / fit["dt"], fit["lag"], fit["gain"], fit["offset"],
                     fit["share"], fit["F"]])
        ctx.metrics[f"kin_{p}_lag"] = fit["lag"]
        ctx.metrics[f"kin_{p}_gain"] = fit["gain"]

        # equation side
        eq = _equation_side(m, ctx, ip, iv, fit)
        data_bad = fit["share"] >= th.kin_min_share and fit["F"] >= th.kin_f
        eq_bad = eq["deviation"] >= th.kin_eq_share
        impact = _impact(m, ctx, ip, iv, fit, th.kin_coef_min_share) if data_bad else []
        big = [r for r in impact if r[4] >= th.kin_coef_change]
        for r in impact:
            impact_rows.append([f"{p}' = {v}"] + r)
        if not (data_bad or eq_bad):
            continue
        flagged = True
        parts = []
        if data_bad:
            bits = []
            if abs(fit["lag"]) >= th.kin_report_lag * fit["dt"]:
                late = "late" if fit["lag"] > 0 else "early"
                bits.append(f"{v} runs {late} by {abs(fit['lag']):.3g} s ({abs(fit['lag']) / fit['dt']:.2g} samples)")
            if abs(fit["gain"] - 1) >= th.kin_report_gain:
                bits.append(f"{v} is scaled by {1 / fit['gain']:.4g} (gain {fit['gain']:.4g} restores {p}' = {v})")
            v_rms = float(np.sqrt(np.mean(np.concatenate([x[:, iv] for x, _ in ctx.pairs]) ** 2)))
            if abs(fit["offset"]) >= th.kin_report_offset * max(v_rms, 1e-300):
                bits.append(f"{v} has a bias of {-fit['offset']:+.3g}")
            what = ("; ".join(bits) if bits else "the mismatch has no simple lag, gain or offset form")
            parts.append(
                f"In the data, {p} does not change by the integral of {v}: {what}. These explain "
                f"{fit['share']:.2%} of the change in {p} (F = {fit['F']:.3g})"
                + (f"; the lag hit the search limit of {th.kin_max_lag} samples, so it may be larger"
                   if fit["at_edge"] else ""))
        if eq_bad:
            parts.append(f"The discovered equation is not {p}' = {v}: {eq['text']}")
            if eq.get("lag_explains") is not None and data_bad and abs(fit["lag"]) > 0:
                parts.append(f"A {v} late by {fit['lag']:.3g} s predicts extra terms equal to {fit['lag']:.3g} "
                             f"times the discovered {v}' equation; that accounts for {eq['lag_explains']:.0%} "
                             f"of them")
        if big:
            parts.append("Correcting " + v + " for this changes the law: " + ", ".join(
                f"{eqn} {t} {c0:+.4g} -> {c1:+.4g}" for eqn, t, c0, c1, _ in big[:4])
                + " (derivative-free refits of the equation's own terms)")
        severe = bool(big) or (data_bad and fit["share"] >= th.kin_error_share) or (
            eq_bad and eq["deviation"] >= th.kin_error_share)
        fix = (f"Check how {v} is measured or filtered: estimate it from {p} with a matched, zero-phase filter, "
               f"or shift it by the reported lag before discovery. Alternatively, fix {p}' = {v} instead of "
               f"discovering it, and interpret the other equations with the reported lag in mind.")
        findings.append(Finding(
            "KINEMATIC_MISMATCH", "error" if severe else "warning", 0, ". ".join(parts) + ".",
            f"{p}' = {v}", fix,
            {"lag": fit["lag"], "lag_samples": fit["lag"] / fit["dt"], "gain": fit["gain"], "offset": fit["offset"],
             "data_share": fit["share"], "F": fit["F"], "equation_deviation": eq["deviation"],
             "extra_terms": eq["extra"], "lag_explains": eq.get("lag_explains"),
             "coefficient_changes": [[e, t, c0, c1, rel] for e, t, c0, c1, rel in impact]}))
    ctx.details["Kinematic pairs (derivative-free)"] = {
        "columns": ["pair", "lag (samples)", "lag (s)", "gain", "offset", "share of change explained", "F"],
        "rows": rows,
        "caption": "Best fit of p' = gain * v(t + lag) + offset on window increments. A positive lag means "
                   "the recorded v runs late. The share is the part of p's change that lag, gain and offset "
                   "explain beyond p' = v.",
    }
    if impact_rows:
        ctx.details["Coefficients after correcting the kinematic mismatch"] = {
            "columns": ["pair", "equation", "term", "as measured", "corrected", "relative change"],
            "rows": impact_rows,
            "caption": "Derivative-free refits of each equation's own terms, before and after replacing v by "
                       "gain * v(t + lag) + offset. Large changes mean the mismatch distorts the law, not just "
                       "the kinematic equation.",
        }
    if len(ctx.pairs) > 1:
        ctx.details["Kinematic fit per trajectory"] = {
            "columns": ["pair", "trajectory", "lag (samples)", "lag (s)", "gain", "offset", "share"],
            "rows": traj_rows,
            "caption": "A lag that changes between trajectories points to non-linear filtering rather than a "
                       "fixed delay.",
        }
    checks.append(Check("KINEMATIC_MISMATCH", 0, "flagged" if flagged else "clean"))


def _equation_side(m: Model, ctx, ip: int, iv: int, fit) -> dict:
    """How far the discovered p' is from v, along the data, and how much a lag explains."""
    vs = sp.Symbol(m.variables[iv])
    X = np.vstack([x for x, _ in ctx.pairs])
    pred = m.predict(X)
    v = X[:, iv]
    dev = pred[:, ip] - v
    scale = max(float(np.linalg.norm(v)), 1e-300)
    out = {"deviation": float(np.linalg.norm(dev) / scale), "extra": {}, "text": ""}
    active = m.active(ip)
    c_v = m.coefficients[ip, m.terms.index(vs)] if vs in m.terms else 0.0
    extra = {term_label(m.terms[j]): float(m.coefficients[ip, j]) for j in active if m.terms[j] != vs}
    out["extra"] = extra
    bits = []
    if c_v == 0:
        bits.append(f"it has no {m.variables[iv]} term")
    elif abs(c_v - 1) > 1e-3:
        bits.append(f"the {m.variables[iv]} coefficient is {c_v:.4g}")
    if extra:
        bits.append("extra terms " + ", ".join(f"{t} ({c:+.4g})" for t, c in extra.items()))
    out["text"] = "; ".join(bits) + f" (they change {m.variables[ip]}' by {out['deviation']:.2%} of {m.variables[iv]})"
    if extra and fit is not None and abs(fit["lag"]) > 0:
        predicted = fit["lag"] * pred[:, iv]          # lag * (discovered v' equation)
        ex = dev - (c_v - 1) * v                      # the extra terms' contribution
        ss = float(np.sum(ex ** 2))
        if ss > 0:
            out["lag_explains"] = float(max(0.0, 1 - np.sum((ex - predicted) ** 2) / ss))
    return out


def _corrected_pairs(ctx, iv, fit):
    out = []
    for x, t in ctx.pairs:
        y = x.copy()
        y[:, iv] = fit["gain"] * CubicSpline(t, x[:, iv])(t + fit["lag"]) + fit["offset"]
        out.append((y, t))
    return out


def _impact(m: Model, ctx, ip: int, iv: int, fit, min_share: float = 0.001) -> list:
    """Coefficient changes in every equation when v is corrected for lag, gain and offset.

    Both fits use the same derivative-free least squares on the equation's own terms, so the
    comparison isolates the correction. The kinematic equation itself is skipped. Terms whose
    contribution is below ``min_share`` of the equation's signal get a change of 0.
    """
    theta = lambda x: m.theta(x)  # noqa: E731
    A0, B0 = N.integral_targets(ctx.pairs, theta, ctx.window)
    A1, B1 = N.integral_targets(_corrected_pairs(ctx, iv, fit), theta, ctx.window)
    rows = []
    for i in range(len(m.variables)):
        if i == ip:
            continue
        act = m.active(i)
        if not act:
            continue
        c0, _ = N.lstsq(A0[:, act], B0[:, i])
        c1, _ = N.lstsq(A1[:, act], B1[:, i])
        signal = max(float(np.linalg.norm(B0[:, i])), 1e-300)
        for k, j in enumerate(act):
            rel = abs(c1[k] - c0[k]) / max(abs(c0[k]), abs(c1[k]), 1e-300)
            # a term that contributes almost nothing to the equation cannot change the law
            weight = max(abs(c0[k]), abs(c1[k])) * float(np.linalg.norm(A0[:, j])) / signal
            if weight < min_share:
                rel = 0.0
            rows.append([f"{m.variables[i]}'", term_label(m.terms[j]), float(c0[k]), float(c1[k]), float(rel)])
    rows.sort(key=lambda r: -r[4])
    return rows
