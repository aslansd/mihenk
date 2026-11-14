"""PHYSICS_VIOLATION: does the discovered equation respect physics you declare?

mihenk cannot know that a free pendulum must lose energy. You can tell it::

    properties={"dissipates": "0.5*omega**2 - 64*cos(theta)"}   # H must not increase
    properties={"conserves":  "0.5*x - log(x) + 0.5*y - log(y)"} # H must stay constant

H is an expression in the state variables. mihenk then compares three rates.

* **Data.** The observed change of H, computed derivative-free from H at the
  two ends of each window (no differentiation of noisy states). Averaged over
  time, this is the rate at which the data lose (or keep) H.
* **Equation, along the data.** ``dH/dt = grad H . f(x)``, using the discovered
  right-hand side ``f`` at the observed states.
* **Equation, along its own simulations.** H along full rollouts from each
  trajectory's first state.

For ``dissipates``:

* an equation that *gains* H on average along the data, or whose simulations
  gain H, is an **error**;
* one that loses H at well under, or well over, the rate the data show is a
  **warning**. For example, the damping was cut in half by a velocity lag.

For ``conserves``:

* an equation whose dH/dt along the data, or whose simulated drift of H, is
  large compared with the size of H's own terms is flagged.

H needs to be only roughly right for ``dissipates``. A small error in its
conservative part (64 vs 63.8) adds a term that averages out over a swing.
For ``conserves``, H must be accurate.
"""

from __future__ import annotations

import numpy as np
import sympy as sp

from . import numerics as N
from .model import Model, parse_term
from .report import Check, Finding

KINDS = ("dissipates", "conserves")


def validate_properties(properties, variables) -> dict:
    """``{"dissipates": "H"}`` / ``["dissipates=H"]`` -> ``{kind: sympy expr}``."""
    if not properties:
        return {}
    if isinstance(properties, (list, tuple)):
        d = {}
        for item in properties:
            if "=" not in item:
                raise ValueError(f"property must look like dissipates=EXPR, got {item!r}")
            k, v = item.split("=", 1)
            d[k.strip()] = v.strip()
        properties = d
    out = {}
    for kind, expr in dict(properties).items():
        if kind not in KINDS:
            raise ValueError(f"unknown property '{kind}'; use one of: {', '.join(KINDS)}")
        e = parse_term(expr, variables) if isinstance(expr, str) else sp.sympify(expr)
        unknown = {str(s) for s in e.free_symbols} - set(variables)
        if unknown:
            raise ValueError(f"{kind}: H uses unknown symbol(s) {', '.join(sorted(unknown))}; "
                             f"H may use only the variables {', '.join(variables)}")
        if not e.free_symbols:
            raise ValueError(f"{kind}: H is a constant")
        out[kind] = e
    return out


def _functions(H, variables):
    syms = [sp.Symbol(v) for v in variables]
    h = sp.lambdify(syms, H, modules="numpy")
    grad = [sp.lambdify(syms, sp.diff(H, s), modules="numpy") for s in syms]

    def H_of(X):
        with np.errstate(all="ignore"):
            v = np.asarray(h(*np.atleast_2d(X).T), float)
        return np.broadcast_to(v, (np.atleast_2d(X).shape[0],)) if v.ndim == 0 else v

    def grad_of(X):
        X = np.atleast_2d(X)
        cols = []
        for g in grad:
            with np.errstate(all="ignore"):
                v = np.asarray(g(*X.T), float)
            cols.append(np.broadcast_to(v, (X.shape[0],)) if v.ndim == 0 else v)
        return np.column_stack(cols)

    return H_of, grad_of


def rates(model: Model, pairs, H, window: int) -> dict:
    """Observed and predicted rates of change of H along the data, plus scales."""
    H_of, grad_of = _functions(H, model.variables)
    dH_obs, T, P_model, P_abs, n = 0.0, 0.0, [], [], 0
    for x, t in pairs:
        h = H_of(x)
        w = min(window, len(t) - 1)
        starts = np.arange(0, len(t) - w, w)               # non-overlapping windows
        dH_obs += float(np.sum(h[starts + w] - h[starts]))
        T += float(np.sum(t[starts + w] - t[starts]))
        g = grad_of(x)
        f = model.predict(x)
        p = np.sum(g * f, axis=1)                          # dH/dt the equation predicts
        P_model.append(p)
        P_abs.append(np.sum(np.abs(g * f), axis=1))        # size of H's terms, for scale
        n += len(t)
    P_model = np.concatenate(P_model)
    P_abs = np.concatenate(P_abs)
    return {"P_data": dH_obs / max(T, 1e-300), "P_model": float(np.mean(P_model)),
            "P_scale": float(np.mean(P_abs)), "P_model_rms": float(np.sqrt(np.mean(P_model ** 2))),
            "H_range": float(np.ptp(np.concatenate([H_of(x) for x, _ in pairs])))}


def natural_scale(model: Model, pairs, H) -> float:
    """How much H changes when each variable moves through its observed range.

    sum_j mean|dH/dx_j| * range(x_j). Unlike H's own range in the data, this
    is not zero when the data conserve H perfectly.
    """
    _, grad_of = _functions(H, model.variables)
    X = np.vstack([x for x, _ in pairs])
    return float(np.sum(np.mean(np.abs(grad_of(X)), axis=0) * np.ptp(X, axis=0))) or 1.0


def rollout_drift(model: Model, pairs, H, bounds) -> dict:
    """Largest gain and largest absolute drift of H along full rollouts, relative to H's natural scale."""
    H_of, _ = _functions(H, model.variables)
    span = natural_scale(model, pairs, H)
    gains, drifts, statuses = [], [], []
    for x, t in pairs:
        y, status = N.rollout(model, x[0], t, bounds)
        statuses.append(status)
        ok = np.all(np.isfinite(y), axis=1)
        if ok.sum() < 2:
            continue
        h = H_of(y[ok])
        gains.append(float((np.max(h) - h[0]) / span))
        drifts.append(float(np.max(np.abs(h - h[0])) / span))
    return {"gain": max(gains, default=float("nan")), "drift": max(drifts, default=float("nan")),
            "statuses": statuses}


def check(ctx, properties, findings, checks):
    m, th = ctx.model, ctx.th
    props = validate_properties(properties, m.variables)
    if not props:
        checks.append(Check("PHYSICS_VIOLATION", 0, "not_assessed",
                            "no physical properties declared (e.g. dissipates=H or conserves=H)"))
        return
    allx = np.vstack([x for x, _ in ctx.pairs])
    lo, hi = allx.min(axis=0), allx.max(axis=0)
    span = np.where(hi > lo, hi - lo, 1.0)
    bounds = (lo - th.bounds_factor * span, hi + th.bounds_factor * span)
    flagged, rows, notes = False, [], []
    for kind, H in props.items():
        r = rates(m, ctx.pairs, H, ctx.window)
        d = rollout_drift(m, ctx.pairs, H, bounds)
        rows.append([kind, str(H), r["P_data"], r["P_model"], r["P_scale"], d["gain"], d["drift"]])
        ctx.metrics[f"phys_{kind}_P_data"] = r["P_data"]
        ctx.metrics[f"phys_{kind}_P_model"] = r["P_model"]
        problems, severe = [], False
        if kind == "dissipates":
            loss = -r["P_data"]
            if loss <= th.phys_data_min * r["P_scale"]:
                notes.append(f"dissipates: the data barely lose H ({r['P_data']:+.3g} per unit time), so the "
                             "dissipation check has little to compare against")
            ref = max(loss, th.phys_data_min * r["P_scale"])
            if r["P_model"] > th.phys_inject * ref:
                severe = True
                problems.append(f"along the observed states the equation *gains* H at {r['P_model']:+.3g} per unit "
                                f"time, while the data lose it at {r['P_data']:+.3g}")
            elif loss > th.phys_data_min * r["P_scale"]:
                ratio = -r["P_model"] / loss
                if ratio < th.phys_ratio_low or ratio > th.phys_ratio_high:
                    problems.append(f"the equation loses H at {ratio:.0%} of the rate the data show "
                                    f"({r['P_model']:+.3g} vs {r['P_data']:+.3g} per unit time)")
            if d["gain"] == d["gain"] and d["gain"] > th.phys_rollout_gain:
                severe = True
                problems.append(f"simulated, the equation raises H by {d['gain']:.0%} of H's scale on the data")
        else:  # conserves
            rel = r["P_model_rms"] / max(r["P_scale"], 1e-300)
            if rel > th.phys_conserve_rel:
                problems.append(f"along the observed states the equation changes H at {rel:.1%} of the size of "
                                f"H's own terms (rms of dH/dt relative to rms of its parts)")
            if d["drift"] == d["drift"] and d["drift"] > th.phys_conserve_drift:
                severe = d["drift"] > 5 * th.phys_conserve_drift
                problems.append(f"simulated, H drifts by {d['drift']:.1%} of H's scale on the data")
            data_drift = abs(r["P_data"]) / max(r["P_scale"], 1e-300)
            if data_drift > th.phys_conserve_rel:
                notes.append(f"conserves: the data themselves change H at {data_drift:.1%} of its terms' size; "
                             "check the declaration or the energy function")
        if problems:
            flagged = True
            findings.append(Finding(
                "PHYSICS_VIOLATION", "error" if severe else "warning", 0,
                f"You declared that the system {kind} H = {H}. " + "; ".join(problems).capitalize() + ".",
                f"{kind} {H}",
                "Look for terms that add energy where the physics removes it, or a damping term that is "
                "missing or biased (a velocity sensor lag lowers apparent damping; see KINEMATIC_MISMATCH). "
                "Constrain the discovery (e.g. a library or sign constraints) rather than reporting this law.",
                {"P_data": r["P_data"], "P_model": r["P_model"], "P_scale": r["P_scale"],
                 "rollout_gain": d["gain"], "rollout_drift": d["drift"]}))
    ctx.details["Declared physical properties"] = {
        "columns": ["property", "H", "data dH/dt", "equation dH/dt", "size of H's terms", "rollout gain of H",
                    "rollout drift of H"],
        "rows": rows,
        "caption": "dH/dt is averaged over time. The data's rate is derivative-free (H at window ends); the "
                   "equation's is grad H . f at the observed states. Rollout gain and drift are relative to "
                   "H's scale on the data: how much H changes when each variable moves through its range.",
    }
    reason = "; ".join(notes)
    checks.append(Check("PHYSICS_VIOLATION", 0, "flagged" if flagged else "clean", reason))
