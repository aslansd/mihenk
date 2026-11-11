"""Turn evidence into findings.

Every detector reports at the highest tier its inputs allow and records a
:class:`~mihenk.report.Check` either way, so a report always says what was
checked, what was clean, and what could not be assessed and why.

Tier 0 (data + equation) never calls the discovery tool: mihenk refits only
the coefficients of the term sets the equation already has, with ordinary
least squares, and uses a derivative-free integral form as its reference.
Tier 1 calls the tool again through a ``refit`` callable to see whether the
*structure* survives resampling and reasonable analysis choices. Tier 2 needs
the true equation and only exists for synthetic data.

Thresholds are defaults calibrated on mihenk's own synthetic cases
(``mihenk validate``), not on experimental data. Treat them as starting points;
every one can be overridden.
"""

from __future__ import annotations

import itertools
from collections import Counter
from dataclasses import asdict, dataclass, fields
from typing import Callable

import numpy as np
import sympy as sp

from . import __version__
from . import kinematics as _kinematics
from . import numerics as N
from .bundle import Bundle
from .model import Model, candidate_library, evaluate_terms, format_rhs, polynomial_library, term_label
from .report import Check, Finding, Report


@dataclass
class Thresholds:
    equifinal_tol: float = 0.05          # an alternative fits "as well" within this fraction of residual
    nullspace_tol: float = 0.01          # normalized singular value below which library terms are dependent
    coef_rel_se: float = 0.1             # relative standard error that counts as 'not pinned down'
    weak_term_f: float = 10.0            # F statistic for dropping a term (integral form) below which it is weak
    deriv_rel_warning: float = 0.3       # relative gap, model vs derivative-free coefficient
    underfit_gain: float = 0.3           # a left-out term that removes this share of the residual
    underfit_f: float = 100.0            # ...with at least this F statistic
    underfit_share: float = 0.01         # ...and this share of the signal itself
    rollout_fit_r2: float = 0.9          # derivative fit counted as "good" above this
    rollout_r2: float = 0.5              # short-horizon rollout R² below this is a trap
    rollout_horizon: float = 0.1         # fraction of each test segment used for the short rollout
    bounds_factor: float = 2.0           # escape = leaving the data's box widened by this many ranges
    stats_std_ratio: float = 0.5         # long rollout std must be within this factor of the data's
    stats_mean_shift: float = 0.5        # long rollout mean shift, in data standard deviations
    n_boot: int = 30                     # tier 1 refits for structural stability
    n_blocks: int = 10                   # blocks per trajectory for the block bootstrap
    stable_freq_warning: float = 0.6     # a found term should reappear in this share of refits
    stable_freq_error: float = 0.3
    new_term_freq: float = 0.5           # a term absent from the equation but this frequent is flagged
    fork_min_share: float = 0.5          # the structure should come back under this share of variants
    gt_coef_rel: float = 0.1             # tier 2: coefficient relative error counted as wrong
    kin_max_lag: int = 10                # KINEMATIC_MISMATCH: lag search range, in samples either way
    kin_min_share: float = 0.001         # ...lag/gain/offset must explain this share of the change in p
    kin_f: float = 10.0                  # ...with at least this F statistic (no-mismatch runs stayed below 0.5)
    kin_eq_share: float = 0.005          # ...or the discovered p' must differ from v by this share of v
    kin_error_share: float = 0.01        # ...an error at this share, a warning below it
    kin_report_lag: float = 0.1          # report a lag from this many samples
    kin_report_gain: float = 0.002       # report a gain this far from 1
    kin_report_offset: float = 0.002     # report an offset from this fraction of v's RMS
    kin_coef_change: float = 0.25        # correcting the mismatch moves a coefficient this much: error
    kin_coef_min_share: float = 0.001    # ...counting only terms carrying this share of the equation

    def update(self, overrides: dict) -> "Thresholds":
        names = {f.name: f.type for f in fields(self)}
        for key, value in overrides.items():
            if key not in names:
                raise KeyError(f"unknown threshold '{key}'; known: {', '.join(names)}")
            setattr(self, key, int(value) if names[key] in (int, "int") else float(value))
        return self


Refit = Callable[[list], Model]   # list of (x, t) pairs -> Model


def scan(bundle: Bundle, model: Model, refit: Refit | None = None, variants: dict | None = None,
         thresholds: Thresholds | None = None, seed: int = 0, kinematics=None) -> Report:
    """Check a discovered equation against the data it came from.

    ``refit`` (tier 1) re-runs the discovery on a list of ``(x, t)`` pairs and
    returns a :class:`Model`; ``variants`` maps names to such callables, each a
    reasonable alternative analysis (another derivative method, threshold, ...).
    The PySINDy adapter builds both from a fitted model.

    ``kinematics`` declares derivative pairs, ``{"theta": "omega"}`` for
    theta' = omega (default: ``bundle.meta["kinematics"]``), and unlocks
    ``KINEMATIC_MISMATCH``.
    """
    th = thresholds or Thresholds()
    if list(model.variables) != list(bundle.variables):
        raise ValueError(f"model variables {model.variables} do not match bundle variables {bundle.variables}")
    ctx = _Context(bundle, model, th)
    findings: list = []
    checks: list = []
    for detector in (_check_equifinality, _check_derivative_artifact, _check_weak_terms, _check_underfit, _check_rollout):
        detector(ctx, findings, checks)
    _kinematics.check(ctx, kinematics if kinematics is not None else bundle.meta.get("kinematics"),
                      findings, checks)
    rng = np.random.default_rng(seed)
    if refit is None:
        for cid in ("UNSTABLE_STRUCTURE",):
            checks.append(Check(cid, 1, "not_assessed", "needs a refit callable (e.g. a PySINDy model)"))
    else:
        _check_stability(ctx, refit, rng, findings, checks)
    if not variants:
        checks.append(Check("FORKING_PATHS", 1, "not_assessed", "needs analysis variants (e.g. a PySINDy model)"))
    else:
        _check_forking_paths(ctx, variants, findings, checks)
    if bundle.true_equations:
        _check_ground_truth(ctx, findings, checks)
    else:
        checks.append(Check("RECOVERY", 2, "not_assessed", "no true equations (synthetic data only)"))

    n_err = sum(f.severity == "error" for f in findings if f.tier < 2)
    n_warn = sum(f.severity == "warning" for f in findings if f.tier < 2)
    if n_err:
        summary = f"{n_err} error(s), {n_warn} warning(s): do not interpret this equation as the law yet"
    elif n_warn:
        summary = f"{n_warn} warning(s): the equation may be real, but parts of it are not yet supported"
    else:
        summary = "no findings: no tested failure explains this equation (that is not proof it is the law)"
    notes = [
        f"held-out data: {ctx.split_desc}",
        f"derivative-free reference: integral form over windows of {ctx.window} samples",
        f"candidate library for alternatives: {len(ctx.library)} terms"
        + ("" if model.library is not None else " (polynomial stand-in; pass the real library for sharper checks)"),
    ]
    return Report(
        name=bundle.name, equation=model.pretty(), summary=summary, metrics=ctx.metrics,
        findings=findings, checks=checks, notes=notes, thresholds=asdict(th),
        version=__version__, source=model.source, details=ctx.details,
    )


# =================================================================== context
class _Context:
    def __init__(self, bundle: Bundle, model: Model, th: Thresholds):
        self.bundle, self.model, self.th = bundle, model, th
        self.pairs = list(zip(bundle.xs, bundle.ts))
        self.train, self.test, self.split_desc = bundle.split()
        self.library = candidate_library(model)
        self.lib_index = {t: i for i, t in enumerate(self.library)}
        nmin = min(len(t) for _, t in self.pairs)
        self.window = N.integral_window(nmin)
        self.metrics: dict = {"variables": len(model.variables), "trajectories": bundle.n_traj,
                              "samples": int(sum(len(t) for _, t in self.pairs)),
                              "active_terms": int(sum(len(model.active(i)) for i in range(len(model.variables))))}
        self.details: dict = {}
        # derivative-free targets on the full library
        theta_lib = lambda x: evaluate_terms(self.library, model.variables, x)  # noqa: E731
        self.A_int, self.B_int = N.integral_targets(self.pairs, theta_lib, self.window)
        # derivative estimates
        self.deriv = {}
        methods = N.derivative_methods(nmin)
        trim = 2
        for name, fn in methods.items():
            X, D = N.stack_derivatives(self.pairs, fn, trim=trim)
            self.deriv[name] = (X, D)
        if bundle.x_dots is not None:
            X = np.vstack([x[trim:-trim] for x in bundle.xs])
            D = np.vstack([d[trim:-trim] for d in bundle.x_dots])
            self.deriv = {"measured": (X, D), **self.deriv}
        self.ref_name = "measured" if "measured" in self.deriv else list(methods)[1]
        Xr, Dr = self.deriv[self.ref_name]
        self.theta_ref = evaluate_terms(self.library, model.variables, Xr)
        self.metrics["derivative_r2"] = N.r2(Dr, model.predict(Xr))
        model_cols = [self.lib_index[t] for t in model.terms]
        pred_int = self.A_int[:, model_cols] @ model.coefficients.T
        self.metrics["integral_r2"] = N.r2(self.B_int, pred_int)

    def cols(self, i: int) -> list:
        """Library columns used by equation ``i``."""
        return [self.lib_index[self.model.terms[j]] for j in self.model.active(i)]


def _say(term) -> str:
    """A term as it reads in a sentence."""
    return "the constant term" if term == 1 else term_label(term)


def _where(model: Model, i: int, term=None) -> str:
    s = f"{model.variables[i]}'"
    return s if term is None else f"{s}, term {term_label(term)}"


# ======================================================= DERIVATIVE_ARTIFACT
def _check_derivative_artifact(ctx: _Context, findings, checks):
    m, th = ctx.model, ctx.th
    rows = []
    flagged = False
    for i in range(len(m.variables)):
        act = m.active(i)
        if not act:
            continue
        cols = ctx.cols(i)
        flips, drifts = [], []
        c_int, _ = N.lstsq(ctx.A_int[:, cols], ctx.B_int[:, i])
        by_method = {}
        for name, (X, D) in ctx.deriv.items():
            th_m = evaluate_terms([ctx.library[c] for c in cols], m.variables, X)
            by_method[name], _ = N.lstsq(th_m, D[:, i])
        for k, j in enumerate(act):
            term = m.terms[j]
            c_model = m.coefficients[i, j]
            ci = c_int[k]
            denom = max(abs(c_model), abs(ci), 1e-12)
            rel = abs(c_model - ci) / denom
            row = [f"{m.variables[i]}'", term_label(term), c_model, ci] + [by_method[n][k] for n in ctx.deriv]
            if cols[k] in ctx.related:
                rows.append(row[:2] + ["not identified (see relations)"] + row[2:])
                continue
            rows.append(row[:2] + [""] + row[2:])
            if np.sign(c_model) != np.sign(ci) and abs(ci) > 0:
                flips.append((term, c_model, ci))
            elif rel > th.deriv_rel_warning:
                drifts.append((term, c_model, ci, rel))
        if flips or drifts:
            flagged = True
            parts = []
            if flips:
                parts.append("Without numerical differentiation the data give "
                             + ", ".join(f"{_say(t)} the opposite sign ({ci:+.4g} vs {cm:+.4g})"
                                         for t, cm, ci in flips)
                             + "; such terms are most likely produced by noise amplified while estimating "
                               "derivatives")
            if drifts:
                parts.append("Coefficients that move by more than "
                             f"{th.deriv_rel_warning:.0%} when estimated derivative-free: "
                             + ", ".join(f"{_say(t)} {cm:+.4g} -> {ci:+.4g}" for t, cm, ci, _ in drifts))
            findings.append(Finding(
                "DERIVATIVE_ARTIFACT", "error" if flips else "warning", 0,
                ". ".join(parts) + ".",
                f"{m.variables[i]}'",
                "Refit with smoothed or weak-form derivatives (e.g. PySINDy's SmoothedFiniteDifference) and "
                "report the derivative-free coefficients next to the discovered ones.",
                {"sign_flips": {str(t): [cm, ci] for t, cm, ci in flips},
                 "drifts": {str(t): [cm, ci, r] for t, cm, ci, r in drifts}}))
    ctx.details["Coefficients by estimation method"] = {
        "columns": ["equation", "term", "note", "as discovered", "derivative-free"] + list(ctx.deriv),
        "rows": rows,
        "caption": "Least-squares refits of each equation's own term set. Agreement means the coefficient "
                   "does not depend on how derivatives were estimated.",
    }
    checks.append(Check("DERIVATIVE_ARTIFACT", 0, "flagged" if flagged else "clean"))


# =============================================================== WEAK_TERM
def _check_weak_terms(ctx: _Context, findings, checks):
    m, th = ctx.model, ctx.th
    A, B = ctx.A_int, ctx.B_int
    n = A.shape[0]
    rows = []
    flagged = False
    for i in range(len(m.variables)):
        act = m.active(i)
        cols = ctx.cols(i)
        if not cols:
            continue
        _, r_full = N.lstsq(A[:, cols], B[:, i])
        weak = []
        dof = max(n - len(cols), 1)
        sigma2 = max(r_full ** 2 / dof, 1e-300)
        for k, j in enumerate(act):
            if (i, cols[k]) in ctx.absorbed or (i in ctx.equifinal_eqs and cols[k] in ctx.related):
                rows.append([f"{m.variables[i]}'", term_label(m.terms[j]), float("nan"), "see EQUIFINALITY"])
                continue
            rest = cols[:k] + cols[k + 1:]
            _, r_drop = N.lstsq(A[:, rest], B[:, i])
            F = (r_drop ** 2 - r_full ** 2) / sigma2
            share = (r_drop - r_full) / max(np.linalg.norm(B[:, i]), 1e-300)
            rows.append([f"{m.variables[i]}'", term_label(m.terms[j]), F, share])
            if F < th.weak_term_f:
                weak.append((m.terms[j], F))
        if weak:
            flagged = True
            findings.append(Finding(
                "WEAK_TERM", "warning", 0,
                "Dropping " + ", ".join(f"{_say(t)} (F = {F:.3g})" for t, F in weak)
                + f" barely changes the fit (F below {th.weak_term_f:g}): the data cannot tell the equation "
                  "with these terms from the equation without them.",
                f"{m.variables[i]}'",
                "Do not interpret these terms without more or less noisy data; a sparser equation fits as well.",
                {"terms": {str(t): F for t, F in weak}}))
    ctx.details["Term strength (derivative-free)"] = {
        "columns": ["equation", "term", "F for dropping it", "share of signal lost"],
        "rows": rows,
        "caption": "Large F means removing the term clearly hurts the fit. Windows overlap by half, so F is "
                   "a screening statistic, not an exact test.",
    }
    checks.append(Check("WEAK_TERM", 0, "flagged" if flagged else "clean"))


# ============================================================ EQUIFINALITY
def _find_relations(ctx: _Context):
    """Combinations of candidate terms that are (nearly) zero on the observed states."""
    lib, th = ctx.library, ctx.th
    Z = N.normalized_columns(ctx.theta_ref)
    _, s, Vt = np.linalg.svd(Z, full_matrices=False)  # unit-norm columns: s is the relation's relative size
    relations = []
    for sv, v in zip(s[::-1], Vt[::-1]):
        if sv >= th.nullspace_tol:
            break
        v = v / np.max(np.abs(v))
        support = [c for c in range(len(lib)) if abs(v[c]) > 0.05]
        relations.append((float(sv), {lib[c]: float(v[c]) for c in support}, set(support)))
    ctx.relations = relations
    ctx.related = set().union(*[r[2] for r in relations]) if relations else set()


def _fmt_rhs(terms, coefs) -> str:
    return format_rhs(terms, coefs)


def _check_equifinality(ctx: _Context, findings, checks):
    m, th = ctx.model, ctx.th
    lib = ctx.library
    A, B = ctx.A_int, ctx.B_int
    _find_relations(ctx)
    ctx.absorbed = set()
    ctx.equifinal_eqs = set()
    flagged = False
    used = set().union(*[set(ctx.cols(i)) for i in range(len(m.variables))])
    rel_rows = []
    for sv, rel, support in ctx.relations:
        text = " ".join(f"{c:+.3g}" + ("" if t == 1 else f"*{term_label(t)}") for t, c in rel.items()) + " ≈ 0"
        rel_rows.append([text, sv, ", ".join(term_label(lib[c]) for c in sorted(support & used)) or "none"])
    if rel_rows:
        ctx.details["Relations among candidate terms on your data"] = {
            "columns": ["relation", "size (smaller = more exact)", "equation terms involved"],
            "rows": rel_rows,
            "caption": "On the observed states these combinations are (nearly) zero, so the terms in each "
                       "relation can stand in for one another.",
        }
    alt_rows = []
    for i in range(len(m.variables)):
        cols = ctx.cols(i)
        if not cols:
            continue
        _, r0 = N.lstsq(A[:, cols], B[:, i])
        signal = max(np.linalg.norm(B[:, i]), 1e-300)
        tol = th.equifinal_tol * max(r0, 1e-6 * signal)
        alternatives = []
        # 1) sparser equivalents: the smallest subset (dropping only related terms) that fits as well
        dropped = []
        must = [c for c in cols if c not in ctx.related]
        optional = [c for c in cols if c in ctx.related]
        found = None
        for size in range(len(must), len(cols)):
            n_opt = size - len(must)
            if n_opt < 0:
                continue
            best = None
            for combo in itertools.islice(itertools.combinations(optional, n_opt), 5000):
                trial = [c for c in cols if c in must or c in combo]
                if not trial:
                    continue
                coef, r = N.lstsq(A[:, trial], B[:, i])
                if r <= r0 + tol and (best is None or r < best[2]):
                    best = (trial, coef, r)
            if best is not None:
                found = best
                break
        if found is not None:
            cur, coef, r = found
            dropped = [c for c in cols if c not in cur]
            alternatives.append(("sparser", [lib[c] for c in cur], coef, r))
            ctx.absorbed |= {(i, c) for c in dropped}
        # 2) one-term swaps that keep the size and the fit
        for k, ck in enumerate(cols):
            if ck in dropped:
                continue
            for cj in range(len(lib)):
                if cj in cols or not (ck in ctx.related or cj in ctx.related):
                    continue
                trial = cols[:k] + [cj] + cols[k + 1:]
                coef, r = N.lstsq(A[:, trial], B[:, i])
                contrib = abs(coef[k]) * np.linalg.norm(A[:, cj]) / signal
                if r <= r0 + tol and contrib > 0.01:
                    alternatives.append(("swap", [lib[c] for c in trial], coef, r))
                    ctx.absorbed.add((i, ck))
        if not alternatives:
            continue
        flagged = True
        ctx.equifinal_eqs.add(i)
        alternatives.sort(key=lambda a: (len(a[1]), a[3]))
        shown = alternatives[:3]
        lines = [f"{m.variables[i]}' = {_fmt_rhs(terms, coef)}" for _, terms, coef, _ in shown]
        for kind, terms, coef, r in alternatives:
            alt_rows.append([f"{m.variables[i]}'", kind, _fmt_rhs(terms, coef), r0, r])
        sparsest = alternatives[0]
        n_s = len(sparsest[1])
        head = (f"A sparser equation fits your data exactly as well ({n_s} term{'s' if n_s != 1 else ''} "
                f"instead of {len(cols)})" if sparsest[0] == "sparser" else "A different equation of the same size fits your "
                "data exactly as well")
        findings.append(Finding(
            "EQUIFINALITY", "error", 0,
            head + ": " + "; or ".join(lines)
            + (f" (and {len(alternatives) - 3} more)" if len(alternatives) > 3 else "")
            + ". On the observed states some candidate terms are combinations of others, so the data cannot "
              "choose between these laws.",
            f"{m.variables[i]}'",
            "Collect data that leaves the current set of states (other initial conditions, transients, "
            "perturbations) so the alternatives make different predictions.",
            {"alternatives": [{"kind": k, "terms": [str(t) for t in ts], "coefficients": list(map(float, c)),
                               "residual": r} for k, ts, c, r in alternatives], "residual_original": r0}))
    if alt_rows:
        ctx.details["Equally good alternative equations"] = {
            "columns": ["equation", "kind", "alternative right-hand side", "residual (yours)", "residual (alt.)"],
            "rows": alt_rows,
            "caption": "Derivative-free least-squares fits on the same data.",
        }
    if not flagged and ctx.relations:
        # a near relation matters only if it leaves the equation's own coefficients loosely determined
        loose = []
        n = A.shape[0]
        for i in range(len(m.variables)):
            cols = ctx.cols(i)
            touched = [k for k, c in enumerate(cols) if c in ctx.related]
            if not touched:
                continue
            coef, r = N.lstsq(A[:, cols], B[:, i])
            sigma2 = r ** 2 / max(n - len(cols), 1)
            try:
                cov = sigma2 * np.linalg.pinv(A[:, cols].T @ A[:, cols])
            except np.linalg.LinAlgError:  # pragma: no cover
                continue
            for k in touched:
                rse = np.sqrt(max(cov[k, k], 0)) / max(abs(coef[k]), 1e-300)
                if rse > th.coef_rel_se:
                    loose.append((i, lib[cols[k]], rse))
        if loose:
            flagged = True
            sv = ctx.relations[0][0]
            findings.append(Finding(
                "EQUIFINALITY", "warning", 0,
                f"Your data lie close to a curve or surface on which some candidate terms nearly cancel "
                f"(relation size {sv:.2g}), and at this noise level that leaves coefficients loosely determined: "
                + ", ".join(f"{m.variables[i]}' {_say(t)} ±{rse:.0%}" for i, t, rse in loose)
                + " (relative standard error, derivative-free fit).",
                "library", "See the relations table; data off this set of states would separate the terms.",
                {"loose": [[m.variables[i], str(t), rse] for i, t, rse in loose]}))
    ctx.metrics["library_relations"] = len(ctx.relations)
    checks.append(Check("EQUIFINALITY", 0, "flagged" if flagged else "clean"))


# ================================================================ UNDERFIT
def _check_underfit(ctx: _Context, findings, checks):
    """Does a term the equation left out explain a large share of what it misses?"""
    m, th = ctx.model, ctx.th
    lib = ctx.library
    deg = max(2, m.max_degree())
    extra = [t for t in polynomial_library(m.variables, deg + 1) if t not in ctx.lib_index]
    if extra:
        theta_x = lambda x: evaluate_terms(extra, m.variables, x)  # noqa: E731
        A_x, _ = N.integral_targets(ctx.pairs, theta_x, ctx.window)
        A_all = np.hstack([ctx.A_int, A_x])
    else:
        A_all = ctx.A_int
    all_terms = lib + extra
    n = A_all.shape[0]
    flagged = False
    rows = []
    for i in range(len(m.variables)):
        cols = ctx.cols(i)
        _, r0 = N.lstsq(A_all[:, cols], ctx.B_int[:, i])
        signal = max(np.linalg.norm(ctx.B_int[:, i]), 1e-300)
        if r0 <= 1e-6 * signal:
            continue  # already exact
        best = None
        for c in range(len(all_terms)):
            if c in cols or (c in ctx.related) or np.ptp(A_all[:, c]) == 0 and all_terms[c] != 1:
                continue
            _, r = N.lstsq(A_all[:, cols + [c]], ctx.B_int[:, i])
            gain = 1 - r / r0
            F = (r0 ** 2 - r ** 2) / max(r ** 2 / max(n - len(cols) - 1, 1), 1e-300)
            if best is None or gain > best[1]:
                best = (c, gain, F)
        if best is None:
            continue
        c, gain, F = best
        in_lib = c < len(lib)
        rows.append([f"{m.variables[i]}'", term_label(all_terms[c]), "yes" if in_lib else "no", gain, F])
        share = (r0 - (1 - gain) * r0) / signal
        rows[-1].insert(4, share)
        if gain >= th.underfit_gain and F >= th.underfit_f and share >= th.underfit_share:
            flagged = True
            where_from = ("It is in the candidate library, so the sparsity step removed it"
                          if in_lib else "It is not in the candidate library at all")
            findings.append(Finding(
                "UNDERFIT", "error" if gain >= 2 * th.underfit_gain else "warning", 0,
                f"Adding {_say(all_terms[c])} cuts this equation's derivative-free residual by {gain:.0%} "
                f"(F = {F:.3g}). {where_from}; the other coefficients may be compensating for it.",
                _where(m, i, all_terms[c]),
                "Refit with a lower threshold or a richer library and check whether the coefficients change."
                if in_lib else "Refit with a library that contains this term (e.g. one degree higher).",
                {"term": str(all_terms[c]), "residual_gain": gain, "F": F, "in_library": in_lib}))
    ctx.details["Best left-out term per equation"] = {
        "columns": ["equation", "term", "in library", "residual reduction", "share of signal", "F"],
        "rows": rows,
        "caption": "For each equation, the single term not in it that most reduces the derivative-free residual "
                   "(searching the library and one polynomial degree beyond it).",
    }
    checks.append(Check("UNDERFIT", 0, "flagged" if flagged else "clean"))


# ============================================================ ROLLOUT_TRAP
def _check_rollout(ctx: _Context, findings, checks):
    m, th = ctx.model, ctx.th
    allx = np.vstack([x for x, _ in ctx.pairs])
    lo, hi = allx.min(axis=0), allx.max(axis=0)
    span = np.where(hi > lo, hi - lo, 1.0)
    bounds = (lo - th.bounds_factor * span, hi + th.bounds_factor * span)
    var_all = allx.var(axis=0)  # rollouts are scored on the scale of the whole data set
    flagged = False
    short_r2, statuses, rows = [], [], []
    for k, (x, t) in enumerate(ctx.test):
        h = max(20, int(th.rollout_horizon * len(t)))
        h = min(h, len(t))
        y, status = N.rollout(m, x[0], t[:h], bounds)
        statuses.append(status)
        ok = np.all(np.isfinite(y), axis=1)
        r = N.r2(x[:h][ok], y[ok], var=var_all) if ok.sum() > 2 else float("nan")
        short_r2.append(r)
        rows.append([f"test segment {k}", "short", h, status, r])
    # long rollouts from the start of each trajectory
    long_status, stats_bad = [], []
    for k, (x, t) in enumerate(ctx.pairs):
        y, status = N.rollout(m, x[0], t, bounds)
        long_status.append(status)
        ok = np.all(np.isfinite(y), axis=1)
        if status == "ok":
            sd_d, sd_m = x.std(axis=0), y.std(axis=0)
            mean_shift = np.abs(y.mean(axis=0) - x.mean(axis=0)) / np.where(sd_d > 0, sd_d, 1)
            ratio = np.where(sd_d > 0, sd_m / np.where(sd_d > 0, sd_d, 1), 1)
            bad = (ratio < th.stats_std_ratio) | (ratio > 1 / th.stats_std_ratio) | (mean_shift > th.stats_mean_shift)
            if bad.any():
                stats_bad.append((k, [m.variables[q] for q in np.where(bad)[0]], ratio, mean_shift))
            rows.append([f"trajectory {k}", "full", len(t), status, N.r2(x[ok], y[ok], var=var_all)])
        else:
            rows.append([f"trajectory {k}", "full", len(t), status, float("nan")])
    ctx.details["Rollouts (integrating the equation)"] = {
        "columns": ["segment", "horizon", "samples", "status", "R²"],
        "rows": rows,
        "caption": "Short rollouts start each held-out segment; full rollouts integrate each trajectory from its "
                   "first sample. For chaotic systems a low full-length R² is expected; the long-run statistics "
                   "check compares spread and mean instead.",
    }
    r2s = [r for r in short_r2 if r == r]
    ctx.metrics["rollout_r2_short"] = float(np.median(r2s)) if r2s else float("nan")
    escaped = [s for s in statuses + long_status if s != "ok"]
    if escaped:
        flagged = True
        findings.append(Finding(
            "ROLLOUT_TRAP", "error", 0,
            f"Integrating the equation does not stay with the data: {len(escaped)} of "
            f"{len(statuses) + len(long_status)} rollouts {', '.join(sorted(set(escaped)))} "
            f"(left the data's range widened by {th.bounds_factor:g}x, or the integrator failed). "
            f"Fitting derivatives well does not make an equation that reproduces the system.",
            "rollout",
            "Check for missing stabilizing terms or a wrong library; consider constrained/stable SINDy "
            "(e.g. TrappingSR3) and always validate by simulation.",
            {"statuses": statuses + long_status}))
    elif r2s and ctx.metrics["derivative_r2"] >= th.rollout_fit_r2 and min(r2s) < th.rollout_r2:
        flagged = True
        findings.append(Finding(
            "ROLLOUT_TRAP", "warning", 0,
            f"The equation fits derivatives well (R² {ctx.metrics['derivative_r2']:.3f}) but a short simulation "
            f"from held-out initial states does not follow the data (worst R² {min(r2s):.3f}).",
            "rollout", "Judge the equation by simulation, not by its derivative fit.",
            {"rollout_r2": short_r2}))
    if stats_bad and not escaped:
        flagged = True
        k, vars_, ratio, shift = stats_bad[0]
        findings.append(Finding(
            "ROLLOUT_TRAP", "warning", 0,
            f"Simulated over the full trajectory, the equation settles into different behaviour from the data "
            f"(variables {', '.join(vars_)}: spread ratio {', '.join(f'{r:.2g}' for r in ratio)}, mean shift "
            f"{', '.join(f'{s:.2g}' for s in shift)} SD). It may decay, grow or orbit where the system does not.",
            f"trajectory {k}", "Compare simulated and observed long-run behaviour before interpreting the equation.",
            {"trajectories": [b[0] for b in stats_bad]}))
    checks.append(Check("ROLLOUT_TRAP", 0, "flagged" if flagged else "clean"))


# ================================================ tier 1: UNSTABLE_STRUCTURE
def _resample(ctx: _Context, rng):
    """Bootstrap: whole trajectories when there are several, else blocks of each trajectory.

    Blocks keep the amount and spread of data the same as the original, so the
    check asks "is the structure stable under resampling?", not "would a shorter
    recording give the same law?".
    """
    pairs = ctx.pairs
    if len(pairs) >= 3:
        idx = rng.integers(0, len(pairs), len(pairs))
        return [pairs[i] for i in idx]
    blocks = []
    for x, t in pairs:
        n_blocks = int(max(2, min(ctx.th.n_blocks, len(t) // 40)))
        edges = np.linspace(0, len(t), n_blocks + 1).astype(int)
        blocks += [(x[a:b], t[a:b]) for a, b in zip(edges[:-1], edges[1:])]
    idx = rng.integers(0, len(blocks), len(blocks))
    return [blocks[i] for i in idx]


def _check_stability(ctx: _Context, refit: Refit, rng, findings, checks):
    m, th = ctx.model, ctx.th
    counts = [Counter() for _ in m.variables]
    exact, ok_runs, errors = 0, 0, []
    original = m.structure()
    for _ in range(th.n_boot):
        try:
            fit = refit(_resample(ctx, rng))
        except Exception as exc:  # noqa: BLE001 - a failing refit is evidence, not a crash
            errors.append(str(exc))
            continue
        ok_runs += 1
        s = fit.structure()
        exact += s == original
        for i, ts in enumerate(s):
            counts[i].update(ts)
    if ok_runs == 0:
        checks.append(Check("UNSTABLE_STRUCTURE", 1, "not_assessed", f"every refit failed: {errors[0][:120]}"))
        return
    flagged = False
    rows = []
    for i, v in enumerate(m.variables):
        orig = original[i]
        shaky, extra = [], []
        for term in sorted(set(counts[i]) | set(orig), key=str):
            freq = counts[i][term] / ok_runs
            rows.append([f"{v}'", term_label(term), "yes" if term in orig else "no", freq])
            if term in orig and freq < th.stable_freq_warning:
                shaky.append((term, freq))
            elif term not in orig and freq >= th.new_term_freq:
                extra.append((term, freq))
        if not (shaky or extra):
            continue
        flagged = True
        parts = []
        if shaky:
            parts.append("Terms in the equation that often disappear: "
                         + ", ".join(f"{_say(t)} ({f:.0%})" for t, f in shaky))
        if extra:
            parts.append("Terms absent from the equation that often appear: "
                         + ", ".join(f"{_say(t)} ({f:.0%})" for t, f in extra))
        sev = "error" if any(f < th.stable_freq_error for _, f in shaky) else "warning"
        findings.append(Finding(
            "UNSTABLE_STRUCTURE", sev, 1,
            f"Refitting on {ok_runs} resamples of your data changes this equation's terms. "
            + ". ".join(parts) + " (share of refits). The structure depends on which part of the data "
            "the discovery saw.",
            f"{v}'",
            "Treat these terms as tentative; check them on independent trajectories, or with more or "
            "less noisy data.",
            {"disappearing": {str(t): f for t, f in shaky}, "appearing": {str(t): f for t, f in extra},
             "refits": ok_runs}))
    ctx.metrics["structure_reproduced"] = exact / ok_runs
    ctx.details["Term frequency across resampled refits"] = {
        "columns": ["equation", "term", "in equation", "share of refits containing it"],
        "rows": rows,
        "caption": f"{ok_runs} refits through the discovery tool on resampled data; "
                   f"the exact structure came back in {exact} of them.",
    }
    reason = f"{len(errors)} of {th.n_boot} refits failed" if errors else ""
    checks.append(Check("UNSTABLE_STRUCTURE", 1, "flagged" if flagged else "clean", reason))


# ===================================================== tier 1: FORKING_PATHS
def _check_forking_paths(ctx: _Context, variants: dict, findings, checks):
    m, th = ctx.model, ctx.th
    original = m.structure()
    rows, same, ok_runs = [], 0, 0
    distinct = Counter()
    for name, fn in variants.items():
        try:
            fit = fn(ctx.pairs)
        except Exception as exc:  # noqa: BLE001
            rows.append([name, "failed", str(exc)[:80]])
            continue
        ok_runs += 1
        s = fit.structure()
        distinct[s] += 1
        same += s == original
        diff = _structure_diff(original, s, m.variables)
        rows.append([name, "same" if s == original else "different", diff or "-"])
    if ok_runs == 0:
        checks.append(Check("FORKING_PATHS", 1, "not_assessed", "every variant failed"))
        return
    share = same / ok_runs
    ctx.metrics["structure_under_variants"] = share
    ctx.metrics["distinct_structures"] = len(distinct)
    ctx.details["Structure under analysis variants"] = {
        "columns": ["variant", "structure", "difference from the equation"],
        "rows": rows,
        "caption": "Each variant is a defensible alternative analysis run on the same data.",
    }
    flagged = share < th.fork_min_share
    if flagged:
        findings.append(Finding(
            "FORKING_PATHS", "warning", 1,
            f"The equation's structure comes back under only {same} of {ok_runs} defensible analysis variants "
            f"({len(distinct)} distinct structures in total). Which law you report depends on analysis choices.",
            "analysis", "Report the variants, or pick the analysis on held-out data before looking at the equation.",
            {"share": share, "distinct": len(distinct)}))
    checks.append(Check("FORKING_PATHS", 1, "flagged" if flagged else "clean"))


def _structure_diff(a, b, variables) -> str:
    parts = []
    for v, sa, sb in zip(variables, a, b):
        plus = sorted(map(term_label, sb - sa))
        minus = sorted(map(term_label, sa - sb))
        if plus or minus:
            parts.append(f"{v}': " + " ".join([f"+{p}" for p in plus] + [f"-{q}" for q in minus]))
    return "; ".join(parts)


# ===================================================== tier 2: ground truth
def _check_ground_truth(ctx: _Context, findings, checks):
    m, th = ctx.model, ctx.th
    truth = Model.from_equations(ctx.bundle.true_equations, variables=m.variables)
    flagged = False
    rows = []
    worst = 0.0
    for i, v in enumerate(m.variables):
        t_coef = {truth.terms[j]: truth.coefficients[i, j] for j in truth.active(i)}
        m_coef = {m.terms[j]: m.coefficients[i, j] for j in m.active(i)}
        for term in sorted(set(t_coef) | set(m_coef), key=str):
            ct, cm = t_coef.get(term, 0.0), m_coef.get(term, 0.0)
            rows.append([f"{v}'", term_label(term), ct, cm])
            if term in t_coef and term not in m_coef:
                flagged = True
                findings.append(Finding("MISSED_TERM", "error", 2, f"{v}': true term {term_label(term)} "
                                        f"({ct:+.4g}) is missing.", _where(m, i, term)))
            elif term in m_coef and term not in t_coef:
                flagged = True
                findings.append(Finding("SPURIOUS_TERM", "error", 2, f"{v}': term {term_label(term)} "
                                        f"({cm:+.4g}) is not in the true equation.", _where(m, i, term)))
            else:
                rel = abs(cm - ct) / abs(ct)
                worst = max(worst, rel)
                if rel > th.gt_coef_rel:
                    flagged = True
                    findings.append(Finding("COEF_ERROR", "warning", 2, f"{v}': {term_label(term)} has "
                                            f"{cm:+.4g}, true {ct:+.4g} ({rel:.0%} off).", _where(m, i, term)))
    ctx.metrics["gt_structure_correct"] = m.structure() == truth.structure()
    ctx.metrics["gt_worst_coef_error"] = worst
    ctx.details["Ground truth comparison"] = {"columns": ["equation", "term", "true", "discovered"], "rows": rows}
    checks.append(Check("RECOVERY", 2, "flagged" if flagged else "clean"))
