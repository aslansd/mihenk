"""Self-validation: does mihenk flag the failures it claims to, and stay quiet otherwise?

Each case simulates a system with a known law, runs a discovery (PySINDy, or a
hand-written equation for cases no fit produces reliably), scans it, and
compares mihenk's tier-0/1 flags with what the case is built to show. Ground
truth (tier 2) is reported alongside so you can see when the discovered law
was actually wrong.

The cases are toy problems. Passing them shows the detectors work as designed,
not that the thresholds are calibrated for your data.
"""

from __future__ import annotations

import time
from dataclasses import dataclass, field

from . import synth
from .model import Model
from .scan import scan


@dataclass
class Case:
    name: str
    system: str
    expect: set                     # flags that must be raised (tier <= 1); empty set = must be clean
    noise: float = 0.0
    threshold: float = 0.1
    degree: int | None = None
    diff: str = "fd"
    n_traj: int = 1
    equation: str | None = None     # a hand-written equation instead of a PySINDy fit
    why: str = ""
    quick: bool = True
    any_of: bool = False            # pass if any expected flag is raised (noisy cases vary by seed)
    forbid: set = field(default_factory=set)   # flags that must not be raised
    sensor: dict | None = None      # imperfect measurement, passed to synth.simulate


CASES = [
    Case("lorenz_clean", "lorenz", set(), threshold=0.5, why="the textbook success; must be clean"),
    Case("lorenz_noise1", "lorenz", set(), noise=0.01, threshold=0.5,
         why="mild noise that SINDy handles; must stay clean"),
    Case("van_der_pol_clean", "van_der_pol", set(), threshold=0.05, why="limit cycle that is not a conic"),
    Case("duffing_clean", "duffing", set(), threshold=0.05, why="cubic system with the right library"),
    Case("linear_clean", "linear_damped", set(), threshold=0.05, why="linear oscillator"),
    Case("lotka_volterra_many", "lotka_volterra", set(), threshold=0.05, n_traj=4,
         why="several orbits separate the polynomial terms"),
    Case("stuart_landau_many", "stuart_landau", set(), threshold=0.05, n_traj=6,
         why="the remedy for the on-cycle case: transients from several starts"),
    Case("stuart_landau_on_cycle", "stuart_landau_on_cycle", {"EQUIFINALITY"}, threshold=0.05,
         why="data confined to x^2 + y^2 = 1: several laws fit exactly"),
    Case("lotka_volterra_one_orbit", "lotka_volterra", set(), threshold=0.05,
         why="one orbit lies near an ellipse, but without noise the coefficients are still pinned down"),
    Case("lotka_volterra_one_orbit_noisy", "lotka_volterra", {"UNSTABLE_STRUCTURE", "EQUIFINALITY"},
         noise=0.02, threshold=0.05, any_of=True,
         why="the same orbit with noise: the near-ellipse lets terms come and go"),
    Case("duffing_small_library", "duffing", {"UNDERFIT"}, threshold=0.05, degree=2,
         why="the cubic term is missing from the library"),
    Case("lorenz_noise5_low_threshold", "lorenz", {"WEAK_TERM", "UNSTABLE_STRUCTURE", "DERIVATIVE_ARTIFACT"},
         noise=0.05, threshold=0.1, any_of=True, why="a spurious constant from differentiated noise"),
    Case("van_der_pol_noise10", "van_der_pol", {"DERIVATIVE_ARTIFACT", "WEAK_TERM"},
         noise=0.1, threshold=0.01, why="heavy noise, dense overfit equation: both must fire", quick=False),
    Case("pendulum_kinematics_clean", "pendulum", set(), threshold=0.02, forbid={"KINEMATIC_MISMATCH"},
         why="declared theta' = omega with an exact sensor: the kinematic check must stay quiet"),
    Case("pendulum_velocity_lag", "pendulum", {"KINEMATIC_MISMATCH"}, threshold=0.02,
         sensor={"omega": {"delay": 7e-4}},
         why="omega recorded 0.7 ms late, as in the Brunton-lab experiment: extra theta' terms, damping halved"),
    Case("pendulum_velocity_gain", "pendulum", {"KINEMATIC_MISMATCH"}, threshold=0.02,
         sensor={"omega": {"gain": 1.02}}, why="omega recorded 2% too large (a scale error)"),
    Case("pendulum_velocity_offset", "pendulum", {"KINEMATIC_MISMATCH"}, threshold=0.02,
         sensor={"omega": {"offset": 0.05}}, why="omega recorded with a bias of 0.05 rad/s"),
    Case("linear_wrong_sign", "linear_damped", {"ROLLOUT_TRAP"},
         equation="x' = 0.1 x + 2 y\ny' = -2 x + 0.1 y",
         why="hand-edited damping sign: close derivatives, but the simulation spirals out"),
]


@dataclass
class CaseResult:
    case: Case
    flagged: set
    errors: set
    gt_correct: bool | None
    passed: bool
    seconds: float
    report: object = field(repr=False, default=None)

    @property
    def false_warning(self) -> bool:
        """A clean case that raised warnings (but no errors)."""
        return not self.case.expect and self.gt_correct is not False and bool(self.flagged) and not self.errors


def run_case(case: Case, seed: int = 0) -> CaseResult:
    t0 = time.time()
    b = synth.simulate(case.system, noise=case.noise, n_traj=case.n_traj, seed=seed, sensor=case.sensor)
    b.meta["name"] = case.name
    refit = variants = None
    if case.equation is not None:
        model = Model.from_text(case.equation, variables=b.variables, source="hand-written")
    else:
        from .adapters import pysindy as A

        degree = case.degree or synth.SYSTEMS[case.system].degree
        sindy, model = A.fit_bundle(b, degree=degree, threshold=case.threshold, diff=case.diff,
                                    fourier=synth.SYSTEMS[case.system].fourier)
        refit = A.make_refit(sindy, b.variables)
        variants = A.make_variants(sindy, b.variables)
    report = scan(b, model, refit=refit, variants=variants, seed=seed)
    flagged = report.flagged(tier_max=1)
    errors = {f.id for f in report.findings if f.tier < 2 and f.severity == "error"}
    gt = report.metrics.get("gt_structure_correct")
    if gt is False and not case.expect:
        passed = bool(flagged)                   # the fit went wrong on a "clean" case: it must be caught
    elif not case.expect:
        passed = not errors                      # a correct law must never get an error
    elif case.any_of:
        passed = bool(case.expect & flagged)
    else:
        passed = case.expect <= flagged
    if gt is False and not flagged:
        passed = False                           # a wrong law that nothing flagged is the worst outcome
    if case.forbid & flagged:
        passed = False
    return CaseResult(case, flagged, errors, gt, passed, time.time() - t0, report)


def run_all(quick: bool = False, seed: int = 0) -> list:
    return [run_case(c, seed) for c in CASES if c.quick or not quick]


def _flags(r) -> str:
    return ", ".join(f + ("" if f in r.errors else "(w)") for f in sorted(r.flagged)) or "(clean)"


def render(results) -> str:
    lines = [f"{'case':<32} {'expected':<40} {'flagged ((w) = warning only)':<48} {'truth':<6} result"]
    for r in results:
        exp = (" | " if r.case.any_of else ", ").join(sorted(r.case.expect)) or "(no errors)"
        truth = "n/a" if r.gt_correct is None else ("right" if r.gt_correct else "wrong")
        lines.append(f"{r.case.name:<32} {exp:<40} {_flags(r):<48} {truth:<6} {'pass' if r.passed else 'FAIL'}")
    n = sum(r.passed for r in results)
    fw = sum(r.false_warning for r in results)
    n_clean = sum(not r.case.expect for r in results)
    lines += ["", f"{n}/{len(results)} cases pass; {fw}/{n_clean} clean cases raised warnings.",
              "'truth' says whether the discovered structure was the true law. 'a | b' means any of them."]
    return "\n".join(lines)


def sweep(seeds, quick: bool = False):
    """Run every case over several seeds and summarize pass and false-warning rates."""
    rows: dict = {}
    for seed in seeds:
        for r in run_all(quick=quick, seed=seed):
            d = rows.setdefault(r.case.name, {"pass": 0, "fw": 0, "wrong": 0, "wrong_flagged": 0, "n": 0,
                                              "clean": not r.case.expect})
            d["n"] += 1
            d["pass"] += r.passed
            d["fw"] += r.false_warning
            if r.gt_correct is False:
                d["wrong"] += 1
                d["wrong_flagged"] += bool(r.flagged)
    lines = [f"{'case':<32} {'passed':>8} {'warned (clean)':>15} {'wrong law, flagged':>20}"]
    for name, d in rows.items():
        fw = f"{d['fw']}/{d['n']}" if d["clean"] else "-"
        wf = f"{d['wrong_flagged']}/{d['wrong']}" if d["wrong"] else "-"
        lines.append(f"{name:<32} {d['pass']:>4}/{d['n']:<3} {fw:>15} {wf:>20}")
    tot = sum(d["pass"] for d in rows.values()), sum(d["n"] for d in rows.values())
    lines += ["", f"{tot[0]}/{tot[1]} case runs pass over seeds {list(seeds)[0]}-{list(seeds)[-1]}."]
    return "\n".join(lines), tot[0] == tot[1]
