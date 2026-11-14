import json
import subprocess
import sys

import numpy as np
import pytest

import mihenk
from mihenk import synth
from mihenk.bundle import Bundle
from mihenk.model import Model, parse_term
from mihenk.scan import Thresholds, scan

from mihenk.adapters import pysindy as A

try:
    import pysindy  # noqa: F401

    HAS_PYSINDY = True
except ImportError:
    HAS_PYSINDY = False

# only the tests that fit with PySINDy need it; the core checks run on numpy/scipy/sympy alone
needs_pysindy = pytest.mark.skipif(not HAS_PYSINDY, reason='needs PySINDy: pip install "mihenk[pysindy]"')


# ----------------------------------------------------------------- model
def test_parse_pysindy_print_output():
    m = Model.from_text("(x)' = -9.999 x + 9.999 y\n(y)' = 27.992 x + -0.999 y + -1.000 x z\n"
                        "(z)' = -2.666 z + 1.000 x y")
    assert m.variables == ["x", "y", "z"]
    assert len(m.active(1)) == 3
    X = np.array([[1.0, 2.0, 3.0]])
    np.testing.assert_allclose(m.predict(X), [[9.999, 27.992 - 1.998 - 3.0, -7.998 + 2.0]], rtol=1e-12)


@pytest.mark.parametrize("lhs", ["x'", "dx/dt", "d(x)/dt", "x_dot", "(x)'", " x ' "])
def test_lhs_forms(lhs):
    m = Model.from_text(f"{lhs} = -x + x^3")
    assert m.pretty() == "x' = -x + x^3"


def test_names_starting_with_d():
    m = Model.from_text("delta' = -delta")
    assert m.variables == ["delta"]


def test_multichar_names_not_split():
    t = parse_term("theta omega", ["theta", "omega"])
    assert str(t) == "omega*theta"


def test_unknown_symbol_rejected():
    with pytest.raises(ValueError, match="unknown symbol"):
        Model.from_equations({"x": "x + q"}, variables=["x"])


def test_model_json_roundtrip(tmp_path):
    m = synth.truth_model("lorenz")
    p = m.save(tmp_path / "m.json")
    m2 = Model.load(p)
    assert m2.structure() == m.structure()
    np.testing.assert_allclose(m2.coefficients, m.coefficients)


def test_bundle_roundtrip_and_validation(tmp_path):
    b = synth.simulate("lotka_volterra", n_traj=2, seed=1)
    p = b.save(tmp_path / "b.npz")
    b2 = Bundle.load(p)
    assert b2.n_traj == 2 and b2.variables == ["x", "y"] and b2.true_equations
    with pytest.raises(ValueError, match="strictly increasing"):
        Bundle([np.zeros((30, 1))], [np.zeros(30)], ["x"])


# --------------------------------------------------------------- adapter
def _fit(system, threshold, degree=None, noise=0.0, n_traj=1, seed=0):
    b = synth.simulate(system, noise=noise, n_traj=n_traj, seed=seed)
    sindy, model = A.fit_bundle(b, degree=degree or synth.SYSTEMS[system].degree, threshold=threshold)
    return b, sindy, model


@needs_pysindy
def test_from_pysindy_matches_pysindy_predictions():
    b, sindy, model = _fit("lorenz", 0.5)
    X = b.xs[0][:100]
    np.testing.assert_allclose(model.predict(X), sindy.predict(X), rtol=1e-8, atol=1e-8)


# ----------------------------------------------------------------- scans
@needs_pysindy
def test_clean_lorenz_is_clean():
    b, sindy, _ = _fit("lorenz", 0.5)
    r = mihenk.scan_pysindy(sindy, b.xs, b.ts, meta={"true_equations": b.true_equations})
    assert not r.flagged(), r.render_text()
    assert r.metrics["gt_structure_correct"] is True
    assert not r.has_errors


@needs_pysindy
def test_on_cycle_equifinality_finds_sparser_law():
    b, sindy, model = _fit("stuart_landau_on_cycle", 0.05)
    r = scan(b, model)
    eq = r.by_id("EQUIFINALITY")
    assert eq and eq[0].severity == "error"
    assert "x' = -y" in eq[0].what


@needs_pysindy
def test_missing_library_term_is_underfit():
    b, sindy, model = _fit("duffing", 0.05, degree=2)
    r = scan(b, model)
    assert any(f.id == "UNDERFIT" and "x^3" in f.where for f in r.findings)


def test_wrong_sign_rollout_trap():
    b = synth.simulate("linear_damped")
    m = Model.from_text("x' = 0.1 x + 2 y\ny' = -2 x + 0.1 y")
    r = scan(b, m)
    assert "ROLLOUT_TRAP" in r.flagged()


@needs_pysindy
def test_noise_spurious_constant_flagged():
    b, sindy, _ = _fit("lorenz", 0.1, noise=0.05)
    r = mihenk.scan_pysindy(sindy, b.xs, b.ts)
    assert {"WEAK_TERM", "UNSTABLE_STRUCTURE"} & r.flagged()


def test_tier1_not_assessed_without_refit():
    b = synth.simulate("lorenz")
    r = scan(b, synth.truth_model("lorenz"))
    status = {c.id: c.status for c in r.checks}
    assert status["UNSTABLE_STRUCTURE"] == "not_assessed"
    assert status["FORKING_PATHS"] == "not_assessed"
    assert status["RECOVERY"] == "clean"


def test_thresholds_override_and_unknown():
    th = Thresholds().update({"weak_term_f": "3", "n_boot": "5"})
    assert th.weak_term_f == 3.0 and th.n_boot == 5
    with pytest.raises(KeyError):
        Thresholds().update({"nope": 1})


def test_report_files(tmp_path):
    b = synth.simulate("stuart_landau_on_cycle")
    r = scan(b, Model.from_text("x' = -y\ny' = x"))
    out = r.save(tmp_path / "rep")
    d = json.loads((out / "report.json").read_text())
    assert d["mihenk_version"] == mihenk.__version__
    html = (out / "report.html").read_text()
    assert "<!DOCTYPE html>" in html and "EQUIFINALITY" in html


# ------------------------------------------------------------------- CLI
def _cli(*args, cwd):
    return subprocess.run([sys.executable, "-m", "mihenk.cli", *args], cwd=cwd, capture_output=True, text=True)


@needs_pysindy
def test_cli_end_to_end(tmp_path):
    assert _cli("synth", "stuart_landau_on_cycle", "-o", "b.npz", cwd=tmp_path).returncode == 0
    r = _cli("fit", "b.npz", "-o", "m.json", "--degree", "3", "--threshold", "0.05", cwd=tmp_path)
    assert r.returncode == 0, r.stderr
    r = _cli("scan", "b.npz", "--model", "m.json", "--out", "rep", cwd=tmp_path)
    assert r.returncode == 1, r.stdout + r.stderr  # errors gate the pipeline
    assert "EQUIFINALITY" in r.stdout
    assert (tmp_path / "rep" / "report.html").exists()


def test_cli_clean_exit_zero(tmp_path):
    _cli("synth", "lorenz", "-o", "b.npz", cwd=tmp_path)
    (tmp_path / "eq.txt").write_text("x' = -10 x + 10 y\ny' = 28 x - y - x z\nz' = -2.6667 z + x y\n")
    r = _cli("scan", "b.npz", "--model", "eq.txt", cwd=tmp_path)
    assert r.returncode == 0, r.stdout + r.stderr


def test_cli_usage_error(tmp_path):
    _cli("synth", "lorenz", "-o", "b.npz", cwd=tmp_path)
    (tmp_path / "eq.txt").write_text("x' = -10 x\n")
    r = _cli("scan", "b.npz", "--model", "eq.txt", cwd=tmp_path)
    assert r.returncode == 2 and "do not match" in r.stderr


def test_rollout_r2_on_resting_segment_is_finite():
    """A held-out segment where the system is at rest must not give a meaningless R² (0.1.0 bug)."""
    t = np.arange(0, 20, 0.01)
    x = np.column_stack([np.exp(-0.5 * t) * np.cos(2 * t), np.exp(-0.5 * t) * (-2 * np.sin(2 * t))])
    x[int(0.8 * len(t)):] = x[int(0.8 * len(t))]          # at rest for the last 20%
    m = Model.from_text("x' = y\ny' = -4 x - y")
    r = scan(Bundle([x], [t], ["x", "y"]), m)
    assert r.metrics["rollout_r2_short"] > -10


# ------------------------------------------------------------ kinematics
from mihenk.kinematics import fit_pair, suggest_pairs, validate_pairs  # noqa: E402


def test_kinematic_pairs_parse_and_reject():
    assert validate_pairs(["theta=omega"], ["theta", "omega"]) == {"theta": "omega"}
    assert validate_pairs({"theta'": "omega"}, ["theta", "omega"]) == {"theta": "omega"}
    with pytest.raises(ValueError, match="not a variable"):
        validate_pairs({"theta": "v"}, ["theta", "omega"])
    with pytest.raises(ValueError, match="itself"):
        validate_pairs({"theta": "theta"}, ["theta", "omega"])


@pytest.mark.parametrize("sensor, field, expected", [
    ({"omega": {"delay": 7e-4}}, "lag", 7e-4),
    ({"omega": {"gain": 1.02}}, "gain", 1 / 1.02),
    ({"omega": {"offset": 0.05}}, "offset", -0.05),
])
def test_fit_pair_recovers_sensor_imperfection(sensor, field, expected):
    b = synth.simulate("pendulum", sensor=sensor)
    f = fit_pair(list(zip(b.xs, b.ts)), 0, 1, window=50)
    assert f[field] == pytest.approx(expected, rel=0.02, abs=2e-5)


def test_kinematic_not_assessed_without_declaration_suggests_pair():
    b = synth.simulate("pendulum")
    b.meta.pop("kinematics")
    r = scan(b, synth.truth_model("pendulum"))
    c = {c.id: c for c in r.checks}["KINEMATIC_MISMATCH"]
    assert c.status == "not_assessed" and "theta=omega" in c.reason


def test_kinematic_clean_with_exact_sensor():
    b = synth.simulate("pendulum")              # declares theta = omega in meta
    r = scan(b, synth.truth_model("pendulum"))
    assert {c.id: c.status for c in r.checks}["KINEMATIC_MISMATCH"] == "clean"


def test_kinematic_lag_flagged_and_explained():
    b = synth.simulate("pendulum", sensor={"omega": {"delay": 7e-4}})
    m = Model.from_text("theta' = omega - 0.0448 sin(theta)\nomega' = -64 sin(theta) - 0.032 omega",
                        variables=b.variables)
    r = scan(b, m)
    f = r.by_id("KINEMATIC_MISMATCH")[0]
    assert f.severity == "error"
    assert f.evidence["lag"] == pytest.approx(7e-4, rel=0.02)
    assert f.evidence["lag_explains"] > 0.9            # the lag predicts the extra theta' term
    damping = [c for c in f.evidence["coefficient_changes"] if c[1] == "omega"][0]
    assert damping[3] == pytest.approx(-0.07, rel=0.05)  # correcting omega restores the damping


def test_kinematics_argument_overrides_meta():
    b = synth.simulate("pendulum", sensor={"omega": {"gain": 1.02}})
    r = scan(b, synth.truth_model("pendulum"), kinematics={})
    assert {c.id: c.status for c in r.checks}["KINEMATIC_MISMATCH"] == "not_assessed"


def test_suggest_pairs():
    assert suggest_pairs(synth.truth_model("pendulum")) == {"theta": "omega"}


def test_cli_kinematic_flag(tmp_path):
    from mihenk.bundle import Bundle as B

    b = synth.simulate("pendulum", sensor={"omega": {"delay": 7e-4}})
    b.meta.pop("kinematics")
    b.save(tmp_path / "p.npz")
    (tmp_path / "eq.txt").write_text("theta' = omega - 0.0448 sin(theta)\nomega' = -64 sin(theta) - 0.032 omega\n")
    r = _cli("scan", "p.npz", "--model", "eq.txt", "--kinematic", "theta=omega", cwd=tmp_path)
    assert r.returncode == 1 and "KINEMATIC_MISMATCH" in r.stdout, r.stdout + r.stderr
    r = _cli("scan", "p.npz", "--model", "eq.txt", cwd=tmp_path)
    assert "KINEMATIC_MISMATCH" not in [ln.split()[1] for ln in r.stdout.splitlines()
                                        if ln.startswith(("ERROR", "WARNING"))]


def test_report_save_formats(tmp_path):
    r = scan(synth.simulate("lorenz"), synth.truth_model("lorenz"))
    out = r.save(tmp_path / "rep", formats=("html", "json"))
    assert sorted(p.name for p in out.iterdir()) == ["report.html", "report.json"]
    with pytest.raises(ValueError, match="unknown report format"):
        r.save(tmp_path / "bad", formats=("pdf",))


# ------------------------------------------------------- 0.3: physics, rollouts, underfit, relations
from mihenk.physics import validate_properties  # noqa: E402

H_PEND = "0.5*omega**2 - 64*cos(theta)"
H_LV = "0.5*x - log(x) + 0.5*y - log(y)"


def _phys(r):
    return {c.id: c.status for c in r.checks}["PHYSICS_VIOLATION"], r.by_id("PHYSICS_VIOLATION")


def test_properties_parse_and_reject():
    assert set(validate_properties({"dissipates": H_PEND}, ["theta", "omega"])) == {"dissipates"}
    assert set(validate_properties([f"conserves={H_LV}"], ["x", "y"])) == {"conserves"}
    with pytest.raises(ValueError, match="unknown property"):
        validate_properties({"heats": "x"}, ["x"])
    with pytest.raises(ValueError, match="unknown symbol"):
        validate_properties({"dissipates": "0.5*omega**2 - g*cos(theta)"}, ["theta", "omega"])


def test_physics_not_assessed_without_declaration():
    r = scan(synth.simulate("pendulum"), synth.truth_model("pendulum"))
    assert _phys(r)[0] == "not_assessed"


def test_dissipation_true_law_clean_injection_error_overdamping_warning():
    b = synth.simulate("pendulum")
    props = {"dissipates": H_PEND}
    assert _phys(scan(b, synth.truth_model("pendulum"), properties=props))[0] == "clean"
    inject = Model.from_text("theta' = omega\nomega' = -64 sin(theta) + 0.03 omega", variables=b.variables)
    status, f = _phys(scan(b, inject, properties=props))
    assert status == "flagged" and f[0].severity == "error" and "gains" in f[0].what
    over = Model.from_text("theta' = omega\nomega' = -64 sin(theta) - 0.21 omega", variables=b.variables)
    status, f = _phys(scan(b, over, properties=props))
    assert status == "flagged" and f[0].severity == "warning"


def test_conservation_true_law_clean_damped_flagged():
    b = synth.simulate("lotka_volterra")
    props = {"conserves": H_LV}
    assert _phys(scan(b, synth.truth_model("lotka_volterra"), properties=props))[0] == "clean"
    damped = Model.from_text("x' = 0.95 x - 0.5 x y\ny' = -y + 0.5 x y", variables=b.variables)
    assert _phys(scan(b, damped, properties=props))[0] == "flagged"


def test_long_rollout_drift_flagged_for_undamped_pendulum():
    b = synth.simulate("pendulum")
    m = Model.from_text("theta' = omega\nomega' = -64 sin(theta)", variables=b.variables)
    r = scan(b, m, kinematics={})
    assert any("drifts away" in f.what for f in r.by_id("ROLLOUT_TRAP"))


def test_twin_rollouts_classify_lorenz_as_chaotic():
    r = scan(synth.simulate("lorenz"), synth.truth_model("lorenz"))
    full = [row for row in r.details["Rollouts (integrating the equation)"]["rows"] if row[1] == "full"]
    assert full[0][5].startswith("chaotic") and full[0][6] > 1000
    assert not any("drifts away" in f.what for f in r.by_id("ROLLOUT_TRAP"))


@needs_pysindy
def test_underfit_finds_degree_five_beyond_cubic_library():
    # theta^5 is the best left-out term for every seed; it passes the warning gate in about
    # 8 of 10 (experiments/physics_synthetic). Seed 1 is one that passes.
    b = synth.simulate("pendulum", seed=1)
    _, m = A.fit_bundle(b, degree=3, threshold=0.02)
    r = scan(b, m, kinematics={})
    rows = r.details["Best left-out term per equation"]["rows"]
    assert max(rows, key=lambda row: row[3])[1] == "theta^5"
    assert [f for f in r.by_id("UNDERFIT") if "theta^5" in f.where]


def test_centred_relation_reports_sin_theta_about_pi():
    b = synth.simulate("pendulum_small", sensor={"theta": {"offset": float(np.pi)}})
    m = Model.from_text("theta' = omega\nomega' = 64 sin(theta)", variables=b.variables,)
    m.library = [parse_term(t, b.variables) for t in ["1", "theta", "omega", "sin(theta)", "cos(theta)"]]
    r = scan(b, m, kinematics={})
    rows = r.details["Relations among candidate terms on your data"]["rows"]
    assert any("sin(theta)" in row[0] and "theta" in row[0].replace("sin(theta)", "") for row in rows)
    assert not any(row[0].replace(" ", "") in ("+1*theta-1≈0", "-1*theta+1≈0") for row in rows)


def test_cli_dissipates_flag(tmp_path):
    synth.simulate("pendulum").save(tmp_path / "p.npz")
    (tmp_path / "eq.txt").write_text("theta' = omega\nomega' = -64 sin(theta) + 0.03 omega\n")
    r = _cli("scan", "p.npz", "--model", "eq.txt", "--dissipates", H_PEND, cwd=tmp_path)
    assert r.returncode == 1 and "PHYSICS_VIOLATION" in r.stdout, r.stdout + r.stderr
