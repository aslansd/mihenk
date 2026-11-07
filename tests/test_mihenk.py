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
