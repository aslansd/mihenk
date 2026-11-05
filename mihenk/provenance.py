"""Optional provenance through daftar (``pip install "mihenk[provenance]"``).

Each tracked scan records the bundle and model files, the thresholds, mihenk's
version and environment, and the headline metrics and findings, so two scans
can be compared later with ``daftar diff``.
"""

from __future__ import annotations

from contextlib import contextmanager
from pathlib import Path


def available() -> bool:
    try:
        import daftar  # noqa: F401
    except ImportError:
        return False
    return True


@contextmanager
def tracked(label: str, params: dict, inputs=(), seed: int | None = None):
    """Yield a daftar run, or ``None`` when daftar is not installed."""
    if not available():
        yield None
        return
    import daftar

    with daftar.track(label, params=params, seed=seed, tags=["mihenk"]) as run:
        for path in inputs:
            run.add_input(path)
        yield run


def log_report(run, report, out_dir=None):
    if run is None:
        return
    for key, value in report.metrics.items():
        if isinstance(value, (int, float)) and not isinstance(value, bool):
            run.log_result(f"metrics.{key}", float(value))
    run.log_result("equation", report.equation)
    run.log_result("findings", sorted(f"{f.severity}:{f.id}:{f.where}" for f in report.findings if f.tier < 2))
    if out_dir is not None:
        for name in ("report.json", "report.html"):
            path = Path(out_dir) / name
            if path.exists():
                run.add_output(path)
