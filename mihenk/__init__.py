"""mihenk: a touchstone for discovered equations.

Before you interpret the law that PySINDy (or any equation-discovery tool)
found, mihenk checks whether it could be an artifact of noise, of how
derivatives were estimated, of which states the data happened to visit, or of
analysis choices.

    import mihenk
    report = mihenk.scan(bundle, model)            # tier 0: data + equation
    report = mihenk.scan_pysindy(sindy, x, t)      # tiers 0-1 with a fitted PySINDy model
"""

__version__ = "0.3.0"

from .bundle import Bundle  # noqa: E402
from .model import Model  # noqa: E402
from .report import Check, Finding, Report  # noqa: E402
from .scan import Thresholds, scan  # noqa: E402


def scan_pysindy(sindy_model, x, t, *, feature_names=None, name="pysindy", val=None,
                 thresholds=None, variants=True, seed=0, meta=None, kinematics=None, properties=None):
    """Scan a fitted PySINDy model on the data it was fitted to (tiers 0 and 1).

    ``kinematics={"theta": "omega"}`` declares theta' = omega and unlocks KINEMATIC_MISMATCH;
    ``properties={"dissipates": "H"}`` (or ``"conserves"``) unlocks PHYSICS_VIOLATION.
    """
    from .adapters import pysindy as _ps

    return _ps.scan_pysindy(sindy_model, x, t, feature_names=feature_names, name=name, val=val,
                            thresholds=thresholds, variants=variants, seed=seed, meta=meta,
                            kinematics=kinematics, properties=properties)


__all__ = ["Bundle", "Model", "Report", "Finding", "Check", "Thresholds", "scan", "scan_pysindy",
           "__version__"]
