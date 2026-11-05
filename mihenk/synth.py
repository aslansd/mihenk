"""Synthetic systems with known equations, for self-validation and for planning.

Each system is integrated from its true equation, then Gaussian noise is added
as a fraction of each variable's standard deviation. The true equation goes
into the bundle's meta, which unlocks mihenk's tier-2 checks.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Callable

import numpy as np
from scipy.integrate import solve_ivp

from .bundle import Bundle
from .model import Model


@dataclass
class System:
    variables: list
    equations: dict
    x0: Callable           # rng -> initial state
    dt: float
    duration: float
    degree: int            # polynomial library degree that contains the truth
    note: str = ""


def _ring(r=1.0):
    def f(rng):
        a = rng.uniform(0, 2 * np.pi)
        return np.array([r * np.cos(a), r * np.sin(a)])
    return f


SYSTEMS = {
    "lorenz": System(
        ["x", "y", "z"],
        {"x": "-10*x + 10*y", "y": "28*x - y - x*z", "z": "-8/3*z + x*y"},
        lambda rng: rng.normal([0, 0, 25], [8, 8, 8]), 0.002, 10.0, 2,
        "chaotic; the standard SINDy demonstration"),
    "lotka_volterra": System(
        ["x", "y"], {"x": "x - 0.5*x*y", "y": "-y + 0.5*x*y"},
        lambda rng: rng.uniform(1.0, 3.5, 2), 0.01, 20.0, 2, "predator-prey cycles"),
    "van_der_pol": System(
        ["x", "y"], {"x": "y", "y": "-x + y - x**2*y"},
        lambda rng: rng.uniform(-3, 3, 2), 0.01, 20.0, 3, "relaxation oscillator with a limit cycle"),
    "duffing": System(
        ["x", "y"], {"x": "y", "y": "-x - 0.1*y - x**3"},
        lambda rng: rng.uniform(-2, 2, 2), 0.01, 20.0, 3, "damped cubic oscillator"),
    "stuart_landau": System(
        ["x", "y"], {"x": "x - y - x**3 - x*y**2", "y": "x + y - x**2*y - y**3"},
        lambda rng: rng.uniform(-1.5, 1.5, 2), 0.01, 20.0, 3, "Hopf normal form; attracting unit circle"),
    "stuart_landau_on_cycle": System(
        ["x", "y"], {"x": "x - y - x**3 - x*y**2", "y": "x + y - x**2*y - y**3"},
        _ring(1.0), 0.01, 20.0, 3,
        "the same system observed only on its limit cycle, where x^2 + y^2 = 1 and terms become interchangeable"),
    "linear_damped": System(
        ["x", "y"], {"x": "-0.1*x + 2*y", "y": "-2*x - 0.1*y"},
        lambda rng: rng.uniform(-2, 2, 2), 0.01, 20.0, 1, "damped linear oscillator"),
}


def truth_model(system: str) -> Model:
    s = SYSTEMS[system]
    return Model.from_equations(s.equations, variables=s.variables, source=f"truth:{system}")


def simulate(system: str, *, noise: float = 0.0, n_traj: int = 1, dt: float | None = None,
             duration: float | None = None, seed: int = 0, x0=None) -> Bundle:
    """Simulate ``system`` and return a :class:`Bundle` with the truth in ``meta``."""
    if system not in SYSTEMS:
        raise KeyError(f"unknown system '{system}'; known: {', '.join(SYSTEMS)}")
    s = SYSTEMS[system]
    rng = np.random.default_rng(seed)
    dt = dt or s.dt
    duration = duration or s.duration
    t = np.arange(0, duration, dt)
    f = truth_model(system).rhs_function()
    xs, clean = [], []
    for k in range(n_traj):
        start = np.asarray(x0[k] if x0 is not None else s.x0(rng), float)
        sol = solve_ivp(f, (t[0], t[-1]), start, t_eval=t, method="LSODA", rtol=1e-10, atol=1e-12)
        if sol.status != 0 or sol.y.shape[1] != len(t):
            raise RuntimeError(f"integration of {system} failed from {start}")
        clean.append(sol.y.T)
    sd = np.vstack(clean).std(axis=0)
    for c in clean:
        xs.append(c + noise * sd * rng.standard_normal(c.shape))
    meta = {"name": f"{system}" + (f"_noise{noise:g}" if noise else ""), "system": system,
            "true_equations": s.equations, "noise": noise, "dt": dt, "seed": seed,
            "library_degree": s.degree}
    return Bundle(xs, [t] * n_traj, list(s.variables), meta=meta)
