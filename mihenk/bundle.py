"""The mihenk bundle: the observations an equation was discovered from, as plain arrays.

A bundle is one ``.npz`` file holding one or more trajectories. Shapes use T
time samples and n state variables; trajectories may differ in length.

Keys
    x             (T, n)   one trajectory, or ``x_0, x_1, ...`` for several
    t             (T,)     sample times (or ``t_0, t_1, ...``); a scalar ``dt`` also works
    x_dot         (T, n)   optional measured derivatives (``x_dot_0, ...``)
    val           (K,)     optional bool: trajectories held out from discovery.
                           With one trajectory, the last 20% is held out instead.
    meta          JSON string: variables, name, units, true_equations (synthetic only), ...

``true_equations`` in meta (``{"x": "-10*x + 10*y", ...}``) unlocks the tier-2
ground-truth checks; only synthetic data has it.
"""

from __future__ import annotations

import json
from dataclasses import dataclass, field
from pathlib import Path

import numpy as np


@dataclass
class Bundle:
    xs: list                       # list of (T_k, n) arrays
    ts: list                       # list of (T_k,) arrays
    variables: list
    x_dots: list | None = None
    val: np.ndarray | None = None  # (K,) bool, held-out trajectories
    meta: dict = field(default_factory=dict)

    def __post_init__(self):
        if isinstance(self.xs, np.ndarray):
            self.xs = [self.xs]
        if isinstance(self.ts, np.ndarray):
            self.ts = [self.ts]
        self.xs = [np.atleast_2d(np.asarray(x, float)) for x in self.xs]
        self.ts = [np.asarray(t, float).ravel() for t in self.ts]
        if len(self.xs) != len(self.ts):
            raise ValueError(f"{len(self.xs)} trajectories but {len(self.ts)} time vectors")
        n = self.xs[0].shape[1]
        for k, (x, t) in enumerate(zip(self.xs, self.ts)):
            if x.shape[1] != n:
                raise ValueError(f"trajectory {k} has {x.shape[1]} variables, expected {n}")
            if x.shape[0] != t.shape[0]:
                raise ValueError(f"trajectory {k}: {x.shape[0]} samples but {t.shape[0]} times")
            if x.shape[0] < 20:
                raise ValueError(f"trajectory {k} has only {x.shape[0]} samples; mihenk needs at least 20")
            if np.any(np.diff(t) <= 0):
                raise ValueError(f"trajectory {k}: times must be strictly increasing")
            if not np.all(np.isfinite(x)):
                raise ValueError(f"trajectory {k} contains NaN or inf")
        if len(self.variables) != n:
            raise ValueError(f"{len(self.variables)} variable names for {n} columns")
        self.variables = [str(v) for v in self.variables]
        if self.x_dots is not None:
            self.x_dots = [np.asarray(d, float) for d in self.x_dots]
        if self.val is not None:
            self.val = np.asarray(self.val, bool).ravel()
            if self.val.shape[0] != len(self.xs):
                raise ValueError("val must have one entry per trajectory")
            if self.val.all():
                raise ValueError("every trajectory is marked held-out; nothing is left for fitting")

    # ------------------------------------------------------------ split
    @property
    def n_traj(self) -> int:
        return len(self.xs)

    @property
    def n_vars(self) -> int:
        return self.xs[0].shape[1]

    @property
    def name(self) -> str:
        return self.meta.get("name", "bundle")

    @property
    def true_equations(self) -> dict | None:
        return self.meta.get("true_equations")

    def split(self, holdout_fraction: float = 0.2):
        """(train, test) lists of (x, t) pairs.

        Held-out trajectories when ``val`` is given (or, with several
        trajectories, the last 20% of them); otherwise the tail of the single
        trajectory.
        """
        pairs = list(zip(self.xs, self.ts))
        if self.val is not None:
            train = [p for p, v in zip(pairs, self.val) if not v]
            test = [p for p, v in zip(pairs, self.val) if v]
            return train, test, "held-out trajectories (val)"
        if self.n_traj >= 3:
            k = max(1, int(round(holdout_fraction * self.n_traj)))
            return pairs[:-k], pairs[-k:], f"last {k} of {self.n_traj} trajectories"
        train, test = [], []
        for x, t in pairs:
            cut = int(round((1 - holdout_fraction) * len(t)))
            train.append((x[:cut], t[:cut]))
            test.append((x[cut:], t[cut:]))
        return train, test, f"last {int(holdout_fraction * 100)}% of each trajectory"

    # ------------------------------------------------------------ io
    def save(self, path) -> Path:
        path = Path(path)
        arrays = {}
        for k, (x, t) in enumerate(zip(self.xs, self.ts)):
            arrays[f"x_{k}"] = x
            arrays[f"t_{k}"] = t
            if self.x_dots is not None:
                arrays[f"x_dot_{k}"] = self.x_dots[k]
        if self.val is not None:
            arrays["val"] = self.val
        meta = dict(self.meta)
        meta["variables"] = self.variables
        arrays["meta"] = np.array(json.dumps(meta))
        np.savez_compressed(path, **arrays)
        return path

    @classmethod
    def load(cls, path) -> "Bundle":
        path = Path(path)
        if path.suffix == ".csv":
            return cls.from_csv(path)
        with np.load(path, allow_pickle=False) as z:
            keys = set(z.files)
            meta = json.loads(str(z["meta"])) if "meta" in keys else {}
            if "x" in keys:
                xs = [z["x"]]
                ts = [_times(z, "t", xs[0].shape[0])]
                dots = [z["x_dot"]] if "x_dot" in keys else None
            else:
                idx = sorted(int(k[2:]) for k in keys if k.startswith("x_") and k[2:].isdigit())
                if not idx:
                    raise ValueError(f"{path}: no 'x' or 'x_0' array")
                xs = [z[f"x_{k}"] for k in idx]
                ts = [_times(z, f"t_{k}", xs[i].shape[0]) for i, k in enumerate(idx)]
                dots = [z[f"x_dot_{k}"] for k in idx] if all(f"x_dot_{k}" in keys for k in idx) else None
            val = z["val"] if "val" in keys else None
        variables = meta.get("variables") or [f"x{i}" for i in range(xs[0].shape[1])]
        meta.setdefault("name", path.stem)
        return cls(xs, ts, variables, x_dots=dots, val=val, meta=meta)

    @classmethod
    def from_csv(cls, path, time_column: str = "t") -> "Bundle":
        """One trajectory per file; a header row names the columns, one of them the time."""
        path = Path(path)
        with open(path) as fh:
            header = [h.strip() for h in fh.readline().split(",")]
        data = np.loadtxt(path, delimiter=",", skiprows=1, ndmin=2)
        if time_column not in header:
            raise ValueError(f"{path}: no '{time_column}' column; columns are {header}")
        ti = header.index(time_column)
        cols = [i for i in range(len(header)) if i != ti]
        return cls([data[:, cols]], [data[:, ti]], [header[i] for i in cols], meta={"name": path.stem})


def _times(z, key, n):
    if key in z.files:
        return z[key]
    if "dt" in z.files:
        return np.arange(n) * float(z["dt"])
    raise ValueError(f"no '{key}' times and no 'dt'")
