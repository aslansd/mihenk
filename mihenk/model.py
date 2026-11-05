"""A discovered equation, as mihenk sees it: variables, terms and coefficients.

mihenk never needs the tool that found the equation. A :class:`Model` is a set
of right-hand sides ``x_i' = sum_j c_ij * theta_j(x)`` where every term
``theta_j`` is a SymPy expression in the state variables. It can be built from

* text, including PySINDy's own ``model.print()`` output::

      (x)' = -9.999 x + 9.999 y
      (y)' = 27.992 x + -0.999 y + -1.000 x z

* a fitted PySINDy model (:func:`mihenk.adapters.pysindy.from_pysindy`);
* SymPy expressions (:meth:`Model.from_sympy`);
* a JSON file written by :meth:`Model.save`.
"""

from __future__ import annotations

import itertools
import json
import re
from dataclasses import dataclass, field
from pathlib import Path

import numpy as np
import sympy as sp
from sympy.parsing.sympy_parser import (
    convert_xor,
    implicit_multiplication,
    parse_expr,
    standard_transformations,
)

_TRANSFORMS = standard_transformations + (convert_xor, implicit_multiplication)
_LHS_DDT = re.compile(r"^\s*d\s*\(?\s*([A-Za-z_][A-Za-z_0-9]*)\s*\)?\s*/\s*dt\s*$")
_LHS = re.compile(r"^\s*\(?\s*([A-Za-z_][A-Za-z_0-9]*?)\s*\)?\s*(?:'|_dot)?\s*$")


def parse_term(text: str, variables) -> sp.Expr:
    """Parse one expression; accepts ``x^2``, ``x y`` (implicit product) and ``sin(1 x)``."""
    local = {v: sp.Symbol(v) for v in variables}
    text = text.replace("·", "*")
    return parse_expr(text, local_dict=local, transformations=_TRANSFORMS)


def split_terms(expr: sp.Expr) -> dict:
    """Split an expanded expression into ``{term: coefficient}`` with numeric coefficients."""
    expr = sp.expand(sp.sympify(expr))
    out: dict = {}
    for term, coef in expr.as_coefficients_dict().items():
        # as_coefficients_dict keeps numbers inside functions; pull out remaining numeric factors
        num, rest = sp.Mul(coef, term).as_coeff_Mul()
        out[rest] = out.get(rest, 0.0) + float(num)
    return {k: v for k, v in out.items() if v != 0.0}


def term_label(term: sp.Expr) -> str:
    return "1" if term == 1 else str(term).replace("**", "^")


def format_rhs(terms, coefs, precision: int = 4) -> str:
    """``x - 0.5 x*y`` style right-hand side."""
    out = ""
    for t, c in zip(terms, coefs):
        c = float(c)
        if c == 0:
            continue
        mag = abs(c)
        label = term_label(t)
        num = f"{mag:.{precision}g}"
        body = num if label == "1" else (label if num == "1" else f"{num} {label}")
        if not out:
            out = ("-" if c < 0 else "") + body
        else:
            out += (" - " if c < 0 else " + ") + body
    return out or "0"


def term_degree(term: sp.Expr, symbols) -> int:
    try:
        return int(sp.Poly(term, *symbols).total_degree())
    except sp.PolynomialError:
        return -1  # not a polynomial term


@dataclass
class Model:
    """Right-hand sides ``x' = Xi @ theta(x)`` over a shared list of terms."""

    variables: list
    terms: list                      # sympy expressions, shared by all equations
    coefficients: np.ndarray         # (n_variables, n_terms)
    library: list | None = None      # candidate terms the discovery searched over, if known
    source: str = ""
    meta: dict = field(default_factory=dict)

    def __post_init__(self):
        self.coefficients = np.asarray(self.coefficients, dtype=float)
        if self.coefficients.shape != (len(self.variables), len(self.terms)):
            raise ValueError(
                f"coefficients have shape {self.coefficients.shape}, expected "
                f"({len(self.variables)}, {len(self.terms)})"
            )
        self._fns = None

    # ----------------------------------------------------------- building
    @classmethod
    def from_equations(cls, equations: dict, variables=None, library=None, source="") -> "Model":
        """``{"x": "-10 x + 10 y", ...}`` (strings or SymPy) in variable order."""
        variables = list(variables or equations.keys())
        missing = [v for v in variables if v not in equations]
        if missing:
            raise ValueError(f"no equation for variable(s): {', '.join(missing)}")
        split = []
        for v in variables:
            rhs = equations[v]
            expr = parse_term(rhs, variables) if isinstance(rhs, str) else sp.sympify(rhs)
            unknown = {str(s) for s in expr.free_symbols} - set(variables)
            if unknown:
                raise ValueError(f"equation for {v} uses unknown symbol(s): {', '.join(sorted(unknown))}")
            split.append(split_terms(expr))
        terms = _ordered_terms(itertools.chain.from_iterable(s.keys() for s in split), variables)
        coefs = np.array([[s.get(t, 0.0) for t in terms] for s in split])
        lib = None
        if library is not None:
            lib = [parse_term(t, variables) if isinstance(t, str) else sp.sympify(t) for t in library]
        return cls(variables, terms, coefs, library=lib, source=source)

    @classmethod
    def from_sympy(cls, equations: dict, **kw) -> "Model":
        return cls.from_equations(equations, **kw)

    @classmethod
    def from_text(cls, text: str, variables=None, source="text") -> "Model":
        """Parse one equation per line: ``x' = ...``, ``(x)' = ...``, ``dx/dt = ...`` or ``x_dot = ...``."""
        eqs = {}
        for raw in text.splitlines():
            line = raw.split("#", 1)[0].strip()
            if not line:
                continue
            if "=" not in line:
                raise ValueError(f"cannot read equation line (no '='): {raw!r}")
            lhs, rhs = line.split("=", 1)
            m = _LHS_DDT.match(lhs) or _LHS.match(lhs)
            if not m:
                raise ValueError(f"cannot read left-hand side {lhs.strip()!r}; write it as x' = ...")
            eqs[m.group(1)] = rhs.strip()
        if not eqs:
            raise ValueError("no equations found")
        return cls.from_equations(eqs, variables=variables or list(eqs), source=source)

    @classmethod
    def load(cls, path) -> "Model":
        path = Path(path)
        text = path.read_text()
        if path.suffix == ".json":
            d = json.loads(text)
            m = cls.from_equations(d["equations"], variables=d["variables"], library=d.get("library"),
                                   source=d.get("source", str(path)))
            m.meta = d.get("meta", {})
            return m
        return cls.from_text(text, source=str(path))

    def save(self, path) -> Path:
        path = Path(path)
        d = {
            "variables": self.variables,
            "equations": self.equations(precision=None),
            "library": None if self.library is None else [str(t) for t in self.library],
            "source": self.source,
            "meta": self.meta,
        }
        path.write_text(json.dumps(d, indent=2))
        return path

    # ----------------------------------------------------------- queries
    @property
    def symbols(self):
        return [sp.Symbol(v) for v in self.variables]

    def active(self, i: int) -> list:
        """Indices of the terms used by equation ``i``."""
        return [j for j in range(len(self.terms)) if self.coefficients[i, j] != 0.0]

    def term_set(self, i: int) -> frozenset:
        return frozenset(self.terms[j] for j in self.active(i))

    def structure(self) -> tuple:
        return tuple(self.term_set(i) for i in range(len(self.variables)))

    def rhs(self, i: int, precision: int | None = 3) -> sp.Expr:
        expr = 0
        for j in self.active(i):
            c = self.coefficients[i, j]
            c = round(float(c), precision) if precision is not None else float(c)
            expr += sp.Float(c) * self.terms[j]
        return expr

    def equations(self, precision: int | None = 3) -> dict:
        out = {}
        for i, v in enumerate(self.variables):
            parts = []
            for j in self.active(i):
                c = self.coefficients[i, j]
                cs = repr(float(c)) if precision is None else f"{c:.{precision}g}"
                t = self.terms[j]
                parts.append(cs if t == 1 else f"{cs}*{t}")
            out[v] = " + ".join(parts) if parts else "0"
        return out

    def pretty(self, precision: int = 4) -> str:
        return "\n".join(
            f"{v}' = " + format_rhs([self.terms[j] for j in self.active(i)],
                                    [self.coefficients[i, j] for j in self.active(i)], precision)
            for i, v in enumerate(self.variables))

    def max_degree(self) -> int:
        degs = [term_degree(t, self.symbols) for t in self.terms]
        return max([d for d in degs if d >= 0], default=1)

    # ----------------------------------------------------------- numerics
    def theta(self, X: np.ndarray, terms=None) -> np.ndarray:
        """Evaluate terms on states ``X`` (T, n) -> (T, n_terms)."""
        return evaluate_terms(self.terms if terms is None else terms, self.variables, X)

    def predict(self, X: np.ndarray) -> np.ndarray:
        return self.theta(X) @ self.coefficients.T

    def rhs_function(self):
        fns = _lambdify(self.terms, self.variables)
        C = self.coefficients

        def f(_t, x):
            vals = np.array([fn(*x) for fn in fns], dtype=float)
            return C @ vals

        return f

    def with_coefficients(self, coefficients) -> "Model":
        return Model(list(self.variables), list(self.terms), np.asarray(coefficients, float),
                     library=self.library, source=self.source, meta=dict(self.meta))


def polynomial_library(variables, degree: int) -> list:
    syms = [sp.Symbol(v) for v in variables]
    out = [sp.Integer(1)]
    for d in range(1, degree + 1):
        for combo in itertools.combinations_with_replacement(syms, d):
            out.append(sp.Mul(*combo))
    return out


def candidate_library(model: Model, degree: int | None = None) -> list:
    """The library the discovery searched over, or a polynomial stand-in plus the model's own terms."""
    if model.library is not None:
        base = list(model.library)
    else:
        base = polynomial_library(model.variables, degree or max(2, model.max_degree()))
    seen = set(base)
    for t in model.terms:
        if t not in seen:
            base.append(t)
            seen.add(t)
    return base


_CACHE: dict = {}


def _lambdify(terms, variables):
    key = (tuple(terms), tuple(variables))
    if key not in _CACHE:
        syms = [sp.Symbol(v) for v in variables]
        _CACHE[key] = [sp.lambdify(syms, t, modules="numpy") for t in terms]
    return _CACHE[key]


def evaluate_terms(terms, variables, X) -> np.ndarray:
    X = np.atleast_2d(np.asarray(X, dtype=float))
    cols = []
    for fn in _lambdify(list(terms), list(variables)):
        with np.errstate(all="ignore"):
            v = np.asarray(fn(*X.T), dtype=float)
        cols.append(np.broadcast_to(v, (X.shape[0],)) if v.ndim == 0 else v)
    return np.column_stack(cols) if cols else np.zeros((X.shape[0], 0))


def _ordered_terms(terms, variables) -> list:
    syms = [sp.Symbol(v) for v in variables]
    uniq = list(dict.fromkeys(terms))

    def key(t):
        d = term_degree(t, syms)
        if d < 0:
            return (99, str(t))
        try:
            mono = sp.Poly(t, *syms).monoms()[0]
        except Exception:  # pragma: no cover
            mono = ()
        return (d, tuple(-m for m in mono))

    return sorted(uniq, key=key)
