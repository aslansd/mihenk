# Changelog

## Unreleased

These changes are on GitHub only. Nothing inside the wheel changed, so there
is no new PyPI release; they ship with the next version.

### Added
- `INSTALL.md`: extras, versions, installing on an unreliable connection
  (including offline installs from downloaded wheels), and editable installs.
- `TESTING.md`: the test suite, the synthetic validation and its pass rules,
  and what to do when a validation run fails.
- `PUBLISHING.md`: release steps, GitHub-only changes, and safe re-uploads.
- `ROADMAP.md`: what 0.1 did and learned, and the plan for 0.2 (experimental
  data, a calibration/test split, `mihenk design`, `mihenk hybrid`).
- A line in the README linking these documents.
- Project URLs in `pyproject.toml`. They reach PyPI with the next release.

### Fixed
- **Tests no longer skip the whole module when PySINDy is missing.** Only the
  6 tests that fit with PySINDy are skipped (`18 passed, 6 skipped`), so the
  core checks are tested on numpy, scipy and sympy alone. Before this, a fresh
  environment without PySINDy reported `1 skipped` and tested nothing.

## 0.1.0 — 2026-11-04

First release.

- Tier-0 checks from data and equation alone: `EQUIFINALITY`, `DERIVATIVE_ARTIFACT`, `WEAK_TERM`, `UNDERFIT`, `ROLLOUT_TRAP`.
- Tier-1 checks through PySINDy refits: `UNSTABLE_STRUCTURE` (bootstrap) and `FORKING_PATHS` (derivative method × threshold).
- Tier-2 ground-truth comparison for synthetic data.
- PySINDy 2.x adapter; equations from PySINDy models, text (including `model.print()` output), SymPy or JSON.
- CLI: `scan`, `fit`, `synth`, `validate` (with `--seeds`), `info`. HTML, JSON, Markdown and text reports.
- Optional daftar provenance (`--daftar`).
- 7 synthetic systems and 14 self-validation cases.
