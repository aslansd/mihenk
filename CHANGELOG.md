# Changelog

## 0.1.0

First release.

- Tier-0 checks from data and equation alone: `EQUIFINALITY`, `DERIVATIVE_ARTIFACT`, `WEAK_TERM`, `UNDERFIT`, `ROLLOUT_TRAP`.
- Tier-1 checks through PySINDy refits: `UNSTABLE_STRUCTURE` (bootstrap) and `FORKING_PATHS` (derivative method × threshold).
- Tier-2 ground-truth comparison for synthetic data.
- PySINDy 2.x adapter; equations from PySINDy models, text (including `model.print()` output), SymPy or JSON.
- CLI: `scan`, `fit`, `synth`, `validate` (with `--seeds`), `info`. HTML, JSON, Markdown and text reports.
- Optional daftar provenance (`--daftar`).
- 7 synthetic systems and 14 self-validation cases.
