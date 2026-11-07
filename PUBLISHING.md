# Publishing a release

## Steps

1. **Bump the version** in `mihenk/__init__.py` (`__version__`). It is the only
   place: `pyproject.toml` reads it dynamically.
2. **Write the changelog entry** in `CHANGELOG.md` before doing anything else.
   Move the `Unreleased` items under the new version and date. If the entry is
   hard to write, the release is not ready.
3. **Run the tests**, with PySINDy installed:
   ```bash
   python -m pytest -q                  # must be 25 passed, none skipped
   ```
4. **Re-run the validation numbers quoted in the README** if any detector or
   threshold changed:
   ```bash
   mihenk validate --seeds 20
   ```
   Paste the new table and recount the summary sentences under it: wrong laws
   flagged, correct-law errors, correct-law warnings. Do not leave old numbers
   next to new code.
5. **Build and check**:
   ```bash
   rm -rf dist build mihenk.egg-info
   python -m build
   twine check dist/*
   ```
6. **Install the wheel in a fresh virtualenv** and run it there. Tests in the
   source tree cannot catch packaging mistakes such as a missing module or a
   wrong entry point. Check twice: first the core alone, then with PySINDy.
   ```bash
   python -m venv /tmp/mihenk-check && /tmp/mihenk-check/bin/pip install dist/mihenk-*.whl
   cd /tmp && /tmp/mihenk-check/bin/mihenk info
   /tmp/mihenk-check/bin/mihenk synth lorenz -o /tmp/l.npz
   printf "x' = -10 x + 10 y\ny' = 28 x - y - x z\nz' = -2.6667 z + x y\n" > /tmp/eq.txt
   /tmp/mihenk-check/bin/mihenk scan /tmp/l.npz --model /tmp/eq.txt   # exit 0, tier 1 not assessed

   /tmp/mihenk-check/bin/pip install pysindy
   /tmp/mihenk-check/bin/mihenk fit /tmp/l.npz -o /tmp/m.json --threshold 0.5
   /tmp/mihenk-check/bin/mihenk scan /tmp/l.npz --model /tmp/m.json  # tier 1 now runs
   ```
7. **Upload, tag, push**:
   ```bash
   twine upload dist/*
   git tag -a v0.1.1 -m "mihenk 0.1.1"
   git push && git push --tags
   ```
8. **Verify from PyPI**. The index can serve the old version for a few minutes:
   ```bash
   pip install --no-cache-dir -U mihenk && mihenk --version
   ```

## Changes that go to GitHub only

Documentation (`*.md`), `tests/` and `examples/` are not inside the wheel.
Changing them needs a commit and a push, but no version bump and no PyPI
upload. Record them under `Unreleased` in the changelog so the next release
carries them.

Two exceptions:

- The **README is the PyPI project page**. A README edit shows on GitHub at
  once, but on PyPI only with the next release.
- **`pyproject.toml` metadata** works the same way: project URLs,
  classifiers and keywords reach PyPI only with the next release.

## Uploading on an unreliable connection

If `twine upload` breaks partway through, it is safe to run it again:

```bash
twine upload --skip-existing dist/*
```

- `--skip-existing` skips files PyPI already accepted, and uploads the rest.
- If PyPI answers `File already exists` for a file you are not sure
  got through, it did. Check the release page before building anything new.
- Never rebuild and re-upload the same version with different contents.
  PyPI refuses it (see below).

## Cautions

- A version number can never be reused on PyPI, even after deleting the
  release. If a release is broken, ship the next patch version.
- Use a project-scoped API token for `mihenk`, not an account-wide one. If you
  keep it in `~/.pypirc`, make that file readable only by you (`chmod 600`).
- TestPyPI needs a separate account and token. For a small pure-Python package,
  the fresh-virtualenv check in step 6 catches the same problems.
- `dist/` is tracked in this repository. After step 5, commit the new files
  in `dist/`, and remove the old ones, so `dist/` and PyPI hold the same
  version.

## Licensing

`LICENSE` holds the full Apache-2.0 text, and `pyproject.toml` declares it as
an SPDX expression (`license = "Apache-2.0"`). mihenk builds on PySINDy (MIT)
through its public API. It ships no PySINDy code. The papers it builds on are
cited in the README.
