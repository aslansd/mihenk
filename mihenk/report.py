"""Findings and reports: what was found, where, how sure, and what was not checked."""

from __future__ import annotations

import html
import json
import textwrap
from dataclasses import asdict, dataclass, field
from pathlib import Path

import numpy as np

SEVERITIES = ("error", "warning", "info")
TIER_NAMES = {0: "data + equation", 1: "refit through the discovery tool", 2: "ground truth"}


@dataclass
class Finding:
    id: str
    severity: str
    tier: int
    what: str
    where: str = ""
    fix: str = ""
    evidence: dict = field(default_factory=dict)


@dataclass
class Check:
    """Outcome of one detector: flagged, clean, or not assessed (with why)."""

    id: str
    tier: int
    status: str  # "flagged" | "clean" | "not_assessed"
    reason: str = ""


@dataclass
class Report:
    name: str
    equation: str
    summary: str
    metrics: dict
    findings: list
    checks: list
    notes: list = field(default_factory=list)
    thresholds: dict = field(default_factory=dict)
    version: str = ""
    source: str = ""
    details: dict = field(default_factory=dict)   # tables for the HTML report

    # ------------------------------------------------------------- queries
    @property
    def has_errors(self) -> bool:
        return any(f.severity == "error" and f.tier < 2 for f in self.findings)

    def flagged(self, tier_max: int = 1) -> set:
        return {c.id for c in self.checks if c.status == "flagged" and c.tier <= tier_max}

    def ground_truth_flags(self) -> set:
        return {c.id for c in self.checks if c.status == "flagged" and c.tier == 2}

    def by_id(self, finding_id: str) -> list:
        return [f for f in self.findings if f.id == finding_id]

    def to_dict(self) -> dict:
        return {
            "mihenk_version": self.version,
            "name": self.name,
            "source": self.source,
            "equation": self.equation,
            "summary": self.summary,
            "metrics": self.metrics,
            "findings": [asdict(f) for f in self.findings],
            "checks": [asdict(c) for c in self.checks],
            "notes": self.notes,
            "thresholds": self.thresholds,
            "details": self.details,
        }

    def __getitem__(self, key):
        return self.to_dict()[key]

    FORMATS = ("json", "md", "txt", "html")

    def save(self, out_dir, formats=FORMATS) -> Path:
        """Write ``report.<format>`` for each of ``formats`` (default: json, md, txt, html)."""
        unknown = set(formats) - set(self.FORMATS)
        if unknown:
            raise ValueError(f"unknown report format(s): {', '.join(sorted(unknown))}; "
                             f"choose from {', '.join(self.FORMATS)}")
        out = Path(out_dir)
        out.mkdir(parents=True, exist_ok=True)
        writers = {
            "json": lambda: json.dumps(self.to_dict(), indent=2, default=_jsonable),
            "md": self.render_markdown,
            "txt": self.render_text,
            "html": self.render_html,
        }
        for fmt in formats:
            (out / f"report.{fmt}").write_text(writers[fmt]())
        return out

    # ------------------------------------------------------------ rendering
    def render_text(self) -> str:
        lines = [f"mihenk {self.version} -- {self.name}", f"  {self.summary}", "", "equation"]
        lines += ["  " + ln for ln in self.equation.splitlines()]
        lines += ["", "metrics"]
        for key, value in self.metrics.items():
            lines.append(f"  {key:<30} {_fmt(value)}")
        lines.append("")
        main = [f for f in self.findings if f.tier < 2]
        gt = [f for f in self.findings if f.tier == 2]
        if not main:
            lines.append("no findings")
        for f in _sorted(main):
            lines.append(f"{f.severity.upper():<9} {f.id}  [tier {f.tier}: {TIER_NAMES[f.tier]}]")
            if f.where:
                lines.append(f"          where: {f.where}")
            lines.append(_wrap(f.what, "          "))
            if f.fix:
                lines.append(_wrap("fix: " + f.fix, "          "))
            lines.append("")
        if gt:
            lines.append("ground truth says")
            for f in _sorted(gt):
                lines.append(_wrap(f"{f.id}: {f.what}", "  "))
            lines.append("")
        na = [c for c in self.checks if c.status == "not_assessed"]
        if na:
            lines.append("not assessed (a check that did not run is not a clean result)")
            for c in na:
                lines.append(f"  {c.id:<24} tier {c.tier}: {c.reason}")
            lines.append("")
        clean = sorted({c.id for c in self.checks if c.status == "clean" and c.tier < 2})
        if clean:
            lines.append("checked and clean: " + ", ".join(clean))
        for note in self.notes:
            lines.append(_wrap("note: " + note, ""))
        return "\n".join(lines).rstrip() + "\n"

    def render_markdown(self) -> str:
        out = [f"# mihenk report: {self.name}", "", f"*mihenk {self.version}*", "", self.summary, "",
               "## Equation", "", "```", self.equation, "```", "", "## Findings", ""]
        main = _sorted([f for f in self.findings if f.tier < 2])
        if not main:
            out.append("No findings.")
        for f in main:
            out.append(f"### {f.severity.upper()} `{f.id}` (tier {f.tier})")
            if f.where:
                out.append(f"*where:* {f.where}  ")
            out.append(f.what)
            if f.fix:
                out.append(f"\n*fix:* {f.fix}")
            out.append("")
        gt = [f for f in self.findings if f.tier == 2]
        if gt:
            out += ["## Ground truth", ""] + [f"- `{f.id}`: {f.what}" for f in _sorted(gt)] + [""]
        out += ["## Checks", "", "| check | tier | status | note |", "|---|---|---|---|"]
        out += [f"| {c.id} | {c.tier} | {c.status} | {c.reason} |" for c in self.checks]
        out += ["", "## Metrics", "", "| metric | value |", "|---|---|"]
        out += [f"| {k} | {_fmt(v)} |" for k, v in self.metrics.items()]
        if self.notes:
            out += ["", "## Notes", ""] + [f"- {n}" for n in self.notes]
        return "\n".join(out) + "\n"

    def render_html(self) -> str:
        return _HTML.format(
            name=_e(self.name), version=_e(self.version), summary=_e(self.summary),
            verdict_class="bad" if self.has_errors else ("warn" if any(
                f.severity == "warning" and f.tier < 2 for f in self.findings) else "ok"),
            equation=_e(self.equation),
            findings=_html_findings(self),
            checks="".join(
                f"<tr><td><code>{_e(c.id)}</code></td><td>{c.tier}</td>"
                f"<td class='st-{c.status}'>{_e(c.status.replace('_', ' '))}</td><td>{_e(c.reason)}</td></tr>"
                for c in self.checks),
            details=_html_details(self.details),
            metrics="".join(f"<tr><td>{_e(k)}</td><td>{_e(_fmt(v))}</td></tr>" for k, v in self.metrics.items()),
            notes="".join(f"<li>{_e(n)}</li>" for n in self.notes),
            thresholds="".join(f"<tr><td><code>{_e(k)}</code></td><td>{_e(_fmt(v))}</td></tr>"
                               for k, v in self.thresholds.items()),
        )


def _html_findings(r: Report) -> str:
    items = []
    for f in _sorted([f for f in r.findings if f.tier < 2]):
        tier = f"tier {f.tier}: {TIER_NAMES[f.tier]}"
        items.append(
            f"<div class='finding {f.severity}'><div class='fh'><span class='sev'>{f.severity}</span>"
            f"<code>{_e(f.id)}</code><span class='tier'>{_e(tier)}</span></div>"
            + (f"<div class='where'>where: {_e(f.where)}</div>" if f.where else "")
            + f"<p>{_e(f.what)}</p>"
            + (f"<p class='fix'><b>fix:</b> {_e(f.fix)}</p>" if f.fix else "")
            + "</div>")
    out = "".join(items) or "<p class='none'>No findings.</p>"
    gt = _sorted([f for f in r.findings if f.tier == 2])
    if gt:
        out += ("<h3>Ground truth (synthetic data only)</h3><ul class='gt'>"
                + "".join(f"<li><code>{_e(f.id)}</code> {_e(f.what)}</li>" for f in gt) + "</ul>")
    elif any(c.id == "RECOVERY" and c.status == "clean" for c in r.checks):
        out += "<h3>Ground truth (synthetic data only)</h3><p>The discovered equation matches the true law.</p>"
    return out


def _html_details(details: dict) -> str:
    out = []
    for title, table in details.items():
        if not table or not isinstance(table, dict) or "rows" not in table:
            continue
        head = "".join(f"<th>{_e(h)}</th>" for h in table["columns"])
        rows = "".join("<tr>" + "".join(f"<td>{_e(_fmt(c))}</td>" for c in row) + "</tr>"
                       for row in table["rows"])
        caption = f"<p class='cap'>{_e(table.get('caption', ''))}</p>" if table.get("caption") else ""
        out.append(f"<h3>{_e(title)}</h3>{caption}<table><tr>{head}</tr>{rows}</table>")
    return "".join(out)


_HTML = """<!DOCTYPE html>
<html lang="en"><head><meta charset="utf-8"><title>mihenk: {name}</title>
<meta name="viewport" content="width=device-width, initial-scale=1">
<style>
:root {{ --ink:#1d2330; --muted:#5b6475; --line:#e2e5ea; --bg:#fafaf8;
  --err:#b42318; --errbg:#fef3f2; --warn:#9a6700; --warnbg:#fffaeb; --info:#175cd3; --infobg:#eff8ff;
  --ok:#067647; --okbg:#ecfdf3; }}
body {{ font: 15px/1.5 system-ui, -apple-system, Segoe UI, sans-serif; color:var(--ink); background:var(--bg);
  max-width: 980px; margin: 2rem auto; padding: 0 1.2rem; }}
h1 {{ font-size: 1.5rem; margin: 0 0 .2rem; }} h2 {{ font-size:1.15rem; margin-top:2rem; border-bottom:1px solid var(--line); padding-bottom:.3rem; }}
h3 {{ font-size:1rem; margin-top:1.4rem; }}
.sub {{ color:var(--muted); margin:0 0 1rem; }}
.verdict {{ padding:.8rem 1rem; border-radius:8px; font-weight:600; }}
.verdict.bad {{ background:var(--errbg); color:var(--err); }} .verdict.warn {{ background:var(--warnbg); color:var(--warn); }}
.verdict.ok {{ background:var(--okbg); color:var(--ok); }}
pre {{ background:#fff; border:1px solid var(--line); border-radius:8px; padding:.8rem 1rem; overflow-x:auto; }}
.finding {{ border-left:4px solid; border-radius:6px; padding:.6rem 1rem; margin:.7rem 0; background:#fff; }}
.finding.error {{ border-color:var(--err); }} .finding.warning {{ border-color:var(--warn); }} .finding.info {{ border-color:var(--info); }}
.fh {{ display:flex; gap:.6rem; align-items:baseline; flex-wrap:wrap; }}
.sev {{ text-transform:uppercase; font-size:.75rem; font-weight:700; letter-spacing:.04em; }}
.error .sev {{ color:var(--err); }} .warning .sev {{ color:var(--warn); }} .info .sev {{ color:var(--info); }}
.tier, .where, .cap {{ color:var(--muted); font-size:.85rem; }}
.finding p {{ margin:.4rem 0; }} ul.gt {{ font-size:.88rem; color:var(--muted); }} .fix {{ color:var(--muted); }}
table {{ border-collapse:collapse; width:100%; background:#fff; font-size:.9rem; }}
td, th {{ border:1px solid var(--line); padding:.35rem .6rem; text-align:left; vertical-align:top; }}
th {{ background:#f3f4f6; }}
.st-flagged {{ color:var(--err); font-weight:600; }} .st-clean {{ color:var(--ok); }} .st-not_assessed {{ color:var(--muted); font-style:italic; }}
code {{ font-family: ui-monospace, SFMono-Regular, Menlo, monospace; font-size:.88em; }}
footer {{ color:var(--muted); font-size:.8rem; margin:3rem 0 1rem; }}
</style></head><body>
<h1>mihenk report: {name}</h1>
<p class="sub">mihenk {version} &middot; is the discovered equation real?</p>
<div class="verdict {verdict_class}">{summary}</div>
<h2>Equation</h2><pre>{equation}</pre>
<h2>Findings</h2>{findings}
<h2>Checks</h2>
<p class="cap">A check that did not run is not a clean result.</p>
<table><tr><th>check</th><th>tier</th><th>status</th><th>note</th></tr>{checks}</table>
<h2>Evidence</h2>{details}
<h2>Metrics</h2><table>{metrics}</table>
<h2>Notes</h2><ul>{notes}</ul>
<h2>Thresholds</h2><table>{thresholds}</table>
<footer>mihenk rules out named ways a discovered equation can be wrong. Passing every check does not prove the equation is the law of your system.</footer>
</body></html>
"""


def _e(s) -> str:
    return html.escape(str(s))


def _sorted(findings):
    order = {s: i for i, s in enumerate(SEVERITIES)}
    return sorted(findings, key=lambda f: (f.tier == 2, order.get(f.severity, 9), f.id))


def _fmt(v) -> str:
    if isinstance(v, (bool, np.bool_)):
        return str(bool(v))
    if isinstance(v, (int, np.integer)):
        return str(int(v))
    if isinstance(v, (float, np.floating)):
        v = float(v)
        if v != v:
            return "n/a"
        return f"{v:.4g}"
    if isinstance(v, (list, tuple)):
        return ", ".join(_fmt(x) for x in v)
    return str(v)


def _wrap(text: str, indent: str) -> str:
    return textwrap.fill(text, width=100, initial_indent=indent, subsequent_indent=indent)


def _jsonable(o):
    if isinstance(o, np.integer):
        return int(o)
    if isinstance(o, np.floating):
        return float(o)
    if isinstance(o, np.ndarray):
        return o.tolist()
    if isinstance(o, (set, frozenset)):
        return sorted(map(str, o))
    return str(o)
