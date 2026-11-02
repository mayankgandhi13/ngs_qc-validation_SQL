"""Rebuild ngs_lab.db, seed it, run every query in queries.sql, and write report.html.

Usage: python run.py
"""

import html
import sqlite3
from datetime import date
from pathlib import Path

import seed_data

ROOT = Path(__file__).resolve().parent
DB_PATH = ROOT / "ngs_lab.db"
SCHEMA_PATH = ROOT / "schema.sql"
QUERIES_PATH = ROOT / "queries.sql"
REPORT_PATH = ROOT / "report.html"


# --------------------------------------------------------------------------
# query file parsing and execution
# --------------------------------------------------------------------------

def parse_queries(text):
    """Split queries.sql into blocks on '-- name: <id> | <title>' lines."""
    blocks, current = [], None
    for line in text.splitlines():
        if line.startswith("-- name:"):
            qid, _, title = line[len("-- name:"):].partition("|")
            current = {"id": qid.strip(), "title": title.strip(), "question": [], "sql": []}
            blocks.append(current)
        elif current is None:
            continue
        elif line.startswith("--") and not current["sql"]:
            current["question"].append(line[2:].strip())
        else:
            current["sql"].append(line)
    for b in blocks:
        b["question"] = " ".join(b["question"])
        b["sql"] = "\n".join(b["sql"]).strip()
    return blocks


def run_query(conn, block):
    cur = conn.execute(block["sql"])
    cols = [d[0] for d in cur.description]
    return cols, cur.fetchall()


# --------------------------------------------------------------------------
# terminal output
# --------------------------------------------------------------------------

def fmt(value):
    if value is None:
        return "-"
    if isinstance(value, float):
        return f"{value:g}" if abs(value) < 1e6 else f"{value:,.0f}"
    return str(value)


def print_table(block, cols, rows):
    print(f"\n== {block['id']} | {block['title']}  ({len(rows)} rows)")
    if not rows:
        print("   (no rows)")
        return
    cells = [[fmt(v) for v in r] for r in rows]
    widths = [max(len(c), *(len(r[i]) for r in cells)) for i, c in enumerate(cols)]
    numeric = [all(isinstance(r[i], (int, float)) or r[i] is None for r in rows) for i in range(len(cols))]

    def line(values):
        return "  ".join(v.rjust(w) if num else v.ljust(w) for v, w, num in zip(values, widths, numeric))

    print(line(cols))
    print("  ".join("-" * w for w in widths))
    for r in cells:
        print(line(r))


# --------------------------------------------------------------------------
# HTML report
# --------------------------------------------------------------------------

CSS = """
:root {
  --bg: #f4f6f8; --surface: #ffffff; --surface-2: #eef2f5; --border: #d5dde4;
  --text: #1f2a33; --muted: #5d6b78; --accent: #3d6a8a; --accent-soft: #dfe9f1;
  --fail: #a63d3d; --fail-bg: #f7e3e1; --row-fail: #fbefed; --pass: #2f6b4f;
  --mono: ui-monospace, SFMono-Regular, Menlo, Consolas, monospace;
}
@media (prefers-color-scheme: dark) {
  :root {
    --bg: #141a1f; --surface: #1c242b; --surface-2: #232d35; --border: #33414c;
    --text: #dbe3ea; --muted: #93a3b1; --accent: #8db4d1; --accent-soft: #24384a;
    --fail: #f0a29a; --fail-bg: #4a2522; --row-fail: #2f2020; --pass: #8fd0ae;
  }
}
* { box-sizing: border-box; }
body {
  margin: 0; background: var(--bg); color: var(--text);
  font: 15px/1.5 system-ui, -apple-system, "Segoe UI", Roboto, sans-serif;
}
header { background: var(--surface); border-bottom: 1px solid var(--border); }
.wrap { max-width: 1180px; margin: 0 auto; padding: 0 16px; }
header .wrap { padding-top: 24px; padding-bottom: 16px; }
h1 { margin: 0 0 4px; font-size: 1.5rem; letter-spacing: -0.01em; }
.sub { margin: 0; color: var(--muted); font-size: 0.9rem; }
.stats { display: grid; grid-template-columns: repeat(auto-fit, minmax(140px, 1fr)); gap: 10px; margin: 18px 0 14px; }
.stat { background: var(--surface-2); border: 1px solid var(--border); border-radius: 8px; padding: 10px 12px; }
.stat b { display: block; font-size: 1.35rem; font-variant-numeric: tabular-nums; }
.stat span { color: var(--muted); font-size: 0.8rem; text-transform: uppercase; letter-spacing: 0.04em; }
.stat.alert b { color: var(--fail); }
nav { display: flex; flex-wrap: wrap; gap: 6px; }
nav a {
  color: var(--accent); background: var(--accent-soft); text-decoration: none;
  padding: 4px 10px; border-radius: 999px; font-size: 0.85rem;
}
nav a:hover { text-decoration: underline; }
main .wrap { padding-top: 8px; padding-bottom: 40px; }
section {
  background: var(--surface); border: 1px solid var(--border); border-radius: 10px;
  padding: 18px 18px 12px; margin-top: 18px; scroll-margin-top: 12px;
}
h2 { margin: 0 0 4px; font-size: 1.1rem; }
h2 code { font-family: var(--mono); color: var(--muted); font-size: 0.85rem; font-weight: normal; margin-left: 6px; }
.question { margin: 0 0 12px; color: var(--muted); }
.table-wrap { overflow-x: auto; -webkit-overflow-scrolling: touch; border: 1px solid var(--border); border-radius: 6px; }
table { border-collapse: collapse; width: 100%; font-size: 0.85rem; }
th, td { padding: 6px 10px; border-bottom: 1px solid var(--border); white-space: nowrap; text-align: left; }
th { background: var(--surface-2); position: sticky; top: 0; font-weight: 600; }
tbody tr:last-child td { border-bottom: 0; }
tbody tr:hover td { background: var(--surface-2); }
td.num { text-align: right; font-variant-numeric: tabular-nums; }
td.null { color: var(--muted); }
td.fail { color: var(--fail); background: var(--fail-bg); font-weight: 600; }
td.pass { color: var(--pass); }
tr.row-fail td { background: var(--row-fail); }
.count { color: var(--muted); font-size: 0.8rem; margin: 8px 0 0; }
details { margin-top: 10px; }
summary { cursor: pointer; color: var(--accent); font-size: 0.85rem; }
pre {
  background: var(--surface-2); border: 1px solid var(--border); border-radius: 6px;
  padding: 12px; overflow-x: auto; font: 0.8rem/1.45 var(--mono); margin: 8px 0 0;
}
footer { color: var(--muted); font-size: 0.8rem; text-align: center; padding: 0 16px 32px; }
"""


def cell_html(value):
    text = html.escape(fmt(value))
    classes = []
    if isinstance(value, (int, float)):
        classes.append("num")
    if value is None:
        classes.append("null")
    if value == "FAIL":
        classes.append("fail")
    elif value == "PASS":
        classes.append("pass")
    cls = f' class="{" ".join(classes)}"' if classes else ""
    return f"<td{cls}>{text}</td>"


def table_html(cols, rows):
    if not rows:
        return '<p class="count">No rows returned.</p>'
    fail_idx = cols.index("pct_fail") if "pct_fail" in cols else None
    head = "".join(f"<th>{html.escape(c)}</th>" for c in cols)
    body = []
    for r in rows:
        row_cls = ' class="row-fail"' if fail_idx is not None and (r[fail_idx] or 0) >= 50 else ""
        body.append(f"<tr{row_cls}>" + "".join(cell_html(v) for v in r) + "</tr>")
    return (f'<div class="table-wrap"><table><thead><tr>{head}</tr></thead>'
            f'<tbody>{"".join(body)}</tbody></table></div>'
            f'<p class="count">{len(rows)} row{"s" if len(rows) != 1 else ""}</p>')


def summary_stats(conn):
    one = lambda sql: conn.execute(sql).fetchone()[0]
    return [
        ("Runs", one("SELECT COUNT(*) FROM runs"), False),
        ("Samples", one("SELECT COUNT(*) FROM samples"), False),
        ("Patients", one("SELECT COUNT(*) FROM patients"), False),
        ("Variant calls", one("SELECT COUNT(*) FROM variants"), False),
        ("QC failures", one("SELECT COUNT(*) FROM v_sample_qc WHERE qc_pass = 0"), True),
        ("Validation reps", one("SELECT COUNT(*) FROM samples WHERE is_validation = 1"), False),
    ]


def write_report(conn, results):
    first, last = conn.execute("SELECT MIN(run_date), MAX(run_date) FROM runs").fetchone()
    stats = "".join(
        f'<div class="stat{" alert" if alert and value else ""}"><b>{value}</b><span>{label}</span></div>'
        for label, value, alert in summary_stats(conn))
    nav = "".join(f'<a href="#{b["id"]}">{html.escape(b["title"])}</a>' for b, _, _ in results)
    sections = "".join(f"""
<section id="{b['id']}">
  <h2>{html.escape(b['title'])}<code>{b['id']}</code></h2>
  <p class="question">{html.escape(b['question'])}</p>
  {table_html(cols, rows)}
  <details><summary>Show SQL</summary><pre><code>{html.escape(b['sql'])}</code></pre></details>
</section>""" for b, cols, rows in results)

    page = f"""<!doctype html>
<html lang="en">
<head>
<meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1">
<title>NGS Lab Report</title>
<style>{CSS}</style>
</head>
<body>
<header><div class="wrap">
  <h1>Clinical NGS Lab: QC &amp; Validation Report</h1>
  <p class="sub">Synthetic data only, no PHI. Runs {first} to {last}. Generated {date.today().isoformat()}.</p>
  <div class="stats">{stats}</div>
  <nav>{nav}</nav>
</div></header>
<main><div class="wrap">{sections}
</div></main>
<footer>clinical-ngs-lab-db: SQLite + Python standard library</footer>
</body>
</html>
"""
    REPORT_PATH.write_text(page, encoding="utf-8")


# --------------------------------------------------------------------------

def main():
    if DB_PATH.exists():
        DB_PATH.unlink()
    conn = sqlite3.connect(DB_PATH)
    conn.execute("PRAGMA foreign_keys = ON")
    conn.executescript(SCHEMA_PATH.read_text(encoding="utf-8"))
    counts = seed_data.seed(conn)
    print(f"Built {DB_PATH.name}: {counts['samples']} samples, {counts['variants']} variant calls, "
          f"{counts['truth']} truth rows")

    results = []
    for block in parse_queries(QUERIES_PATH.read_text(encoding="utf-8")):
        cols, rows = run_query(conn, block)
        print_table(block, cols, rows)
        results.append((block, cols, rows))

    write_report(conn, results)
    conn.close()
    print(f"\nWrote {REPORT_PATH.name} ({len(results)} sections)")


if __name__ == "__main__":
    main()
