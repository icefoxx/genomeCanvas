"""Render docs/PROFILESEG_REPORT.md as a self-contained HTML page (no external assets).

    python benchmarks/profileseg_report_html.py [docs/PROFILESEG_REPORT.md] [docs/PROFILESEG_REPORT.html]
"""

from __future__ import annotations

import sys
from pathlib import Path

import markdown

CSS = """
:root { --bg:#ffffff; --fg:#1f2328; --muted:#59636e; --line:#d1d9e0; --head:#f6f8fa; --accent:#0969da; --code:#eff1f3; }
@media (prefers-color-scheme: dark) { :root:not([data-theme="light"]) {
  --bg:#0d1117; --fg:#e6edf3; --muted:#9198a1; --line:#3d444d; --head:#151b23; --accent:#4493f8; --code:#262c36; } }
:root[data-theme="dark"] { --bg:#0d1117; --fg:#e6edf3; --muted:#9198a1; --line:#3d444d; --head:#151b23; --accent:#4493f8; --code:#262c36; }
* { box-sizing: border-box; }
body { margin:0; background:var(--bg); color:var(--fg);
  font: 15px/1.6 -apple-system, BlinkMacSystemFont, "Segoe UI", Helvetica, Arial, sans-serif; }
main { max-width: 1040px; margin: 0 auto; padding: 32px 16px 64px; }
h1 { font-size: 1.9rem; margin: 0 0 .4rem; }
h2 { font-size: 1.35rem; margin: 2.2rem 0 .6rem; padding-bottom: .3rem; border-bottom: 1px solid var(--line); }
h3 { font-size: 1.1rem; margin: 1.6rem 0 .5rem; }
p, li { max-width: 72ch; }
em { color: var(--muted); }
code { background: var(--code); padding: .1em .35em; border-radius: 4px; font-size: .9em; }
a { color: var(--accent); }
.table { overflow-x: auto; margin: .8rem 0 1.2rem; border: 1px solid var(--line); border-radius: 6px; }
table { border-collapse: collapse; width: 100%; font-size: .9rem; }
th, td { text-align: left; vertical-align: top; padding: 7px 10px; border-bottom: 1px solid var(--line); }
th { background: var(--head); font-weight: 600; white-space: nowrap; }
tr:last-child td { border-bottom: none; }
td strong { color: var(--fg); }
"""


def render(md_path: Path, html_path: Path) -> None:
    text = md_path.read_text(encoding="utf-8")
    body = markdown.markdown(text, extensions=["tables", "fenced_code"])
    heading = next((ln[2:].strip() for ln in text.splitlines() if ln.startswith("# ")), "ProfileSeg Evaluation")
    title = heading if len(heading) <= 40 or " — " not in heading else heading.split(" — ")[-1].strip().capitalize()
    body = body.replace("<table>", '<div class="table"><table>').replace("</table>", "</table></div>")
    html = f"""<!doctype html>
<html lang="en">
<head>
<meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1">
<title>{title}</title>
<style>{CSS}</style>
</head>
<body><main>
{body}
</main></body>
</html>
"""
    html_path.write_text(html, encoding="utf-8")


if __name__ == "__main__":
    src = Path(sys.argv[1]) if len(sys.argv) > 1 else Path("docs/PROFILESEG_REPORT.md")
    dst = Path(sys.argv[2]) if len(sys.argv) > 2 else src.with_suffix(".html")
    render(src, dst)
    print(dst)
