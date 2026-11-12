# GenomeCanvas plotting skills

Agent skills for every GenomeCanvas Studio plot. Each folder holds a `SKILL.md` (front matter +
instructions) that Claude Code, or the Studio's future agent mode, can load. Each skill shows how to
call the plot and which parameters matter.

| Skill | Plots |
|---|---|
| [genomecanvas](genomecanvas/SKILL.md) | **start here**: survey the data, recommend figures, agree a plan with the user, run, review, report |
| [genomecanvas-plot](genomecanvas-plot/SKILL.md) | the runner (`gcplot.py`): request format, local vs server mode, shared style/export settings |
| [genomecanvas-tracks](genomecanvas-tracks/SKILL.md) | stacked tracks, overlaid signals, gene models |
| [genomecanvas-enrichment](genomecanvas-enrichment/SKILL.md) | heatmap, profile, profile + heatmap, clustered heatmap |
| [genomecanvas-relationships](genomecanvas-relationships/SKILL.md) | correlation, PCA, PairGrid over regions / genome bins / windows; Sample QC for large collections |
| [genomecanvas-differential](genomecanvas-differential/SKILL.md) | volcano, peak-center comparison |
| [genomecanvas-transform](genomecanvas-transform/SKILL.md) | difference, ratio, log₂ ratio, local background |

## Quick start

```bash
python3 skills/genomecanvas-plot/scripts/gcplot.py survey path/to/data --write /tmp/plan   # what is here, what to plot
python3 skills/genomecanvas-plot/scripts/gcplot.py list
python3 skills/genomecanvas-plot/scripts/gcplot.py run skills/genomecanvas-plot/examples/relationships-pca.json --out /tmp/pca
# same request against a running Studio (python3 web/server.py --port 8765)
python3 skills/genomecanvas-plot/scripts/gcplot.py run skills/genomecanvas-plot/examples/relationships-pca.json \
    --server http://127.0.0.1:8765 --project demo --out /tmp/pca_server
```

`examples/` has one request per plot on the bundled demo data (`tests/fixtures/trackplot_demo`).

## Two ways to call

* **Local**: `gcplot.py` imports `genomecanvas.studio.engine.run_plot` from `src/` and renders in
  process on local paths. Needs the studio extras (`pip install -e ".[studio]"`).
* **Agent / web server**: the same request goes to `POST /api/projects/<id>/jobs`. The server then
  applies authorisation, quotas, the job queue and worker isolation. An agent embedded in the web
  server should use this path rather than call the engine directly.

Both modes apply the Studio UI defaults (`ui_defaults` in `gcplot.py`, mirroring `effectiveSettings`
in `web/app.js`), so a request renders the same figure in either mode. If a UI default changes,
update `ui_defaults` too.

To make Claude Code load these automatically in this repo, link them into `.claude/skills/`:

```bash
mkdir -p .claude/skills && for d in skills/genomecanvas skills/genomecanvas-*; do ln -s "../../$d" .claude/skills/; done
```
