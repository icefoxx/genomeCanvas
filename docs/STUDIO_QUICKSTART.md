# GenomeCanvas Studio: quick start

Copy-paste recipes for running the studio and its command-line tools. For production
deployment (HTTPS, SSO, sizing, backups), see [STUDIO_DEPLOY.md](STUDIO_DEPLOY.md).

## 0. Install

```bash
git clone https://github.com/<you>/genomeCanvas.git && cd genomeCanvas
python3 -m venv .venv && source .venv/bin/activate      # Python ≥ 3.11
pip install -e ".[studio]"                              # engine + pyBigWig, scipy, scikit-learn, plotly
gcanvas-studio --help
```

## 1. Personal use: one person, own data, no sign-in

```bash
# Built-in public demo data only
python3 web/server.py --port 8765

# Your own folder of BigWig / bedGraph / BED / GTF files (read in place, no upload limit)
python3 web/server.py --port 8765 --data-dir ~/data/chipseq_project
```

Open <http://127.0.0.1:8765>. The project menu shows **Local · chipseq_project** and the demo
projects. Pick a plot type, tick files, then press **Plot** (⌘/Ctrl + Enter).

Tips:
* Drop extra files anywhere on the page to upload them.
* **Sample relationships** → *Summarise signal over* → **Genome-wide bins** gives
  multiBigwigSummary-style correlation/PCA without a BED file.
* PCA opens in an **Interactive** view: hover for details, drag labels. *Download* also offers an
  interactive HTML file.
* Every figure ships with `methods.txt`. Use *Download → Copy methods paragraph*.

## 2. Small multi-user setup on one machine (password accounts)

```bash
mkdir -p ~/studio && cd ~/studio
gcanvas-studio example-config > studio.toml
# edit studio.toml: auth = "local", data_root = "/path/to/deliveries", workspace = "/path/to/studio-data",
#                   public_url = "http://<this-host>:8765", host = "0.0.0.0" (only on a trusted network)

gcanvas-studio --config studio.toml user add you@core.edu --admin     # prints an invitation link
gcanvas-studio --config studio.toml serve
```

Open the invitation link, choose a password, and you land in the studio. The **Admin console**
is at `/admin`. Over the internet, put it behind HTTPS ([STUDIO_DEPLOY.md](STUDIO_DEPLOY.md)).

## 3. Giving a client a project

Prepare the delivery folder (labels, groups and colours come from `samples.tsv`):

```text
/path/to/deliveries/smith/P2026-017/
├── samples.tsv      # see deploy/example-project/samples.tsv
├── project.yaml     # optional: name, genome, defaults, templates
├── bigwig/*.bw
└── peaks/*.narrowPeak
```

Then, from the admin console or the CLI:

```bash
gcs() { gcanvas-studio --config studio.toml "$@"; }   # works in bash and zsh
gcs project add P2026-017 --name "Smith lab · H3K27ac" --dir smith/P2026-017 --genome mm10 --expires 2027-04-30
gcs project member P2026-017 jane@uni.edu --role editor     # new address → account + invitation link
gcs project member P2026-017 pi@uni.edu --role viewer       # can open/download saved figures
gcs project scan P2026-017                                  # rerun after changing delivered files
gcs project list
```

## 4. Everyday admin commands

```bash
gcs user list
gcs user invite jane@uni.edu          # new link (also a password reset)
gcs user disable former@uni.edu
gcs project archive P2026-017         # access ends, nothing is deleted
gcs project activate P2026-017
gcs project scan                      # all active projects + shared annotations
gcs cleanup                           # apply retention rules now (also runs hourly)
```

## 5. Docker (production)

```bash
cd deploy
cp studio.toml studio.local.toml     # edit public_url, auth, paths
mkdir certs                          # fullchain.pem + privkey.pem
docker compose up -d --build
docker compose exec studio gcanvas-studio --config /etc/genomecanvas/studio.toml user add you@core.edu --admin
```

Single sign-on: `auth = "proxy"` in `studio.local.toml`, mount `nginx-sso.conf`, fill in
`oauth2-proxy.cfg`, then run `docker compose --profile sso up -d`.

## 6. The figure engine without the web UI

The same engine powers the original command line:

```bash
gcanvas init figure.yaml                   # write a starter figure specification
gcanvas validate figure.yaml
gcanvas render figure.yaml --output figure.pdf
```

It can also be used from Python, for example to compute a genome-wide 10 kb correlation table:

```python
from genomecanvas.pipeline.genome_bins import genome_summary
summary, report = genome_summary(
    {"WT": "WT.bw", "KO": "KO.bw"}, size=10_000, skip_zeros=True, blacklist="mm10-blacklist.bed",
)
print(report.windows_kept, summary.values.shape)    # bins × samples
```

## 7. Configuration via environment variables

Every `studio.toml` key can be set as `GC_STUDIO_<KEY>`, which is handy for containers:

```bash
GC_STUDIO_AUTH=local GC_STUDIO_WORKSPACE=/srv/genomecanvas GC_STUDIO_WORKERS=4 \
  gcanvas-studio serve --host 127.0.0.1 --port 8765
```
