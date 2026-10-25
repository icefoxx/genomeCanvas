# Deploying GenomeCanvas Studio for a sequencing core

GenomeCanvas Studio is a web server that lets clients make publication-ready figures from
the data a core facility delivered to them: genome tracks, TSS/peak heatmaps and profiles,
sample correlation and PCA (on regions or genome-wide bins), and signal transformations.

This guide covers the multi-user deployment. For a personal, single-user copy run
`python3 web/server.py --data-dir /path/to/bigwigs` and open http://127.0.0.1:8765.

## 1. How it fits together

```
client browser ──HTTPS──▶ nginx ──▶ (oauth2-proxy, optional SSO)
                            │
                            ▼
                    gcanvas-studio serve  (HTTP API + static UI, port 8765)
                     │            │
               SQLite db     worker processes × N  (one figure each, timeout + memory cap)
                     │            │
   /srv/genomecanvas (rw)    /data/deliverables (ro)   /data/annotations (ro)
   uploads, results, cache   one folder per project    shared GTFs per genome
```

* **Project = isolation boundary.** A client sees only projects they are a member of
  (plus the public demo projects, if enabled). Every file, job, result and saved figure is
  checked against project membership on every request.
* **Roles:** administrators (core staff) see everything; project *editors* can plot, upload
  and save figures; *viewers* (for example a PI or collaborator) can open and download saved figures.
* **Deliverables are read in place** from a read-only mount. Clients may upload extra files
  into a per-project quota; uploads expire after `upload_retention_days`.
* **Figures run in separate worker processes** with a per-job timeout and an address-space
  limit, so one heavy request cannot take the service down. Region matrices and summaries are
  cached, so restyling a figure (colours, fonts, size) is near-instant.

## 2. Quick start with Docker Compose

Requirements: a Linux VM with Docker, read access to the delivery storage, a DNS name and a
TLS certificate.

```sh
git clone <this repository> genomecanvas && cd genomecanvas/deploy
cp studio.toml studio.local.toml        # edit public_url, paths, auth, workers
mkdir certs && cp /path/to/fullchain.pem /path/to/privkey.pem certs/
# edit the two host paths under services.studio.volumes in docker-compose.yml
docker compose up -d --build
docker compose exec studio gcanvas-studio --config /etc/genomecanvas/studio.toml \
    user add you@core.example.edu --admin          # prints an invitation link
```

Open the invitation link, set a password, and you are in the **Admin console** (`/admin`).

### Single sign-on instead of passwords

1. Register an OIDC application with your identity provider (Entra ID/Azure AD, Google,
   Okta, Keycloak, CILogon, ...). The redirect URL is `https://<host>/oauth2/callback`.
2. `cp oauth2-proxy.cfg.example oauth2-proxy.cfg` and fill in the issuer, client ID/secret
   and a cookie secret.
3. In `studio.local.toml` set `auth = "proxy"`. In `docker-compose.yml` mount
   `./nginx-sso.conf` instead of `./nginx.conf`.
4. `docker compose --profile sso up -d`

The first sign-in creates an account with **no project access**; an administrator then adds
the person to projects. Identity headers are accepted only from `trusted_proxies`, and nginx
overwrites any identity header a client sends.

### Without containers (systemd)

```sh
python3 -m venv /opt/genomecanvas/venv
/opt/genomecanvas/venv/bin/pip install -e "/opt/genomecanvas/src[studio]"
gcanvas-studio example-config > /etc/genomecanvas/studio.toml   # then edit
cp deploy/genomecanvas-studio.service /etc/systemd/system/
systemctl enable --now genomecanvas-studio
```

Put nginx in front using `deploy/nginx.conf`, with `proxy_pass http://127.0.0.1:8765`.

## 3. Preparing project data

Each project points at one folder under `data_root`. The **manifest** `samples.tsv` decides
what clients see and how it is labelled; without a manifest, every supported file in the
folder is listed using its file name.

```
/data/core/deliverables/smith/P2026-017/
├── samples.tsv          # optional but recommended
├── project.yaml         # optional: name, genome, default settings, templates
├── bigwig/WT_rep1.bw ...
└── peaks/consensus_peaks.narrowPeak
```

`samples.tsv` is tab-separated. Lines starting with `#` are ignored:

| column | required | meaning |
|---|---|---|
| `file` | yes | path relative to the project folder |
| `label` | | name shown in figures (default: file name) |
| `group` | | condition/replicate group, used for volcano plots and PCA markers |
| `color` | | sample colour, e.g. `#1F6F68` (default: palette) |
| `cmap` | | per-sample heatmap colour scale (Matplotlib name or `#fff,#08306b`) |
| `role` | | `signal`, `regions` or `genes` (default from the extension) |
| `track` | | `chromhmm` to draw a BED as chromatin states |

Supported formats: BigWig (`.bw`, `.bigwig`), bedGraph (`.bedgraph`, `.bdg`, `.bg`),
BED / narrowPeak / broadPeak, and GTF. Compressed `.gz` text files must be decompressed.
Rows that point to missing files are reported in the scan result rather than silently dropped.

`project.yaml` example: see `deploy/example-project/project.yaml`. Templates defined there
appear under **Settings › Templates** for that project and are refreshed on every scan.

**Shared annotations** (`shared_root`) are organised as `<genome>/<file>`, for example
`mm10/gencode.vM25.basic.gtf`. They are offered to every project whose genome matches.

## 4. Day-to-day administration

The admin console (`/admin`) covers projects, membership, invitations, templates, usage and
the activity log. Everything is also available from the command line, which suits LIMS
automation:

```sh
gcs() { gcanvas-studio --config /etc/genomecanvas/studio.toml "$@"; }   # bash and zsh
gcs project add P2026-017 --name "Smith lab · H3K27ac" --dir smith/P2026-017 \
    --genome mm10 --expires 2027-04-30 --quota-gb 20
gcs project member P2026-017 jane.smith@uni.edu --role editor   # new users get an invite link
gcs project member P2026-017 pi@uni.edu --role viewer
gcs project scan P2026-017        # after adding or changing delivered files
gcs project scan                  # every active project + shared annotations
gcs project archive P2026-017     # access ends; data and saved figures are kept
gcs user invite jane.smith@uni.edu  # new link = password reset
gcs user disable former@uni.edu
gcs cleanup                       # apply retention now (also runs hourly)
```

When a project's `expires` date passes, members lose access automatically, without
deleting anything.

## 5. Configuration reference

Settings are read from the TOML file (`[studio]` table). Any key can be overridden with
an environment variable `GC_STUDIO_<KEY>`, or on the command line for `serve`
(`--host`, `--port`, `--workers`).

| key | default | notes |
|---|---|---|
| `host`, `port` | `127.0.0.1`, `8765` | keep on localhost, or on the container network, behind nginx |
| `public_url` | — | base URL for invitation links; also accepted as Origin |
| `secure_cookies` | `false` | set `true` behind HTTPS |
| `auth` | `none` | `none` (single local user), `local` (accounts), `proxy` (SSO header) |
| `proxy_user_header` / `proxy_name_header` | `X-Auth-Request-Email` / `-User` | headers set by oauth2-proxy |
| `trusted_proxies` | `127.0.0.1`, `::1` | IPs or CIDRs; identity and `X-Forwarded-For` are only trusted from these |
| `proxy_auto_create` | `true` | create accounts on first SSO login (no project access) |
| `session_hours` | `12` | password-session lifetime |
| `data_root` | — | read-only root of project folders |
| `shared_root` | — | read-only `<genome>/<annotation>` folders |
| `workspace` | system temp | database, uploads, results, cache (back this up) |
| `static_dir` | repository `web/` | UI files |
| `demo_project` | `true` | public demo projects for every user |
| `workers` | `2` | parallel figures (CPU cores) |
| `job_timeout` | `300` | seconds per figure |
| `worker_memory_mb` | `4096` | per-worker address-space cap (Linux) |
| `max_jobs_per_user`, `max_queue` | `3`, `64` | fairness and back-pressure |
| `max_upload_mb` | `2048` | per file; match nginx `client_max_body_size` |
| `default_quota_gb` | `20` | per-project upload quota |
| `job_retention_days` | `14` | results deleted afterwards unless saved as a figure |
| `upload_retention_days` | `30` | client uploads deleted afterwards |
| `cache_max_gb` | `20` | matrix/summary cache, least-recently-used eviction |

## 6. Sizing

The numbers below were measured on Apple-silicon, single-threaded per figure, and include
figure rendering.

| task (3 samples, mouse genome) | time | peak memory |
|---|---|---|
| Correlation/PCA on 10 k promoters | ~2.5 s | ~100 MB |
| TSS heatmap, 10 k rows, ±3 kb, 50 bp | ~8 s (restyle ~4 s) | ~300 MB |
| Genome-wide 10 kb bins, dense 400 MB BigWigs | ~23 s (cached reruns ~1–2 s) | ~350 MB |
| Genome-wide 1 kb bins, same data | ~26 s | ~1.3 GB |

As a starting point: `workers` = number of cores minus one, `worker_memory_mb` = 3–4 GB, and
container memory ≥ `workers × worker_memory_mb + 2 GB`. Fine genome-wide bins (≤ 1 kb) on many
samples are the most memory-hungry request; the studio refuses requests whose table would
exceed its limits and suggests larger bins.

## 7. Backups, upgrades, monitoring

* **Back up** `workspace/studio.db` (SQLite; safe to copy with `sqlite3 studio.db ".backup out.db"`)
  and, if saved figures matter, `workspace/jobs/`. Uploads and the cache can be regenerated or are
  short-lived.
* **Upgrade:** `git pull && docker compose up -d --build`. The schema is created idempotently on
  start-up. Jobs that were running when the server stopped are marked as interrupted.
* **Health:** `GET /api/health` (used by the container health check). The activity log (admin
  console, or the `audit` table) records sign-ins, uploads, jobs, figure saves and admin actions,
  with the client IP taken from `X-Forwarded-For` when the request comes through a trusted proxy.

## 8. Security notes

* Always serve over HTTPS (`secure_cookies = true`). The server itself speaks plain HTTP and
  is meant to sit behind nginx.
* Sessions are HttpOnly, SameSite=Lax cookies. State-changing requests require a custom header
  plus a same-origin check (CSRF). Pages ship with a strict Content-Security-Policy.
* File paths never come from the browser: deliverables are registered by scans confined to
  `data_root`. Uploads are written under the workspace with sanitised names and are
  format-checked before they can be used.
* Sign-in throttling counts failures per account+IP (10 per 10 min) and per IP (60).
* The *Advanced* JSON panel only accepts known style fields and an allow-list of Matplotlib
  rcParams. Anything that could reach the file system or external programs, such as
  `text.usetex`, is refused.
* Result HTML files (interactive PCA) are always served as downloads, never rendered inside
  the studio's origin.
* Interactive plots use the Plotly bundle shipped in the `plotly` Python package and served at
  `/vendor/plotly.min.js`, so no third-party CDN is needed under the CSP.

## 9. Scientific caveats to share with clients

* Differential "volcano" plots use a normal-approximation test and are **exploratory**. Use
  DESeq2/edgeR/DiffBind for inferential claims. Peak-centre p-values are unadjusted.
* Signal is plotted as delivered. Normalisation (CPM, RPGC, spike-in) happens upstream and
  should be stated in the methods. Every result includes `methods.txt` and `provenance.json`
  (input hashes, software versions, parameters).
