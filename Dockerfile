# GenomeCanvas Studio: multi-user web server for publication-ready genomic figures.
# Build:  docker build -t genomecanvas-studio .
# Run:    see deploy/docker-compose.yml and docs/STUDIO_DEPLOY.md
FROM python:3.12-slim

ENV PYTHONDONTWRITEBYTECODE=1 \
    PYTHONUNBUFFERED=1 \
    PIP_NO_CACHE_DIR=1 \
    MPLBACKEND=Agg

# pyBigWig compiles against libcurl/zlib; fonts give Matplotlib Arial-metric (Liberation) and DejaVu faces.
RUN apt-get update \
 && apt-get install -y --no-install-recommends \
      build-essential libcurl4-openssl-dev zlib1g-dev \
      fonts-dejavu-core fonts-liberation2 curl tini \
 && rm -rf /var/lib/apt/lists/*

WORKDIR /app
COPY pyproject.toml README.md LICENSE ./
COPY src ./src
COPY web ./web
# Small public demo datasets (the large reference_demo fixtures are excluded by .dockerignore).
COPY tests/fixtures/trackplot_demo ./tests/fixtures/trackplot_demo
COPY tests/fixtures/seqplots_demo ./tests/fixtures/seqplots_demo

RUN pip install -e ".[studio]" \
 && useradd --system --uid 10001 --home-dir /srv/genomecanvas --shell /usr/sbin/nologin studio \
 && mkdir -p /srv/genomecanvas /etc/genomecanvas \
 && chown studio:studio /srv/genomecanvas

USER studio
ENV GC_STUDIO_WORKSPACE=/srv/genomecanvas \
    GC_STUDIO_STATIC_DIR=/app/web \
    MPLCONFIGDIR=/srv/genomecanvas/mpl
VOLUME ["/srv/genomecanvas"]
EXPOSE 8765

HEALTHCHECK --interval=30s --timeout=5s --start-period=20s \
  CMD curl -fsS http://127.0.0.1:8765/api/health || exit 1

ENTRYPOINT ["tini", "--", "gcanvas-studio"]
CMD ["--config", "/etc/genomecanvas/studio.toml", "serve", "--host", "0.0.0.0"]
