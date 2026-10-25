"""Single-user launcher kept for compatibility: python3 web/server.py [--port N] [--data-dir PATH].

For multi-user deployments use `gcanvas-studio serve --config studio.toml` (see docs/STUDIO_DEPLOY.md).
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from genomecanvas.studio.config import load_config  # noqa: E402
from genomecanvas.studio.db import Database  # noqa: E402
from genomecanvas.studio.projects import Projects  # noqa: E402
from genomecanvas.studio.server import serve  # noqa: E402


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--port", type=int, default=8765)
    parser.add_argument("--config", help="TOML configuration (enables multi-user options)")
    parser.add_argument("--workers", type=int)
    parser.add_argument(
        "--data-dir",
        type=Path,
        help="local folder to expose as project 'local' (files are read in place)",
    )
    args = parser.parse_args()
    overrides = {"port": args.port, "workers": args.workers}
    if args.data_dir:
        overrides["data_root"] = str(args.data_dir.expanduser().resolve())
    config = load_config(args.config, overrides)
    if args.data_dir:
        config.workspace_path.mkdir(parents=True, exist_ok=True)
        projects = Projects(Database(config.db_path), config)
        if not projects.db.one("SELECT id FROM projects WHERE id = 'local'"):
            projects.create("local", name=f"Local · {args.data_dir.name}", data_dir=".")
        else:
            projects.db.update(
                "projects", "id", "local", data_dir=".", name=f"Local · {args.data_dir.name}"
            )
    serve(config)


if __name__ == "__main__":
    main()
