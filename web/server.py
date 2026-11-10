"""Single-user launcher kept for compatibility: python3 web/server.py [--port N] [--data-dir PATH].

Without --data-dir (or $GENOMECANVAS_DATA_DIR) it opens the GSE189561 mESC demo, downloading it
(about 3.7 GB, see web/demo_data.py) on first start; pass --no-download-data to skip the download.

For multi-user deployments use `gcanvas-studio serve --config studio.toml` (see docs/STUDIO_DEPLOY.md).
"""

from __future__ import annotations

import argparse
import os
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

import demo_data  # noqa: E402

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
        default=os.environ.get("GENOMECANVAS_DATA_DIR") or None,
        help="local folder to expose as project 'local' (files are read in place; "
        "default: $GENOMECANVAS_DATA_DIR, else the GSE189561 demo)",
    )
    parser.add_argument(
        "--download-data",
        action=argparse.BooleanOptionalAction,
        default=True,
        help="fetch missing GSE189561 demo files into the demo folder before starting "
        "(never into a folder of your own data)",
    )
    args = parser.parse_args()
    if args.data_dir is None and not args.config:
        args.data_dir = demo_data.DEFAULT_DIR
    if args.data_dir:
        folder = args.data_dir.expanduser().resolve()
        demo = folder == demo_data.DEFAULT_DIR.resolve() or demo_data.is_demo(folder)
        if demo and args.download_data:
            if not demo_data.ensure(folder):
                print("Demo data incomplete; starting with the files present.", file=sys.stderr)
        elif not folder.is_dir():
            if demo:
                print(f"No demo data at {folder}; run without --no-download-data to fetch it.")
                args.data_dir = None
            else:
                parser.error(f"--data-dir {folder} does not exist")
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
