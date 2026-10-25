"""`gcanvas-studio` command: run the server and administer users, projects and storage."""

from __future__ import annotations

import argparse
import getpass
import json
import sys

from .config import example_config, load_config
from .db import Database
from .projects import Projects


def _context(args):
    overrides = {
        k: getattr(args, k)
        for k in ("workspace", "data_root", "shared_root", "auth")
        if getattr(args, k, None)
    }
    config = load_config(args.config, overrides)
    config.workspace_path.mkdir(parents=True, exist_ok=True)
    db = Database(config.db_path)
    return config, db, Projects(db, config)


def _invite(config, db, user):
    if config.auth != "local":
        return f"(users sign in via {config.auth} authentication)"
    token = db.create_invite(user["id"])
    base = config.public_url.rstrip("/") or f"http://{config.host}:{config.port}"
    return f"{base}/invite?token={token}"


def main(argv=None):
    parser = argparse.ArgumentParser(prog="gcanvas-studio", description=__doc__)
    parser.add_argument("--config", help="TOML configuration file")
    parser.add_argument("--workspace")
    parser.add_argument("--data-root", dest="data_root")
    parser.add_argument("--shared-root", dest="shared_root")
    parser.add_argument("--auth", choices=["none", "local", "proxy"])
    sub = parser.add_subparsers(dest="command", required=True)

    serve = sub.add_parser("serve", help="run the web server")
    serve.add_argument("--host")
    serve.add_argument("--port", type=int)
    serve.add_argument("--workers", type=int)

    sub.add_parser("example-config", help="print a commented configuration file")

    user = sub.add_parser("user", help="manage accounts")
    user_sub = user.add_subparsers(dest="action", required=True)
    add = user_sub.add_parser("add")
    add.add_argument("email")
    add.add_argument("--name")
    add.add_argument("--admin", action="store_true")
    add.add_argument(
        "--password", action="store_true", help="set a password now instead of an invite link"
    )
    user_sub.add_parser("list")
    invite = user_sub.add_parser("invite", help="new invitation / password-reset link")
    invite.add_argument("email")
    for action in ("disable", "enable"):
        user_sub.add_parser(action).add_argument("email")

    project = sub.add_parser("project", help="manage projects")
    project_sub = project.add_subparsers(dest="action", required=True)
    padd = project_sub.add_parser("add")
    padd.add_argument("id")
    padd.add_argument("--name")
    padd.add_argument("--dir", default="", help="folder relative to data_root")
    padd.add_argument("--genome", default="")
    padd.add_argument("--description", default="")
    padd.add_argument("--expires", help="YYYY-MM-DD; access ends after this date")
    padd.add_argument("--quota-gb", type=float)
    project_sub.add_parser("list")
    scan = project_sub.add_parser("scan")
    scan.add_argument(
        "id", nargs="?", help="omit to rescan every active project and shared annotations"
    )
    member = project_sub.add_parser("member")
    member.add_argument("id")
    member.add_argument("email")
    member.add_argument("--role", choices=["editor", "viewer"], default="editor")
    member.add_argument("--remove", action="store_true")
    for action in ("archive", "activate"):
        project_sub.add_parser(action).add_argument("id")

    sub.add_parser("cleanup", help="apply retention rules now")
    args = parser.parse_args(argv)

    if args.command == "example-config":
        print(example_config())
        return 0
    if args.command == "serve":
        from .server import serve as run

        overrides = {
            k: getattr(args, k)
            for k in ("workspace", "data_root", "shared_root", "auth", "host", "port", "workers")
            if getattr(args, k, None)
        }
        run(load_config(args.config, overrides))
        return 0

    config, db, projects = _context(args)
    if args.command == "user":
        if args.action == "add":
            password = None
            if args.password:
                password = getpass.getpass("Password (≥10 characters): ")
            created = db.create_user(
                args.email, name=args.name, password=password, is_admin=args.admin
            )
            print(f"created {created['email']}{' (admin)' if args.admin else ''}")
            if not password:
                print("invitation link (valid 7 days):", _invite(config, db, created))
        elif args.action == "list":
            for u in db.all("SELECT * FROM users ORDER BY email"):
                flags = ("admin " if u["is_admin"] else "") + ("" if u["active"] else "disabled")
                print(f"{u['email']:40s} {u['name'] or '':24s} {flags}")
        elif args.action == "invite":
            target = db.user_by_email(args.email) or sys.exit(f"no user {args.email}")
            print(_invite(config, db, target))
        else:
            target = db.user_by_email(args.email) or sys.exit(f"no user {args.email}")
            db.update("users", "id", target["id"], active=int(args.action == "enable"))
            if args.action == "disable":
                db.run("DELETE FROM sessions WHERE user_id = ?", target["id"])
            print(f"{args.action}d {target['email']}")
    elif args.command == "project":
        if args.action == "add":
            created = projects.create(
                args.id,
                name=args.name,
                data_dir=args.dir,
                genome=args.genome,
                description=args.description,
                expires=args.expires,
                quota_gb=args.quota_gb,
            )
            print(f"created project {created['id']}")
            if args.dir:
                print(json.dumps(projects.scan(args.id), indent=2, default=str))
        elif args.action == "list":
            for p in db.all("SELECT * FROM projects ORDER BY created DESC"):
                n = db.one("SELECT COUNT(*) AS n FROM members WHERE project_id = ?", p["id"])["n"]
                print(
                    f"{p['id']:24s} {p['status']:9s} {p['genome'] or '-':8s} members={n:<3d} dir={p['data_dir'] or '-'}  expires={p['expires'] or '-'}"
                )
        elif args.action == "scan":
            ids = (
                [args.id]
                if args.id
                else [p["id"] for p in db.all("SELECT id FROM projects WHERE status = 'active'")]
            )
            for pid in ids:
                print(pid, json.dumps(projects.scan(pid), default=str))
            if not args.id:
                print("shared", json.dumps(projects.scan_shared()))
        elif args.action == "member":
            target = db.user_by_email(args.email)
            if target is None:
                if args.remove:
                    sys.exit(f"no user {args.email}")
                target = db.create_user(args.email)
                print(f"created {target['email']}; invitation link:", _invite(config, db, target))
            projects.set_member(args.id, target["id"], None if args.remove else args.role)
            print(("removed " if args.remove else f"{args.role}: ") + target["email"])
        else:
            db.update(
                "projects",
                "id",
                args.id,
                status="archived" if args.action == "archive" else "active",
            )
            print(f"{args.id}: {'archived' if args.action == 'archive' else 'active'}")
    elif args.command == "cleanup":
        from .jobs import JobRunner

        print(json.dumps(JobRunner(db, config).cleanup()))
    return 0


if __name__ == "__main__":
    sys.exit(main())
