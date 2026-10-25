"""SQLite persistence for users, projects, files, jobs, saved figures, templates and audit."""

from __future__ import annotations

import hashlib
import json
import secrets
import sqlite3
import time
import uuid
from contextlib import contextmanager
from pathlib import Path

SCHEMA = """
PRAGMA journal_mode=WAL;
CREATE TABLE IF NOT EXISTS users (
  id TEXT PRIMARY KEY, email TEXT UNIQUE NOT NULL, name TEXT, password_hash TEXT,
  is_admin INTEGER DEFAULT 0, active INTEGER DEFAULT 1, created REAL, last_login REAL
);
CREATE TABLE IF NOT EXISTS sessions (
  token_hash TEXT PRIMARY KEY, user_id TEXT NOT NULL, created REAL, expires REAL
);
CREATE TABLE IF NOT EXISTS invites (
  token_hash TEXT PRIMARY KEY, user_id TEXT NOT NULL, expires REAL
);
CREATE TABLE IF NOT EXISTS projects (
  id TEXT PRIMARY KEY, name TEXT, description TEXT, data_dir TEXT, genome TEXT,
  status TEXT DEFAULT 'active', expires TEXT, quota_bytes INTEGER, defaults TEXT,
  created REAL
);
CREATE TABLE IF NOT EXISTS members (
  project_id TEXT, user_id TEXT, role TEXT DEFAULT 'editor', PRIMARY KEY (project_id, user_id)
);
CREATE TABLE IF NOT EXISTS files (
  id TEXT PRIMARY KEY, project_id TEXT, path TEXT NOT NULL, name TEXT, role TEXT, label TEXT,
  grp TEXT, color TEXT, cmap TEXT, track TEXT, source TEXT, size INTEGER, mtime REAL,
  sort INTEGER DEFAULT 0, info TEXT, uploaded_by TEXT, created REAL
);
CREATE INDEX IF NOT EXISTS files_project ON files(project_id);
CREATE TABLE IF NOT EXISTS jobs (
  id TEXT PRIMARY KEY, project_id TEXT, user_id TEXT, status TEXT, payload TEXT, result TEXT,
  error TEXT, created REAL, started REAL, finished REAL, elapsed REAL, bytes INTEGER DEFAULT 0
);
CREATE INDEX IF NOT EXISTS jobs_project ON jobs(project_id, created);
CREATE TABLE IF NOT EXISTS figures (
  id TEXT PRIMARY KEY, project_id TEXT, name TEXT, notes TEXT, config TEXT, job_id TEXT,
  created_by TEXT, created REAL, updated REAL
);
CREATE TABLE IF NOT EXISTS templates (
  id TEXT PRIMARY KEY, project_id TEXT, name TEXT, description TEXT, config TEXT,
  created_by TEXT, created REAL
);
CREATE TABLE IF NOT EXISTS audit (
  id INTEGER PRIMARY KEY AUTOINCREMENT, ts REAL, user_id TEXT, project_id TEXT, action TEXT,
  detail TEXT, ip TEXT
);
CREATE INDEX IF NOT EXISTS audit_ts ON audit(ts);
"""

JSON_COLUMNS = {"info", "payload", "result", "config", "defaults"}


def token_hash(token: str) -> str:
    return hashlib.sha256(token.encode()).hexdigest()


def hash_password(password: str) -> str:
    if len(password) < 10:
        raise ValueError("Passwords need at least 10 characters")
    salt = secrets.token_hex(16)
    digest = hashlib.pbkdf2_hmac("sha256", password.encode(), bytes.fromhex(salt), 310_000).hex()
    return f"pbkdf2_sha256$310000${salt}${digest}"


def check_password(password: str, stored: str | None) -> bool:
    if not stored:
        return False
    try:
        _, rounds, salt, digest = stored.split("$")
    except ValueError:
        return False
    candidate = hashlib.pbkdf2_hmac(
        "sha256", password.encode(), bytes.fromhex(salt), int(rounds)
    ).hex()
    return secrets.compare_digest(candidate, digest)


def new_id() -> str:
    return uuid.uuid4().hex


class Database:
    def __init__(self, path: Path):
        self.path = Path(path)
        self.path.parent.mkdir(parents=True, exist_ok=True)
        with self.connect() as db:
            db.executescript(SCHEMA)

    @contextmanager
    def connect(self):
        db = sqlite3.connect(self.path, timeout=30)
        db.row_factory = sqlite3.Row
        db.execute("PRAGMA foreign_keys=ON")
        try:
            yield db
            db.commit()
        finally:
            db.close()

    @staticmethod
    def _row(row):
        if row is None:
            return None
        out = dict(row)
        for key in JSON_COLUMNS & out.keys():
            if out[key]:
                out[key] = json.loads(out[key])
        return out

    def one(self, sql, *args):
        with self.connect() as db:
            return self._row(db.execute(sql, args).fetchone())

    def all(self, sql, *args):
        with self.connect() as db:
            return [self._row(r) for r in db.execute(sql, args).fetchall()]

    def run(self, sql, *args):
        with self.connect() as db:
            return db.execute(sql, args).rowcount

    def insert(self, table, **values):
        for key in JSON_COLUMNS & values.keys():
            if values[key] is not None and not isinstance(values[key], str):
                values[key] = json.dumps(values[key])
        cols = ", ".join(values)
        marks = ", ".join("?" for _ in values)
        with self.connect() as db:
            db.execute(f"INSERT INTO {table} ({cols}) VALUES ({marks})", tuple(values.values()))
        return values

    def update(self, table, key, key_value, **values):
        if not values:
            return 0
        for col in JSON_COLUMNS & values.keys():
            if values[col] is not None and not isinstance(values[col], str):
                values[col] = json.dumps(values[col])
        sets = ", ".join(f"{col} = ?" for col in values)
        return self.run(f"UPDATE {table} SET {sets} WHERE {key} = ?", *values.values(), key_value)

    # ------------------------------------------------------------ users & sessions
    def user_by_email(self, email):
        return self.one("SELECT * FROM users WHERE email = ? COLLATE NOCASE", email.strip())

    def create_user(self, email, name=None, password=None, is_admin=False):
        email = email.strip().lower()
        if "@" not in email and email != "local":
            raise ValueError("Use an e-mail address as the user name")
        if self.user_by_email(email):
            raise ValueError(f"User {email} already exists")
        user = dict(
            id=new_id(),
            email=email,
            name=name or email.split("@")[0],
            is_admin=int(is_admin),
            active=1,
            created=time.time(),
            password_hash=hash_password(password) if password else None,
        )
        self.insert("users", **user)
        return user

    def create_session(self, user_id, hours):
        token = secrets.token_urlsafe(32)
        now = time.time()
        self.insert(
            "sessions",
            token_hash=token_hash(token),
            user_id=user_id,
            created=now,
            expires=now + hours * 3600,
        )
        self.run("UPDATE users SET last_login = ? WHERE id = ?", now, user_id)
        self.run("DELETE FROM sessions WHERE expires < ?", now)
        return token

    def session_user(self, token):
        if not token:
            return None
        return self.one(
            "SELECT u.* FROM sessions s JOIN users u ON u.id = s.user_id "
            "WHERE s.token_hash = ? AND s.expires > ? AND u.active = 1",
            token_hash(token),
            time.time(),
        )

    def create_invite(self, user_id, days=7):
        token = secrets.token_urlsafe(24)
        self.run("DELETE FROM invites WHERE user_id = ?", user_id)
        self.insert(
            "invites",
            token_hash=token_hash(token),
            user_id=user_id,
            expires=time.time() + days * 86400,
        )
        return token

    def redeem_invite(self, token, password):
        invite = self.one(
            "SELECT * FROM invites WHERE token_hash = ? AND expires > ?",
            token_hash(token),
            time.time(),
        )
        if not invite:
            raise ValueError(
                "This invitation link is invalid or has expired; ask the core for a new one"
            )
        self.run(
            "UPDATE users SET password_hash = ?, active = 1 WHERE id = ?",
            hash_password(password),
            invite["user_id"],
        )
        self.run("DELETE FROM invites WHERE user_id = ?", invite["user_id"])
        self.run("DELETE FROM sessions WHERE user_id = ?", invite["user_id"])
        return self.one("SELECT * FROM users WHERE id = ?", invite["user_id"])

    # ------------------------------------------------------------ projects & membership
    def project_role(self, user, project_id):
        """'admin', 'editor', 'viewer' or None."""
        if user is None:
            return None
        if user["is_admin"]:
            return "admin"
        row = self.one(
            "SELECT role FROM members WHERE project_id = ? AND user_id = ?", project_id, user["id"]
        )
        return row["role"] if row else None

    def projects_for(self, user):
        if user["is_admin"]:
            rows = self.all("SELECT *, 'admin' AS role FROM projects ORDER BY created DESC")
        else:
            rows = self.all(
                "SELECT p.*, m.role FROM projects p JOIN members m ON m.project_id = p.id "
                "WHERE m.user_id = ? ORDER BY p.created DESC",
                user["id"],
            )
        return rows

    def audit(self, user, action, project_id=None, detail=None, ip=None):
        self.insert(
            "audit",
            ts=time.time(),
            user_id=user["id"] if user else None,
            project_id=project_id,
            action=action,
            detail=json.dumps(detail) if isinstance(detail, (dict, list)) else detail,
            ip=ip,
        )
