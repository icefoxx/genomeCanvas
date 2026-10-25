"""Plot job execution: persistent worker processes with timeouts, plus housekeeping."""

from __future__ import annotations

import multiprocessing as mp
import shutil
import threading
import time
from collections import deque
from multiprocessing.connection import wait
from pathlib import Path

from .db import Database


def friendly(exc: BaseException) -> str:
    if isinstance(exc, MemoryError):
        return "The figure needed more memory than allowed; use fewer regions or a larger bin size"
    if isinstance(exc, (ValueError, FileNotFoundError, KeyError)):
        return str(exc).strip("'\"")
    return f"{type(exc).__name__}: {exc}"


def _worker_main(conn, memory_mb: int, cache_dir: str) -> None:
    """Long-lived worker: imports the engine once, then renders jobs sent over the pipe."""
    try:
        import resource

        limit = memory_mb * 1024 * 1024
        resource.setrlimit(resource.RLIMIT_AS, (limit, limit))
    except (ImportError, ValueError, OSError):
        pass  # not enforceable on this platform (e.g. macOS); the timeout still applies
    import os

    os.environ.setdefault("MPLCONFIGDIR", str(Path(cache_dir).parent / "mpl"))
    from .engine import run_plot

    while True:
        message = conn.recv()
        if message is None:
            return
        job_id, payload, records, destination = message
        began = time.time()
        try:
            warnings, counts = run_plot(payload, records, destination, cache_dir)
            conn.send(
                dict(
                    status="complete", warnings=warnings, counts=counts, elapsed=time.time() - began
                )
            )
        except BaseException as exc:  # noqa: BLE001 - everything is reported to the client
            conn.send(dict(status="error", error=friendly(exc), elapsed=time.time() - began))


class _Slot:
    def __init__(self, context, memory_mb, cache_dir):
        self.context, self.memory_mb, self.cache_dir = context, memory_mb, cache_dir
        self.job = None
        self.started = 0.0
        self.start()

    def start(self):
        self.conn, child = self.context.Pipe()
        self.process = self.context.Process(
            target=_worker_main, args=(child, self.memory_mb, self.cache_dir), daemon=True
        )
        self.process.start()

    def restart(self):
        try:
            self.process.kill()
            self.process.join(5)
        except Exception:  # noqa: BLE001
            pass
        self.start()


class JobRunner:
    def __init__(self, db: Database, config):
        self.db = db
        self.config = config
        self.root = config.workspace_path / "jobs"
        self.cache_dir = config.workspace_path / "cache"
        self.root.mkdir(parents=True, exist_ok=True)
        self.cache_dir.mkdir(parents=True, exist_ok=True)
        self.queue: deque = deque()
        self.lock = threading.Condition()
        self.slots: list[_Slot] = []
        self.stopping = False
        # Jobs that were in flight when the server stopped will never finish.
        self.db.run(
            "UPDATE jobs SET status = 'error', error = 'Interrupted by a server restart; plot again' "
            "WHERE status IN ('queued', 'running')"
        )

    def start(self):
        context = mp.get_context("spawn")
        self.slots = [
            _Slot(context, self.config.worker_memory_mb, str(self.cache_dir))
            for _ in range(max(1, self.config.workers))
        ]
        threading.Thread(target=self._loop, name="gc-jobs", daemon=True).start()
        threading.Thread(target=self._housekeeping, name="gc-cleanup", daemon=True).start()

    def stop(self):
        self.stopping = True
        for slot in self.slots:
            try:
                slot.conn.send(None)
            except Exception:  # noqa: BLE001
                pass

    def job_dir(self, project_id, job_id) -> Path:
        return self.root / project_id / job_id

    def submit(self, job_id, project_id, payload, records):
        destination = self.job_dir(project_id, job_id)
        destination.mkdir(parents=True)
        with self.lock:
            self.queue.append((job_id, payload, records, str(destination)))
            self.lock.notify()

    def position(self, job_id):
        with self.lock:
            for index, item in enumerate(self.queue):
                if item[0] == job_id:
                    return index + 1
        return 0

    # ------------------------------------------------------------ scheduling loop
    def _finish(self, job_id, result):
        row = self.db.one("SELECT project_id FROM jobs WHERE id = ?", job_id)
        directory = self.job_dir(row["project_id"], job_id) if row else None
        size = (
            sum(p.stat().st_size for p in directory.iterdir() if p.is_file())
            if directory and directory.is_dir()
            else 0
        )
        self.db.update(
            "jobs",
            "id",
            job_id,
            status=result["status"],
            error=result.get("error"),
            result={k: v for k, v in result.items() if k not in ("status", "error")},
            finished=time.time(),
            elapsed=round(result.get("elapsed", 0), 2),
            bytes=size,
        )

    def _loop(self):
        while not self.stopping:
            with self.lock:
                for slot in self.slots:
                    if slot.job is None and self.queue:
                        item = self.queue.popleft()
                        slot.job, slot.started = item[0], time.time()
                        self.db.update(
                            "jobs", "id", item[0], status="running", started=slot.started
                        )
                        try:
                            slot.conn.send(item)
                        except (BrokenPipeError, OSError):
                            slot.restart()
                            slot.conn.send(item)
            busy = [s for s in self.slots if s.job]
            if not busy:
                with self.lock:
                    self.lock.wait(timeout=1.0)
                continue
            ready = wait([s.conn for s in busy] + [s.process.sentinel for s in busy], timeout=0.25)
            now = time.time()
            for slot in busy:
                if slot.conn in ready:
                    try:
                        result = slot.conn.recv()
                    except (EOFError, OSError):
                        result = dict(
                            status="error",
                            error="The plotting worker stopped unexpectedly (out of memory?)",
                        )
                        slot.restart()
                    self._finish(slot.job, result)
                    slot.job = None
                elif not slot.process.is_alive():
                    self._finish(
                        slot.job,
                        dict(
                            status="error",
                            error="The plotting worker crashed (likely out of memory)",
                        ),
                    )
                    slot.job = None
                    slot.restart()
                elif now - slot.started > self.config.job_timeout:
                    self._finish(
                        slot.job,
                        dict(
                            status="error",
                            elapsed=now - slot.started,
                            error=f"Stopped after {self.config.job_timeout}s; use fewer regions, a larger bin size or a smaller window",
                        ),
                    )
                    slot.job = None
                    slot.restart()

    # ------------------------------------------------------------ housekeeping
    def _housekeeping(self):
        while not self.stopping:
            try:
                self.cleanup()
            except Exception as exc:  # noqa: BLE001 - never let cleanup kill the server
                print(f"cleanup failed: {exc}", flush=True)
            time.sleep(3600)

    def cleanup(self, now=None):
        """Apply retention: old jobs (unless saved as figures), old uploads, cache size."""
        now = now or time.time()
        report = dict(jobs=0, uploads=0, cache=0)
        keep = {
            r["job_id"] for r in self.db.all("SELECT job_id FROM figures WHERE job_id IS NOT NULL")
        }
        cutoff = now - self.config.job_retention_days * 86400
        for job in self.db.all(
            "SELECT id, project_id FROM jobs WHERE created < ? AND status IN ('complete','error')",
            cutoff,
        ):
            if job["id"] in keep:
                continue
            shutil.rmtree(self.job_dir(job["project_id"], job["id"]), ignore_errors=True)
            self.db.run("DELETE FROM jobs WHERE id = ?", job["id"])
            report["jobs"] += 1
        cutoff = now - self.config.upload_retention_days * 86400
        uploads = self.config.workspace_path / "uploads"
        for row in self.db.all(
            "SELECT * FROM files WHERE source = 'upload' AND created < ?", cutoff
        ):
            folder = Path(row["path"]).parent
            if folder.parent.parent == uploads:
                shutil.rmtree(folder, ignore_errors=True)
            self.db.run("DELETE FROM files WHERE id = ?", row["id"])
            report["uploads"] += 1
        entries = sorted(
            (p for p in self.cache_dir.iterdir() if p.is_file()), key=lambda p: p.stat().st_mtime
        )
        total = sum(p.stat().st_size for p in entries)
        limit = self.config.cache_max_gb * 1024**3
        while entries and total > limit:
            victim = entries.pop(0)
            total -= victim.stat().st_size
            victim.unlink(missing_ok=True)
            report["cache"] += 1
        self.db.run("DELETE FROM sessions WHERE expires < ?", now)
        self.db.run("DELETE FROM invites WHERE expires < ?", now)
        self.db.run("DELETE FROM audit WHERE ts < ?", now - 365 * 86400)
        return report
