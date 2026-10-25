"""GenomeCanvas Studio: engine, projects/manifests, multi-user isolation, jobs and retention."""

from __future__ import annotations

import http.cookiejar
import json
import shutil
import struct
import threading
import time
import urllib.error
import urllib.request
from http.server import ThreadingHTTPServer
from pathlib import Path

import pytest

pytest.importorskip("pyBigWig")

from genomecanvas.studio import engine  # noqa: E402
from genomecanvas.studio.config import load_config  # noqa: E402
from genomecanvas.studio.db import Database  # noqa: E402
from genomecanvas.studio.projects import DEMOS, Projects  # noqa: E402
from genomecanvas.studio.server import Studio, make_handler  # noqa: E402

FIX = Path(__file__).resolve().parent / "fixtures" / "trackplot_demo"
H1 = DEMOS["demo"]["defaults"]


# ---------------------------------------------------------------- fixtures


@pytest.fixture()
def data_root(tmp_path):
    """A core-facility delivery folder with a manifest and project.yaml."""
    root = tmp_path / "deliveries"
    project = root / "smith" / "P1"
    (project / "bw").mkdir(parents=True)
    for name in ("H1_k4me3", "H1_k27ac", "H1_Ctcf", "H1_Oct4"):
        shutil.copy(FIX / "bw" / f"{name}.bw", project / "bw")
    shutil.copy(FIX / "narrowpeak" / "H1_chromHMM.bed", project / "states.bed")
    shutil.copy(FIX / "narrowpeak" / "H1_k4me1.bed", project / "peaks.bed")
    shutil.copy(FIX / "narrowpeak" / "H1_Ctcf.bed", project / "ctcf.bed")
    shutil.copy(FIX / "chr6_refGene.gtf", project / "genes.gtf")
    (project / "samples.tsv").write_text(
        "file\tlabel\tgroup\tcolor\n"
        "bw/H1_k4me3.bw\tWT-1\tWT\t#123456\n"
        "bw/H1_k27ac.bw\tWT-2\tWT\t\n"
        "bw/H1_Ctcf.bw\tKO-1\tKO\t\n"
        "bw/H1_Oct4.bw\tKO-2\tKO\t\n"
        "states.bed\tChromatin states\tchromHMM\t\n"
        "peaks.bed\tConsensus peaks\tpeaks\t\n"
        "ctcf.bed\tCTCF sites\tpeaks\t\n"
        "genes.gtf\trefGene\tannotation\t\n"
        "missing.bw\tGhost\tWT\t\n"
    )
    (project / "project.yaml").write_text(
        "name: Smith lab ChIP-seq\ngenome: hg19\n"
        "defaults: {chrom: chr6, start: 31125000, end: 31145000, flank: 1000, binsize: 50}\n"
        "templates:\n  - name: TSS heatmap\n    config: {category: 1, plot: 0, settings: {sort: descend}}\n"
    )
    other = root / "jones" / "P2"
    other.mkdir(parents=True)
    shutil.copy(FIX / "bw" / "H1_H2az.bw", other)
    shutil.copy(FIX / "narrowpeak" / "H1_k4me1.bed", other / "peaks.bed")
    return root


def make_config(tmp_path, data_root, **extra):
    return load_config(
        None,
        dict(
            workspace=str(tmp_path / "ws"),
            data_root=str(data_root),
            auth="local",
            workers=1,
            demo_project=False,
            port=0,
            **extra,
        ),
    )


@pytest.fixture()
def studio(tmp_path, data_root):
    config = make_config(tmp_path, data_root)
    app = Studio(config)
    app.start()
    db, projects = app.db, app.projects
    db.create_user("core@uni.edu", password="core-password", is_admin=True)
    projects.create("P1", data_dir="smith/P1")
    projects.create("P2", data_dir="jones/P2", quota_gb=1e-8)  # ~10 bytes
    projects.scan("P1")
    projects.scan("P2")
    users = {}
    for email, project, role in [
        ("alice@lab.org", "P1", "editor"),
        ("bob@lab.org", "P1", "viewer"),
        ("carol@other.org", "P2", "editor"),
    ]:
        users[email] = db.create_user(email, password=f"{email}-pw-123")
        projects.set_member(project, users[email]["id"], role)
    server = ThreadingHTTPServer(("127.0.0.1", 0), make_handler(app))
    threading.Thread(target=server.serve_forever, daemon=True).start()
    yield app, f"http://127.0.0.1:{server.server_address[1]}"
    server.shutdown()
    app.jobs.stop()


class Client:
    def __init__(self, base):
        self.base = base
        self.opener = urllib.request.build_opener(
            urllib.request.HTTPCookieProcessor(http.cookiejar.CookieJar())
        )

    def call(self, method, path, body=None, headers=None, raw=None):
        data = raw if raw is not None else (json.dumps(body).encode() if body is not None else None)
        request = urllib.request.Request(
            self.base + path,
            data=data,
            method=method,
            headers={
                "X-GenomeCanvas": "studio",
                "Content-Type": "application/json",
                **(headers or {}),
            },
        )
        try:
            with self.opener.open(request) as response:
                payload = response.read()
                return response.status, (
                    json.loads(payload)
                    if response.headers.get_content_type() == "application/json"
                    else payload
                )
        except urllib.error.HTTPError as exc:
            payload = exc.read()
            try:
                return exc.code, json.loads(payload)
            except json.JSONDecodeError:
                return exc.code, payload

    def login(self, email, password):
        status, body = self.call("POST", "/api/login", {"email": email, "password": password})
        assert status == 200, body
        return self

    def wait(self, job_id, timeout=120):
        deadline = time.time() + timeout
        while time.time() < deadline:
            status, job = self.call("GET", f"/api/jobs/{job_id}")
            assert status == 200, job
            if job["status"] in ("complete", "error"):
                return job
            time.sleep(0.3)
        raise AssertionError("job did not finish")


# ---------------------------------------------------------------- engine


def demo_records(tmp_path):
    config = load_config(None, dict(workspace=str(tmp_path / "demo-ws")))
    projects = Projects(Database(config.db_path), config)
    projects.ensure_demos()
    project = projects.db.one("SELECT * FROM projects WHERE id = 'demo'")
    return projects.resolve(project, [dict(id=f["id"]) for f in projects.files(project)])


@pytest.mark.parametrize(("category", "plot"), [(0, 0), (1, 2), (2, 0), (4, 2)])
def test_engine_renders_every_format(tmp_path, category, plot):
    records = demo_records(tmp_path)
    out = tmp_path / "out"
    out.mkdir()
    engine.run_plot(
        dict(category=category, plot=plot, settings={**H1, "dpi": 72}, palette="npg"), records, out
    )
    for name in ("figure.png", "figure.svg", "figure.pdf", "methods.txt", "provenance.json"):
        assert (out / name).stat().st_size > 0
    provenance = json.loads((out / "provenance.json").read_text())
    assert "matplotlib" in provenance["software"]


def test_exact_size_export(tmp_path):
    records = demo_records(tmp_path)
    settings = {
        **H1,
        "dpi": 100,
        "exact_size": True,
        "size_unit": "mm",
        "fig_width": 180,
        "fig_height": 90,
    }
    engine.run_plot(dict(category=0, plot=0, settings=settings), records, tmp_path)
    width, height = struct.unpack(">II", (tmp_path / "figure.png").read_bytes()[16:24])
    assert abs(width - 180 / 25.4 * 100) <= 1 and abs(height - 90 / 25.4 * 100) <= 1


def test_matrix_cache_reused_and_recoloured(tmp_path):
    records = demo_records(tmp_path)
    cache = tmp_path / "cache"
    for index, color in enumerate(("#111111", "#ee0000")):
        out = tmp_path / f"run{index}"
        out.mkdir()
        recs = [dict(r, color=color) if r["role"] == "signal" else r for r in records]
        engine.run_plot(dict(category=1, plot=1, settings={**H1, "dpi": 72}), recs, out, cache)
    assert len(list(cache.glob("matrix-*.npz"))) == 1
    assert "#ee0000" in (tmp_path / "run1" / "provenance.json").read_text()


def test_advanced_overrides_and_rc_allowlist(tmp_path):
    records = demo_records(tmp_path)
    engine.run_plot(
        dict(
            category=1,
            plot=1,
            settings=H1,
            advanced={"profile": {"figsize": [4, 3]}, "rc": {"axes.linewidth": 2}},
        ),
        records,
        tmp_path,
    )
    with pytest.raises(ValueError, match="Unknown heatmap option"):
        engine.run_plot(
            dict(category=1, plot=0, settings=H1, advanced={"heatmap": {"nope": 1}}),
            records,
            tmp_path,
        )
    with pytest.raises(ValueError, match="cannot be changed"):
        engine.run_plot(
            dict(category=1, plot=1, settings=H1, advanced={"rc": {"text.usetex": True}}),
            records,
            tmp_path,
        )


def test_helpers():
    assert engine.parse_highlights("chr6:31,136,000-31,138,000 #ff0000 0.3", "#000") == [
        dict(start=31136000, end=31138000, color="#ff0000", alpha=0.3)
    ]
    assert engine.palette_colors({"custom_palette": "#111111, 222222, #333"}) == [
        "#111111",
        "#222222",
        "#333333",
    ]
    assert (
        engine.inspect_file(FIX / "narrowpeak" / "H1_chromHMM.bed", "regions")["kind"] == "chromHMM"
    )
    assert engine.gene_search([FIX / "chr6_refGene.gtf"], "POU5F1")[0]["chrom"] == "chr6"


# ---------------------------------------------------------------- projects


def test_manifest_drives_labels_groups_colours(studio):
    app, _ = studio
    project = app.db.one("SELECT * FROM projects WHERE id = 'P1'")
    files = {f["label"]: f for f in app.projects.files(project)}
    assert files["WT-1"]["grp"] == "WT" and files["WT-1"]["color"] == "#123456"
    assert "Ghost" not in files  # missing manifest rows are reported, not registered
    assert project["name"] == "Smith lab ChIP-seq" and project["defaults"]["chrom"] == "chr6"
    assert app.db.one("SELECT name FROM templates WHERE project_id = 'P1'")["name"] == "TSS heatmap"
    assert any("missing.bw" in w for w in app.projects.scan("P1")["warnings"])


def test_project_paths_cannot_escape_data_root(studio):
    app, _ = studio
    with pytest.raises(ValueError, match="under data_root"):
        app.projects.create("P3", data_dir="../../")


# ---------------------------------------------------------------- HTTP / multi-user


def test_authentication_required(studio):
    _, base = studio
    anonymous = Client(base)
    assert anonymous.call("GET", "/api/me")[0] == 401
    status, _ = anonymous.call(
        "POST", "/api/login", {"email": "alice@lab.org", "password": "wrong"}
    )
    assert status == 401
    request = urllib.request.Request(base + "/", method="GET")
    opener = urllib.request.build_opener(
        type(
            "NoRedirect",
            (urllib.request.HTTPRedirectHandler,),
            {"redirect_request": lambda *a, **k: None},
        )
    )
    with pytest.raises(urllib.error.HTTPError) as exc:
        opener.open(request)
    assert exc.value.code == 303 and exc.value.headers["Location"] == "/login"


def test_csrf_header_and_origin(studio):
    _, base = studio
    alice = Client(base).login("alice@lab.org", "alice@lab.org-pw-123")
    assert (
        alice.call("POST", "/api/projects/P1/jobs", {}, headers={"X-GenomeCanvas": "nope"})[0]
        == 403
    )
    assert (
        alice.call("POST", "/api/projects/P1/jobs", {}, headers={"Origin": "https://evil.example"})[
            0
        ]
        == 403
    )


def test_isolation_between_projects_and_roles(studio):
    _, base = studio
    alice = Client(base).login("alice@lab.org", "alice@lab.org-pw-123")
    bob = Client(base).login("bob@lab.org", "bob@lab.org-pw-123")
    carol = Client(base).login("carol@other.org", "carol@other.org-pw-123")

    assert [p["id"] for p in alice.call("GET", "/api/me")[1]["projects"]] == ["P1"]
    assert carol.call("GET", "/api/projects/P1")[0] == 404
    status, project = alice.call("GET", "/api/projects/P1")
    files = project["files"]
    signals = [f for f in files if f["role"] == "signal"]
    regions = [f for f in files if f["role"] == "regions"]
    assert {f["label"] for f in signals} == {"WT-1", "WT-2", "KO-1", "KO-2"}

    body = dict(
        category=1, plot=1, settings=H1, files=[{"id": f["id"]} for f in signals + regions[:1]]
    )
    status, job = alice.call("POST", "/api/projects/P1/jobs", body)
    assert status == 202, job
    result = alice.wait(job["id"])
    assert result["status"] == "complete", result

    # Carol cannot use P1 file ids from her own project, nor see Alice's job or results.
    assert carol.call("POST", "/api/projects/P2/jobs", body)[0] == 400
    assert carol.call("GET", f"/api/jobs/{job['id']}")[0] == 404
    assert carol.call("GET", f"/results/{job['id']}/figure.png")[0] == 404

    # Bob (viewer) can read but not plot or save.
    assert bob.call("POST", "/api/projects/P1/jobs", body)[0] == 403
    status, saved = alice.call(
        "POST",
        "/api/projects/P1/figures",
        {"name": "Fig 1", "config": {"category": 1}, "job_id": job["id"]},
    )
    assert status == 201
    figures = bob.call("GET", "/api/projects/P1/figures")[1]["figures"]
    assert figures[0]["name"] == "Fig 1" and figures[0]["thumb"]
    assert bob.call("GET", figures[0]["thumb"])[0] == 200
    assert bob.call("DELETE", f"/api/figures/{saved['figure']['id']}")[0] == 403
    assert carol.call("GET", "/api/projects/P1/figures")[0] == 404


def test_uploads_quota_and_ownership(studio):
    _, base = studio
    alice = Client(base).login("alice@lab.org", "alice@lab.org-pw-123")
    carol = Client(base).login("carol@other.org", "carol@other.org-pw-123")
    content = (FIX / "narrowpeak" / "H1_Ctcf.bed").read_bytes()
    status, uploaded = alice.call(
        "POST",
        "/api/projects/P1/uploads",
        raw=content,
        headers={"X-Filename": "extra.bed", "Content-Type": "application/octet-stream"},
    )
    assert status == 201 and uploaded["info"]["rows"] == 5
    assert (
        alice.call(
            "POST", "/api/projects/P1/uploads", raw=b"x", headers={"X-Filename": "../../evil.sh"}
        )[0]
        == 400
    )
    # P2 has a tiny quota.
    assert (
        carol.call(
            "POST", "/api/projects/P2/uploads", raw=content, headers={"X-Filename": "a.bed"}
        )[0]
        == 413
    )
    assert carol.call("DELETE", f"/api/projects/P2/files/{uploaded['id']}")[0] == 404
    assert alice.call("DELETE", f"/api/projects/P1/files/{uploaded['id']}")[0] == 200


def test_admin_invites_and_archiving(studio):
    app, base = studio
    core = Client(base).login("core@uni.edu", "core-password")
    status, created = core.call(
        "POST", "/api/admin/projects/P1/members", {"email": "dan@lab.org", "role": "editor"}
    )
    assert status == 200 and "/invite?token=" in created["invite_url"]
    token = created["invite_url"].split("token=")[1]
    dan = Client(base)
    assert dan.call("POST", "/api/invite", {"token": token, "password": "dan-secret-1"})[0] == 200
    assert [p["id"] for p in dan.call("GET", "/api/me")[1]["projects"]] == ["P1"]
    assert dan.call("POST", "/api/invite", {"token": token, "password": "again-secret"})[0] == 400

    alice = Client(base).login("alice@lab.org", "alice@lab.org-pw-123")
    assert alice.call("GET", "/api/admin/overview")[0] == 403
    assert core.call("PATCH", "/api/admin/projects/P1", {"status": "archived"})[0] == 200
    assert alice.call("GET", "/api/projects/P1")[0] == 404
    assert (
        core.call("PATCH", "/api/admin/projects/P1", {"status": "active", "expires": "2000-01-01"})[
            0
        ]
        == 200
    )
    assert alice.call("GET", "/api/projects/P1")[0] == 404
    overview = core.call("GET", "/api/admin/overview")[1]
    assert {p["id"] for p in overview["projects"]} == {"P1", "P2"}
    actions = {a["action"] for a in core.call("GET", "/api/admin/audit")[1]["audit"]}
    assert {"login", "member.set", "invite.accepted", "project.updated"} <= actions


# ---------------------------------------------------------------- jobs & retention


def test_job_timeout_kills_worker(tmp_path, data_root):
    config = make_config(tmp_path, data_root, job_timeout=1)
    app = Studio(config)
    app.start()
    try:
        app.projects.create("P1", data_dir="smith/P1")
        app.projects.scan("P1")
        project = app.db.one("SELECT * FROM projects WHERE id = 'P1'")
        records = app.projects.resolve(
            project, [dict(id=f["id"]) for f in app.projects.files(project) if f["role"] != "genes"]
        )
        payload = dict(category=2, plot=2, settings=H1)  # PairGrid takes several seconds
        app.db.insert(
            "jobs",
            id="a" * 32,
            project_id="P1",
            user_id="u",
            status="queued",
            payload=payload,
            created=time.time(),
        )
        app.jobs.submit("a" * 32, "P1", payload, records)
        deadline = time.time() + 60
        while time.time() < deadline:
            job = app.db.one("SELECT * FROM jobs WHERE id = ?", "a" * 32)
            if job["status"] in ("complete", "error"):
                break
            time.sleep(0.2)
        assert job["status"] == "error" and "Stopped after 1s" in job["error"]
    finally:
        app.jobs.stop()


def test_cleanup_respects_saved_figures(tmp_path, data_root):
    config = make_config(tmp_path, data_root, job_retention_days=1)
    app = Studio(config)
    old = time.time() - 3 * 86400
    for job_id in ("1" * 32, "2" * 32):
        app.db.insert(
            "jobs", id=job_id, project_id="P1", user_id="u", status="complete", created=old
        )
        app.jobs.job_dir("P1", job_id).mkdir(parents=True)
    app.db.insert(
        "figures", id="f", project_id="P1", name="keep", job_id="1" * 32, created=old, updated=old
    )
    report = app.jobs.cleanup()
    assert report["jobs"] == 1
    assert (
        app.jobs.job_dir("P1", "1" * 32).is_dir() and not app.jobs.job_dir("P1", "2" * 32).exists()
    )


# ---------------------------------------------------------------- genome bins / sliding windows


def test_genome_bins_match_pybigwig_and_windows_match_bins():
    import numpy as np
    import pyBigWig

    from genomecanvas.analysis.summary import extract_summary
    from genomecanvas.pipeline.genome_bins import genome_summary, region_summary, select_chromosomes

    sources = {n: str(FIX / "bw" / f"H1_{n}.bw") for n in ("k4me3", "Ctcf")}
    locus = ("chr6", 31125000, 31145000)
    handle = pyBigWig.open(sources["k4me3"])
    for stat in ("mean", "max"):
        summary, _ = genome_summary(sources, size=1000, region=locus, stat=stat, skip_zeros=False)
        expected = np.nan_to_num(
            np.array(handle.stats(*locus, type=stat, nBins=20, exact=True), dtype=float)
        )
        assert np.allclose(summary.values[:, 0], expected)
    windows, _ = genome_summary(
        sources, size=2000, step=500, region=locus, stat="mean0", skip_zeros=False
    )
    fine, _ = genome_summary(sources, size=500, region=locus, stat="mean0", skip_zeros=False)
    assert np.allclose(
        windows.values[:, 1], np.convolve(fine.values[:, 1], np.ones(4) / 4, "valid")
    )
    assert (
        windows.starts[1] - windows.starts[0] == 500 and windows.ends[0] - windows.starts[0] == 2000
    )

    bed = FIX / "narrowpeak" / "H1_chromHMM.bed"
    old = extract_summary(sources, bed, engine="pybigwig")
    new, _ = region_summary(sources, [bed], stat="sum")
    assert np.allclose(old.values, new.values)

    assert select_chromosomes(
        {"chr1": 1, "chr10": 1, "chr2": 1, "chrM": 1, "chrUn_x": 1, "chr1_random": 1, "chrX": 1}
    ) == ["chr1", "chr2", "chr10", "chrX"]
    with pytest.raises(ValueError, match="multiple of the step"):
        genome_summary(sources, size=1000, step=300, region=locus)


def test_genome_bins_blacklist_and_zero_filter():
    from genomecanvas.pipeline.genome_bins import genome_summary

    sources = {n: str(FIX / "bw" / f"H1_{n}.bw") for n in ("k4me3", "k27ac")}
    full, report = genome_summary(sources, size=1000, chroms=["chr6"], skip_zeros=True)
    assert report.windows_total > report.windows_kept == len(full.regions) > 0
    blocked, report = genome_summary(
        sources, size=1000, chroms=["6"], blacklist=FIX / "narrowpeak" / "H1_Ctcf.bed"
    )
    assert report.blacklisted > 0 and len(blocked.regions) < len(full.regions)


@pytest.mark.parametrize("mode", ["bins", "windows"])
def test_relationship_modes_via_api(studio, mode):
    _, base = studio
    alice = Client(base).login("alice@lab.org", "alice@lab.org-pw-123")
    files = alice.call("GET", "/api/projects/P1")[1]["files"]
    signals = [{"id": f["id"]} for f in files if f["role"] == "signal"]
    blacklist = next(f["id"] for f in files if f["label"] == "CTCF sites")
    settings = {
        **H1,
        "rel_mode": mode,
        "rel_bin": 1000,
        "rel_step": 250,
        "rel_scope": "chroms",
        "rel_chroms": "chr6",
        "rel_blacklist": blacklist,
        "dpi": 72,
    }
    for plot in (0, 1):
        status, job = alice.call(
            "POST",
            "/api/projects/P1/jobs",
            dict(category=2, plot=plot, settings=settings, files=signals),
        )
        assert status == 202, job
        result = alice.wait(job["id"])
        assert result["status"] == "complete", result
    names = {f["name"] for f in result["files"]}
    assert {"summary.tsv", "methods.txt"} <= names
    methods = alice.call("GET", f"/results/{job['id']}/methods.txt")[1].decode()
    assert (
        "sliding 1 kb windows with a 250 bp step"
        if mode == "windows"
        else "non-overlapping 1 kb genomic bins"
    ) in methods
    assert "CTCF sites" in methods
    header = alice.call("GET", f"/results/{job['id']}/summary.tsv")[1].decode().splitlines()[0]
    assert header.startswith("chrom\tstart\tend\tWT-1")


def test_dense_and_sparse_read_paths_agree():
    import numpy as np

    from genomecanvas.pipeline import genome_bins as gb

    path = str(FIX / "bw" / "H1_k27ac.bw")
    for stat in gb.STATS:
        results = []
        for mode in ("intervals", "values"):
            reader = gb.SignalReader(path, read_mode=mode)
            total, cover, peak = gb.bin_chromosome(
                reader, "chr6", 31_120_000, 31_150_000, 500, stat
            )
            results.append(gb._finish(stat, total, cover, peak, 500))
            reader.close()
        assert np.allclose(*results), stat


def test_reference_point_matrix_orients_minus_strand(tmp_path):
    """Minus-strand rows are reversed and mirrored around the anchor base (deepTools/bwtool)."""
    import numpy as np
    import pyBigWig

    from genomecanvas.pipeline.matrix import compute_region_matrix
    from genomecanvas.pipeline.metadata import SampleMeta

    path = str(FIX / "bw" / "H1_k27ac.bw")
    tss = 31132113
    bed = tmp_path / "pair.bed"
    bed.write_text(f"chr6\t{tss}\t{tss + 1}\tplus\t0\t+\nchr6\t{tss}\t{tss + 1}\tminus\t0\t-\n")
    handle = pyBigWig.open(path)

    def stats(start, end):
        return np.nan_to_num(
            np.array(
                handle.stats("chr6", start, end, type="mean", nBins=40, exact=True), dtype=float
            )
        )

    sample = [SampleMeta(path, "k27ac", "c", "m", "Blues", "auto")]
    for anchor in ("TSS", "center"):
        m = compute_region_matrix(
            sample,
            [(str(bed), "x")],
            flank=2000,
            binsize=100,
            ref_point=anchor,
            skip_zeros=False,
            engine="pybigwig",
        )
        assert np.allclose(m.values[0], stats(tss - 2000, tss + 2000)), anchor
        assert np.allclose(m.values[1], stats(tss - 1999, tss + 2001)[::-1]), anchor


@pytest.mark.parametrize("threshold", [1e9, 0.0], ids=["interval-path", "per-base-path"])
def test_window_matrix_matches_pybigwig_and_handles_edges(monkeypatch, threshold):
    import numpy as np
    import pyBigWig

    from genomecanvas.pipeline import genome_bins
    from genomecanvas.pipeline.genome_bins import window_matrix

    monkeypatch.setattr(genome_bins, "DENSE_BYTES_PER_BASE", threshold)  # force each read path

    path = str(FIX / "bw" / "H1_k4me3.bw")
    handle = pyBigWig.open(path)
    length = handle.chroms()["chr6"]
    windows = [("chr6", 31125000 + 700 * i, 31125000 + 700 * i + 3000) for i in range(25)]
    windows += [("chr6", length - 1000, length + 2000), ("chr6", -500, 2500), ("chrNope", 0, 3000)]
    got = window_matrix(path, windows, 30, "mean")
    for row, (chrom, start, end) in enumerate(windows[:25]):
        expected = np.nan_to_num(
            np.array(handle.stats(chrom, start, end, type="mean", nBins=30, exact=True), dtype=float)
        )
        assert np.allclose(got[row], expected)
    edge = got[25]
    assert np.isnan(edge[10:]).all() and np.isfinite(edge[:10]).all()  # past the chromosome end
    assert np.isnan(got[26][:5]).all() and np.isnan(got[27]).all()


def test_pca_outputs_interactive_data_and_repelled_labels(tmp_path):
    records = [r for r in demo_records(tmp_path) if r["role"] == "signal"]
    settings = {
        **H1,
        "dpi": 72,
        "rel_mode": "bins",
        "rel_bin": 500,
        "rel_scope": "chroms",
        "rel_chroms": "chr6",
    }
    engine.run_plot(dict(category=2, plot=1, title="PCA", settings=settings), records, tmp_path)
    data = json.loads((tmp_path / "pca.json").read_text())
    v = data["variance"]  # PC1-3 of 7 samples: decreasing, below 100 % in total
    assert len(data["points"]) == 7 and v[0] >= v[1] >= v[2] and sum(v) <= 100 + 1e-6
    assert data["use_groups"] and {p["symbol"] for p in data["points"]} == {"circle", "square"}
    assert (tmp_path / "pca_scores.tsv").read_text().startswith("sample\tgroup\tfile\tPC1")
    # OCT4 and NANOG nearly coincide: their labels must have been pushed apart.
    pts = {p["name"]: p for p in data["points"]}
    oct4, nanog = pts["OCT4"], pts["NANOG"]
    assert abs(oct4["pcs"][0] - nanog["pcs"][0]) < 2
    assert oct4["label_offset"] != nanog["label_offset"]
    pytest.importorskip("plotly")
    assert "plotly" in (tmp_path / "pca_interactive.html").read_text()[:5000].lower()


def test_login_throttling_counts_failures_per_account(studio):
    app, base = studio
    attacker = Client(base)
    for _ in range(10):
        assert (
            attacker.call(
                "POST", "/api/login", {"email": "alice@lab.org", "password": "nope-nope"}
            )[0]
            == 401
        )
    assert (
        attacker.call(
            "POST", "/api/login", {"email": "alice@lab.org", "password": "alice@lab.org-pw-123"}
        )[0]
        == 429
    )
    # A different account from the same address is unaffected, and successes are not counted.
    for _ in range(12):
        Client(base).login("bob@lab.org", "bob@lab.org-pw-123")


def test_trusted_proxy_cidr_and_forwarded_ip(tmp_path, data_root):
    app = Studio(make_config(tmp_path, data_root, trusted_proxies=["10.0.0.0/8", "127.0.0.1"]))
    assert app.trusted("10.2.3.4") and app.trusted("127.0.0.1")
    assert not app.trusted("192.168.1.5") and not app.trusted("not-an-ip")
