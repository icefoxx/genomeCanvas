"""Download the GEO GSE189561 mouse ESC demo (about 3.7 GB) used as Studio's default local project.

    python web/demo_data.py [DEST] [--verify]

SOX2 and CTCF ChIP-seq with MACS2 peaks, H3K4me3 / H3K27ac / H3K4me1 and input from asynchronous
mESCs (mm10, PMID 37225980), one replicate each. Files already present at the expected size are
skipped (--verify also checks their SHA-256). A samples.tsv manifest and project.yaml are written
when missing, and the peak files are decompressed because Studio reads plain narrowPeak.
"""

from __future__ import annotations

import argparse
import gzip
import hashlib
import os
import shutil
import sys
import urllib.request
from pathlib import Path

DEFAULT_DIR = Path.home() / ".genomecanvas" / "data" / "GSE189561_mESC"
GEO = "https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM5703nnn/{gsm}/suppl/{name}"
MARKER = "# GSE189561"

# name -> (size, sha256, label, group)
FILES: dict[str, tuple[int, str, str, str]] = {
    "GSM5703952_asynchronousmESCs_SOX2_rep2.ChIP.bw": (
        568995356,
        "c0ec460c6dead8ec68c475d823ecbe81828b347e0c1a68f6851b79f0cc9ee224",
        "SOX2",
        "TF",
    ),
    "GSM5703956_asynchronousmESCs_CTCF_rep1.ChIP.bw": (
        1047180870,
        "de4acb97c5defe0b64533423fc8af57a01735d86dfe8ae60d3e7a34d4f9c6fe3",
        "CTCF",
        "TF",
    ),
    "GSM5703916_asynchronousmESCs_H3K4me3_rep3.ChIP.bw": (
        425958078,
        "c139bd40d9cf6f59903da582f8314a6c1430ab75319f0dd3daad4c741e10fdb8",
        "H3K4me3",
        "Histone",
    ),
    "GSM5703918_asynchronousmESCs_H3K27ac_rep3.ChIP.bw": (
        422083054,
        "dff19ff2e40d9e13814f0fa106b6045c66a1306a30caa08ed3d23ab6c6502a06",
        "H3K27ac",
        "Histone",
    ),
    "GSM5703908_asynchronousmESCs_H3K4me1_rep3.ChIP.bw": (
        588617646,
        "e97a103df0fbc8b6f6d17d49f2b6f0a51d59a0ac23dc5ca7ea276dcadb423ee8",
        "H3K4me1",
        "Histone",
    ),
    "GSM5703911_asynchronousmESCs_input_rep7.ChIP.bw": (
        615689764,
        "24134e0f5ed0ada2f7cfd911c686e798e8cad86771ad7202ffc865231cfcadf0",
        "Input",
        "Input",
    ),
    "GSM5703952_asynchronousmESCs_SOX2_rep2.ChIP.macs2.narrowPeak.gz": (
        1536184,
        "795d432143d2692cea17e05d33d6680feaad7e7b018fe10c36b764b0d24b7996",
        "SOX2 peaks",
        "SOX2 peaks",
    ),
    "GSM5703956_asynchronousmESCs_CTCF_rep1.ChIP.macs2.narrowPeak.gz": (
        1847520,
        "35166b2d80ce4807d3c9c8ad81edf77daccbc3c74c3281c6cf5a6232d8123fb3",
        "CTCF peaks",
        "CTCF peaks",
    ),
}

PROJECT_YAML = """\
# Read on every scan; see genomeCanvas docs/STUDIO_DEPLOY.md.
name: GSE189561 · mESC SOX2 / CTCF / histone marks
genome: mm10
description: GEO GSE189561 (PMID 37225980), asynchronous mESC, RPM bigWigs + MACS2 peaks
defaults:
  chrom: chr3
  start: 34630000
  end: 34670000
  assembly: mm10
  flank: 1000
  binsize: 50
  alignment: center
  xlab: distance from peak center
  colorbar_label: RPM
  numerator: TF
  denominator: Histone
  title: mESC · Sox2 locus
"""


def is_demo(folder: Path) -> bool:
    """True when `folder` holds (part of) this demo, judged by its manifest."""
    manifest = folder / "samples.tsv"
    return manifest.is_file() and manifest.read_text(encoding="utf-8").startswith(MARKER)


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1 << 20), b""):
            digest.update(chunk)
    return digest.hexdigest()


def fetch(folder: Path, name: str, size: int, expected: str, verify: bool) -> bool:
    target = folder / name
    if (
        target.is_file()
        and target.stat().st_size == size
        and (not verify or sha256(target) == expected)
    ):
        print(f"ok        {name}")
        return True
    partial = target.with_name(target.name + ".part")
    url = GEO.format(gsm=name.split("_")[0], name=name)
    print(f"download  {name} ({size / 1e6:.0f} MB) <- {url}", flush=True)
    try:
        with urllib.request.urlopen(url, timeout=60) as response, partial.open("wb") as out:
            shutil.copyfileobj(response, out, 1 << 20)
    except OSError as exc:
        print(f"  failed: {exc}", file=sys.stderr)
        partial.unlink(missing_ok=True)
        return False
    if sha256(partial) != expected:
        print("  failed: SHA-256 mismatch", file=sys.stderr)
        partial.unlink(missing_ok=True)
        return False
    partial.replace(target)
    return True


def ensure(folder: Path, *, verify: bool = False) -> bool:
    """Download missing demo files into `folder` and write its manifest; True when complete."""
    folder.mkdir(parents=True, exist_ok=True)
    failed = [
        name
        for name, (size, expected, _label, _group) in FILES.items()
        if not fetch(folder, name, size, expected, verify)
    ]
    for name in FILES:
        gz = folder / name
        plain = gz.with_suffix("")
        if name.endswith(".gz") and gz.is_file() and not plain.is_file():
            with gzip.open(gz, "rb") as src, plain.open("wb") as out:
                shutil.copyfileobj(src, out)
    manifest = folder / "samples.tsv"
    if not manifest.exists():
        rows = [
            f"{MARKER} · asynchronous mESC (mm10), one replicate per mark",
            "file\tlabel\tgroup",
        ]
        rows += [
            f"{name.removesuffix('.gz')}\t{label}\t{group}"
            for name, (_size, _sha, label, group) in FILES.items()
        ]
        manifest.write_text("\n".join(rows) + "\n", encoding="utf-8")
    if not (folder / "project.yaml").exists():
        (folder / "project.yaml").write_text(PROJECT_YAML, encoding="utf-8")
    if failed:
        print(f"could not download: {', '.join(failed)}", file=sys.stderr)
    return not failed


def main() -> int:
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    parser.add_argument(
        "dest",
        nargs="?",
        type=Path,
        default=os.environ.get("GENOMECANVAS_DATA_DIR") or DEFAULT_DIR,
        help=f"target folder (default: $GENOMECANVAS_DATA_DIR or {DEFAULT_DIR})",
    )
    parser.add_argument("--verify", action="store_true", help="also check SHA-256 of present files")
    args = parser.parse_args()
    return 0 if ensure(args.dest.expanduser().resolve(), verify=args.verify) else 1


if __name__ == "__main__":
    sys.exit(main())
