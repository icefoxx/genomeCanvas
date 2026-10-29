"""Download the large public reference_demo inputs (about 104 MB) that are not stored in git.

    python tests/fixtures/reference_demo/download.py

Files already present with the expected SHA-256 are skipped. Sources: SOURCES.md in this directory.
"""

from __future__ import annotations

import hashlib
import shutil
import sys
import urllib.request
from pathlib import Path

HERE = Path(__file__).resolve().parent
GEO = "https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM5746nnn/{gsm}/suppl/{name}"
ENSEMBL = "{scheme}://ftp.ensembl.org/pub/release-102/gtf/mus_musculus/{name}"

# name -> (sha256, candidate URLs tried in order)
FILES: dict[str, tuple[str, list[str]]] = {
    "GSM5746910_MS_2cell_Input_rep1.bigWig": (
        "b1aef1189a574a67ec3dcd57f76987d2896e26ff0f746bce8a92711396eb6eed",
        [GEO.format(gsm="GSM5746910", name="GSM5746910_MS_2cell_Input_rep1.bigWig")],
    ),
    "GSM5746911_MS_2cell_Input_rep2.bigWig": (
        "c2afd64dece622d569f32540c9bb3858fc94fe7a184fa9b2036adbfa8c17b774",
        [GEO.format(gsm="GSM5746911", name="GSM5746911_MS_2cell_Input_rep2.bigWig")],
    ),
    "GSM5746912_MS_2cell_IP_rep1.bigWig": (
        "e7e9b0b25c69d3c55219e861fbc5f040f58fe6cc1125c8334d7a23b769f550f6",
        [GEO.format(gsm="GSM5746912", name="GSM5746912_MS_2cell_IP_rep1.bigWig")],
    ),
    "GSM5746913_MS_2cell_IP_rep2.bigWig": (
        "d023286e8af272a15144a3c312c1a9a29ef8578c2c426d5606f5aa828df26289",
        [GEO.format(gsm="GSM5746913", name="GSM5746913_MS_2cell_IP_rep2.bigWig")],
    ),
    "Mus_musculus.GRCm38.102.gtf.gz": (
        "8321415404aaf788c7da79774488ff227ac006d09a57ce6c616573a510338f64",
        [
            ENSEMBL.format(scheme=scheme, name="Mus_musculus.GRCm38.102.gtf.gz")
            for scheme in ("https", "http")  # the checksum guards the plain-HTTP fallback
        ],
    ),
}


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1 << 20), b""):
            digest.update(chunk)
    return digest.hexdigest()


def fetch(name: str, expected: str, urls: list[str]) -> bool:
    target = HERE / name
    if target.is_file() and sha256(target) == expected:
        print(f"ok        {name}")
        return True
    partial = target.with_name(target.name + ".part")
    for url in urls:
        print(f"download  {name} <- {url}", flush=True)
        try:
            with urllib.request.urlopen(url, timeout=60) as response, partial.open("wb") as out:
                shutil.copyfileobj(response, out, 1 << 20)
        except OSError as exc:
            print(f"  failed: {exc}", file=sys.stderr)
            continue
        if sha256(partial) == expected:
            partial.replace(target)
            return True
        print("  failed: SHA-256 mismatch", file=sys.stderr)
    partial.unlink(missing_ok=True)
    return False


def main() -> int:
    failed = [name for name, (expected, urls) in FILES.items() if not fetch(name, expected, urls)]
    if failed:
        print(f"could not download: {', '.join(failed)}", file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())
