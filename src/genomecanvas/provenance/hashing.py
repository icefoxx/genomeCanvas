"""Content hashing helpers."""

from hashlib import sha256
from pathlib import Path


def hash_file(path: str | Path) -> str:
    """Return the SHA256 digest of a local file."""
    digest = sha256()
    with Path(path).open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()
