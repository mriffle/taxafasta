"""Taxonomy download and caching logic (§6)."""

from __future__ import annotations

import io
import sys
import tarfile
from datetime import datetime, timezone
from pathlib import Path
from urllib.error import URLError
from urllib.request import urlopen

TAXDUMP_URL = "https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump.tar.gz"
_REQUIRED_FILES = ("nodes.dmp", "merged.dmp")
_OPTIONAL_FILES = ("names.dmp", "delnodes.dmp")
_TIMESTAMP_FILE = ".taxafasta_download_timestamp"


def _default_cache_dir() -> Path:
    return Path.home() / ".taxafasta"


def _is_cache_valid(cache_dir: Path) -> bool:
    """Check whether the cache directory has the required taxonomy files."""
    for fname in _REQUIRED_FILES:
        if not (cache_dir / fname).exists():
            return False
    return True


def _write_timestamp(cache_dir: Path) -> None:
    ts = datetime.now(timezone.utc).isoformat()
    (cache_dir / _TIMESTAMP_FILE).write_text(ts, encoding="utf-8")


def read_timestamp(cache_dir: Path) -> str | None:
    """Read the download timestamp from cache, or None."""
    ts_path = cache_dir / _TIMESTAMP_FILE
    if ts_path.exists():
        return ts_path.read_text(encoding="utf-8").strip()
    return None


def _fetch_url(url: str) -> bytes:
    """Fetch URL content, using requests if available, else urllib."""
    try:
        import requests  # type: ignore[import-untyped]

        resp = requests.get(url, timeout=120)
        resp.raise_for_status()
        return resp.content  # type: ignore[no-any-return]
    except ImportError:
        pass
    except Exception as exc:
        raise OSError(str(exc)) from exc
    with urlopen(url, timeout=120) as resp:  # noqa: S310
        return resp.read()  # type: ignore[no-any-return]


def download_taxdump(cache_dir: Path | None = None) -> Path:
    """Download and extract taxdump.tar.gz into the cache directory.

    Parameters
    ----------
    cache_dir : directory to store extracted files; defaults to ~/.taxafasta/

    Returns
    -------
    Path to the cache directory containing extracted .dmp files.
    """
    if cache_dir is None:
        cache_dir = _default_cache_dir()
    cache_dir = Path(cache_dir)

    if _is_cache_valid(cache_dir):
        return cache_dir

    cache_dir.mkdir(parents=True, exist_ok=True)
    print(f"Downloading taxonomy data from {TAXDUMP_URL} ...", file=sys.stderr)

    try:
        data = _fetch_url(TAXDUMP_URL)
    except (URLError, OSError) as exc:
        print(f"Error downloading taxonomy data: {exc}", file=sys.stderr)
        raise SystemExit(1) from exc

    try:
        with tarfile.open(fileobj=io.BytesIO(data), mode="r:gz") as tar:
            wanted = set(_REQUIRED_FILES + _OPTIONAL_FILES)
            for member in tar.getmembers():
                if member.name in wanted:
                    tar.extract(member, path=cache_dir, filter="data")
    except (tarfile.TarError, OSError) as exc:
        print(f"Error extracting taxonomy archive: {exc}", file=sys.stderr)
        raise SystemExit(1) from exc

    if not _is_cache_valid(cache_dir):
        print(
            "Error: Downloaded archive does not contain required files "
            f"({', '.join(_REQUIRED_FILES)}).",
            file=sys.stderr,
        )
        raise SystemExit(1)

    _write_timestamp(cache_dir)
    print(f"Taxonomy data cached in {cache_dir}", file=sys.stderr)
    return cache_dir


def ensure_taxdump(
    taxdump_dir: Path | None,
    cache_dir: Path | None,
) -> Path:
    """Resolve taxonomy dump directory.

    If taxdump_dir is provided and valid, use it.
    Otherwise download/use cache.
    """
    if taxdump_dir is not None:
        taxdump_dir = Path(taxdump_dir)
        if not taxdump_dir.is_dir():
            print(
                f"Error: Taxdump directory does not exist: {taxdump_dir}",
                file=sys.stderr,
            )
            raise SystemExit(1)
        for fname in _REQUIRED_FILES:
            if not (taxdump_dir / fname).exists():
                print(f"Error: {fname} not found in {taxdump_dir}", file=sys.stderr)
                raise SystemExit(1)
        return taxdump_dir

    return download_taxdump(cache_dir)
