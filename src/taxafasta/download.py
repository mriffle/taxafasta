"""Taxonomy and UniProt FASTA download and caching logic (§6)."""

from __future__ import annotations

import gzip
import io
import sys
import tarfile
import time
from datetime import datetime, timezone
from pathlib import Path
from typing import IO, Any
from urllib.error import URLError
from urllib.request import Request, urlopen

TAXDUMP_URL = "https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump.tar.gz"

UNIPROT_TREMBL_URL = (
    "https://ftp.uniprot.org/pub/databases/uniprot/"
    "current_release/knowledgebase/complete/"
    "uniprot_trembl.fasta.gz"
)
UNIPROT_SWISSPROT_URL = (
    "https://ftp.uniprot.org/pub/databases/uniprot/"
    "current_release/knowledgebase/complete/"
    "uniprot_sprot.fasta.gz"
)
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


_MAX_RETRIES = 5
_INITIAL_BACKOFF = 5.0
_BACKOFF_FACTOR = 2.0
_CONNECT_TIMEOUT = 120
_READ_CHUNK = 65536


def _open_raw_stream(
    url: str,
    byte_offset: int = 0,
) -> Any:
    """Open a raw byte stream from *url*, optionally resuming at *byte_offset*."""
    try:
        import requests

        headers = {}
        if byte_offset > 0:
            headers["Range"] = f"bytes={byte_offset}-"
        resp = requests.get(
            url,
            stream=True,
            timeout=_CONNECT_TIMEOUT,
            headers=headers,
        )
        resp.raise_for_status()
        raw = resp.raw
        raw.decode_content = False
        return raw
    except ImportError:
        pass

    req = Request(url)  # noqa: S310
    if byte_offset > 0:
        req.add_header("Range", f"bytes={byte_offset}-")
    return urlopen(req, timeout=_CONNECT_TIMEOUT)  # noqa: S310


class ResilientByteStream(io.RawIOBase):
    """A byte stream that transparently reconnects on transient errors.

    Tracks the number of compressed bytes consumed and, on failure,
    reopens the HTTP connection with a ``Range`` header so the gzip
    decompressor above it sees an uninterrupted stream.
    """

    def __init__(
        self,
        url: str,
        *,
        label: str = "",
        max_retries: int = _MAX_RETRIES,
        initial_backoff: float = _INITIAL_BACKOFF,
    ) -> None:
        super().__init__()
        self._url = url
        self._label = label
        self._max_retries = max_retries
        self._initial_backoff = initial_backoff
        self._bytes_read = 0
        self._inner: Any = _open_raw_stream(url)

    # -- io.RawIOBase interface --

    def readable(self) -> bool:
        return True

    def readinto(self, b: bytearray | memoryview) -> int:  # type: ignore[override]
        """Read up to len(b) bytes into *b*, retrying on network errors."""
        backoff = self._initial_backoff
        for attempt in range(self._max_retries + 1):
            try:
                if hasattr(self._inner, "readinto"):
                    n = self._inner.readinto(b)
                else:
                    n = self._inner_read_into(b)
                if n is None:
                    n = 0
                self._bytes_read += int(n)
                return int(n)
            except (OSError, URLError) as exc:
                if attempt >= self._max_retries:
                    raise
                print(
                    f"\nNetwork error on {self._label} at byte "
                    f"{self._bytes_read:,}: {exc}\n"
                    f"Retrying in {backoff:.0f}s "
                    f"(attempt {attempt + 2}/{self._max_retries + 1})...",
                    file=sys.stderr,
                )
                time.sleep(backoff)
                backoff *= _BACKOFF_FACTOR
                self._reconnect()
        return 0  # unreachable, satisfies mypy

    def _inner_read_into(self, b: bytearray | memoryview) -> int:
        """Fallback when the inner stream lacks readinto."""
        data = self._inner.read(len(b))
        if not data:
            return 0
        n = len(data)
        b[:n] = data
        return n

    def _reconnect(self) -> None:
        """Close the current connection and reopen with a Range header."""
        try:
            self._inner.close()
        except Exception:  # noqa: BLE001
            pass
        self._inner = _open_raw_stream(self._url, self._bytes_read)
        print(
            f"Resumed {self._label} from byte {self._bytes_read:,}",
            file=sys.stderr,
        )

    def close(self) -> None:
        try:
            self._inner.close()
        except Exception:  # noqa: BLE001
            pass
        super().close()


def open_uniprot_stream(
    url: str,
    *,
    label: str = "UniProt FASTA",
    max_retries: int = _MAX_RETRIES,
) -> IO[str]:
    """Open a streaming, gzip-decompressed text stream from a URL.

    The FASTA data is decompressed on the fly and never saved to disk.
    On transient network errors the connection is automatically retried
    with HTTP Range headers so the decompressor sees an uninterrupted
    byte stream.

    Parameters
    ----------
    url : HTTPS URL to a gzipped FASTA file.
    label : human-readable name for error messages.
    max_retries : maximum reconnection attempts before giving up.

    Returns
    -------
    A text-mode readable stream of decompressed FASTA data.
    """
    print(
        f"Streaming {label} from {url} ...",
        file=sys.stderr,
    )
    try:
        raw = ResilientByteStream(
            url,
            label=label,
            max_retries=max_retries,
        )
    except (URLError, OSError) as exc:
        print(
            f"Error downloading {label}: {exc}",
            file=sys.stderr,
        )
        raise SystemExit(1) from exc

    buffered = io.BufferedReader(raw, buffer_size=_READ_CHUNK)
    decompressed = gzip.open(buffered, mode="rb")
    return io.TextIOWrapper(
        decompressed,
        encoding="utf-8",
        errors="replace",
    )
