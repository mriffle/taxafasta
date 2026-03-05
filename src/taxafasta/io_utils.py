"""Smart file opening with gzip detection, buffering, and isal fallback."""

from __future__ import annotations

import gzip
import io
import types
from pathlib import Path
from typing import IO

# Buffer sizes for high-throughput I/O (§4.3)
READ_BUFFER_SIZE = 16 * 1024 * 1024  # 16 MB
WRITE_BUFFER_SIZE = 16 * 1024 * 1024  # 16 MB

# Gzip magic bytes
_GZIP_MAGIC = b"\x1f\x8b"

# Try to import isal for faster gzip (§10)
_igzip: types.ModuleType | None
try:
    import isal.igzip as _igzip

    _HAS_ISAL = True
except ImportError:
    _igzip = None
    _HAS_ISAL = False


def has_isal() -> bool:
    """Return True if python-isal is available."""
    return _HAS_ISAL


def is_gzip_file(path: str | Path) -> bool:
    """Detect gzip format by magic bytes, not file extension (§3.1)."""
    with open(path, "rb") as fh:
        return fh.read(2) == _GZIP_MAGIC


def open_input(path: str | Path) -> IO[str]:
    """Open an input file transparently handling gzip.

    Uses isal.igzip when available for ~5-10x faster decompression.
    Falls back to stdlib gzip. Returns a text-mode stream.
    """
    path = Path(path)
    if is_gzip_file(path):
        if _HAS_ISAL and _igzip is not None:
            raw = _igzip.open(path, mode="rb")
        else:
            raw = gzip.open(path, mode="rb")
        return io.TextIOWrapper(raw, encoding="utf-8", errors="replace")
    return open(
        path,
        buffering=READ_BUFFER_SIZE,
        encoding="utf-8",
        errors="replace",
    )


def open_output(
    path: str | Path,
    *,
    use_gzip: bool = True,
) -> tuple[IO[str], Path]:
    """Open an output file, optionally gzip-compressed.

    Parameters
    ----------
    path : path to the output file
    use_gzip : if True (default), compress with gzip and append .gz if needed

    Returns
    -------
    (text_stream, resolved_path)
        resolved_path has .gz appended when gzip is enabled.
    """
    path = Path(path)
    if use_gzip:
        if not path.name.endswith(".gz"):
            path = path.with_name(path.name + ".gz")
        if _HAS_ISAL and _igzip is not None:
            raw = _igzip.open(path, mode="wb")
        else:
            raw = gzip.open(path, mode="wb")
        stream: IO[str] = io.TextIOWrapper(
            io.BufferedWriter(
                raw,
                buffer_size=WRITE_BUFFER_SIZE,
            ),
            encoding="utf-8",
            newline="",
        )
        return stream, path

    stream = open(
        path,
        "w",
        buffering=WRITE_BUFFER_SIZE,
        encoding="utf-8",
        newline="",
    )
    return stream, path
