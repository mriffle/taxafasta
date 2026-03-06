"""Unit tests for the download/caching module (§11.4)."""

from __future__ import annotations

import gzip
import io
import tarfile
from pathlib import Path
from unittest.mock import patch

import pytest

from taxafasta.download import (
    _is_cache_valid,
    download_taxdump,
    ensure_taxdump,
    read_timestamp,
)


def _make_taxdump_tar_gz(nodes_content: bytes, merged_content: bytes) -> bytes:
    """Create an in-memory taxdump.tar.gz with given file contents."""
    buf = io.BytesIO()
    with tarfile.open(fileobj=buf, mode="w:gz") as tar:
        for name, data in [
            ("nodes.dmp", nodes_content),
            ("merged.dmp", merged_content),
        ]:
            info = tarfile.TarInfo(name=name)
            info.size = len(data)
            tar.addfile(info, io.BytesIO(data))
    return buf.getvalue()


# --- Cache validation ---


def test_is_cache_valid_true(tiny_taxdump_dir: Path) -> None:
    assert _is_cache_valid(tiny_taxdump_dir) is True


def test_is_cache_valid_false(tmp_path: Path) -> None:
    assert _is_cache_valid(tmp_path) is False


# --- Cache hit ---


def test_download_taxdump_cache_hit(tiny_taxdump_dir: Path) -> None:
    """When cache is valid, no download should occur."""
    result = download_taxdump(tiny_taxdump_dir)
    assert result == tiny_taxdump_dir


# --- Successful download (mocked) ---


def test_download_taxdump_success(tmp_path: Path) -> None:
    nodes = b"1\t|\t1\t|\tno rank\t|\n2\t|\t1\t|\tsuperkingdom\t|\n"
    merged = b"50\t|\t7\t|\n"
    archive = _make_taxdump_tar_gz(nodes, merged)

    with patch("taxafasta.download._fetch_url", return_value=archive):
        result = download_taxdump(tmp_path / "cache")

    assert (result / "nodes.dmp").exists()
    assert (result / "merged.dmp").exists()
    assert read_timestamp(result) is not None


# --- Network failure ---


def test_download_taxdump_network_error(tmp_path: Path) -> None:
    side_effect = OSError("Network unreachable")
    with patch("taxafasta.download._fetch_url", side_effect=side_effect):
        with pytest.raises(SystemExit):
            download_taxdump(tmp_path / "cache")


# --- Corrupt archive ---


def test_download_taxdump_corrupt(tmp_path: Path) -> None:
    with patch("taxafasta.download._fetch_url", return_value=b"this is not a tar.gz"):
        with pytest.raises(SystemExit):
            download_taxdump(tmp_path / "cache")


# --- ensure_taxdump ---


def test_ensure_taxdump_with_dir(tiny_taxdump_dir: Path) -> None:
    result = ensure_taxdump(tiny_taxdump_dir, None)
    assert result == tiny_taxdump_dir


def test_ensure_taxdump_missing_dir(tmp_path: Path) -> None:
    with pytest.raises(SystemExit):
        ensure_taxdump(tmp_path / "nonexistent", None)


def test_ensure_taxdump_missing_files(tmp_path: Path) -> None:
    empty_dir = tmp_path / "empty"
    empty_dir.mkdir()
    with pytest.raises(SystemExit):
        ensure_taxdump(empty_dir, None)


# --- Timestamp ---


def test_read_timestamp_none(tmp_path: Path) -> None:
    assert read_timestamp(tmp_path) is None


# --- File permissions ---


def test_cache_dir_created(tmp_path: Path) -> None:
    """download_taxdump should create the cache directory if it doesn't exist."""
    nodes = b"1\t|\t1\t|\tno rank\t|\n"
    merged = b""
    archive = _make_taxdump_tar_gz(nodes, merged)

    cache = tmp_path / "sub" / "deep" / "cache"
    with patch("taxafasta.download._fetch_url", return_value=archive):
        download_taxdump(cache)

    assert cache.exists()
    assert cache.is_dir()


# --- ResilientByteStream ---

# Shared compressed test data: a gzipped FASTA with many lines so that
# mid-stream failures happen after partial delivery of real content.
_TEST_FASTA_LINES = [
    f">sp|P{i:05d}|PROT_{i} OS=Org OX={i} PE=1 SV=1\n{'ACGT' * 20}\n"
    for i in range(500)
]
_TEST_FASTA_TEXT = "".join(_TEST_FASTA_LINES)
_TEST_FASTA_GZ: bytes = gzip.compress(_TEST_FASTA_TEXT.encode())


class _OffsetByteStream:
    """Simulate an HTTP response that serves bytes from a given offset.

    Optionally raises *fail_error* after *fail_after_bytes* bytes have
    been read from this particular stream instance, simulating a
    mid-transfer network failure.
    """

    def __init__(
        self,
        data: bytes,
        offset: int = 0,
        fail_after_bytes: int | None = None,
        fail_error: type[Exception] = BrokenPipeError,
    ) -> None:
        self._data = data[offset:]
        self._pos = 0
        self._fail_after = fail_after_bytes
        self._fail_error = fail_error
        self.closed = False

    def read(self, n: int = -1) -> bytes:
        if self._fail_after is not None and self._pos >= self._fail_after:
            raise self._fail_error("simulated network error")
        if n < 0:
            chunk = self._data[self._pos :]
        else:
            chunk = self._data[self._pos : self._pos + n]
        self._pos += len(chunk)
        return chunk

    def close(self) -> None:
        self.closed = True


def test_resilient_stream_normal_read() -> None:
    """No errors → data is delivered byte-for-byte."""
    from taxafasta.download import ResilientByteStream

    data = b"hello world"
    with patch(
        "taxafasta.download._open_raw_stream",
        return_value=_OffsetByteStream(data),
    ):
        stream = ResilientByteStream(
            "http://example.com/test.gz",
            label="test",
        )
    buf = bytearray(1024)
    n = stream.readinto(buf)
    assert buf[:n] == data
    stream.close()


def test_resilient_stream_byte_offset_on_reconnect() -> None:
    """Reconnect must request the exact byte offset already consumed."""
    from taxafasta.download import ResilientByteStream

    data = b"0123456789ABCDEF" * 10  # 160 bytes
    recorded_offsets: list[int] = []

    def _mock_open(url: str, byte_offset: int = 0) -> _OffsetByteStream:
        recorded_offsets.append(byte_offset)
        if byte_offset == 0:
            # Fail after delivering some bytes
            return _OffsetByteStream(data, offset=0, fail_after_bytes=50)
        return _OffsetByteStream(data, offset=byte_offset)

    with (
        patch("taxafasta.download._open_raw_stream", side_effect=_mock_open),
        patch("taxafasta.download.time") as mock_time,
    ):
        mock_time.sleep = lambda _: None
        stream = ResilientByteStream(
            "http://x", label="test", max_retries=3, initial_backoff=0.0
        )
        result = b""
        while True:
            buf = bytearray(4)  # small reads to exercise the loop
            n = stream.readinto(buf)
            if n == 0:
                break
            result += bytes(buf[:n])

    assert result == data
    assert recorded_offsets[0] == 0
    # The reconnect offset must be > 0 and match the bytes delivered
    assert len(recorded_offsets) == 2
    assert recorded_offsets[1] > 0
    # The two segments must concatenate to the full data
    first_part = data[: recorded_offsets[1]]
    second_part = data[recorded_offsets[1] :]
    assert result == first_part + second_part
    stream.close()


def test_resilient_stream_mid_stream_failure() -> None:
    """Failure after partial data: no bytes lost, no bytes duplicated."""
    from taxafasta.download import ResilientByteStream

    data = b"A" * 100 + b"B" * 100 + b"C" * 100  # 300 bytes

    # Fail after 150 bytes (mid-stream)
    def _mock_open(url: str, byte_offset: int = 0) -> _OffsetByteStream:
        if byte_offset == 0:
            return _OffsetByteStream(data, offset=0, fail_after_bytes=150)
        return _OffsetByteStream(data, offset=byte_offset)

    with (
        patch("taxafasta.download._open_raw_stream", side_effect=_mock_open),
        patch("taxafasta.download.time") as mock_time,
    ):
        mock_time.sleep = lambda _: None
        stream = ResilientByteStream(
            "http://x", label="test", max_retries=3, initial_backoff=0.0
        )
        result = b""
        while True:
            buf = bytearray(32)
            n = stream.readinto(buf)
            if n == 0:
                break
            result += bytes(buf[:n])

    assert result == data, (
        f"Data mismatch: got {len(result)} bytes, expected {len(data)}"
    )
    stream.close()


def test_resilient_stream_multiple_consecutive_failures() -> None:
    """Recover after 3 consecutive failures at different offsets."""
    from taxafasta.download import ResilientByteStream

    data = b"X" * 500
    call_idx = {"n": 0}

    def _mock_open(url: str, byte_offset: int = 0) -> _OffsetByteStream:
        call_idx["n"] += 1
        # Calls 1-3 each fail after delivering ~50 bytes
        if call_idx["n"] <= 3:
            return _OffsetByteStream(data, offset=byte_offset, fail_after_bytes=50)
        # Call 4 succeeds
        return _OffsetByteStream(data, offset=byte_offset)

    with (
        patch("taxafasta.download._open_raw_stream", side_effect=_mock_open),
        patch("taxafasta.download.time") as mock_time,
    ):
        mock_time.sleep = lambda _: None
        stream = ResilientByteStream(
            "http://x", label="test", max_retries=5, initial_backoff=0.0
        )
        result = b""
        while True:
            buf = bytearray(64)
            n = stream.readinto(buf)
            if n == 0:
                break
            result += bytes(buf[:n])

    assert result == data
    assert call_idx["n"] == 4  # 3 failures + 1 success
    stream.close()


def test_resilient_stream_connection_reset_error() -> None:
    """ConnectionResetError is caught and retried, not just BrokenPipeError."""
    from taxafasta.download import ResilientByteStream

    data = b"PAYLOAD"
    call_count = {"n": 0}

    def _mock_open(url: str, byte_offset: int = 0) -> _OffsetByteStream:
        call_count["n"] += 1
        if call_count["n"] == 1:
            return _OffsetByteStream(
                data,
                offset=0,
                fail_after_bytes=0,
                fail_error=ConnectionResetError,
            )
        return _OffsetByteStream(data, offset=byte_offset)

    with (
        patch("taxafasta.download._open_raw_stream", side_effect=_mock_open),
        patch("taxafasta.download.time") as mock_time,
    ):
        mock_time.sleep = lambda _: None
        stream = ResilientByteStream(
            "http://x", label="test", max_retries=3, initial_backoff=0.0
        )
        buf = bytearray(1024)
        n = stream.readinto(buf)

    assert buf[:n] == data
    assert call_count["n"] == 2  # 1 failure + 1 success
    stream.close()


def test_resilient_stream_exhausted_retries() -> None:
    """After max_retries failures the exception propagates."""
    from taxafasta.download import ResilientByteStream

    class _AlwaysFail:
        def read(self, n: int = -1) -> bytes:
            raise BrokenPipeError("Broken pipe")

        def close(self) -> None:
            pass

    with (
        patch(
            "taxafasta.download._open_raw_stream",
            return_value=_AlwaysFail(),
        ),
        patch("taxafasta.download.time") as mock_time,
    ):
        mock_time.sleep = lambda _: None
        stream = ResilientByteStream(
            "http://x", label="test", max_retries=2, initial_backoff=0.0
        )
        buf = bytearray(1024)
        with pytest.raises(BrokenPipeError):
            stream.readinto(buf)
    stream.close()


def test_resilient_stream_bytes_read_tracking() -> None:
    """_bytes_read reflects the total compressed bytes consumed."""
    from taxafasta.download import ResilientByteStream

    data = b"0123456789" * 10  # 100 bytes
    with patch(
        "taxafasta.download._open_raw_stream",
        return_value=_OffsetByteStream(data),
    ):
        stream = ResilientByteStream("http://x", label="test", max_retries=1)
    total = 0
    while True:
        buf = bytearray(7)  # odd size to test partial reads
        n = stream.readinto(buf)
        if n == 0:
            break
        total += n
    assert total == 100
    assert stream._bytes_read == 100  # noqa: SLF001
    stream.close()


def test_full_stack_gzip_integrity_through_reconnect() -> None:
    """End-to-end: gzip → BufferedReader → ResilientByteStream.

    This is the gold-standard data-fidelity test.  We compress a
    multi-line FASTA, simulate a mid-stream failure at a byte boundary
    inside the compressed data, reconnect, and verify every single
    decompressed character matches the original.
    """
    from taxafasta.download import ResilientByteStream

    gz_data = _TEST_FASTA_GZ
    # Fail roughly 40% through the compressed data
    fail_point = len(gz_data) * 2 // 5

    def _mock_open(url: str, byte_offset: int = 0) -> _OffsetByteStream:
        if byte_offset == 0:
            return _OffsetByteStream(gz_data, offset=0, fail_after_bytes=fail_point)
        return _OffsetByteStream(gz_data, offset=byte_offset)

    with (
        patch("taxafasta.download._open_raw_stream", side_effect=_mock_open),
        patch("taxafasta.download.time") as mock_time,
    ):
        mock_time.sleep = lambda _: None
        raw = ResilientByteStream(
            "http://x", label="test", max_retries=3, initial_backoff=0.0
        )
        buffered = io.BufferedReader(raw, buffer_size=256)
        decompressed = gzip.open(buffered, mode="rb")
        text_stream = io.TextIOWrapper(decompressed, encoding="utf-8", errors="replace")

        recovered_text = text_stream.read()

    # Verify complete text match (byte-for-byte fidelity)
    assert recovered_text == _TEST_FASTA_TEXT, (
        f"Text mismatch: got {len(recovered_text)} chars, "
        f"expected {len(_TEST_FASTA_TEXT)} chars"
    )
    # Also verify line-by-line to give better diagnostics on failure
    recovered_lines = recovered_text.splitlines(keepends=True)
    expected_lines = _TEST_FASTA_TEXT.splitlines(keepends=True)
    assert len(recovered_lines) == len(expected_lines)
    for i, (got, expected) in enumerate(zip(recovered_lines, expected_lines)):
        assert got == expected, (
            f"Line {i} differs:\n  got:      {got!r}\n  expected: {expected!r}"
        )


def test_full_stack_gzip_integrity_multi_failure() -> None:
    """End-to-end with 3 failures at different points in compressed data."""
    from taxafasta.download import ResilientByteStream

    gz_data = _TEST_FASTA_GZ
    fail_points = [
        len(gz_data) // 5,
        len(gz_data) * 2 // 5,
        len(gz_data) * 3 // 5,
    ]
    call_idx = {"n": 0}

    def _mock_open(url: str, byte_offset: int = 0) -> _OffsetByteStream:
        idx = call_idx["n"]
        call_idx["n"] += 1
        if idx < len(fail_points):
            # Each connection fails after delivering some bytes
            budget = fail_points[idx] - byte_offset
            if budget > 0:
                return _OffsetByteStream(
                    gz_data,
                    offset=byte_offset,
                    fail_after_bytes=budget,
                )
        return _OffsetByteStream(gz_data, offset=byte_offset)

    with (
        patch("taxafasta.download._open_raw_stream", side_effect=_mock_open),
        patch("taxafasta.download.time") as mock_time,
    ):
        mock_time.sleep = lambda _: None
        raw = ResilientByteStream(
            "http://x", label="test", max_retries=5, initial_backoff=0.0
        )
        buffered = io.BufferedReader(raw, buffer_size=256)
        decompressed = gzip.open(buffered, mode="rb")
        text_stream = io.TextIOWrapper(decompressed, encoding="utf-8", errors="replace")
        recovered_text = text_stream.read()

    assert recovered_text == _TEST_FASTA_TEXT, (
        f"Full text mismatch: got {len(recovered_text)} chars, "
        f"expected {len(_TEST_FASTA_TEXT)} chars"
    )
