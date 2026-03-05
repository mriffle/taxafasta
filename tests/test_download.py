"""Unit tests for the download/caching module (§11.4)."""

from __future__ import annotations

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
