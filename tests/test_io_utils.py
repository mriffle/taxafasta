"""Unit tests for the I/O utilities module (§11.3)."""

from __future__ import annotations

import gzip
from pathlib import Path

from taxafasta.io_utils import has_isal, is_gzip_file, open_input, open_output


# --- Gzip detection ---

def test_is_gzip_file_true(sample_fasta_gz_path: Path) -> None:
    assert is_gzip_file(sample_fasta_gz_path) is True


def test_is_gzip_file_false(sample_fasta_path: Path) -> None:
    assert is_gzip_file(sample_fasta_path) is False


def test_is_gzip_regardless_of_extension(tmp_path: Path) -> None:
    """A gzip file without .gz extension should still be detected."""
    gz_data = gzip.compress(b"hello world\n")
    no_ext = tmp_path / "data.txt"
    no_ext.write_bytes(gz_data)
    assert is_gzip_file(no_ext) is True


def test_is_not_gzip_with_gz_extension(tmp_path: Path) -> None:
    """A plain file with .gz extension should NOT be detected as gzip."""
    fake = tmp_path / "fake.fasta.gz"
    fake.write_bytes(b"not actually gzip\n")
    assert is_gzip_file(fake) is False


# --- Transparent opening ---

def test_open_input_plain(sample_fasta_path: Path) -> None:
    with open_input(sample_fasta_path) as fh:
        content = fh.read()
    assert ">sp|P12345|ARATH_ECOLI" in content


def test_open_input_gzip(sample_fasta_gz_path: Path) -> None:
    with open_input(sample_fasta_gz_path) as fh:
        content = fh.read()
    assert ">sp|P12345|ARATH_ECOLI" in content


def test_open_input_equivalence(sample_fasta_path: Path, sample_fasta_gz_path: Path) -> None:
    """Plain and gzip inputs should produce the same content."""
    with open_input(sample_fasta_path) as fh:
        plain = fh.read()
    with open_input(sample_fasta_gz_path) as fh:
        gzipped = fh.read()
    assert plain == gzipped


# --- Output ---

def test_open_output_gzip_default(tmp_path: Path) -> None:
    out_path = tmp_path / "out.fasta"
    stream, resolved = open_output(out_path, use_gzip=True)
    stream.write(">test\nACGT\n")
    stream.close()
    assert resolved.name == "out.fasta.gz"
    assert resolved.exists()
    # Verify it's valid gzip
    with gzip.open(resolved, "rt") as fh:
        assert ">test" in fh.read()


def test_open_output_gzip_already_has_gz(tmp_path: Path) -> None:
    out_path = tmp_path / "out.fasta.gz"
    stream, resolved = open_output(out_path, use_gzip=True)
    stream.write(">test\nACGT\n")
    stream.close()
    assert resolved.name == "out.fasta.gz"


def test_open_output_no_gzip(tmp_path: Path) -> None:
    out_path = tmp_path / "out.fasta"
    stream, resolved = open_output(out_path, use_gzip=False)
    stream.write(">test\nACGT\n")
    stream.close()
    assert resolved.name == "out.fasta"
    assert resolved.exists()
    content = resolved.read_text()
    assert ">test" in content


# --- isal fallback ---

def test_has_isal_returns_bool() -> None:
    """has_isal() should return a bool regardless of installation status."""
    assert isinstance(has_isal(), bool)
