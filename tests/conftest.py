"""Shared fixtures for taxafasta tests."""

from __future__ import annotations

from pathlib import Path

import pytest

DATA_DIR = Path(__file__).parent / "data"


@pytest.fixture
def data_dir() -> Path:
    """Path to the tests/data directory."""
    return DATA_DIR


@pytest.fixture
def tiny_nodes_path(data_dir: Path) -> Path:
    return data_dir / "tiny_nodes.dmp"


@pytest.fixture
def tiny_merged_path(data_dir: Path) -> Path:
    return data_dir / "tiny_merged.dmp"


@pytest.fixture
def tiny_names_path(data_dir: Path) -> Path:
    return data_dir / "tiny_names.dmp"


@pytest.fixture
def sample_fasta_path(data_dir: Path) -> Path:
    return data_dir / "sample.fasta"


@pytest.fixture
def sample_fasta_gz_path(data_dir: Path) -> Path:
    return data_dir / "sample.fasta.gz"


@pytest.fixture
def tiny_taxdump_dir(data_dir: Path) -> Path:
    """Subdirectory with tiny taxonomy files named nodes.dmp, merged.dmp, names.dmp."""
    return data_dir / "tiny_taxdump"
