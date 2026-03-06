"""Smoke tests that hit real UniProt endpoints.

These are marked ``@pytest.mark.smoke`` and require network access.
Run them explicitly with::

    pytest -m smoke

They are **not** included in the default test suite.
"""

from __future__ import annotations

import pytest

from taxafasta.download import (
    UNIPROT_SWISSPROT_URL,
    UNIPROT_TREMBL_URL,
    open_uniprot_stream,
)


@pytest.mark.smoke
def test_stream_swissprot_first_entries() -> None:
    """Download the first few entries from Swiss-Prot and verify FASTA."""
    stream = open_uniprot_stream(
        UNIPROT_SWISSPROT_URL,
        label="Swiss-Prot smoke test",
    )
    lines_read = 0
    headers_seen = 0
    try:
        for line in stream:
            lines_read += 1
            if line.startswith(">"):
                headers_seen += 1
                # Verify it looks like a UniProt header
                assert "OX=" in line or "OS=" in line
            if headers_seen >= 5:
                break
    finally:
        stream.close()

    assert headers_seen >= 5
    assert lines_read > headers_seen  # sequences too


@pytest.mark.smoke
def test_stream_trembl_first_entries() -> None:
    """Download the first few entries from TrEMBL and verify FASTA."""
    stream = open_uniprot_stream(
        UNIPROT_TREMBL_URL,
        label="TrEMBL smoke test",
    )
    lines_read = 0
    headers_seen = 0
    try:
        for line in stream:
            lines_read += 1
            if line.startswith(">"):
                headers_seen += 1
                assert "OX=" in line or "OS=" in line
            if headers_seen >= 5:
                break
    finally:
        stream.close()

    assert headers_seen >= 5
    assert lines_read > headers_seen
