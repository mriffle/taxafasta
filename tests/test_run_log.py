"""Unit tests for the run log module (§11.5)."""

from __future__ import annotations

from pathlib import Path

from taxafasta.fasta import FilterStats
from taxafasta.run_log import resolve_log_path, write_log

# --- Log file naming ---


def test_resolve_log_path_basic(tmp_path: Path) -> None:
    out = tmp_path / "bacteria.fasta"
    log = resolve_log_path(out)
    assert log == tmp_path / "bacteria.fasta.log"


def test_resolve_log_path_gz_stripped(tmp_path: Path) -> None:
    """Log for bacteria.fasta.gz should be bacteria.fasta.log."""
    out = tmp_path / "bacteria.fasta.gz"
    log = resolve_log_path(out)
    assert log == tmp_path / "bacteria.fasta.log"


def test_resolve_log_path_collision(tmp_path: Path) -> None:
    out = tmp_path / "bacteria.fasta"
    log0 = tmp_path / "bacteria.fasta.log"
    log0.write_text("existing log")
    log = resolve_log_path(out)
    assert log == tmp_path / "bacteria.fasta.log1"


def test_resolve_log_path_multiple_collisions(tmp_path: Path) -> None:
    out = tmp_path / "bacteria.fasta"
    (tmp_path / "bacteria.fasta.log").write_text("log")
    (tmp_path / "bacteria.fasta.log1").write_text("log")
    log = resolve_log_path(out)
    assert log == tmp_path / "bacteria.fasta.log2"


# --- write_log ---


def _make_stats(
    total: int = 100,
    included: int = 60,
    excluded: int = 38,
    no_ox: int = 1,
    unknown_count: int = 1,
) -> FilterStats:
    stats = FilterStats()
    stats.total = total
    stats.included = included
    stats.excluded = excluded
    stats.no_ox = no_ox
    stats.unknown_taxid_count = unknown_count
    if no_ox > 0:
        stats.unparseable_headers.append(">gnl|custom|seq1 hypothetical protein")
    if unknown_count > 0:
        stats.unknown_taxids[99999999] = unknown_count
    return stats


def test_write_log_contains_metadata(tmp_path: Path) -> None:
    log_path = tmp_path / "test.fasta.log"
    stats = _make_stats()
    write_log(
        log_path,
        command_line="taxafasta -i input.fasta -t 2 -o test.fasta",
        input_path=Path("input.fasta"),
        output_path=Path("test.fasta.gz"),
        include_taxids=[2],
        exclude_taxids=None,
        use_merged=True,
        taxdump_source="/data/taxdump",
        allowed_set_size=1000,
        stats=stats,
        elapsed_seconds=10.5,
        names={2: "Bacteria"},
    )
    content = log_path.read_text()
    assert "taxafasta v" in content
    assert "Python" in content
    assert "Command: taxafasta -i input.fasta -t 2 -o test.fasta" in content
    assert "Input: input.fasta" in content
    assert "Output: test.fasta.gz (gzip enabled)" in content
    assert "2 (Bacteria)" in content
    assert "Merged ID resolution: enabled" in content
    assert "/data/taxdump" in content
    assert "1,000" in content  # allowed set size


def test_write_log_contains_warnings(tmp_path: Path) -> None:
    log_path = tmp_path / "test.fasta.log"
    stats = _make_stats(no_ox=1, unknown_count=14)
    stats.unknown_taxids = {99999999: 14}
    write_log(
        log_path,
        command_line="taxafasta -i in.fasta -t 2 -o out.fasta",
        input_path=Path("in.fasta"),
        output_path=Path("out.fasta.gz"),
        include_taxids=[2],
        exclude_taxids=None,
        use_merged=True,
        taxdump_source="/data/taxdump",
        allowed_set_size=1000,
        stats=stats,
        elapsed_seconds=5.0,
    )
    content = log_path.read_text()
    assert "WARNINGS" in content
    assert "[x1] Header could not be parsed" in content
    assert "[x14] Taxonomy ID 99999999" in content


def test_write_log_no_warnings(tmp_path: Path) -> None:
    log_path = tmp_path / "test.fasta.log"
    stats = FilterStats()
    stats.total = 50
    stats.included = 50
    stats.excluded = 0
    write_log(
        log_path,
        command_line="taxafasta -i in.fasta -t 2 -o out.fasta",
        input_path=Path("in.fasta"),
        output_path=Path("out.fasta.gz"),
        include_taxids=[2],
        exclude_taxids=None,
        use_merged=True,
        taxdump_source="/data/taxdump",
        allowed_set_size=1000,
        stats=stats,
        elapsed_seconds=1.0,
    )
    content = log_path.read_text()
    assert "WARNINGS" not in content


def test_write_log_summary_stats(tmp_path: Path) -> None:
    log_path = tmp_path / "test.fasta.log"
    stats = _make_stats(total=100, included=60, excluded=38, no_ox=1, unknown_count=1)
    write_log(
        log_path,
        command_line="taxafasta -i in.fasta -t 2 -o out.fasta",
        input_path=Path("in.fasta"),
        output_path=Path("out.fasta"),
        include_taxids=[2],
        exclude_taxids=None,
        use_merged=False,
        taxdump_source="/data",
        allowed_set_size=500,
        stats=stats,
        elapsed_seconds=2.0,
    )
    content = log_path.read_text()
    assert "SUMMARY:" in content
    assert "100" in content
    assert "60" in content
    assert "Merged ID resolution: disabled" in content


def test_write_log_exclude_taxids(tmp_path: Path) -> None:
    log_path = tmp_path / "test.fasta.log"
    stats = FilterStats()
    stats.total = 10
    stats.included = 5
    stats.excluded = 5
    write_log(
        log_path,
        command_line="taxafasta -i in.fasta -t 2759 -e 40674 -o out.fasta",
        input_path=Path("in.fasta"),
        output_path=Path("out.fasta.gz"),
        include_taxids=[2759],
        exclude_taxids=[40674],
        use_merged=True,
        taxdump_source="/data",
        allowed_set_size=100,
        stats=stats,
        elapsed_seconds=1.0,
        names={2759: "Eukaryota", 40674: "Mammalia"},
    )
    content = log_path.read_text()
    assert "2759 (Eukaryota)" in content
    assert "40674 (Mammalia)" in content


def test_write_log_uncompressed_output(tmp_path: Path) -> None:
    log_path = tmp_path / "test.fasta.log"
    stats = FilterStats()
    stats.total = 10
    stats.included = 10
    write_log(
        log_path,
        command_line="taxafasta -i in.fasta -t 2 -o out.fasta --no-gzip",
        input_path=Path("in.fasta"),
        output_path=Path("out.fasta"),
        include_taxids=[2],
        exclude_taxids=None,
        use_merged=True,
        taxdump_source="/data",
        allowed_set_size=100,
        stats=stats,
        elapsed_seconds=0.5,
    )
    content = log_path.read_text()
    assert "uncompressed" in content
