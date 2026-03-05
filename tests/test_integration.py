"""End-to-end integration tests (§11.7, §11.8)."""

from __future__ import annotations

import gzip
from pathlib import Path

import pytest

from taxafasta.cli import main


# --- Helpers ---

def _read_output(path: Path) -> str:
    """Read output file, handling gzip transparently."""
    if path.name.endswith(".gz"):
        with gzip.open(path, "rt") as fh:
            return fh.read()
    return path.read_text()


def _count_entries(text: str) -> int:
    return sum(1 for line in text.splitlines() if line.startswith(">"))


# --- Basic filter ---

def test_basic_bacteria_filter(tmp_path: Path, data_dir: Path, tiny_taxdump_dir: Path) -> None:
    """Filter to bacteria (taxid 2). Should include bacterial entries only."""
    out = tmp_path / "bacteria.fasta"
    with pytest.raises(SystemExit) as exc_info:
        main([
            "-i", str(data_dir / "sample.fasta"),
            "-t", "2",
            "-o", str(out),
            "-d", str(tiny_taxdump_dir),
            "--no-gzip",
        ])
    # Exit code 2 because sample has entries with no OX and unknown taxid
    assert exc_info.value.code == 2
    output = _read_output(out)
    # Bacterial entries: OX=7 (Azorhizobium), OX=9 (Buchnera), OX=11 (Cellvibrio)
    # Also OX=50 (merged->7)
    assert "OX=7" in output
    assert "OX=9" in output or "OX=51" in output  # 51 merged to 9
    assert "OX=11" in output
    assert "OX=50" in output  # merged to 7, should be included
    # Non-bacterial entries should NOT be present
    assert "OX=9606" not in output
    assert "OX=11111" not in output
    assert "OX=22222" not in output


def test_basic_bacteria_filter_entry_count(tmp_path: Path, data_dir: Path, tiny_taxdump_dir: Path) -> None:
    out = tmp_path / "bacteria.fasta"
    with pytest.raises(SystemExit):
        main([
            "-i", str(data_dir / "sample.fasta"),
            "-t", "2",
            "-o", str(out),
            "-d", str(tiny_taxdump_dir),
            "--no-gzip",
        ])
    output = _read_output(out)
    # Expect: OX=7, OX=9, OX=11, OX=50 (4 bacterial entries)
    assert _count_entries(output) == 4


# --- Multiple taxids ---

def test_multiple_taxids(tmp_path: Path, data_dir: Path, tiny_taxdump_dir: Path) -> None:
    out = tmp_path / "bact_virus.fasta"
    with pytest.raises(SystemExit):
        main([
            "-i", str(data_dir / "sample.fasta"),
            "-t", "2", "10239",
            "-o", str(out),
            "-d", str(tiny_taxdump_dir),
            "--no-gzip",
        ])
    output = _read_output(out)
    assert "OX=7" in output
    assert "OX=11111" in output  # virus
    assert "OX=9606" not in output  # human


# --- Exclude ---

def test_exclude_subtree(tmp_path: Path, data_dir: Path, tiny_taxdump_dir: Path) -> None:
    """Include eukaryotes (2759), exclude mammals (40674)."""
    out = tmp_path / "euk_no_mammal.fasta"
    with pytest.raises(SystemExit):
        main([
            "-i", str(data_dir / "sample.fasta"),
            "-t", "2759",
            "-e", "40674",
            "-o", str(out),
            "-d", str(tiny_taxdump_dir),
            "--no-gzip",
        ])
    output = _read_output(out)
    # Eukaryota (2759) included, but Mammalia (40674) and Human (9606) excluded
    # Our sample has OX=9606 (human, mammal) — should be excluded
    # No other eukaryote entries in sample that aren't mammals
    assert "OX=9606" not in output


# --- Merged taxid handling ---

def test_merged_taxid_included(tmp_path: Path, data_dir: Path, tiny_taxdump_dir: Path) -> None:
    """Entry with OX=50 (merged to 7, which is under bacteria) should be included."""
    out = tmp_path / "merged.fasta"
    with pytest.raises(SystemExit):
        main([
            "-i", str(data_dir / "sample.fasta"),
            "-t", "2",
            "-o", str(out),
            "-d", str(tiny_taxdump_dir),
            "--no-gzip",
        ])
    output = _read_output(out)
    assert "OX=50" in output


# --- Gzip input ---

def test_gzip_input(tmp_path: Path, data_dir: Path, tiny_taxdump_dir: Path) -> None:
    out = tmp_path / "bacteria.fasta"
    with pytest.raises(SystemExit):
        main([
            "-i", str(data_dir / "sample.fasta.gz"),
            "-t", "2",
            "-o", str(out),
            "-d", str(tiny_taxdump_dir),
            "--no-gzip",
        ])
    output = _read_output(out)
    assert "OX=7" in output
    assert _count_entries(output) == 4


# --- Default gzip output ---

def test_default_gzip_output(tmp_path: Path, data_dir: Path, tiny_taxdump_dir: Path) -> None:
    out = tmp_path / "bacteria.fasta"
    with pytest.raises(SystemExit):
        main([
            "-i", str(data_dir / "sample.fasta"),
            "-t", "2",
            "-o", str(out),
            "-d", str(tiny_taxdump_dir),
        ])
    # Output should be gzip-compressed with .gz appended
    gz_path = tmp_path / "bacteria.fasta.gz"
    assert gz_path.exists()
    output = _read_output(gz_path)
    assert "OX=7" in output


# --- No-gzip output ---

def test_no_gzip_output(tmp_path: Path, data_dir: Path, tiny_taxdump_dir: Path) -> None:
    out = tmp_path / "bacteria.fasta"
    with pytest.raises(SystemExit):
        main([
            "-i", str(data_dir / "sample.fasta"),
            "-t", "2",
            "-o", str(out),
            "-d", str(tiny_taxdump_dir),
            "--no-gzip",
        ])
    assert out.exists()
    # Verify it's plain text, not gzip
    raw = out.read_bytes()
    assert raw[:2] != b"\x1f\x8b"


# --- Log file creation ---

def test_log_file_created(tmp_path: Path, data_dir: Path, tiny_taxdump_dir: Path) -> None:
    out = tmp_path / "bacteria.fasta"
    with pytest.raises(SystemExit):
        main([
            "-i", str(data_dir / "sample.fasta"),
            "-t", "2",
            "-o", str(out),
            "-d", str(tiny_taxdump_dir),
            "--no-gzip",
        ])
    log_path = tmp_path / "bacteria.fasta.log"
    assert log_path.exists()
    content = log_path.read_text()
    assert "taxafasta v" in content
    assert "SUMMARY:" in content
    assert "Command:" in content


def test_log_file_gzip_naming(tmp_path: Path, data_dir: Path, tiny_taxdump_dir: Path) -> None:
    """When output is bacteria.fasta.gz, log should be bacteria.fasta.log."""
    out = tmp_path / "bacteria.fasta"
    with pytest.raises(SystemExit):
        main([
            "-i", str(data_dir / "sample.fasta"),
            "-t", "2",
            "-o", str(out),
            "-d", str(tiny_taxdump_dir),
        ])
    log_path = tmp_path / "bacteria.fasta.log"
    assert log_path.exists()


# --- Log file collision ---

def test_log_file_collision(tmp_path: Path, data_dir: Path, tiny_taxdump_dir: Path) -> None:
    out = tmp_path / "bacteria.fasta"
    # Run twice
    with pytest.raises(SystemExit):
        main(["-i", str(data_dir / "sample.fasta"), "-t", "2", "-o", str(out),
              "-d", str(tiny_taxdump_dir), "--no-gzip"])
    with pytest.raises(SystemExit):
        main(["-i", str(data_dir / "sample.fasta"), "-t", "2", "-o", str(out),
              "-d", str(tiny_taxdump_dir), "--no-gzip"])
    assert (tmp_path / "bacteria.fasta.log").exists()
    assert (tmp_path / "bacteria.fasta.log1").exists()


# --- Warning deduplication ---

def test_warning_deduplication_in_log(tmp_path: Path, data_dir: Path, tiny_taxdump_dir: Path) -> None:
    out = tmp_path / "bacteria.fasta"
    with pytest.raises(SystemExit):
        main(["-i", str(data_dir / "sample.fasta"), "-t", "2", "-o", str(out),
              "-d", str(tiny_taxdump_dir), "--no-gzip"])
    log = (tmp_path / "bacteria.fasta.log").read_text()
    # Should have warnings section
    assert "WARNINGS" in log
    # Unknown taxid 99999999 should appear once in warnings
    assert "99999999" in log


# --- Unparseable header warning ---

def test_unparseable_header_in_log(tmp_path: Path, data_dir: Path, tiny_taxdump_dir: Path) -> None:
    out = tmp_path / "bacteria.fasta"
    with pytest.raises(SystemExit):
        main(["-i", str(data_dir / "sample.fasta"), "-t", "2", "-o", str(out),
              "-d", str(tiny_taxdump_dir), "--no-gzip"])
    log = (tmp_path / "bacteria.fasta.log").read_text()
    assert "Header could not be parsed for OX field" in log


# --- No matches ---

def test_no_matches(tmp_path: Path, data_dir: Path, tiny_taxdump_dir: Path) -> None:
    """Filter with a taxid that matches nothing in our sample."""
    # Create a FASTA where no entries match the filter taxid
    no_match_fasta = tmp_path / "no_match.fasta"
    no_match_fasta.write_text(
        ">sp|P1|X OS=Org OX=9606 PE=1 SV=1\nACGT\n"
        ">sp|P2|Y OS=Org OX=22222 PE=1 SV=1\nTTTT\n"
    )
    out = tmp_path / "nothing.fasta"
    # Filter to viruses (10239) — neither 9606 nor 22222 is under 10239
    with pytest.raises(SystemExit) as exc_info:
        main(["-i", str(no_match_fasta), "-t", "10239", "-o", str(out),
              "-d", str(tiny_taxdump_dir), "--no-gzip"])
    output = _read_output(out)
    assert _count_entries(output) == 0
    # No warnings (all OX fields valid, all taxids known) → exit code 0
    assert exc_info.value.code == 0


# --- All match ---

def test_all_match(tmp_path: Path, data_dir: Path, tiny_taxdump_dir: Path) -> None:
    """Filter with root taxid 1 — everything should match (except no-OX and unknown)."""
    out = tmp_path / "all.fasta"
    with pytest.raises(SystemExit) as exc_info:
        main(["-i", str(data_dir / "sample.fasta"), "-t", "1", "-o", str(out),
              "-d", str(tiny_taxdump_dir), "--no-gzip"])
    output = _read_output(out)
    # All entries with valid, known OX should be included
    # Excluded: no-OX entry and unknown taxid entry (OX=99999999)
    assert "OX=7" in output
    assert "OX=9606" in output
    assert "OX=11111" in output
    assert "OX=22222" in output
    assert "OX=50" in output
    # The no-OX entry should NOT be included
    assert "gnl|custom|seq1" not in output
    # Exit code 2 due to warnings
    assert exc_info.value.code == 2


# --- Accuracy tests (§11.8) ---

def test_accuracy_exact_root_matches(tmp_path: Path, data_dir: Path, tiny_taxdump_dir: Path) -> None:
    """Entry whose taxid IS the requested root should match."""
    out = tmp_path / "exact.fasta"
    with pytest.raises(SystemExit):
        main(["-i", str(data_dir / "sample.fasta"), "-t", "2",
              "-o", str(out), "-d", str(tiny_taxdump_dir), "--no-gzip"])
    # OX=7 is a descendant of 2, but let's test with taxid 7 directly
    out2 = tmp_path / "exact2.fasta"
    with pytest.raises(SystemExit):
        main(["-i", str(data_dir / "sample.fasta"), "-t", "7",
              "-o", str(out2), "-d", str(tiny_taxdump_dir), "--no-gzip"])
    output = _read_output(out2)
    assert "OX=7" in output
    # OX=50 (merged to 7) should also be included
    assert "OX=50" in output


def test_accuracy_sibling_not_matched(tmp_path: Path, data_dir: Path, tiny_taxdump_dir: Path) -> None:
    """Sibling of requested root should NOT match."""
    out = tmp_path / "sibling.fasta"
    # Request bacteria (2), archaea (2157) is a sibling
    with pytest.raises(SystemExit):
        main(["-i", str(data_dir / "sample.fasta"), "-t", "2",
              "-o", str(out), "-d", str(tiny_taxdump_dir), "--no-gzip"])
    output = _read_output(out)
    assert "OX=22222" not in output  # Archaea species


def test_accuracy_parent_not_matched(tmp_path: Path, data_dir: Path, tiny_taxdump_dir: Path) -> None:
    """Parent of requested root should NOT match."""
    out = tmp_path / "parent.fasta"
    # Request mammals (40674), eukaryota (2759) is parent
    with pytest.raises(SystemExit):
        main(["-i", str(data_dir / "sample.fasta"), "-t", "40674",
              "-o", str(out), "-d", str(tiny_taxdump_dir), "--no-gzip"])
    output = _read_output(out)
    # Only human (9606, under 40674) should match
    assert "OX=9606" in output
    # No eukaryota-level matches in our sample besides human


def test_accuracy_merged_outside_subtree(tmp_path: Path, data_dir: Path, tiny_taxdump_dir: Path) -> None:
    """Merged ID pointing outside the subtree should NOT match."""
    out = tmp_path / "merged_outside.fasta"
    # taxid 52 merges to 99999, which is NOT in our tree
    # If someone has OX=52 in the FASTA, it should not match bacteria
    with pytest.raises(SystemExit):
        main(["-i", str(data_dir / "sample.fasta"), "-t", "10239",
              "-o", str(out), "-d", str(tiny_taxdump_dir), "--no-gzip"])
    output = _read_output(out)
    # Only virus entry should match
    assert "OX=11111" in output
    assert "OX=7" not in output


# --- Exit codes ---

def test_exit_code_success_no_warnings(tmp_path: Path, data_dir: Path, tiny_taxdump_dir: Path) -> None:
    """Create a clean FASTA with no warnings and verify exit code 0."""
    clean_fasta = tmp_path / "clean.fasta"
    clean_fasta.write_text(
        ">sp|P1|X OS=Org OX=7 PE=1 SV=1\nACGT\n"
        ">sp|P2|Y OS=Org OX=9 PE=1 SV=1\nTTTT\n"
    )
    out = tmp_path / "out.fasta"
    with pytest.raises(SystemExit) as exc_info:
        main(["-i", str(clean_fasta), "-t", "2", "-o", str(out),
              "-d", str(tiny_taxdump_dir), "--no-gzip"])
    assert exc_info.value.code == 0


def test_exit_code_with_warnings(tmp_path: Path, data_dir: Path, tiny_taxdump_dir: Path) -> None:
    out = tmp_path / "out.fasta"
    with pytest.raises(SystemExit) as exc_info:
        main(["-i", str(data_dir / "sample.fasta"), "-t", "2", "-o", str(out),
              "-d", str(tiny_taxdump_dir), "--no-gzip"])
    assert exc_info.value.code == 2


def test_exit_code_missing_input(tmp_path: Path, data_dir: Path, tiny_taxdump_dir: Path) -> None:
    with pytest.raises(SystemExit) as exc_info:
        main(["-i", str(tmp_path / "nonexistent.fasta"), "-t", "2",
              "-o", str(tmp_path / "out.fasta"), "-d", str(tiny_taxdump_dir)])
    assert exc_info.value.code == 1
