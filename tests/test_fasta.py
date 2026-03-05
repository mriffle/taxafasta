"""Unit tests for the FASTA module (§11.2)."""

from __future__ import annotations

import io

from taxafasta.fasta import FilterStats, extract_ox, filter_fasta

# --- OX field extraction ---


def test_extract_ox_swissprot() -> None:
    header = (
        ">sp|P12345|ARATH_ECOLI Some protein"
        " OS=Azorhizobium caulinodans OX=7 GN=gene1 PE=1 SV=1"
    )
    assert extract_ox(header) == 7


def test_extract_ox_trembl() -> None:
    header = (
        ">tr|Q99999|HUMAN_PROT Human protein OS=Homo sapiens OX=9606 GN=TP53 PE=1 SV=2"
    )
    assert extract_ox(header) == 9606


def test_extract_ox_no_gn_field() -> None:
    header = (
        ">sp|P99998|VIRUS_PROT Viral protein OS=Tobacco mosaic virus OX=11111 PE=1 SV=1"
    )
    assert extract_ox(header) == 11111


def test_extract_ox_missing() -> None:
    header = ">gnl|custom|seq1 hypothetical protein no OX field here"
    assert extract_ox(header) is None


def test_extract_ox_malformed_abc() -> None:
    header = ">sp|P00001|FAKE Some protein OS=Org OX=abc PE=1 SV=1"
    assert extract_ox(header) is None


def test_extract_ox_malformed_empty() -> None:
    header = ">sp|P00001|FAKE Some protein OS=Org OX= PE=1 SV=1"
    assert extract_ox(header) is None


def test_extract_ox_large_id() -> None:
    header = ">tr|Z00001|UNKN Some protein OS=Unknown OX=99999999 PE=1 SV=1"
    assert extract_ox(header) == 99999999


# --- FilterStats ---


def test_filter_stats_no_ox_warning(capsys: object) -> None:
    stats = FilterStats()
    stats.record_no_ox(">bad header line 1")
    stats.record_no_ox(">bad header line 2")
    assert stats.no_ox == 2
    assert len(stats.unparseable_headers) == 1  # Only first recorded


def test_filter_stats_unknown_taxid(capsys: object) -> None:
    stats = FilterStats()
    stats.record_unknown_taxid(99999)
    stats.record_unknown_taxid(99999)
    stats.record_unknown_taxid(88888)
    assert stats.unknown_taxid_count == 3
    assert stats.unknown_taxids[99999] == 2
    assert stats.unknown_taxids[88888] == 1


# --- filter_fasta ---


def _make_fasta(*entries: tuple[str, str]) -> str:
    """Helper: build a FASTA string from (header, sequence) tuples."""
    lines = []
    for header, seq in entries:
        lines.append(header + "\n")
        if seq:
            lines.append(seq + "\n")
    return "".join(lines)


def test_filter_fasta_basic_include() -> None:
    fasta = _make_fasta(
        (">sp|P1|X OS=Org OX=7 PE=1 SV=1", "ACGT"),
        (">sp|P2|Y OS=Org OX=9606 PE=1 SV=1", "TTTT"),
    )
    allowed = {7}
    all_known = {7, 9606}
    inp = io.StringIO(fasta)
    out = io.StringIO()
    stats = filter_fasta(inp, out, allowed, all_known)
    output = out.getvalue()
    assert ">sp|P1|X" in output
    assert "ACGT" in output
    assert ">sp|P2|Y" not in output
    assert stats.included == 1
    assert stats.excluded == 1


def test_filter_fasta_multi_line_sequence() -> None:
    header1 = ">sp|P1|X OS=Org OX=7 PE=1 SV=1\n"
    seq1 = "ACGT\nTTTT\nGGGG\n"
    header2 = ">sp|P2|Y OS=Org OX=9606 PE=1 SV=1\n"
    seq2 = "AAAA\nCCCC\n"
    fasta = header1 + seq1 + header2 + seq2
    allowed = {7}
    all_known = {7, 9606}
    inp = io.StringIO(fasta)
    out = io.StringIO()
    stats = filter_fasta(inp, out, allowed, all_known)
    assert out.getvalue() == header1 + seq1
    assert stats.included == 1
    assert stats.excluded == 1


def test_filter_fasta_empty_sequence() -> None:
    fasta = ">sp|P1|X OS=Org OX=7 PE=1 SV=1\n>sp|P2|Y OS=Org OX=9606 PE=1 SV=1\nAAAA\n"
    allowed = {7}
    all_known = {7, 9606}
    inp = io.StringIO(fasta)
    out = io.StringIO()
    stats = filter_fasta(inp, out, allowed, all_known)
    output = out.getvalue()
    assert ">sp|P1|X" in output
    assert "AAAA" not in output
    assert stats.included == 1
    assert stats.excluded == 1


def test_filter_fasta_missing_ox() -> None:
    fasta = ">gnl|custom|seq1 no OX field\nACGT\n"
    allowed = {7}
    all_known = {7}
    inp = io.StringIO(fasta)
    out = io.StringIO()
    stats = filter_fasta(inp, out, allowed, all_known)
    assert out.getvalue() == ""
    assert stats.no_ox == 1
    assert stats.total == 1


def test_filter_fasta_unknown_taxid() -> None:
    fasta = ">sp|P1|X OS=Org OX=99999999 PE=1 SV=1\nACGT\n"
    allowed = {7}
    all_known = {7}
    inp = io.StringIO(fasta)
    out = io.StringIO()
    stats = filter_fasta(inp, out, allowed, all_known)
    assert out.getvalue() == ""
    assert stats.unknown_taxid_count == 1
    assert 99999999 in stats.unknown_taxids


def test_filter_fasta_entry_boundary_correctness() -> None:
    """Interleaved matching and non-matching entries should not cross-contaminate."""
    fasta = _make_fasta(
        (">sp|P1|A OS=Org OX=7 PE=1 SV=1", "AAAA"),
        (">sp|P2|B OS=Org OX=9606 PE=1 SV=1", "BBBB"),
        (">sp|P3|C OS=Org OX=9 PE=1 SV=1", "CCCC"),
        (">sp|P4|D OS=Org OX=2759 PE=1 SV=1", "DDDD"),
    )
    allowed = {7, 9}
    all_known = {7, 9, 9606, 2759}
    inp = io.StringIO(fasta)
    out = io.StringIO()
    stats = filter_fasta(inp, out, allowed, all_known)
    output = out.getvalue()
    assert "AAAA" in output
    assert "BBBB" not in output
    assert "CCCC" in output
    assert "DDDD" not in output
    assert stats.included == 2
    assert stats.excluded == 2


def test_filter_fasta_output_fidelity() -> None:
    """Output should be byte-for-byte identical to manual extraction."""
    line1 = ">sp|P1|A OS=Org OX=7 PE=1 SV=1\n"
    seq1 = "ACGTACGT\n"
    line2 = ">sp|P2|B OS=Org OX=9606 PE=1 SV=1\n"
    seq2 = "TTTTTTTT\n"
    fasta = line1 + seq1 + line2 + seq2
    allowed = {7}
    all_known = {7, 9606}
    inp = io.StringIO(fasta)
    out = io.StringIO()
    filter_fasta(inp, out, allowed, all_known)
    assert out.getvalue() == line1 + seq1


def test_filter_fasta_output_fidelity_interleaved_multiline() -> None:
    """Included entries should preserve exact sequence content and formatting."""
    header1 = ">sp|P1|A OS=Org OX=7 PE=1 SV=1\n"
    seq1 = "AAAA\nBBBB \nCCCC\n"
    header2 = ">sp|P2|B OS=Org OX=9606 PE=1 SV=1\n"
    seq2 = "XXXX\nYYYY\n"
    header3 = ">sp|P3|C OS=Org OX=9 PE=1 SV=1\n"
    seq3 = "DDDD\nEEEE\n"
    fasta = header1 + seq1 + header2 + seq2 + header3 + seq3

    allowed = {7, 9}
    all_known = {7, 9, 9606}
    inp = io.StringIO(fasta)
    out = io.StringIO()
    stats = filter_fasta(inp, out, allowed, all_known)

    assert out.getvalue() == header1 + seq1 + header3 + seq3
    assert stats.included == 2
    assert stats.excluded == 1


def test_filter_fasta_output_fidelity_preserves_crlf_line_endings() -> None:
    """When input uses CRLF, output should preserve CRLF exactly."""
    header1 = ">sp|P1|A OS=Org OX=7 PE=1 SV=1\r\n"
    seq1 = "ACGT\r\nTTTT\r\n"
    header2 = ">sp|P2|B OS=Org OX=9606 PE=1 SV=1\r\n"
    seq2 = "GGGG\r\n"
    fasta = header1 + seq1 + header2 + seq2

    allowed = {7}
    all_known = {7, 9606}
    inp = io.StringIO(fasta)
    out = io.StringIO()
    filter_fasta(inp, out, allowed, all_known)

    assert out.getvalue() == header1 + seq1


def test_filter_fasta_large_header() -> None:
    """Very long header lines should be handled correctly."""
    long_desc = "A" * 5000
    header = f">sp|P1|X {long_desc} OS=Org OX=7 PE=1 SV=1"
    fasta = header + "\nACGT\n"
    allowed = {7}
    all_known = {7}
    inp = io.StringIO(fasta)
    out = io.StringIO()
    stats = filter_fasta(inp, out, allowed, all_known)
    assert stats.included == 1
    assert header in out.getvalue()
