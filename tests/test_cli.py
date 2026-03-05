"""Unit tests for the CLI module (§11.6)."""

from __future__ import annotations

import pytest

from taxafasta import __version__
from taxafasta.cli import build_parser

# --- Argument parsing ---


def test_parser_required_args() -> None:
    parser = build_parser()
    with pytest.raises(SystemExit):
        parser.parse_args([])  # no args


def test_parser_missing_input() -> None:
    parser = build_parser()
    with pytest.raises(SystemExit):
        parser.parse_args(["-t", "2", "-o", "out.fasta"])


def test_parser_missing_taxid() -> None:
    parser = build_parser()
    with pytest.raises(SystemExit):
        parser.parse_args(["-i", "in.fasta", "-o", "out.fasta"])


def test_parser_missing_output() -> None:
    parser = build_parser()
    with pytest.raises(SystemExit):
        parser.parse_args(["-i", "in.fasta", "-t", "2"])


def test_parser_basic() -> None:
    parser = build_parser()
    args = parser.parse_args(["-i", "in.fasta", "-t", "2", "-o", "out.fasta"])
    assert str(args.input) == "in.fasta"
    assert args.taxid == [2]
    assert str(args.output) == "out.fasta"
    assert args.no_gzip is False
    assert args.verbose is False


def test_parser_multiple_taxids() -> None:
    parser = build_parser()
    args = parser.parse_args(
        ["-i", "in.fasta", "-t", "2", "-t", "10239", "-o", "out.fasta"],
    )
    assert args.taxid == [2, 10239]


def test_parser_exclude() -> None:
    parser = build_parser()
    args = parser.parse_args(
        ["-i", "in.fasta", "-t", "2759", "-e", "40674", "-o", "out.fasta"],
    )
    assert args.exclude == [40674]


def test_parser_multiple_excludes() -> None:
    parser = build_parser()
    args = parser.parse_args(
        ["-i", "in.fasta", "-t", "1", "-e", "40674", "-e", "10239", "-o", "out.fasta"],
    )
    assert args.exclude == [40674, 10239]


def test_parser_no_gzip() -> None:
    parser = build_parser()
    args = parser.parse_args(
        ["-i", "in.fasta", "-t", "2", "-o", "out.fasta", "--no-gzip"],
    )
    assert args.no_gzip is True


def test_parser_no_merge() -> None:
    parser = build_parser()
    args = parser.parse_args(
        ["-i", "in.fasta", "-t", "2", "-o", "out.fasta", "--no-merge"],
    )
    assert args.no_merge is True


def test_parser_verbose() -> None:
    parser = build_parser()
    args = parser.parse_args(["-i", "in.fasta", "-t", "2", "-o", "out.fasta", "-v"])
    assert args.verbose is True


# --- Version ---


def test_version(capsys: pytest.CaptureFixture[str]) -> None:
    parser = build_parser()
    with pytest.raises(SystemExit) as exc_info:
        parser.parse_args(["--version"])
    assert exc_info.value.code == 0
    captured = capsys.readouterr()
    assert __version__ in captured.out


# --- Help ---


def test_help(capsys: pytest.CaptureFixture[str]) -> None:
    parser = build_parser()
    with pytest.raises(SystemExit) as exc_info:
        parser.parse_args(["--help"])
    assert exc_info.value.code == 0
    captured = capsys.readouterr()
    assert "taxafasta" in captured.out
    assert "--input" in captured.out
    assert "--taxid" in captured.out
