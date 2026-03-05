"""Argument parsing and main orchestration (§5)."""

from __future__ import annotations

import argparse
import sys
import time
from pathlib import Path

from taxafasta import __version__
from taxafasta.download import ensure_taxdump, read_timestamp
from taxafasta.fasta import filter_fasta
from taxafasta.io_utils import open_input, open_output
from taxafasta.run_log import resolve_log_path, write_log
from taxafasta.taxonomy import build_allowed_set


def build_parser() -> argparse.ArgumentParser:
    """Build the CLI argument parser (§5.1)."""
    parser = argparse.ArgumentParser(
        prog="taxafasta",
        description="Filter UniProt FASTA files by NCBI taxonomy.",
    )
    parser.add_argument(
        "--input",
        "-i",
        required=True,
        type=Path,
        help="Path to input UniProt FASTA file (plain or gzipped).",
    )
    parser.add_argument(
        "--taxid",
        "-t",
        required=True,
        action="append",
        type=int,
        help=(
            "NCBI taxonomy ID to include (all descendants included)."
            " Repeat for multiple IDs."
        ),
    )
    parser.add_argument(
        "--output",
        "-o",
        required=True,
        type=Path,
        help="Path to output FASTA file. Gzip-compressed by default.",
    )
    parser.add_argument(
        "--no-gzip",
        action="store_true",
        default=False,
        help="Disable gzip compression of output.",
    )
    parser.add_argument(
        "--taxdump",
        "-d",
        type=Path,
        default=None,
        help="Path to an already-extracted taxdump directory.",
    )
    parser.add_argument(
        "--cache-dir",
        type=Path,
        default=None,
        help="Directory for caching downloaded taxonomy files.",
    )
    parser.add_argument(
        "--exclude",
        "-e",
        action="append",
        type=int,
        default=None,
        help=(
            "Taxonomy ID to exclude (with descendants),"
            " applied after inclusion. Repeat for multiple IDs."
        ),
    )
    parser.add_argument(
        "--no-merge",
        action="store_true",
        default=False,
        help="Disable merged taxonomy ID resolution.",
    )
    parser.add_argument(
        "--verbose",
        "-v",
        action="store_true",
        default=False,
        help="Periodic progress updates to stderr.",
    )
    parser.add_argument(
        "--version",
        action="version",
        version=f"taxafasta {__version__}",
    )
    return parser


def main(argv: list[str] | None = None) -> None:
    """Entry point for the taxafasta CLI."""
    parser = build_parser()
    args = parser.parse_args(argv)

    # Record the full command line for the log
    if argv is None:
        command_line = "taxafasta " + " ".join(sys.argv[1:])
    else:
        command_line = "taxafasta " + " ".join(argv)

    # Validate input file
    input_path: Path = args.input
    if not input_path.exists():
        print(f"Error: Input file does not exist: {input_path}", file=sys.stderr)
        raise SystemExit(1)
    if not input_path.is_file():
        print(f"Error: Input path is not a file: {input_path}", file=sys.stderr)
        raise SystemExit(1)

    use_gzip = not args.no_gzip
    use_merged = not args.no_merge

    # Resolve taxonomy dump
    taxdump_dir = ensure_taxdump(args.taxdump, args.cache_dir)

    # Build taxonomy source description for the log
    if args.taxdump is not None:
        taxdump_source = str(args.taxdump)
    else:
        ts = read_timestamp(taxdump_dir)
        taxdump_source = f"{taxdump_dir} (downloaded {ts})" if ts else str(taxdump_dir)

    # Build allowed set
    print("Building taxonomy allowed set...", file=sys.stderr)
    allowed_taxids, parent_of, merged_to, names = build_allowed_set(
        taxdump_dir,
        args.taxid,
        args.exclude,
        use_merged=use_merged,
    )
    print(f"Allowed taxid set size: {len(allowed_taxids):,}", file=sys.stderr)

    # Build set of all known taxids for unknown-taxid detection
    all_known: set[int] = set(parent_of.keys()) | set(merged_to.keys())

    # Open I/O streams
    try:
        in_stream = open_input(input_path)
    except OSError as exc:
        print(f"Error opening input file: {exc}", file=sys.stderr)
        raise SystemExit(1) from exc

    try:
        out_stream, resolved_output = open_output(args.output, use_gzip=use_gzip)
    except OSError as exc:
        print(f"Error opening output file: {exc}", file=sys.stderr)
        raise SystemExit(1) from exc

    # Filter
    print("Filtering FASTA...", file=sys.stderr)
    start = time.monotonic()
    try:
        stats = filter_fasta(
            in_stream,
            out_stream,
            allowed_taxids,
            all_known,
            verbose=args.verbose,
        )
    finally:
        out_stream.close()
        in_stream.close()
    elapsed = time.monotonic() - start

    # Write log file
    log_path = resolve_log_path(resolved_output)
    write_log(
        log_path,
        command_line=command_line,
        input_path=input_path,
        output_path=resolved_output,
        include_taxids=args.taxid,
        exclude_taxids=args.exclude,
        use_merged=use_merged,
        taxdump_source=taxdump_source,
        allowed_set_size=len(allowed_taxids),
        stats=stats,
        elapsed_seconds=elapsed,
        names=names,
    )

    # Final summary to stderr
    print(
        f"Done. {stats.included:,} / {stats.total:,} entries "
        f"written to {resolved_output}",
        file=sys.stderr,
    )

    # Exit code (§7)
    has_warnings = stats.no_ox > 0 or stats.unknown_taxid_count > 0
    raise SystemExit(2 if has_warnings else 0)
