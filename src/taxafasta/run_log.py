"""Log file generation and warning accumulation (§5.4, §5.5)."""

from __future__ import annotations

import platform
import sys
from datetime import datetime, timezone
from pathlib import Path

from taxafasta import __version__
from taxafasta.fasta import FilterStats


def resolve_log_path(output_path: Path) -> Path:
    """Determine the log file path, handling collisions (§5.4).

    If the output is ``bacteria.fasta.gz``, the log is ``bacteria.fasta.log``.
    If ``bacteria.fasta.log`` exists, try ``.log1``, ``.log2``, etc.
    """
    # Strip .gz suffix for log naming
    base = output_path
    if base.name.endswith(".gz"):
        base = base.with_name(base.name[:-3])

    log_path = base.with_suffix(base.suffix + ".log")
    if not log_path.exists():
        return log_path

    idx = 1
    while True:
        candidate = log_path.with_name(log_path.name + str(idx))
        if not candidate.exists():
            return candidate
        idx += 1


def write_log(
    log_path: Path,
    *,
    command_line: str,
    input_path: Path,
    output_path: Path,
    include_taxids: list[int],
    exclude_taxids: list[int] | None,
    use_merged: bool,
    taxdump_source: str,
    allowed_set_size: int,
    stats: FilterStats,
    elapsed_seconds: float,
    names: dict[int, str] | None = None,
) -> None:
    """Write the run log file (§5.4)."""
    timestamp = datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")
    py_version = platform.python_version()
    plat = f"{platform.system()} {platform.machine()}"

    rate = stats.total / elapsed_seconds if elapsed_seconds > 0 else 0
    h, rem = divmod(int(elapsed_seconds), 3600)
    m, s = divmod(rem, 60)
    elapsed_str = f"{h:02d}:{m:02d}:{s:02d}"

    lines: list[str] = []
    lines.append(
        f"taxafasta v{__version__} | Python {py_version} | {plat} | {timestamp}"
    )
    lines.append(f"Command: {command_line}")
    lines.append(f"Input: {input_path}")

    gz_note = "gzip enabled" if output_path.name.endswith(".gz") else "uncompressed"
    lines.append(f"Output: {output_path} ({gz_note})")

    # Included taxids with optional names
    inc_parts: list[str] = []
    for tid in include_taxids:
        if names and tid in names:
            inc_parts.append(f"{tid} ({names[tid]})")
        else:
            inc_parts.append(str(tid))
    lines.append(f"Included taxids: {', '.join(inc_parts)}")

    # Excluded taxids
    if exclude_taxids:
        exc_parts: list[str] = []
        for tid in exclude_taxids:
            if names and tid in names:
                exc_parts.append(f"{tid} ({names[tid]})")
            else:
                exc_parts.append(str(tid))
        lines.append(f"Excluded taxids: {', '.join(exc_parts)}")
    else:
        lines.append("Excluded taxids: (none)")

    lines.append(f"Merged ID resolution: {'enabled' if use_merged else 'disabled'}")
    lines.append(f"Taxdump source: {taxdump_source}")
    lines.append(f"Allowed taxid set size: {allowed_set_size:,}")
    lines.append("")

    # Warnings section (§5.5)
    total_warnings = len(stats.unparseable_headers) + len(stats.unknown_taxids)
    if total_warnings > 0:
        lines.append(f"WARNINGS ({total_warnings} unique):")
        for hdr in stats.unparseable_headers:
            lines.append(
                f"  [x{stats.no_ox}] Header could not be parsed for OX field: {hdr}"
            )
        for taxid, count in sorted(stats.unknown_taxids.items()):
            lines.append(
                f"  [x{count}] Taxonomy ID {taxid} not found in NCBI taxonomy data"
            )
        lines.append("")

    # Summary
    lines.append("SUMMARY:")
    lines.append(f"  Total entries processed:  {stats.total:>15,}")
    lines.append(f"  Entries included:         {stats.included:>15,}")
    lines.append(f"  Entries excluded:         {stats.excluded:>15,}")
    lines.append(f"  Entries skipped (no OX):  {stats.no_ox:>15,}")
    lines.append(f"  Entries with unknown taxid:{stats.unknown_taxid_count:>14,}")
    lines.append(f"  Elapsed time:             {elapsed_str:>15}")
    lines.append(f"  Processing rate:          ~{rate:>13,.0f} entries/sec")

    log_path.write_text("\n".join(lines) + "\n", encoding="utf-8")
    print(f"Log written to {log_path}", file=sys.stderr)
