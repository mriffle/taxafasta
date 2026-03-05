"""Streaming FASTA reader/writer with OX field extraction (§4.2)."""

from __future__ import annotations

import re
import sys
import time
from typing import IO

# Compiled regex for OX= field extraction (§4.2 step 2a)
_OX_PATTERN = re.compile(r"OX=(\d+)")


def extract_ox(header: str) -> int | None:
    """Extract the OX= taxonomy ID from a UniProt FASTA header.

    Returns the integer taxonomy ID, or None if no valid OX= field is found.
    """
    match = _OX_PATTERN.search(header)
    if match is None:
        return None
    return int(match.group(1))


class FilterStats:
    """Accumulates statistics and warnings during FASTA filtering."""

    __slots__ = (
        "total",
        "included",
        "excluded",
        "no_ox",
        "unknown_taxid_count",
        "unknown_taxids",
        "unparseable_headers",
        "_first_no_ox_warned",
    )

    def __init__(self) -> None:
        self.total: int = 0
        self.included: int = 0
        self.excluded: int = 0
        self.no_ox: int = 0
        self.unknown_taxid_count: int = 0
        self.unknown_taxids: dict[int, int] = {}  # taxid -> count
        self.unparseable_headers: list[str] = []  # first occurrence(s)
        self._first_no_ox_warned: bool = False

    def record_no_ox(self, header: str) -> None:
        """Record a header with no valid OX= field."""
        self.no_ox += 1
        if not self._first_no_ox_warned:
            self._first_no_ox_warned = True
            self.unparseable_headers.append(header.rstrip())
            print(
                f"WARNING: Header could not be parsed for OX field: {header.rstrip()}",
                file=sys.stderr,
            )

    def record_unknown_taxid(self, taxid: int) -> None:
        """Record an entry with a valid OX= value not in taxonomy data."""
        self.unknown_taxid_count += 1
        if taxid not in self.unknown_taxids:
            self.unknown_taxids[taxid] = 0
            print(
                f"WARNING: Taxonomy ID {taxid} not found in NCBI taxonomy data",
                file=sys.stderr,
            )
        self.unknown_taxids[taxid] += 1


def filter_fasta(
    input_stream: IO[str],
    output_stream: IO[str],
    allowed_taxids: set[int],
    all_known_taxids: set[int],
    *,
    verbose: bool = False,
    progress_interval: int = 1_000_000,
) -> FilterStats:
    """Stream-filter a FASTA file, writing only entries whose OX= taxid is allowed.

    Parameters
    ----------
    input_stream : text-mode readable stream of FASTA data
    output_stream : text-mode writable stream for filtered output
    allowed_taxids : pre-computed set of taxonomy IDs that pass the filter
    all_known_taxids : full set of taxonomy IDs in nodes.dmp + merged.dmp
    verbose : if True, emit periodic progress to stderr
    progress_interval : number of entries between progress reports

    Returns
    -------
    FilterStats with counts and accumulated warnings
    """
    stats = FilterStats()
    include_current = False
    start_time = time.monotonic()

    for line in input_stream:
        if line.startswith(">"):
            # New FASTA entry header
            stats.total += 1
            taxid = extract_ox(line)

            if taxid is None:
                stats.record_no_ox(line)
                include_current = False
            elif taxid in allowed_taxids:
                include_current = True
                stats.included += 1
            else:
                include_current = False
                stats.excluded += 1
                if taxid not in all_known_taxids:
                    stats.record_unknown_taxid(taxid)

            if include_current:
                output_stream.write(line)

            # Verbose progress reporting (§5.3)
            if verbose and stats.total % progress_interval == 0:
                elapsed = time.monotonic() - start_time
                rate = stats.total / elapsed if elapsed > 0 else 0
                h, rem = divmod(int(elapsed), 3600)
                m, s = divmod(rem, 60)
                print(
                    f"[{h:02d}:{m:02d}:{s:02d}] "
                    f"Processed {stats.total:,} entries | "
                    f"Included: {stats.included:,} | "
                    f"Excluded: {stats.total - stats.included:,} | "
                    f"Rate: ~{rate:,.0f} entries/sec",
                    file=sys.stderr,
                )
        else:
            # Sequence line — pass through only if current entry is included
            if include_current:
                output_stream.write(line)

    return stats
