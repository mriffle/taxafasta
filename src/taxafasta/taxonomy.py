"""Parsing NCBI taxonomy dump files and building descendant sets (§4.1)."""

from __future__ import annotations

import sys
from collections import defaultdict, deque
from pathlib import Path
from typing import Optional


def parse_nodes(path: Path) -> dict[int, int]:
    """Parse ``nodes.dmp`` into a {tax_id: parent_tax_id} dictionary.

    The file uses ``\\t|\\t`` as the column delimiter.
    Only the first two columns are needed.
    """
    parent_of: dict[int, int] = {}
    with open(path, "r", encoding="utf-8") as fh:
        for line in fh:
            parts = line.split("\t|\t", 2)
            if len(parts) < 2:
                continue
            tax_id = int(parts[0].strip())
            parent_id = int(parts[1].strip())
            parent_of[tax_id] = parent_id
    return parent_of


def parse_merged(path: Path) -> dict[int, int]:
    """Parse ``merged.dmp`` into a {old_tax_id: new_tax_id} dictionary.

    Handles chain resolution: if old→mid→new, resolves old→new.
    """
    raw: dict[int, int] = {}
    with open(path, "r", encoding="utf-8") as fh:
        for line in fh:
            parts = line.split("\t|\t", 1)
            if not parts:
                continue
            # Last column ends with \t|\n
            old_id = int(parts[0].strip())
            new_id = int(parts[1].rstrip().rstrip("|").strip()) if len(parts) > 1 else old_id
            raw[old_id] = new_id

    # Resolve chains
    merged_to: dict[int, int] = {}
    for old_id, new_id in raw.items():
        seen: set[int] = {old_id}
        current = new_id
        while current in raw and current not in seen:
            seen.add(current)
            current = raw[current]
        merged_to[old_id] = current
    return merged_to


def parse_names(path: Path) -> dict[int, str]:
    """Parse ``names.dmp`` to extract scientific names.

    Only rows with name class ``scientific name`` are included.
    Returns {tax_id: scientific_name}.
    """
    names: dict[int, str] = {}
    with open(path, "r", encoding="utf-8") as fh:
        for line in fh:
            parts = line.split("\t|\t")
            if len(parts) < 4:
                continue
            name_class = parts[3].rstrip().rstrip("|").strip()
            if name_class == "scientific name":
                tax_id = int(parts[0].strip())
                name = parts[1].strip()
                names[tax_id] = name
    return names


def build_children_index(parent_of: dict[int, int]) -> dict[int, set[int]]:
    """Invert parent_of to create children_of: {tax_id: {child_ids}}."""
    children_of: dict[int, set[int]] = defaultdict(set)
    for tax_id, parent_id in parent_of.items():
        if tax_id != parent_id:  # Skip root's self-reference
            children_of[parent_id].add(tax_id)
    return dict(children_of)


def collect_descendants(
    roots: list[int],
    children_of: dict[int, set[int]],
) -> set[int]:
    """BFS from each root to collect all descendant taxonomy IDs.

    The root IDs themselves are included in the result set.
    """
    allowed: set[int] = set()
    queue: deque[int] = deque(roots)
    while queue:
        node = queue.popleft()
        if node in allowed:
            continue
        allowed.add(node)
        for child in children_of.get(node, ()):
            if child not in allowed:
                queue.append(child)
    return allowed


def expand_with_merged(
    allowed_taxids: set[int],
    merged_to: dict[int, int],
) -> set[int]:
    """Expand allowed_taxids with old merged IDs that map into the set (§4.1 step 5)."""
    expanded = set(allowed_taxids)
    for old_id, new_id in merged_to.items():
        if new_id in expanded:
            expanded.add(old_id)
    return expanded


def build_allowed_set(
    taxdump_dir: Path,
    include_taxids: list[int],
    exclude_taxids: Optional[list[int]] = None,
    *,
    use_merged: bool = True,
) -> tuple[set[int], dict[int, int], dict[int, int], dict[int, str]]:
    """Build the complete allowed taxonomy ID set.

    Parameters
    ----------
    taxdump_dir : directory containing nodes.dmp, merged.dmp, etc.
    include_taxids : root taxonomy IDs to include (with all descendants)
    exclude_taxids : root taxonomy IDs to exclude (with all descendants)
    use_merged : whether to resolve merged IDs

    Returns
    -------
    (allowed_taxids, parent_of, merged_to, names)
    """
    nodes_path = taxdump_dir / "nodes.dmp"
    merged_path = taxdump_dir / "merged.dmp"
    names_path = taxdump_dir / "names.dmp"

    if not nodes_path.exists():
        print(f"Error: {nodes_path} not found.", file=sys.stderr)
        raise SystemExit(1)
    if not merged_path.exists() and use_merged:
        print(f"Error: {merged_path} not found.", file=sys.stderr)
        raise SystemExit(1)

    parent_of = parse_nodes(nodes_path)
    merged_to = parse_merged(merged_path) if use_merged and merged_path.exists() else {}
    names = parse_names(names_path) if names_path.exists() else {}

    # Resolve user-supplied taxids through merged mapping
    resolved_includes: list[int] = []
    for tid in include_taxids:
        resolved = merged_to.get(tid, tid) if use_merged else tid
        if resolved not in parent_of:
            print(
                f"Error: Taxonomy ID {tid} not found in nodes.dmp or merged.dmp.",
                file=sys.stderr,
            )
            raise SystemExit(1)
        resolved_includes.append(resolved)

    children_of = build_children_index(parent_of)
    allowed = collect_descendants(resolved_includes, children_of)

    # Apply exclusions
    if exclude_taxids:
        resolved_excludes: list[int] = []
        for tid in exclude_taxids:
            resolved = merged_to.get(tid, tid) if use_merged else tid
            if resolved not in parent_of:
                print(
                    f"Error: Exclude taxonomy ID {tid} not found in nodes.dmp or merged.dmp.",
                    file=sys.stderr,
                )
                raise SystemExit(1)
            resolved_excludes.append(resolved)
        excluded = collect_descendants(resolved_excludes, children_of)
        allowed -= excluded

    # Expand with merged IDs
    if use_merged and merged_to:
        allowed = expand_with_merged(allowed, merged_to)

    return allowed, parent_of, merged_to, names
