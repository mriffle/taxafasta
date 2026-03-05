"""Unit tests for the taxonomy module (§11.1)."""

from __future__ import annotations

from pathlib import Path

from taxafasta.taxonomy import (
    build_allowed_set,
    build_children_index,
    collect_descendants,
    expand_with_merged,
    parse_merged,
    parse_names,
    parse_nodes,
)

# --- Parsing nodes.dmp ---


def test_parse_nodes_basic(tiny_nodes_path: Path) -> None:
    parent_of = parse_nodes(tiny_nodes_path)
    assert parent_of[1] == 1  # root is its own parent
    assert parent_of[2] == 131567  # Bacteria -> cellular organisms
    assert parent_of[131567] == 1
    assert parent_of[9606] == 40674  # Human -> Mammalia


def test_parse_nodes_all_entries(tiny_nodes_path: Path) -> None:
    parent_of = parse_nodes(tiny_nodes_path)
    # Our fixture has 18 entries
    assert len(parent_of) == 18


# --- Parsing merged.dmp ---


def test_parse_merged_basic(tiny_merged_path: Path) -> None:
    merged = parse_merged(tiny_merged_path)
    assert merged[50] == 7
    assert merged[51] == 9


def test_parse_merged_chain_resolution(tiny_merged_path: Path) -> None:
    merged = parse_merged(tiny_merged_path)
    # 52 -> 99999, and 99999 is not in merged, so it stays 99999
    assert merged[52] == 99999


# --- Parsing names.dmp ---


def test_parse_names(tiny_names_path: Path) -> None:
    names = parse_names(tiny_names_path)
    assert names[2] == "Bacteria"
    assert names[9606] == "Homo sapiens"
    assert names[10239] == "Viruses"


# --- Children index ---


def test_build_children_index(tiny_nodes_path: Path) -> None:
    parent_of = parse_nodes(tiny_nodes_path)
    children_of = build_children_index(parent_of)
    # Root (1) has children: 131567, 10239
    assert 131567 in children_of[1]
    assert 10239 in children_of[1]
    # Bacteria (2) has children: 32199, 335928, 1706371
    assert 32199 in children_of[2]
    assert 335928 in children_of[2]
    assert 1706371 in children_of[2]
    # Root should NOT be in its own children (self-reference skipped)
    assert 1 not in children_of.get(1, set())


# --- Descendant set construction ---


def test_collect_descendants_single_root(tiny_nodes_path: Path) -> None:
    parent_of = parse_nodes(tiny_nodes_path)
    children_of = build_children_index(parent_of)
    # Bacteria (2) and all descendants
    desc = collect_descendants([2], children_of)
    assert 2 in desc
    assert 7 in desc  # leaf under 6 under 335928 under 2
    assert 9 in desc  # leaf under 32199 under 2
    assert 11 in desc  # leaf under 10 under 1706371 under 2
    assert 9606 not in desc  # Human is eukaryote
    assert 10239 not in desc  # Viruses


def test_collect_descendants_multiple_roots(tiny_nodes_path: Path) -> None:
    parent_of = parse_nodes(tiny_nodes_path)
    children_of = build_children_index(parent_of)
    desc = collect_descendants([2, 10239], children_of)
    assert 7 in desc
    assert 11111 in desc  # virus
    assert 9606 not in desc


def test_collect_descendants_leaf_node(tiny_nodes_path: Path) -> None:
    parent_of = parse_nodes(tiny_nodes_path)
    children_of = build_children_index(parent_of)
    desc = collect_descendants([9606], children_of)
    assert desc == {9606}


def test_collect_descendants_root_returns_everything(tiny_nodes_path: Path) -> None:
    parent_of = parse_nodes(tiny_nodes_path)
    children_of = build_children_index(parent_of)
    desc = collect_descendants([1], children_of)
    assert desc == set(parent_of.keys())


# --- Merged ID expansion ---


def test_expand_with_merged_included(tiny_merged_path: Path) -> None:
    merged = parse_merged(tiny_merged_path)
    allowed = {7, 9}  # taxid 50->7 and 51->9 should be added
    expanded = expand_with_merged(allowed, merged)
    assert 50 in expanded
    assert 51 in expanded


def test_expand_with_merged_excluded() -> None:
    merged = {100: 200}
    allowed = {7, 9}
    expanded = expand_with_merged(allowed, merged)
    assert 100 not in expanded


# --- Exclude logic ---


def test_exclude_subtree(tiny_taxdump_dir: Path) -> None:
    # Include eukaryotes (2759), exclude mammals (40674)
    allowed, _, _, _ = build_allowed_set(
        tiny_taxdump_dir,
        [2759],
        [40674],
    )
    assert 2759 in allowed
    assert 40674 not in allowed
    assert 9606 not in allowed  # Human is mammal


def test_exclude_multiple_subtrees(tiny_taxdump_dir: Path) -> None:
    # Include bacteria (2), exclude two sub-branches: 335928 (family containing 6,7)
    # and 1706371 (family containing 10,11). Only 32199 branch (containing 9) should remain.
    allowed, _, _, _ = build_allowed_set(
        tiny_taxdump_dir,
        [2],
        [335928, 1706371],
    )
    assert 2 in allowed
    assert 32199 in allowed
    assert 9 in allowed
    # Excluded branches
    assert 335928 not in allowed
    assert 6 not in allowed
    assert 7 not in allowed
    assert 1706371 not in allowed
    assert 10 not in allowed
    assert 11 not in allowed


# --- Unknown taxid handling ---


def test_unknown_user_taxid(tiny_taxdump_dir: Path) -> None:
    """User-supplied taxid not in tree should cause SystemExit."""
    import pytest

    with pytest.raises(SystemExit):
        build_allowed_set(tiny_taxdump_dir, [77777777])


# --- build_allowed_set integration ---


def test_build_allowed_set_bacteria(tiny_taxdump_dir: Path) -> None:
    allowed, parent_of, merged_to, names = build_allowed_set(tiny_taxdump_dir, [2])
    # Should include bacteria subtree
    assert 2 in allowed
    assert 7 in allowed
    assert 9 in allowed
    # Merged IDs mapping into bacteria should be included
    assert 50 in allowed  # 50 -> 7
    assert 51 in allowed  # 51 -> 9
    # Non-bacteria
    assert 9606 not in allowed
    assert 10239 not in allowed
    # names should be populated
    assert names[2] == "Bacteria"


def test_build_allowed_set_no_merge(tiny_taxdump_dir: Path) -> None:
    allowed, _, _, _ = build_allowed_set(
        tiny_taxdump_dir,
        [2],
        use_merged=False,
    )
    assert 2 in allowed
    assert 7 in allowed
    # Merged IDs should NOT be present
    assert 50 not in allowed
    assert 51 not in allowed
