# `taxafasta` — UniProt FASTA Taxonomy Filter

## Software Specification v1.0

---

## 1. Purpose

`taxafasta` is a command-line tool that filters large UniProt protein FASTA files by NCBI taxonomy. Given one or more NCBI taxonomy IDs, it outputs only those FASTA entries whose organism belongs to the specified taxon or any of its descendant taxa. For example, supplying taxonomy ID `2` (Bacteria) produces a bacteria-only subset of the input file.

---

## 2. Scope & Scale

The tool is designed for datasets at the scale of UniProt TrEMBL, which as of early 2025 contains approximately 250 million sequence entries. The compressed FASTA file is roughly 60–70 GB (gzipped), decompressing to several hundred gigabytes. Performance and memory efficiency are first-order design concerns.

The NCBI taxonomy hierarchy (`nodes.dmp`) contains approximately 2.5 million nodes. This is small enough to fit entirely in memory as a flat dictionary.

---

## 3. Input Formats

### 3.1 UniProt FASTA

The input is a standard UniProt FASTA file (Swiss-Prot, TrEMBL, or a concatenation of both). Each entry consists of a header line beginning with `>` followed by one or more sequence lines. The header format is:

```
>db|UniqueIdentifier|EntryName ProteinName OS=OrganismName OX=OrganismIdentifier [GN=GeneName] PE=ProteinExistence SV=SequenceVersion
```

Where `db` is `sp` (Swiss-Prot) or `tr` (TrEMBL). The field `OX=` contains the NCBI taxonomy ID as an integer. This is the field used for filtering.

The input file may be:
- An uncompressed `.fasta` or `.fa` file
- A gzip-compressed `.fasta.gz` or `.fa.gz` file

Detection should be based on magic bytes (`\x1f\x8b`), not file extension, to handle edge cases.

### 3.2 NCBI Taxonomy Dump Files

The program uses the standard NCBI `taxdump.tar.gz` archive, downloadable from `https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump.tar.gz`. The relevant files within this archive are:

- **`nodes.dmp`**: Maps each taxonomy ID to its parent taxonomy ID. Columns are delimited by `\t|\t`. Only the first two columns are needed: `tax_id` and `parent_tax_id`. The root node (taxid `1`) is its own parent.

- **`merged.dmp`**: Maps old (deprecated/merged) taxonomy IDs to their current replacements. This is critical because UniProt FASTA files may contain taxonomy IDs that have been merged since the UniProt release was generated.

- **`names.dmp`** (optional, for human-readable output): Maps taxonomy IDs to scientific names. Only rows where the name class is `scientific name` are relevant.

- **`delnodes.dmp`** (optional, for diagnostics): Lists taxonomy IDs that have been deleted entirely.

---

## 4. Core Algorithm

### 4.1 Pre-computation: Build the Allowed TaxID Set

This is the performance-critical insight. Rather than traversing the taxonomy tree for every FASTA entry, the program should pre-compute a flat `set[int]` of all taxonomy IDs that pass the filter. This reduces the per-entry filtering operation to an O(1) set membership test.

**Steps:**

1. **Parse `nodes.dmp`** into a dictionary: `parent_of: dict[int, int]` mapping each `tax_id` to its `parent_tax_id`. This requires approximately 2.5 million entries and will consume roughly 100–200 MB of memory.

2. **Parse `merged.dmp`** into a dictionary: `merged_to: dict[int, int]` mapping old taxonomy IDs to their current replacements.

3. **Build a children index**: Invert `parent_of` to create `children_of: dict[int, set[int]]` mapping each taxonomy ID to its immediate children. This is a single O(n) pass.

4. **Collect all descendant IDs for each user-supplied root taxid**: Starting from each user-supplied taxonomy ID, perform a breadth-first (or depth-first) traversal of `children_of` to collect every descendant. Union the results across all user-supplied root IDs. The result is a single `set[int]` called `allowed_taxids`.

5. **Expand the set with merged IDs**: For every entry in `merged_to`, if the *new* (target) taxonomy ID is in `allowed_taxids`, also add the *old* (source) taxonomy ID to the set. This ensures that FASTA entries annotated with deprecated taxonomy IDs are still correctly matched.

The resulting `allowed_taxids` set may contain millions of entries (e.g., Bacteria alone has over a million descendant taxids), but Python `set` membership testing remains O(1) and the memory footprint is manageable (tens of MB for a set of integers).

### 4.2 Streaming Filter of FASTA

The FASTA file is streamed line-by-line (never loaded entirely into memory). The algorithm is:

1. Open the input file. If gzipped, open via Python's `gzip.open()` (or `isal.igzip.open()`).
2. Read line-by-line. When a header line (starting with `>`) is encountered:
   a. Extract the `OX=` value from the header using a compiled regular expression: `re.compile(r'OX=(\d+)')`.
   b. If extraction fails (no match), record a warning (see §5.5), increment the unparseable counter, and set `include_current_entry = False`. Entries without valid taxonomy annotation are never included.
   c. If extraction succeeds, convert to `int` and test membership in `allowed_taxids`.
   d. If the taxid is not in `allowed_taxids` and also not in the full taxonomy tree, record a per-ID warning (see §5.5) and increment the unknown-taxid counter.
   e. Set `include_current_entry` accordingly.
3. For sequence lines (non-header), write them to output only if `include_current_entry` is `True`.
4. Write matching header lines and their sequence lines to the output stream.

### 4.3 Buffered I/O

For maximum throughput on files of this size:
- Use a large read buffer (e.g., 8–16 MB) when opening the input file.
- Use a large write buffer on the output stream.
- When reading gzipped input, `gzip.open()` with an explicit buffer size should be used.
- Consider supporting `pigz` / `igzip` (via the `isal` Python package) as a faster alternative to Python's built-in `gzip` for decompression, as `isal.igzip` can be 5–10x faster than `gzip`.

---

## 5. Command-Line Interface

The tool should use Python's `argparse` (or `click`) for CLI parsing. The command name is `taxafasta`.

### 5.1 Arguments

| Argument | Short | Required | Description |
|---|---|---|---|
| `--input` | `-i` | Yes | Path to input UniProt FASTA file (plain or gzipped) |
| `--taxid` | `-t` | Yes | One or more NCBI taxonomy IDs to include (space-separated or repeated). All descendants are automatically included. |
| `--output` | `-o` | Yes | Path to output FASTA file. Output is gzip-compressed by default (see `--no-gzip`). If path does not end in `.gz`, the suffix `.gz` is appended automatically unless `--no-gzip` is set. |
| `--no-gzip` | | No | Disable gzip compression of output. Write uncompressed FASTA. |
| `--taxdump` | `-d` | No | Path to an already-extracted taxdump directory containing `nodes.dmp`, `merged.dmp`, etc. If omitted, the program downloads and caches `taxdump.tar.gz`. |
| `--cache-dir` | | No | Directory for caching downloaded taxonomy files. Defaults to `~/.taxafasta/`. |
| `--exclude` | `-e` | No | One or more NCBI taxonomy IDs to exclude (with descendants), applied after inclusion. |
| `--no-merge` | | No | Flag to disable merged taxonomy ID resolution. |
| `--verbose` | `-v` | No | Periodic progress updates to stderr (see §5.3 Progress Reporting). |
| `--version` | | No | Print version and exit. |

### 5.2 Example Usage

```bash
# Filter to bacteria only (output is gzipped by default)
taxafasta -i uniprot_trembl.fasta.gz -t 2 -o bacteria.fasta
# produces bacteria.fasta.gz and bacteria.fasta.log

# Filter to bacteria and viruses
taxafasta -i uniprot_trembl.fasta.gz -t 2 10239 -o bacteria_and_viruses.fasta

# Uncompressed output
taxafasta -i uniprot_trembl.fasta.gz -t 9606 -o human.fasta --no-gzip

# Use pre-downloaded taxonomy files
taxafasta -i uniprot_trembl.fasta.gz -t 2 --taxdump /data/taxdump/ -o output.fasta

# Include eukaryotes but exclude mammals
taxafasta -i uniprot_trembl.fasta.gz -t 2759 -e 40674 -o euk_no_mammals.fasta

# Verbose progress reporting
taxafasta -i uniprot_trembl.fasta.gz -t 2 -o bacteria.fasta -v
```

### 5.3 Progress Reporting

When `--verbose` is set, the tool prints periodic progress updates to stderr. Rather than attempting a progress bar (which is unreliable for gzipped input where total uncompressed size is unknown), the tool emits a status line every N entries (e.g., every 1 million entries). Each status line includes:

- Total entries processed so far
- Entries included so far
- Entries excluded so far
- Elapsed wall-clock time
- Approximate processing rate (entries/second)

Example:

```
[00:01:32] Processed 10,000,000 entries | Included: 4,231,017 | Excluded: 5,768,983 | Rate: ~108,696 entries/sec
[00:03:05] Processed 20,000,000 entries | Included: 8,490,112 | Excluded: 11,509,888 | Rate: ~107,527 entries/sec
```

These lines are printed to stderr so they do not interfere with FASTA output in any scenario.

### 5.4 Log File

Every run produces a log file alongside the output file. The log file name is derived from the output file path with a `.log` extension. If the output file is `bacteria.fasta` (or `bacteria.fasta.gz` after auto-compression), the log file is `bacteria.fasta.log`. If a log file with that name already exists, a numeric suffix is appended: `.log1`, `.log2`, etc.

The log file contains:

1. **Run metadata**: taxafasta version, Python version, platform, timestamp (UTC ISO 8601).
2. **Command line**: The exact command-line invocation, including all arguments.
3. **Parameters**: All resolved parameters (input path, output path, included taxids, excluded taxids, taxdump path, whether merged ID resolution was enabled, etc.).
4. **Taxonomy version**: The source of the taxonomy data — either the path to the user-supplied taxdump directory, or the download URL and the modification date of the downloaded `taxdump.tar.gz`. If `names.dmp` is available, the log also records the human-readable names of the user-supplied taxonomy IDs (e.g., `2 = Bacteria`, `10239 = Viruses`).
5. **Allowed set size**: The total number of taxonomy IDs in the allowed set after inclusion, exclusion, and merged ID expansion.
6. **Warnings**: All warnings generated during the run (see §5.5).
7. **Summary statistics**: Total entries processed, entries included, entries excluded, entries skipped due to missing/malformed OX field, entries with unrecognized taxonomy IDs, elapsed time, processing rate.

Example log file:

```
taxafasta v1.0.0 | Python 3.12.1 | Linux x86_64 | 2025-03-05T14:22:07Z
Command: taxafasta -i uniprot_trembl.fasta.gz -t 2 10239 -o bacteria_viruses.fasta -v
Input: uniprot_trembl.fasta.gz
Output: bacteria_viruses.fasta.gz (gzip enabled)
Included taxids: 2 (Bacteria), 10239 (Viruses)
Excluded taxids: (none)
Merged ID resolution: enabled
Taxdump source: /home/user/.taxafasta/taxdump/ (downloaded 2025-03-01)
Allowed taxid set size: 1,847,293

WARNINGS (3 unique):
  [x1] Header could not be parsed for OX field: >gnl|custom|seq1 hypothetical protein
  [x14] Taxonomy ID 99999999 not found in NCBI taxonomy data
  [x2] Taxonomy ID 88888888 not found in NCBI taxonomy data

SUMMARY:
  Total entries processed:  251,600,768
  Entries included:         142,309,441
  Entries excluded:         109,291,310
  Entries skipped (no OX):            1
  Entries with unknown taxid:        16
  Elapsed time:             00:47:12
  Processing rate:          ~88,832 entries/sec
```

### 5.5 Warning Behavior

Warnings are printed to stderr during the run and also recorded in the log file. Warning deduplication rules:

- **Unparseable header** (no valid `OX=` field): The first occurrence prints the full header line to stderr. Subsequent occurrences of unparseable headers are counted silently and summarized in the log. These entries are always excluded from output.
- **Unknown taxonomy ID** (valid `OX=` integer, but not found in `nodes.dmp` or `merged.dmp`): A warning is printed to stderr the first time a given unknown taxonomy ID is encountered. Subsequent entries with the same unknown ID are counted silently. All counts are summarized per-ID in the log. These entries are always excluded from output.

This approach avoids flooding stderr with millions of repeated warnings while still surfacing every distinct problem.

---

## 6. Taxonomy Download & Caching

When `--taxdump` is not provided:

1. Check `--cache-dir` (default `~/.taxafasta/`) for an existing `taxdump.tar.gz` and its extraction.
2. If not found, download from `https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump.tar.gz` using Python's `urllib` or the `requests` library.
3. Extract `nodes.dmp`, `merged.dmp`, and optionally `names.dmp` into the cache directory.
4. Store a timestamp or checksum file alongside the extracted data so users can verify freshness.
5. Print a message to stderr noting the download.

The program should not re-download on every run. A `--update-taxonomy` flag could force a fresh download.

---

## 7. Error Handling

The program must handle the following error conditions gracefully, with clear error messages to stderr and appropriate nonzero exit codes:

- Input file does not exist or is not readable.
- Input file is neither valid FASTA nor valid gzip.
- A user-supplied taxonomy ID is not found in `nodes.dmp` or `merged.dmp` (fatal error — the user made a mistake in their request).
- A FASTA entry lacks an `OX=` field or has a malformed one (warning to stderr on first occurrence, counted, always excluded from output). See §5.5.
- A FASTA entry has a valid `OX=` integer that does not appear in the taxonomy data (warning to stderr on first occurrence of each unique ID, counted, always excluded from output). See §5.5.
- `taxdump.tar.gz` download fails (network error, HTTP error, corrupt archive).
- `nodes.dmp` or `merged.dmp` is malformed or missing from the taxdump directory.
- Output file cannot be created or written to.
- Disk space or memory exhaustion during processing.

Exit codes:
- `0`: Success
- `1`: Fatal error (missing input, invalid arguments)
- `2`: Partial success (some entries had warnings, e.g., missing OX fields)

---

## 8. Performance Targets

Given the scale of UniProt TrEMBL (250M+ entries, hundreds of GB uncompressed):

- **Taxonomy loading**: The entire `nodes.dmp` and `merged.dmp` parse and descendant-set construction should complete in under 30 seconds on commodity hardware.
- **FASTA streaming throughput**: The tool should be I/O-bound, not CPU-bound. With `isal`-accelerated gzip decompression, throughput should approach the raw disk/network read speed. A target of 200+ MB/s of uncompressed data processed is reasonable.
- **Memory usage**: The taxonomy data structures and the allowed-set together should consume under 1 GB of RAM. The FASTA stream itself adds negligible memory overhead (single line buffer).

---

## 9. Project Structure

See §16 for the complete project structure including CI/CD workflow files.

---

## 10. Dependencies

### Runtime

| Package | Purpose | Required? |
|---|---|---|
| Python ≥ 3.10 | Language runtime | Yes |
| `isal` (python-isal) | Fast gzip decompression/compression (drop-in replacement for `gzip`) | Recommended (falls back to stdlib `gzip`) |
| `requests` | HTTP download of taxdump (only if auto-download used) | Optional |

The tool should have minimal required dependencies. The `isal` and `requests` packages should be optional extras installable via `pip install taxafasta[fast]` and `pip install taxafasta[download]` respectively, with graceful fallback when absent.

### Development / Testing

| Package | Purpose |
|---|---|
| `pytest` | Test framework |
| `pytest-cov` | Coverage reporting |
| `pytest-benchmark` | Performance regression tests |
| `mypy` | Static type checking |
| `ruff` | Linting and formatting |

---

## 11. Testing Strategy

Testing is paramount. Every module must have thorough unit tests. Accuracy and fidelity are the top priorities.

### 11.1 Taxonomy Module Tests (`test_taxonomy.py`)

- **Parsing `nodes.dmp`**: Verify correct extraction of `tax_id` → `parent_tax_id` from fixture data that replicates the real `\t|\t` delimited format. Test with edge cases: the root node (taxid 1 is its own parent), single-child chains, wide branching.
- **Parsing `merged.dmp`**: Verify old→new mapping. Test with entries where the target has itself been merged (chain resolution).
- **Children index construction**: Given a known tree, verify that `children_of` is correctly inverted from `parent_of`.
- **Descendant set construction**: For a fixture taxonomy tree, verify that requesting a subtree root returns exactly the expected set of descendants. Test with multiple roots, overlapping subtrees (should union correctly), leaf nodes (should return just themselves), and the root node (should return everything).
- **Merged ID expansion**: Verify that old taxonomy IDs mapping into the allowed set are included. Verify that old IDs mapping outside the set are excluded.
- **Exclude logic**: Verify that exclusion taxids and their descendants are correctly removed from the allowed set, including when inclusion and exclusion overlap.
- **Unknown taxid handling**: Verify appropriate warning when a user-supplied taxid does not exist in the tree.

### 11.2 FASTA Module Tests (`test_fasta.py`)

- **OX field extraction**: Test the regex against a variety of real-world UniProt header formats — Swiss-Prot entries, TrEMBL entries, entries with and without `GN=` fields, entries where `OX=` appears at different positions in the header. Verify correct integer extraction.
- **Missing OX field**: Test behavior when a header line lacks `OX=`. Verify the entry is excluded, a warning is generated on the first occurrence, and subsequent occurrences are counted silently.
- **Malformed OX field**: Test `OX=abc`, `OX=`, `OX=123abc` — verify graceful handling, entry excluded, warning generated.
- **Unknown taxonomy ID**: Test an entry with a valid `OX=` integer that does not exist in the taxonomy data. Verify the entry is excluded, a warning is emitted for the first occurrence of that ID, and subsequent entries with the same unknown ID produce no additional stderr output (only increment the count).
- **Multi-line sequences**: Verify that all sequence lines following a matching header are included, and all sequence lines following a non-matching header are excluded.
- **Entry boundary correctness**: Use a fixture FASTA with interleaved matching and non-matching entries. Verify the output contains exactly the correct entries with no cross-contamination of sequence lines between entries.
- **Empty sequence**: Test entries with a header but no sequence lines.
- **Large headers**: Test with very long header lines (thousands of characters).
- **Output fidelity**: Verify that the output FASTA is byte-for-byte identical to what a manual extraction would produce — headers and sequences are not modified, line endings are preserved, no trailing whitespace is added or removed.

### 11.3 I/O Module Tests (`test_io_utils.py`)

- **Gzip detection**: Test that gzip files are correctly identified by magic bytes regardless of file extension. Test that plain text files are not misidentified.
- **Transparent opening**: Verify that the same API works for both gzip and plain files, returning equivalent content.
- **isal fallback**: Verify that when `isal` is not installed, the module falls back to `gzip` without error.
- **Gzip output default**: Verify that when no `--no-gzip` flag is set, the output is valid gzip and `.gz` is appended to the output path.
- **No-gzip output**: Verify that `--no-gzip` produces a plain text file without `.gz` suffix appended.
- **Buffer sizes**: Verify that custom buffer sizes are applied (may require mocking).

### 11.4 Download Module Tests (`test_download.py`)

- **Successful download** (mocked): Mock HTTP responses to verify correct download, extraction, and caching behavior.
- **Cache hit**: Verify that a second call with existing cache does not re-download.
- **Network failure**: Mock connection errors and verify clean error handling with appropriate messages.
- **Corrupt archive**: Provide a truncated/corrupt tar.gz fixture and verify error handling.
- **File permissions**: Verify the cache directory is created with appropriate permissions.

### 11.5 Log Module Tests (`test_run_log.py`)

- **Log file naming**: Verify the log file is created alongside the output file with `.log` extension. Verify that if `bacteria.fasta.log` exists, the next run produces `bacteria.fasta.log1`, then `.log2`, etc.
- **Log file naming with gzip**: Verify that when the output is `bacteria.fasta.gz`, the log file is `bacteria.fasta.log` (not `bacteria.fasta.gz.log`).
- **Run metadata**: Verify the log contains the taxafasta version, Python version, platform, and a valid UTC ISO 8601 timestamp.
- **Command line recording**: Verify the log contains the exact command-line string used to invoke the tool.
- **Parameter recording**: Verify all resolved parameters appear in the log (input path, output path, taxids, excluded taxids, taxdump source, merge setting).
- **Taxonomy version recording**: Verify the taxdump source path or download URL and date appear in the log. Verify human-readable taxid names are included when `names.dmp` is available.
- **Warning accumulation and deduplication**: Verify that warnings are accumulated with counts. Verify that the log records each unique warning with its occurrence count (e.g., `[x14] Taxonomy ID 99999999 not found`).
- **Summary statistics**: Verify the log contains correct counts for total, included, excluded, skipped (no OX), unknown taxid entries, elapsed time, and processing rate.
- **Empty warnings**: Verify the log omits the warnings section entirely when no warnings were generated.

### 11.6 CLI Tests (`test_cli.py`)

- **Argument parsing**: Test all argument combinations — required args missing, invalid values, conflicting options.
- **Help output**: Verify `--help` produces expected output.
- **Version output**: Verify `--version` prints the version string.
- **Exit codes**: Verify correct exit codes for success, failure, and partial-success scenarios.

### 11.7 Integration Tests (`test_integration.py`)

These are end-to-end tests that invoke the CLI on fixture files and verify output.

- **Basic filter**: A small FASTA with entries from bacteria (taxid 2 descendants), archaea, and eukaryotes. Filter with `-t 2`, verify output contains only bacterial entries.
- **Multiple taxids**: Filter with `-t 2 10239`, verify output contains bacteria and viruses but nothing else.
- **Exclude**: Filter with `-t 2759 -e 40674`, verify eukaryotes are included but mammals are excluded.
- **Merged taxid handling**: A FASTA entry annotated with an old (merged) taxonomy ID that maps into the included subtree. Verify it is included.
- **Gzip input**: Same test as basic filter but with gzip-compressed input. Verify identical output.
- **Default gzip output**: Verify that output is gzip-compressed by default. Verify `.gz` suffix is appended to the output path if not already present.
- **No-gzip output**: Verify `--no-gzip` produces uncompressed output without `.gz` suffix.
- **Log file creation**: Verify a `.log` file is created alongside every output file. Verify its contents include all required sections (metadata, command, parameters, taxonomy version, summary).
- **Log file collision**: Run twice with the same output path. Verify the second run creates `.log1` (not overwriting `.log`).
- **Warning deduplication**: A FASTA with multiple entries having the same unknown taxonomy ID. Verify stderr shows the warning only once, and the log file records the count.
- **Unparseable header warning**: A FASTA with an entry missing `OX=`. Verify stderr shows the warning once with the full header, the entry is excluded, and the log records the count.
- **No matches**: Verify behavior when no entries match the filter (empty output, exit code 0, log reflects this).
- **All match**: Verify behavior when every entry matches.

### 11.8 Accuracy Tests

A dedicated set of tests using carefully constructed fixture data that exercises every branch of the algorithm:

- An entry whose taxid is exactly the requested root taxid (should match).
- An entry whose taxid is a leaf descendant many levels deep (should match).
- An entry whose taxid is a sibling of the requested root (should NOT match).
- An entry whose taxid is the parent of the requested root (should NOT match).
- An entry whose taxid is a merged ID pointing into the subtree (should match).
- An entry whose taxid is a merged ID pointing outside the subtree (should NOT match).
- An entry whose taxid is in `delnodes.dmp` (should NOT match, warn).
- An entry with no OX field (should be excluded, warning generated).

### 11.9 Performance / Regression Tests

Using `pytest-benchmark`:

- Benchmark taxonomy loading and descendant-set construction with a realistically sized fixture (or the actual `nodes.dmp` if available in CI).
- Benchmark FASTA throughput on a generated fixture of 100K entries.
- Track these benchmarks over time to detect performance regressions.

---

## 12. Naming Conventions

Consistent naming across all distribution channels:

| Context | Name | Notes |
|---|---|---|
| PyPI package | `taxafasta` | Lowercase, no hyphens or underscores. Installed via `pip install taxafasta`. |
| Python import | `taxafasta` | `import taxafasta` or `from taxafasta import ...` |
| CLI command | `taxafasta` | Console script entry point. |
| GitHub repository | `taxafasta` | e.g., `github.com/<org>/taxafasta` |
| Docker image (GHCR) | `ghcr.io/<org>/taxafasta` | Lowercase, matches the GitHub repository name. |
| Docker image tags | `<version>`, `<major>.<minor>`, `<major>`, `latest` | Semantic version tags derived from the GitHub release tag (e.g., release `v1.2.3` produces image tags `1.2.3`, `1.2`, `1`, and `latest`). |

The `<org>` placeholder refers to the GitHub user or organization that owns the repository.

---

## 13. Distribution

### 13.1 PyPI (pip)

The project uses `pyproject.toml` with a PEP 621 compliant build configuration (e.g., `hatchling` or `setuptools` as the build backend). A console script entry point maps `taxafasta` to the CLI:

```toml
[project]
name = "taxafasta"
dynamic = ["version"]
description = "Filter UniProt FASTA files by NCBI taxonomy"
readme = "README.md"
license = {text = "MIT"}
requires-python = ">=3.10"
classifiers = [
    "Development Status :: 4 - Beta",
    "Environment :: Console",
    "Intended Audience :: Science/Research",
    "Topic :: Scientific/Engineering :: Bio-Informatics",
    "License :: OSI Approved :: Apache License",
    "Programming Language :: Python :: 3",
]

[project.scripts]
taxafasta = "taxafasta.cli:main"

[project.optional-dependencies]
fast = ["isal>=1.0"]
download = ["requests>=2.28"]
all = ["isal>=1.0", "requests>=2.28"]
dev = ["pytest", "pytest-cov", "pytest-benchmark", "mypy", "ruff"]
```

### 13.2 Version Management

The package version is the single source of truth, defined in `src/taxafasta/__init__.py` as `__version__`. The `pyproject.toml` reads it dynamically. Git release tags must follow the format `v<major>.<minor>.<patch>` (e.g., `v1.2.3`). The CI workflows derive all Docker image tags and PyPI version strings from this Git tag.

### 13.3 Docker

A multi-stage Dockerfile:

1. **Build stage**: Based on `python:3.12-slim`. Install build dependencies, copy source, build the wheel.
2. **Test stage**: Install the wheel and dev dependencies, run the full test suite. This stage is used in CI but not shipped.
3. **Runtime stage**: Based on `python:3.12-slim`. Install only the wheel and runtime dependencies (`isal`, `requests`). The entrypoint is the `taxafasta` command.

The Dockerfile should accept a `VERSION` build argument so the image can be labeled with the version at build time:

```dockerfile
ARG VERSION=dev
LABEL org.opencontainers.image.version="${VERSION}"
LABEL org.opencontainers.image.source="https://github.com/<org>/taxafasta"
LABEL org.opencontainers.image.description="Filter UniProt FASTA files by NCBI taxonomy"
```

Example usage:

```bash
# Pull a specific version
docker pull ghcr.io/<org>/taxafasta:1.2.3

# Run it
docker run --rm -v /data:/data ghcr.io/<org>/taxafasta:1.2.3 \
  -i /data/uniprot_trembl.fasta.gz -t 2 -o /data/bacteria.fasta
# produces /data/bacteria.fasta.gz and /data/bacteria.fasta.log

# Use the latest release
docker run --rm -v /data:/data ghcr.io/<org>/taxafasta:latest \
  -i /data/uniprot_trembl.fasta.gz -t 2 10239 -o /data/bacteria_viruses.fasta
```

---

## 14. GitHub Actions CI/CD

The project includes two GitHub Actions workflow files in `.github/workflows/`.

### 14.1 `ci.yml` — Continuous Integration (runs on every push and PR)

**Trigger**: Every push to any branch and every pull request.

**Jobs**:

1. **`lint`**: Runs `ruff check` and `ruff format --check` on the codebase. Fails the build on any lint or formatting violation.

2. **`typecheck`**: Runs `mypy src/` with strict mode. Fails on type errors.

3. **`test`**: Runs the full test suite across a matrix of Python versions (3.10, 3.11, 3.12). Steps:
   - Check out the repository.
   - Set up the Python version.
   - Install the package in editable mode with dev dependencies: `pip install -e ".[all,dev]"`.
   - Run `pytest --cov=taxafasta --cov-report=xml --cov-report=term tests/`.
   - Upload the coverage report as an artifact (and optionally to a coverage service).
   - Fail the build if coverage drops below a configured threshold (e.g., 90%).

All three jobs run in parallel. The full workflow should complete in under 5 minutes for the unit/integration test suite.

### 14.2 `release.yml` — Release & Publish (runs on GitHub Release creation)

**Trigger**: `on: release: types: [published]`. This fires when a maintainer creates a GitHub Release with a tag matching `v*.*.*` (e.g., `v1.2.3`).

**Jobs**:

1. **`test`**: Identical to the `ci.yml` test job. The full test suite runs again against the exact release commit to ensure nothing was missed. The remaining jobs depend on this one passing.

2. **`publish-pypi`**: Builds the sdist and wheel, then publishes to PyPI using `pypa/gh-action-pypi-publish`. Steps:
   - Check out the repository.
   - Set up Python.
   - Install build tools: `pip install build`.
   - Build: `python -m build`.
   - Publish using the trusted publisher mechanism (OIDC, no API token needed if configured in PyPI).

3. **`publish-docker`**: Builds and pushes the Docker image to GitHub Container Registry (GHCR). Steps:
   - Check out the repository.
   - Set up Docker Buildx (`docker/setup-buildx-action@v3`).
   - Log in to GHCR (`docker/login-action@v3` with `registry: ghcr.io`, using `github.actor` and `secrets.GITHUB_TOKEN`).
   - Extract metadata and generate tags (`docker/metadata-action@v5`):
     ```yaml
     images: ghcr.io/${{ github.repository }}
     tags: |
       type=semver,pattern={{version}}
       type=semver,pattern={{major}}.{{minor}}
       type=semver,pattern={{major}}
       type=raw,value=latest
     ```
     For a release tagged `v1.2.3`, this produces four image tags: `1.2.3`, `1.2`, `1`, and `latest`.
   - Build and push (`docker/build-push-action@v6`):
     ```yaml
     context: .
     push: true
     tags: ${{ steps.meta.outputs.tags }}
     labels: ${{ steps.meta.outputs.labels }}
     build-args: VERSION=${{ github.ref_name }}
     ```

**Required repository settings**:
- The repository's `GITHUB_TOKEN` must have `packages: write` permission (set in the workflow's `permissions` block).
- For PyPI trusted publishing, the repository must be registered as a trusted publisher on pypi.org under the `taxafasta` project.

### 14.3 Workflow File Locations

```
.github/
└── workflows/
    ├── ci.yml        # Lint + typecheck + test on push/PR
    └── release.yml   # Test + publish to PyPI + publish Docker image on release
```

### 14.4 CI Tests for CI Configuration (`test_ci.py`, optional)

While not strictly unit tests, the project may include a smoke test that verifies:
- The `Dockerfile` builds successfully (can be run locally via `docker build --target test .`).
- The built image's entrypoint works (`docker run --rm <image> --version` prints the version).

---

## 15. README

The repository must include a `README.md` at the root. It serves as the primary documentation and is also displayed on PyPI and the GitHub repository page. The README should contain the following sections:

### 15.1 Header & Badges

A one-line project name and tagline, followed by CI/CD status badges:

```
# taxafasta

Filter UniProt protein FASTA files by NCBI taxonomy.

[![CI](https://github.com/<org>/taxafasta/actions/workflows/ci.yml/badge.svg)](...)
[![PyPI version](https://img.shields.io/pypi/v/taxafasta)](https://pypi.org/project/taxafasta/)
[![Docker](https://img.shields.io/badge/ghcr.io-taxafasta-blue)](https://ghcr.io/<org>/taxafasta)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)
```

### 15.2 Overview

A concise (3–5 sentence) description of what the tool does and why it exists. Key points to cover:
- Filters large UniProt FASTA files (Swiss-Prot and/or TrEMBL) to include only proteins from specified NCBI taxonomy subtrees.
- Designed for files at the scale of full UniProt TrEMBL (250M+ entries).
- Uses the NCBI taxonomy hierarchy to include all descendants of specified taxonomy IDs.
- Handles merged/deprecated taxonomy IDs automatically.

### 15.3 How It Works

A brief (1–2 paragraph) description of the architecture, aimed at someone who wants to understand the tool's approach without reading the full specification. Cover:
- The tool parses NCBI taxonomy dump files (`nodes.dmp`, `merged.dmp`) to build a parent-child tree in memory.
- From user-supplied taxonomy IDs, it pre-computes a flat set of all allowed taxonomy IDs (the specified IDs plus all their descendants). This makes filtering a single O(1) set-membership check per FASTA entry.
- The FASTA file is streamed line-by-line and never loaded into memory. Each entry's `OX=` field is extracted and checked against the pre-computed set. Matching entries are written to the gzip-compressed output.
- A log file is generated for every run recording parameters, taxonomy version, warnings, and summary statistics.

### 15.4 Installation

Show both pip and Docker installation methods:

```
## Installation

### pip

```bash
pip install taxafasta

# With recommended performance dependencies:
pip install taxafasta[all]
```

### Docker

```bash
docker pull ghcr.io/<org>/taxafasta:latest
```
```

### 15.5 Quick Start

A minimal working example that someone can copy-paste:

```
## Quick Start

# Download NCBI taxonomy (automatic on first run, or provide manually)
taxafasta -i uniprot_trembl.fasta.gz -t 2 -o bacteria.fasta

# This produces:
#   bacteria.fasta.gz   — gzip-compressed FASTA with only bacterial proteins
#   bacteria.fasta.log  — run log with parameters, warnings, and statistics
```

### 15.6 Usage Examples

A more comprehensive set of examples showing common use cases:

```
## Usage

### Filter to a single taxonomic group
taxafasta -i uniprot_trembl.fasta.gz -t 2 -o bacteria.fasta

### Filter to multiple groups (bacteria + viruses)
taxafasta -i uniprot_trembl.fasta.gz -t 2 10239 -o bacteria_viruses.fasta

### Exclude a subtree (eukaryotes minus mammals)
taxafasta -i uniprot_trembl.fasta.gz -t 2759 -e 40674 -o euk_no_mammals.fasta

### Use pre-downloaded taxonomy files
taxafasta -i uniprot_trembl.fasta.gz -t 2 --taxdump /path/to/taxdump/ -o bacteria.fasta

### Uncompressed output
taxafasta -i uniprot_trembl.fasta.gz -t 9606 -o human.fasta --no-gzip

### Verbose progress
taxafasta -i uniprot_trembl.fasta.gz -t 2 -o bacteria.fasta -v
```

### 15.7 Docker Usage

Show equivalent Docker examples with volume mounting:

```
## Docker Usage

docker run --rm -v /data:/data ghcr.io/<org>/taxafasta:1.2.3 \
  -i /data/uniprot_trembl.fasta.gz -t 2 -o /data/bacteria.fasta

# With pre-downloaded taxonomy
docker run --rm \
  -v /data:/data \
  -v /taxonomy:/taxonomy:ro \
  ghcr.io/<org>/taxafasta:latest \
  -i /data/uniprot_trembl.fasta.gz -t 2 --taxdump /taxonomy -o /data/bacteria.fasta
```

### 15.8 Common Taxonomy IDs

A reference table of commonly used taxonomy IDs for convenience:

| Taxonomy ID | Name |
|---|---|
| 2 | Bacteria |
| 2157 | Archaea |
| 2759 | Eukaryota |
| 10239 | Viruses |
| 9606 | Homo sapiens |
| 7742 | Vertebrata |
| 40674 | Mammalia |
| 33208 | Metazoa |
| 3193 | Embryophyta (land plants) |
| 4751 | Fungi |

### 15.9 NCBI Taxonomy Data

Brief explanation of how the tool obtains taxonomy data:
- By default, the tool automatically downloads and caches `taxdump.tar.gz` from NCBI's FTP server on first run.
- Users can supply pre-downloaded taxonomy files with `--taxdump`.
- Link to the NCBI taxonomy FTP: `https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/`

### 15.10 Development

Brief instructions for contributors:

```
## Development

git clone https://github.com/<org>/taxafasta.git
cd taxafasta
pip install -e ".[all,dev]"

# Run tests
pytest

# Run with coverage
pytest --cov=taxafasta

# Lint and format
ruff check src/ tests/
ruff format src/ tests/

# Type check
mypy src/
```

### 15.11 License

State the license (MIT recommended for a bioinformatics utility).

---

## 16. Project Structure (updated)

```
taxafasta/
├── .github/
│   └── workflows/
│       ├── ci.yml                # Lint + typecheck + test on push/PR
│       └── release.yml           # Test + publish PyPI + publish Docker on release
├── pyproject.toml                # Project metadata, build config, dependencies
├── Dockerfile                    # Multi-stage build (build → test → runtime)
├── README.md                     # Project documentation (see §15)
├── LICENSE                       # MIT license
├── src/
│   └── taxafasta/
│       ├── __init__.py           # Package init, __version__
│       ├── __main__.py           # Entry point: `python -m taxafasta`
│       ├── cli.py                # Argument parsing and main orchestration
│       ├── taxonomy.py           # Parsing nodes.dmp, merged.dmp; building descendant sets
│       ├── fasta.py              # Streaming FASTA reader/writer with OX extraction
│       ├── download.py           # Taxonomy download and caching logic
│       ├── io_utils.py           # Smart file opening (gzip detection, buffering, isal)
│       └── run_log.py            # Log file generation and warning accumulation
└── tests/
    ├── conftest.py               # Shared fixtures (small taxonomy trees, sample FASTA)
    ├── test_taxonomy.py          # Unit tests for taxonomy module
    ├── test_fasta.py             # Unit tests for FASTA parsing and filtering
    ├── test_download.py          # Unit tests for download/caching (mocked network)
    ├── test_io_utils.py          # Unit tests for file detection and opening
    ├── test_run_log.py           # Unit tests for log file generation and warning dedup
    ├── test_cli.py               # Integration tests for CLI argument handling
    ├── test_integration.py       # End-to-end tests with real-format fixture files
    └── data/                     # Small fixture files for testing
        ├── tiny_nodes.dmp
        ├── tiny_merged.dmp
        ├── tiny_names.dmp
        ├── sample.fasta
        └── sample.fasta.gz
```

---

## 17. Design Decisions & Rationale

### Why a pre-computed set instead of per-entry tree traversal?

Per-entry tree traversal would require, for each of the 250M+ FASTA entries, walking up the taxonomy tree from the entry's taxid to the root to check if any ancestor matches a requested taxid. Even with the tree in memory, this would be O(depth) per entry (tree depth can be 30+ levels), resulting in billions of dictionary lookups.

The pre-computed set approach does O(n) work once (where n is the number of nodes in the taxonomy, ~2.5M) and then O(1) work per FASTA entry. For 250M entries, this is dramatically faster. The set construction takes seconds; the savings during filtering are enormous.

### Why not use BioPython's SeqIO?

BioPython's `SeqIO.parse()` is convenient but creates full `SeqRecord` objects for every entry, involving string parsing, object allocation, and attribute assignment that is entirely unnecessary for this use case. We only need to read the header line, extract the OX value, and pass through the raw text. A line-by-line streaming approach avoids all overhead and is substantially faster at this scale.

### Why support merged taxonomy IDs?

UniProt releases and NCBI taxonomy releases are not synchronized. A UniProt FASTA file may contain taxonomy IDs that have since been merged into other IDs in a newer taxonomy dump. Without merged ID resolution, entries annotated with old IDs would be silently dropped even though they belong to the correct taxonomic group. This is a correctness concern, not just a convenience.

### Why `isal` over `gzip`?

Python's built-in `gzip` module is pure Python (wrapping zlib) and is significantly slower than the ISA-L library, which uses hardware-accelerated compression/decompression. On large files like TrEMBL, this can mean the difference between hours and tens of minutes. The `python-isal` package provides a drop-in replacement API.

---

## 18. Future Considerations

These are explicitly out of scope for v1.0 but may be considered later:

- **Taxonomy name lookup**: Accept organism names (e.g., "Bacteria") in addition to numeric taxids, resolved via `names.dmp`.
- **Parallel processing**: Shard the input file and process chunks in parallel (complex for gzipped input without seekable index).
- **FASTA indexing**: Create an index of taxid→byte-offset for repeated queries on the same file.
- **UniRef support**: Handle UniRef FASTA headers, which have a different format.
- **Reporting**: Generate a report of taxonomy distribution in the input and output files.
- **Streaming download**: Accept a URL as input and stream the FASTA directly from UniProt's FTP server.

---

## 19. Open Questions

1. **Compressed output format**: Should the tool support writing `.zst` (Zstandard) compressed output in addition to gzip? Zstandard offers better compression ratios and speed, but gzip is more universally supported in bioinformatics toolchains.

2. **`--exclude` in v1.0**: The `--exclude` flag adds implementation and testing complexity. Should it be included in v1.0 or deferred to a later release?

## 20. Resolved Decisions

- **Entries without OX fields**: Always excluded. No configuration flag. A warning is printed on the first occurrence. (§5.5)
- **Output compression**: Gzip-compressed by default, controllable via `--no-gzip`. (§5.1)
- **Progress reporting**: Periodic entry-count status lines to stderr when `--verbose` is set. No progress bar. (§5.3)
- **Taxonomy version pinning**: Recorded in a per-run log file alongside the output. (§5.4)
- **Warning deduplication**: Each distinct warning type/ID is printed to stderr once. Counts are accumulated and summarized in the log file. (§5.5)
