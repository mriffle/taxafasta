# taxafasta

Filter UniProt protein FASTA files by NCBI taxonomy.

[![CI](https://github.com/mriffle/taxafasta/actions/workflows/ci.yml/badge.svg)](https://github.com/mriffle/taxafasta/actions/workflows/ci.yml)
[![License](https://img.shields.io/badge/License-Apache_2.0-blue.svg)](LICENSE)

## Overview

`taxafasta` filters large UniProt FASTA files (Swiss-Prot and/or TrEMBL) to include only proteins from specified NCBI taxonomy subtrees. It is designed for files at the scale of full UniProt TrEMBL (250M+ entries, hundreds of GB). The tool uses the NCBI taxonomy hierarchy to automatically include all descendants of specified taxonomy IDs and handles merged/deprecated taxonomy IDs transparently.

## How It Works

The tool parses NCBI taxonomy dump files (`nodes.dmp`, `merged.dmp`) to build a parent→child tree in memory. From user-supplied taxonomy IDs, it pre-computes a flat set of all allowed taxonomy IDs (the specified IDs plus all their descendants), reducing per-entry filtering to an O(1) set-membership check.

The FASTA file is streamed line-by-line and never loaded into memory. Each entry's `OX=` field is extracted and checked against the pre-computed set. Matching entries are written to the (gzip-compressed by default) output. A log file is generated for every run recording parameters, taxonomy version, warnings, and summary statistics.

## Requirements

- **Python 3.10 or newer** (Or Docker)

## Installation

### pip

```bash
pip install taxafasta

# With recommended performance dependencies:
pip install taxafasta[all]
```

> **Troubleshooting:** If you see an error like:
> ```
> ERROR: Could not find a version that satisfies the requirement taxafasta (from versions: none)
> ERROR: No matching distribution found for taxafasta
> ```
> Your Python version is likely too old. Verify with `python --version` — taxafasta requires Python 3.10+.

### Docker

```bash
docker pull ghcr.io/mriffle/taxafasta:latest
```

## Quick Start

```bash
# Download NCBI taxonomy (automatic on first run, or provide manually)
taxafasta -i uniprot_trembl.fasta.gz -t 2 -o bacteria.fasta

# This produces:
#   bacteria.fasta.gz   — gzip-compressed FASTA with only bacterial proteins
#   bacteria.fasta.log  — run log with parameters, warnings, and statistics
```

## Usage

### Filter to a single taxonomic group

```bash
taxafasta -i uniprot_trembl.fasta.gz -t 2 -o bacteria.fasta
```

### Filter to multiple groups (bacteria + viruses)

```bash
taxafasta -i uniprot_trembl.fasta.gz -t 2 10239 -o bacteria_viruses.fasta
```

### Exclude a subtree (eukaryotes minus mammals)

```bash
taxafasta -i uniprot_trembl.fasta.gz -t 2759 -e 40674 -o euk_no_mammals.fasta
```

### Use pre-downloaded taxonomy files

```bash
taxafasta -i uniprot_trembl.fasta.gz -t 2 --taxdump /path/to/taxdump/ -o bacteria.fasta
```

### Uncompressed output

```bash
taxafasta -i uniprot_trembl.fasta.gz -t 9606 -o human.fasta --no-gzip
```

### Verbose progress

```bash
taxafasta -i uniprot_trembl.fasta.gz -t 2 -o bacteria.fasta -v
```

## Docker Usage

```bash
docker run --rm -v /data:/data ghcr.io/mriffle/taxafasta:latest \
  -i /data/uniprot_trembl.fasta.gz -t 2 -o /data/bacteria.fasta

# With pre-downloaded taxonomy
docker run --rm \
  -v /data:/data \
  -v /taxonomy:/taxonomy:ro \
  ghcr.io/mriffle/taxafasta:latest \
  -i /data/uniprot_trembl.fasta.gz -t 2 --taxdump /taxonomy -o /data/bacteria.fasta
```

## Common Taxonomy IDs

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

## NCBI Taxonomy Data

By default, the tool automatically downloads and caches `taxdump.tar.gz` from NCBI's FTP server on first run. Users can supply pre-downloaded taxonomy files with `--taxdump`. See: https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/

## Development

```bash
git clone https://github.com/mriffle/taxafasta.git
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

## License

Apache 2.0 — see [LICENSE](LICENSE) for details.
