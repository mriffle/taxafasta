"""taxafasta — Filter UniProt FASTA files by NCBI taxonomy."""

try:
    from taxafasta._version import __version__
except ModuleNotFoundError:  # pragma: no cover — editable install fallback
    __version__ = "0.0.0"

__all__ = ["__version__"]
