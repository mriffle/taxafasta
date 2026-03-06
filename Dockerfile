# Multi-stage Dockerfile for taxafasta (§14.3)

ARG PYTHON_VERSION=3.12

# ---- Build stage ----
FROM python:${PYTHON_VERSION}-slim AS build

ARG VERSION=0.0.0
ENV SETUPTOOLS_SCM_PRETEND_VERSION=${VERSION}

WORKDIR /app
COPY pyproject.toml README.md LICENSE ./
COPY src/ src/

RUN pip install --no-cache-dir build hatch-vcs && \
    python -m build --wheel --outdir /app/dist

# ---- Test stage (used in CI, not shipped) ----
FROM python:${PYTHON_VERSION}-slim AS test

WORKDIR /app
COPY --from=build /app/dist/*.whl /app/dist/

RUN pip install --no-cache-dir "$(ls /app/dist/*.whl)[all]" && \
    pip install --no-cache-dir pytest pytest-cov pytest-benchmark mypy ruff

COPY tests/ tests/
COPY src/ src/
COPY pyproject.toml .

RUN pytest tests/

# ---- Runtime stage ----
FROM python:${PYTHON_VERSION}-slim AS runtime

ARG VERSION=dev
LABEL org.opencontainers.image.version="${VERSION}"
LABEL org.opencontainers.image.source="https://github.com/mriffle/taxafasta"
LABEL org.opencontainers.image.description="Filter UniProt FASTA files by NCBI taxonomy"

COPY --from=build /app/dist/*.whl /tmp/

RUN pip install --no-cache-dir "$(ls /tmp/*.whl)[all]" && \
    rm -rf /tmp/*.whl

ENTRYPOINT ["taxafasta"]
