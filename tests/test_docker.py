"""Docker image build and smoke tests (§15.4).

These tests require Docker to be installed and accessible.
They are skipped automatically if Docker is not available.
"""

from __future__ import annotations

import subprocess

import pytest

# Detect whether Docker is available
_docker_available = False
try:
    result = subprocess.run(
        ["docker", "info"],
        capture_output=True,
        timeout=10,
    )
    _docker_available = result.returncode == 0
except (FileNotFoundError, subprocess.TimeoutExpired):
    pass

pytestmark = pytest.mark.skipif(
    not _docker_available,
    reason="Docker not available",
)

IMAGE_TAG = "taxafasta-test:latest"


@pytest.fixture(scope="module")
def docker_image() -> str:
    """Build the Docker runtime image once per test module."""
    result = subprocess.run(
        [
            "docker",
            "build",
            "--build-arg",
            "VERSION=0.0.0",
            "--target",
            "runtime",
            "-t",
            IMAGE_TAG,
            ".",
        ],
        capture_output=True,
        text=True,
        timeout=300,
    )
    assert result.returncode == 0, f"Docker build failed:\n{result.stderr}"
    return IMAGE_TAG


def test_docker_version(docker_image: str) -> None:
    """Verify the Docker image runs and prints a version string."""
    result = subprocess.run(
        ["docker", "run", "--rm", docker_image, "--version"],
        capture_output=True,
        text=True,
        timeout=30,
    )
    assert result.returncode == 0, f"--version failed:\n{result.stderr}"
    assert "taxafasta" in result.stdout.lower()


def test_docker_help(docker_image: str) -> None:
    """Verify the Docker image --help works."""
    result = subprocess.run(
        ["docker", "run", "--rm", docker_image, "--help"],
        capture_output=True,
        text=True,
        timeout=30,
    )
    assert result.returncode == 0, f"--help failed:\n{result.stderr}"
    assert "--input" in result.stdout
    assert "--taxid" in result.stdout


def test_docker_filter(
    docker_image: str,
    tmp_path_factory: pytest.TempPathFactory,
) -> None:
    """Run an actual filter through the Docker image."""
    host_dir = tmp_path_factory.mktemp("docker_filter")

    # Create a tiny FASTA
    fasta = host_dir / "input.fasta"
    fasta.write_text(
        ">sp|P1|X OS=Org OX=7 PE=1 SV=1\nACGT\n"
        ">sp|P2|Y OS=Org OX=9606 PE=1 SV=1\nTTTT\n"
    )

    # Create a tiny taxdump
    taxdump = host_dir / "taxdump"
    taxdump.mkdir()
    (taxdump / "nodes.dmp").write_text(
        "1\t|\t1\t|\tno rank\t|\n"
        "131567\t|\t1\t|\tno rank\t|\n"
        "2\t|\t131567\t|\tsuperkingdom\t|\n"
        "6\t|\t2\t|\tgenus\t|\n"
        "7\t|\t6\t|\tspecies\t|\n"
        "2759\t|\t131567\t|\tsuperkingdom\t|\n"
        "40674\t|\t2759\t|\tclass\t|\n"
        "9606\t|\t40674\t|\tspecies\t|\n"
    )
    (taxdump / "merged.dmp").write_text("")

    result = subprocess.run(
        [
            "docker",
            "run",
            "--rm",
            "-v",
            f"{host_dir}:/data",
            docker_image,
            "-i",
            "/data/input.fasta",
            "-t",
            "2",
            "-o",
            "/data/output.fasta",
            "-d",
            "/data/taxdump",
            "--no-gzip",
        ],
        capture_output=True,
        text=True,
        timeout=30,
    )
    # Exit code 0 (no warnings in our clean input)
    assert result.returncode == 0, (
        f"Filter failed (rc={result.returncode}):\n{result.stderr}"
    )

    output = (host_dir / "output.fasta").read_text()
    assert "OX=7" in output
    assert "OX=9606" not in output

    # Log file should exist
    assert (host_dir / "output.fasta.log").exists()
