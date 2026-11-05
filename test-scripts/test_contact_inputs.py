#!/usr/bin/env python3
"""Small end-to-end fixtures for pairs, pairs.gz, and BAM APA input."""

from __future__ import annotations

import argparse
import gzip
import math
import shutil
import subprocess
import tempfile
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
DEFAULT_BINARY = ROOT / "apa-c-hic-main" / "apa_hic"


def read_matrix(path: Path) -> list[list[float]]:
    return [[float(value) for value in line.split()] for line in path.read_text().splitlines()]


def assert_matrix(actual: list[list[float]], expected: list[list[float]], label: str) -> None:
    if len(actual) != len(expected) or any(len(a) != len(e) for a, e in zip(actual, expected)):
        raise AssertionError(f"{label}: wrong matrix shape: {actual}")
    for row, (actual_row, expected_row) in enumerate(zip(actual, expected)):
        for col, (actual_value, expected_value) in enumerate(zip(actual_row, expected_row)):
            if not math.isclose(actual_value, expected_value, rel_tol=1e-5, abs_tol=1e-6):
                raise AssertionError(
                    f"{label}: cell ({row}, {col}) is {actual_value}, expected {expected_value}"
                )


def run_apa(binary: Path, contact_file: Path, loops: Path, output: Path) -> None:
    subprocess.run(
        [
            str(binary),
            "--window",
            "1",
            "--min-dist",
            "0",
            "--threads",
            "2",
            "--res",
            "10",
            str(contact_file),
            str(loops),
            str(output),
        ],
        check=True,
    )


def write_bam(samtools: str, directory: Path) -> Path:
    sam = directory / "contacts.sam"
    sam.write_text(
        "@HD\tVN:1.6\tSO:queryname\n"
        "@SQ\tSN:chr1\tLN:100\n"
        "r1a\t65\tchr1\t12\t60\t1M\t=\t32\t0\tA\tI\n"
        "r1a\t129\tchr1\t32\t60\t1M\t=\t12\t0\tA\tI\n"
        "r1b\t65\tchr1\t12\t60\t1M\t=\t32\t0\tA\tI\n"
        "r1b\t129\tchr1\t32\t60\t1M\t=\t12\t0\tA\tI\n"
        "r2\t65\tchr1\t22\t60\t1M\t=\t32\t0\tA\tI\n"
        "r2\t129\tchr1\t32\t60\t1M\t=\t22\t0\tA\tI\n"
        "r3\t65\tchr1\t12\t60\t1M\t=\t42\t0\tA\tI\n"
        "r3\t129\tchr1\t42\t60\t1M\t=\t12\t0\tA\tI\n"
        "r4\t65\tchr1\t2\t60\t1M\t=\t2\t0\tA\tI\n"
        "r4\t129\tchr1\t2\t60\t1M\t=\t2\t0\tA\tI\n"
    )
    bam = directory / "contacts.bam"
    subprocess.run([samtools, "view", "-b", "-o", str(bam), str(sam)], check=True)
    return bam


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--binary", type=Path, default=DEFAULT_BINARY)
    args = parser.parse_args()
    if not args.binary.exists():
        raise FileNotFoundError(f"Build apa_hic first: {args.binary}")

    expected = [
        [0.0, 0.0, 0.0],
        [0.0, 0.49382716, 0.74074074],
        [0.0, 0.74074074, 0.0],
    ]

    with tempfile.TemporaryDirectory(prefix="apa_contact_inputs_") as temp:
        directory = Path(temp)
        loops = directory / "loops.bedpe"
        loops.write_text("chr1\t10\t20\tchr1\t30\t40\n")

        pairs_text = (
            "## pairs format v1.0\n"
            "#chromsize: chr1 100\n"
            "#columns: readID chrom1 pos1 chrom2 pos2 strand1 strand2 pair_type\n"
            "r1a chr1 12 chr1 32 + - UU\n"
            "r1b chr1 12 chr1 32 + - UU\n"
            "r2 chr1 22 chr1 32 + - UU\n"
            "r3 chr1 12 chr1 42 + - UU\n"
            "r4 chr1 2 chr1 2 + - UU\n"
        )
        pairs = directory / "contacts.pairs"
        pairs.write_text(pairs_text)
        pairs_output = directory / "pairs.txt"
        run_apa(args.binary, pairs, loops, pairs_output)
        assert_matrix(read_matrix(pairs_output), expected, "pairs")

        pairs_gz = directory / "contacts.pairs.gz"
        with gzip.open(pairs_gz, "wt") as handle:
            handle.write(pairs_text)
        pairs_gz_output = directory / "pairs_gz.txt"
        run_apa(args.binary, pairs_gz, loops, pairs_gz_output)
        assert_matrix(read_matrix(pairs_gz_output), expected, "pairs.gz")

        reordered = directory / "reordered.pairs"
        reordered.write_text(
            "## pairs format v1.0\n"
            "#columns: pos2 chrom2 readID pos1 chrom1\n"
            "32 chr1 r1a 12 chr1\n"
            "32 chr1 r1b 12 chr1\n"
            "32 chr1 r2 22 chr1\n"
            "42 chr1 r3 12 chr1\n"
            "2 chr1 r4 2 chr1\n"
        )
        reordered_output = directory / "reordered.txt"
        run_apa(args.binary, reordered, loops, reordered_output)
        assert_matrix(read_matrix(reordered_output), expected, "reordered header/inferred sizes")

        samtools = shutil.which("samtools")
        if samtools is None:
            print("PASS: pairs variants (samtools unavailable; BAM fixture skipped)")
            return 0

        bam = write_bam(samtools, directory)
        bam_output = directory / "bam.txt"
        run_apa(args.binary, bam, loops, bam_output)
        assert_matrix(read_matrix(bam_output), expected, "BAM")

    print("PASS: pairs variants and BAM produce the expected aggregate-normalized APA")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
