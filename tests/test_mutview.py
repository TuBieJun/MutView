import os
import shutil
import subprocess
import sys
from pathlib import Path

import pytest

REFERENCE_REGION = (
    "TCTCCATAGTACCCCCAAATCAATAAATCATCAAGTCTTATTCTACCTTCCAAAGAGCCTTACATATGTTCCTTTATTTT"
    "CATCTGTAACACCACTATTCCTGTCTAAGCCTACCTATGTCATTTTTGGAAGAGAATATAGTCACCTATGCGACCTTCCC"
    "ACTTAAAATCCTACTATTTACGCTTCAGTAAAAGAAAAAAAATTTTTAATCTAAGTATGTAATTCTTTTGCTGAAGACAC"
    "TTCACTTGCTTCTGTGCCCTTAAACTGGTATGTTATCATGGTATAGTAGGCCATCCAAGAC"
)
REF_CHROM = "1"
REF_START = 55149
REF_END = 55449


def build_small_reference(ref_fasta, samtools_path):
    ref_index = ref_fasta.parent / (ref_fasta.name + ".fai")
    if ref_fasta.exists() and ref_index.exists():
        return

    full_length = REF_END
    prefix = "N" * (REF_START - 1)
    reference_sequence = prefix + REFERENCE_REGION
    if len(reference_sequence) != full_length:
        raise ValueError("Generated reference sequence length is incorrect")

    ref_lines = [">1"]
    for i in range(0, len(reference_sequence), 60):
        ref_lines.append(reference_sequence[i : i + 60])
    ref_fasta.write_text("\n".join(ref_lines) + "\n")

    completed = subprocess.run(
        [samtools_path, "faidx", str(ref_fasta)],
        capture_output=True,
        text=True,
        check=False,
    )
    if completed.returncode != 0:
        raise RuntimeError(
            f"samtools faidx failed for generated reference fasta:\n"
            f"stdout:\n{completed.stdout}\n"
            f"stderr:\n{completed.stderr}"
        )


def test_mutview_output_matches_fixture(tmp_path):
    repo_root = Path(__file__).resolve().parents[1]
    test_dir = repo_root / "tests"
    expected_html = test_dir / "test_out.html"
    bam_file = test_dir / "test.bam"
    ref_fasta = test_dir / "ref.fa"
    samtools_path = shutil.which("samtools")

    if not samtools_path:
        pytest.skip("samtools is not available in PATH")
    if not bam_file.exists():
        pytest.skip(f"fixture BAM file not found: {bam_file}")
    if not expected_html.exists():
        pytest.skip(f"fixture HTML file not found: {expected_html}")

    build_small_reference(ref_fasta, samtools_path)

    output_html = tmp_path / "test_out.generated.html"
    env = os.environ.copy()
    env["PYTHONPATH"] = str(repo_root / "src")

    completed = subprocess.run(
        [
            sys.executable,
            "-m",
            "mutview.process",
            "-c",
            "1",
            "-p",
            "55299",
            "-b",
            str(bam_file),
            "-r",
            str(ref_fasta),
            "-s",
            "-S",
            samtools_path,
            "-f",
            "Menlo",
            "-o",
            str(output_html),
        ],
        cwd=str(repo_root),
        env=env,
        capture_output=True,
        text=True,
    )

    assert completed.returncode == 0, (
        f"mutview exit status: {completed.returncode}\n"
        f"stdout:\n{completed.stdout}\n"
        f"stderr:\n{completed.stderr}"
    )
    assert output_html.exists(), "Expected output HTML file was not generated"

    expected_bytes = expected_html.read_bytes()
    actual_bytes = output_html.read_bytes()
    if expected_bytes != actual_bytes:
        expected_lines = expected_html.read_text(errors="replace").splitlines()
        actual_lines = output_html.read_text(errors="replace").splitlines()
        diff = "\n".join(
            list(
                __import__("difflib").unified_diff(
                    expected_lines,
                    actual_lines,
                    fromfile="expected",
                    tofile="actual",
                    lineterm="",
                )
            )[:50]
        )
        pytest.fail("Generated HTML does not match expected fixture.\n" + diff)
