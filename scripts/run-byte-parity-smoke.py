#!/usr/bin/env python3
"""Run a small byte-for-byte C-vs-Rust samtools parity smoke suite.

This is intentionally dev-facing. The upstream Perl harness remains the CI
authority for the promoted fixture groups; this helper gives Phase 5 a quick
direct executable-vs-executable diff loop for commands that produce stable text
without harness normalization.
"""

from __future__ import annotations

import argparse
import difflib
import os
import subprocess
import tempfile
from dataclasses import dataclass
from pathlib import Path


@dataclass(frozen=True)
class Case:
    name: str
    args: tuple[str, ...]
    files: tuple[tuple[str, bytes], ...] = ()
    side_files: tuple[tuple[str, bytes], ...] = ()
    setup_external: tuple[tuple[str, ...], ...] = ()
    setup: tuple[tuple[str, ...], ...] = ()
    output_files: tuple[str, ...] = ()
    rendered_output_files: tuple[tuple[str, tuple[str, ...]], ...] = ()
    c_args: tuple[str, ...] | None = None
    rust_args: tuple[str, ...] | None = None
    compare_stderr: bool = True
    env: tuple[tuple[str, str], ...] = ()
    strip_lines_containing: tuple[str, ...] = ()


def repo_root() -> Path:
    return Path(__file__).resolve().parent.parent


def default_cases(root: Path) -> list[Case]:
    test = root / "samtools" / "test"
    dat = test / "dat"
    bedcov = test / "bedcov"
    checksum = test / "checksum"
    ampliconclip = test / "ampliconclip"
    ampliconstats = test / "ampliconstats"
    addrprg = test / "addrprg"
    collate = test / "collate"
    consensus = test / "consensus"
    cram_size = test / "cram_size"
    fixmate = test / "fixmate"
    large_pos = test / "large_pos"
    markdup = test / "markdup"
    reset = test / "reset"
    stat = test / "stat"
    quickcheck = test / "quickcheck"
    reference_setup = (
        (
            "view",
            "-e",
            "pos<1000||pos>1200",
            "-O",
            "cram,embed_ref=1",
            "-T",
            str(dat / "mpileup.ref.fa"),
            "-o",
            "{tmp}/reference-embed.cram",
            str(dat / "mpileup.1.sam"),
        ),
        ("index", "{tmp}/reference-embed.cram"),
    )
    mpileup_setup = (
        (
            "view",
            "-b",
            "-o",
            "{tmp}/mpileup.1.bam",
            str(dat / "mpileup.1.sam"),
        ),
        (
            "view",
            "-b",
            "-o",
            "{tmp}/mpileup.2.bam",
            str(dat / "mpileup.2.sam"),
        ),
        (
            "view",
            "-b",
            "-o",
            "{tmp}/mpileup.3.bam",
            str(dat / "mpileup.3.sam"),
        ),
        ("index", "{tmp}/mpileup.1.bam"),
        ("index", "{tmp}/mpileup.2.bam"),
        ("index", "{tmp}/mpileup.3.bam"),
    )
    view_filter_sam = (
        b"@HD\tVN:1.6\n"
        b"@SQ\tSN:chr1\tLN:100\n"
        b"@RG\tID:rg1\tLB:lib1\n"
        b"@RG\tID:rg2\tLB:lib2\n"
        b"r1\t0\tchr1\t1\t60\t4M\t*\t0\t0\tACGT\t!!!!\tRG:Z:rg1\tXX:i:7\n"
        b"r2\t0\tchr1\t2\t10\t4M\t*\t0\t0\tTGCA\t####\tRG:Z:rg2\tXX:i:8\n"
        b"r3\t4\t*\t0\t0\t*\t*\t0\t0\tNN\t!!\tXX:i:7\n"
    )
    view_filter_ref = b">chr1\n" + b"N" * 100 + b"\n"
    view_binary_aux_sam = (
        b"@HD\tVN:1.6\n"
        b"@SQ\tSN:ref\tLN:12\n"
        b"@RG\tID:rg1\n"
        b"r1\t0\tref\t1\t20\t2M\t*\t0\t0\tAC\t!!"
        b"\tRG:Z:rg1\tNM:i:0\tXX:Z:drop\n"
    )
    view_binary_aux_ref = b">ref\nACGTACGTACGT\n"
    return [
        Case(
            "view-sam-text",
            ("view", "--no-PG", "-h", str(dat / "view.001.sam")),
        ),
        Case(
            "view-bam-header-text",
            ("view", "--no-PG", "-h", str(dat / "test_input_1_a.bam")),
        ),
        Case(
            "view-bam-count",
            ("view", "-c", str(dat / "test_input_1_a.bam")),
        ),
        Case(
            "view-bam-count-mapq",
            ("view", "-c", "-q", "20", str(dat / "test_input_1_a.bam")),
        ),
        Case(
            "view-bam-count-required-flag",
            ("view", "-c", "-f", "0x2", str(dat / "test_input_1_a.bam")),
        ),
        Case(
            "view-bam-count-filtering-flag",
            ("view", "-c", "-F", "0x4", str(dat / "test_input_1_a.bam")),
        ),
        Case(
            "view-bam-region-text",
            ("view", "--no-PG", str(dat / "test_input_1_a.bam"), "ref1:1-30"),
        ),
        Case(
            "view-expr-mapq-bam-text",
            ("view", "--no-PG", "-e", "mapq>=20", str(dat / "test_input_1_a.bam")),
        ),
        Case(
            "view-expr-mapq-bam-count",
            ("view", "-c", "-e", "mapq>=20", str(dat / "test_input_1_a.bam")),
        ),
        Case(
            "view-expr-flag-name-bam-text",
            (
                "view",
                "--no-PG",
                "-e",
                "flag.proper_pair",
                str(dat / "test_input_1_a.bam"),
            ),
        ),
        Case(
            "view-expr-endpos-bam-text",
            ("view", "--no-PG", "-e", "endpos>=40", str(dat / "test_input_1_a.bam")),
        ),
        Case(
            "view-expr-rlen-bam-text",
            ("view", "--no-PG", "-e", "rlen>=10", str(dat / "test_input_1_a.bam")),
        ),
        Case(
            "view-expr-sclen-bam-text",
            ("view", "--no-PG", "-e", "sclen==0", str(dat / "test_input_1_a.bam")),
        ),
        Case(
            "view-expr-hclen-bam-text",
            ("view", "--no-PG", "-e", "hclen==0", str(dat / "test_input_1_a.bam")),
        ),
        Case(
            "view-expr-mpos-bam-text",
            ("view", "--no-PG", "-e", "mpos>0", str(dat / "test_input_1_a.bam")),
        ),
        Case(
            "view-expr-pnext-bam-text",
            ("view", "--no-PG", "-e", "pnext>0", str(dat / "test_input_1_a.bam")),
        ),
        Case(
            "view-expr-tlen-bam-text",
            ("view", "--no-PG", "-e", "tlen!=0", str(dat / "test_input_1_a.bam")),
        ),
        Case(
            "view-expr-rnext-bam-text",
            ("view", "--no-PG", "-e", 'rnext=="ref1"', str(dat / "test_input_1_a.bam")),
        ),
        Case(
            "view-expr-mrname-bam-text",
            ("view", "--no-PG", "-e", 'mrname=="ref1"', str(dat / "test_input_1_a.bam")),
        ),
        Case(
            "view-expr-refid-bam-text",
            ("view", "--no-PG", "-e", "refid>=0", str(dat / "test_input_1_a.bam")),
        ),
        Case(
            "view-expr-mrefid-bam-text",
            ("view", "--no-PG", "-e", "mrefid>=0", str(dat / "test_input_1_a.bam")),
        ),
        Case(
            "view-expr-ncigar-bam-text",
            ("view", "--no-PG", "-e", "ncigar>0", str(dat / "test_input_1_a.bam")),
        ),
        Case(
            "view-expr-refid-sam-text",
            ("view", "--no-PG", "-e", "refid>=0", str(dat / "view.001.sam")),
        ),
        Case(
            "view-expr-mrefid-sam-text",
            ("view", "--no-PG", "-e", "mrefid>=0", str(dat / "view.001.sam")),
        ),
        Case(
            "view-expr-rnext-sam-text",
            ("view", "--no-PG", "-e", 'rnext=="*"', str(dat / "view.001.sam")),
        ),
        Case(
            "view-expr-ncigar-sam-text",
            ("view", "--no-PG", "-e", "ncigar>0", str(dat / "view.001.sam")),
        ),
        Case(
            "view-expr-region-bam-text",
            (
                "view",
                "--no-PG",
                "-e",
                "mapq>=20",
                str(dat / "test_input_1_a.bam"),
                "ref1:1-60",
            ),
        ),
        Case(
            "view-qname-file-text",
            ("view", "--no-PG", "-N", "{tmp}/qnames.txt", "{tmp}/view-filters.sam"),
            files=(
                ("view-filters.sam", view_filter_sam),
                ("qnames.txt", b"r1\nr3\n"),
            ),
        ),
        Case(
            "view-qname-file-negated-text",
            ("view", "--no-PG", "-N", "^{tmp}/qnames.txt", "{tmp}/view-filters.sam"),
            files=(
                ("view-filters.sam", view_filter_sam),
                ("qnames.txt", b"r1\nr3\n"),
            ),
        ),
        Case(
            "view-read-group-text",
            ("view", "--no-PG", "-r", "rg1", "{tmp}/view-filters.sam"),
            files=(("view-filters.sam", view_filter_sam),),
        ),
        Case(
            "view-read-group-file-text",
            ("view", "--no-PG", "-R", "{tmp}/rgs.txt", "{tmp}/view-filters.sam"),
            files=(
                ("view-filters.sam", view_filter_sam),
                ("rgs.txt", b"rg2\n"),
            ),
        ),
        Case(
            "view-read-group-exclude-no-rg-text",
            ("view", "--no-PG", "-r", "rg1", "-n", "{tmp}/view-filters.sam"),
            files=(("view-filters.sam", view_filter_sam),),
        ),
        Case(
            "view-library-text",
            ("view", "--no-PG", "-l", "lib2", "{tmp}/view-filters.sam"),
            files=(("view-filters.sam", view_filter_sam),),
        ),
        Case(
            "view-aux-tag-text",
            ("view", "--no-PG", "-d", "XX:7", "{tmp}/view-filters.sam"),
            files=(("view-filters.sam", view_filter_sam),),
        ),
        Case(
            "view-aux-tag-file-text",
            ("view", "--no-PG", "-D", "XX:{tmp}/aux-values.txt", "{tmp}/view-filters.sam"),
            files=(
                ("view-filters.sam", view_filter_sam),
                ("aux-values.txt", b"7\n"),
            ),
        ),
        Case(
            "view-count-save-counts-json",
            (
                "view",
                "-c",
                "--save-counts",
                "{side_tmp}/view-counts.json",
                "-f",
                "0x2",
                str(dat / "test_input_1_a.bam"),
            ),
            output_files=("view-counts.json",),
        ),
        Case(
            "view-cram-no-ref-expr-mapq-save-counts-json",
            (
                "view",
                "-c",
                "--save-counts",
                "{side_tmp}/view-cram-expr-mapq-counts.json",
                "-o",
                "{side_tmp}/view-cram-expr-mapq-count.txt",
                "-e",
                "mapq>=20",
                str(dat / "test_input_1_a.cram"),
            ),
            output_files=(
                "view-cram-expr-mapq-counts.json",
                "view-cram-expr-mapq-count.txt",
            ),
        ),
        Case(
            "view-cram-no-ref-expr-flag-save-counts-json",
            (
                "view",
                "-c",
                "--save-counts",
                "{side_tmp}/view-cram-expr-flag-counts.json",
                "-o",
                "{side_tmp}/view-cram-expr-flag-count.txt",
                "-e",
                "flag.proper_pair",
                str(dat / "test_input_1_a.cram"),
            ),
            output_files=(
                "view-cram-expr-flag-counts.json",
                "view-cram-expr-flag-count.txt",
            ),
        ),
        Case(
            "view-cram-no-ref-read-group-save-counts-json",
            (
                "view",
                "-c",
                "--save-counts",
                "{side_tmp}/view-cram-rg-counts.json",
                "-o",
                "{side_tmp}/view-cram-rg-count.txt",
                "-r",
                "rg1",
                "-n",
                "{tmp}/view-filters.cram",
            ),
            files=(
                ("view-filters.sam", view_filter_sam),
                ("view-filters.fa", view_filter_ref),
            ),
            setup=(
                (
                    "view",
                    "--no-PG",
                    "-C",
                    "-T",
                    "{tmp}/view-filters.fa",
                    "-o",
                    "{tmp}/view-filters.cram",
                    "{tmp}/view-filters.sam",
                ),
            ),
            output_files=("view-cram-rg-counts.json", "view-cram-rg-count.txt"),
        ),
        Case(
            "view-cram-no-ref-library-save-counts-json",
            (
                "view",
                "-c",
                "--save-counts",
                "{side_tmp}/view-cram-library-counts.json",
                "-o",
                "{side_tmp}/view-cram-library-count.txt",
                "-l",
                "lib2",
                "{tmp}/view-filters.cram",
            ),
            files=(
                ("view-filters.sam", view_filter_sam),
                ("view-filters.fa", view_filter_ref),
            ),
            setup=(
                (
                    "view",
                    "--no-PG",
                    "-C",
                    "-T",
                    "{tmp}/view-filters.fa",
                    "-o",
                    "{tmp}/view-filters.cram",
                    "{tmp}/view-filters.sam",
                ),
            ),
            output_files=("view-cram-library-counts.json", "view-cram-library-count.txt"),
        ),
        Case(
            "view-cram-no-ref-aux-save-counts-json",
            (
                "view",
                "-c",
                "--save-counts",
                "{side_tmp}/view-cram-aux-counts.json",
                "-o",
                "{side_tmp}/view-cram-aux-count.txt",
                "-d",
                "XX:8",
                "{tmp}/view-filters.cram",
            ),
            files=(
                ("view-filters.sam", view_filter_sam),
                ("view-filters.fa", view_filter_ref),
            ),
            setup=(
                (
                    "view",
                    "--no-PG",
                    "-C",
                    "-T",
                    "{tmp}/view-filters.fa",
                    "-o",
                    "{tmp}/view-filters.cram",
                    "{tmp}/view-filters.sam",
                ),
            ),
            output_files=("view-cram-aux-counts.json", "view-cram-aux-count.txt"),
        ),
        Case(
            "view-bam-input-bam-output-remove-tag",
            (
                "view",
                "--no-PG",
                "-b",
                "-x",
                "NM",
                "-o",
                "{side_tmp}/view-aux-out.bam",
                "{tmp}/view-aux-in.bam",
            ),
            files=(("view-aux.sam", view_binary_aux_sam),),
            setup=(
                (
                    "view",
                    "--no-PG",
                    "-b",
                    "-o",
                    "{tmp}/view-aux-in.bam",
                    "{tmp}/view-aux.sam",
                ),
            ),
            rendered_output_files=(
                ("view-aux-out.bam", ("view", "--no-PG", "{output}")),
            ),
        ),
        Case(
            "view-cram-input-cram-output-remove-tag",
            (
                "view",
                "--no-PG",
                "-C",
                "-T",
                "{tmp}/view-aux.fa",
                "-x",
                "XX",
                "-o",
                "{side_tmp}/view-aux-out.cram",
                "{tmp}/view-aux-in.cram",
            ),
            files=(
                ("view-aux.sam", view_binary_aux_sam),
                ("view-aux.fa", view_binary_aux_ref),
            ),
            setup=(
                (
                    "view",
                    "--no-PG",
                    "-C",
                    "-T",
                    "{tmp}/view-aux.fa",
                    "-o",
                    "{tmp}/view-aux-in.cram",
                    "{tmp}/view-aux.sam",
                ),
            ),
            rendered_output_files=(
                (
                    "view-aux-out.cram",
                    ("view", "--no-PG", "-T", "{tmp}/view-aux.fa", "{output}"),
                ),
            ),
        ),
        Case(
            "sort-sam-text",
            (
                "sort",
                "--no-PG",
                "-O",
                "sam",
                "-o",
                "-",
                str(dat / "sort_name_input_1.sam"),
            ),
        ),
        Case(
            "flagstat-bam-text",
            ("flagstat", str(dat / "test_input_1_a.bam")),
        ),
        Case(
            "flagstat-bam-json",
            ("flagstat", "-O", "json", str(dat / "test_input_1_a.bam")),
        ),
        Case(
            "flagstat-bam-tsv",
            ("flagstat", "-O", "tsv", str(dat / "test_input_1_a.bam")),
        ),
        Case(
            "flagstat-cram-no-reference-text",
            ("flagstat", str(dat / "test_input_1_a.cram")),
        ),
        Case(
            "idxstats-bam-text",
            ("idxstats", str(dat / "test_input_1_a.bam")),
        ),
        Case(
            "idxstats-cram-no-reference-text",
            ("idxstats", str(dat / "test_input_1_a.cram")),
        ),
        Case(
            "stats-bam-text",
            ("stats", str(dat / "test_input_1_a.bam")),
        ),
        Case(
            "stats-remove-dups-text",
            ("stats", "-d", str(dat / "test_input_1_a.bam")),
        ),
        Case(
            "stats-required-flag-text",
            ("stats", "-f", "0x2", str(dat / "test_input_1_a.bam")),
        ),
        Case(
            "stats-filtering-flag-text",
            ("stats", "-F", "0x400", str(dat / "test_input_1_a.bam")),
        ),
        Case(
            "stats-insert-size-zero-text",
            ("stats", "-i", "0", str(dat / "test_input_1_a.bam")),
        ),
        Case(
            "stats-coverage-bins-text",
            ("stats", "-c", "1,20,1", str(dat / "test_input_1_a.bam")),
        ),
        Case(
            "stats-trim-quality-text",
            ("stats", "-q", "20", str(dat / "test_input_1_a.bam")),
        ),
        Case(
            "stats-read-length-text",
            ("stats", "-l", "10", str(dat / "test_input_1_a.bam")),
        ),
        Case(
            "stats-most-inserts-text",
            ("stats", "-m", "0.5", str(dat / "test_input_1_a.bam")),
        ),
        Case(
            "stats-target-file-text",
            ("stats", "-t", str(stat / "11.stats.targets"), str(stat / "11_target.bam")),
        ),
        Case(
            "stats-cov-threshold-text",
            (
                "stats",
                "-g",
                "1",
                "-t",
                str(stat / "11.stats.targets"),
                str(stat / "11_target.bam"),
            ),
        ),
        Case(
            "stats-reference-text",
            (
                "stats",
                "-r",
                str(stat / "test.fa"),
                str(stat / "1_map_cigar.sam"),
            ),
        ),
        Case(
            "stats-ref-stats-text",
            ("stats", "--ref-stats", str(stat / "11_target.sam")),
        ),
        Case(
            "stats-reference-ref-stats-text",
            (
                "stats",
                "-r",
                str(stat / "test1.fa"),
                "--ref-stats",
                str(stat / "11_target.sam"),
            ),
        ),
        Case(
            "stats-region-text",
            (
                "stats",
                "-r",
                str(stat / "test1.fa"),
                str(stat / "11_target.bam"),
                "alpha:10-20",
            ),
        ),
        Case(
            "index-bai-output",
            (
                "index",
                "-o",
                "{side_tmp}/test_input_1_a.bam.bai",
                str(dat / "test_input_1_a.bam"),
            ),
            output_files=("test_input_1_a.bam.bai",),
        ),
        Case(
            "checksum-bam-text",
            ("checksum", str(checksum / "chk1.bam")),
        ),
        Case(
            "checksum-output-file",
            ("checksum", "-o", "{side_tmp}/checksum.txt", str(checksum / "chk1.bam")),
            output_files=("checksum.txt",),
        ),
        Case(
            "faidx-text",
            ("faidx", str(dat / "view.001.fa"), "ref1:1-12"),
        ),
        Case(
            "faidx-index-output",
            ("faidx", "{side_tmp}/ref.fa"),
            side_files=(("ref.fa", (dat / "view.001.fa").read_bytes()),),
            output_files=("ref.fa.fai",),
        ),
        Case(
            "fqidx-text",
            ("fqidx", "{tmp}/reads.fq", "r1:1-4"),
            files=(("reads.fq", b"@r1\nACGTACGT\n+\nabcdefgh\n"),),
        ),
        Case(
            "fqidx-index-output",
            ("fqidx", "{side_tmp}/reads.fq"),
            side_files=(
                ("reads.fq", b"@r1\nACGTACGT\n+\nabcdefgh\n@r2\nTTAA\n+\n!!!!\n"),
            ),
            output_files=("reads.fq.fai",),
        ),
        Case(
            "faidx-zero-warning",
            ("faidx", str(dat / "view.001.fa"), "ref1:10000000-10000005"),
        ),
        Case(
            "faidx-truncated-warning",
            ("faidx", str(dat / "view.001.fa"), "ref1:50-80"),
        ),
        Case(
            "faidx-continue-missing-warning",
            (
                "faidx",
                "--continue",
                str(dat / "view.001.fa"),
                "ref1:1-12",
                "missing",
                "ref1:50-80",
            ),
        ),
        Case(
            "fqidx-zero-warning",
            ("fqidx", "{tmp}/reads.fq", "r1:99-100"),
            files=(("reads.fq", b"@r1\nACGTACGT\n+\nabcdefgh\n"),),
        ),
        Case(
            "fqidx-truncated-warning",
            ("fqidx", "{tmp}/reads.fq", "r1:5-99"),
            files=(("reads.fq", b"@r1\nACGTACGT\n+\nabcdefgh\n"),),
        ),
        Case(
            "fqidx-continue-missing-warning",
            (
                "fqidx",
                "--continue",
                "{tmp}/reads.fq",
                "r1:1-4",
                "missing",
                "r1:5-99",
            ),
            files=(("reads.fq", b"@r1\nACGTACGT\n+\nabcdefgh\n"),),
        ),
        Case(
            "bedcov-bam-text",
            (
                "bedcov",
                str(bedcov / "bedcov.bed"),
                str(bedcov / "bedcov.bam"),
            ),
        ),
        Case(
            "bedcov-count-header-text",
            (
                "bedcov",
                "-c",
                "-H",
                str(bedcov / "bedcov.bed"),
                str(bedcov / "bedcov.bam"),
            ),
        ),
        Case(
            "bedcov-skip-deletions-text",
            (
                "bedcov",
                "-j",
                str(bedcov / "bedcov.bed"),
                str(bedcov / "bedcov.bam"),
            ),
        ),
        Case(
            "bedcov-flag-mask-text",
            (
                "bedcov",
                "-g512",
                "-G2048",
                str(bedcov / "bedcov.bed"),
                str(bedcov / "bedcov.bam"),
            ),
        ),
        Case(
            "coverage-bam-text",
            ("coverage", str(dat / "test_input_1_a.bam")),
        ),
        Case(
            "coverage-output-file",
            ("coverage", "-o", "{side_tmp}/coverage.txt", str(dat / "test_input_1_a.bam")),
            output_files=("coverage.txt",),
        ),
        Case(
            "coverage-ascii-histogram-text",
            ("coverage", "-m", "-A", "-w", "20", "{tmp}/coverage-hist.sam"),
            files=(
                (
                    "coverage-hist.sam",
                    b"@HD\tVN:1.6\tSO:coordinate\n"
                    b"@SQ\tSN:chr1\tLN:80\n"
                    b"r1\t0\tchr1\t1\t60\t40M\t*\t0\t0\t"
                    b"ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT\t"
                    b"!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n"
                    b"r2\t0\tchr1\t1\t60\t40M\t*\t0\t0\t"
                    b"ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT\t"
                    b"!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n",
                ),
            ),
        ),
        Case(
            "coverage-ascii-depth-plot-text",
            ("coverage", "-D", "-A", "-w", "20", "{tmp}/coverage-depth.sam"),
            files=(
                (
                    "coverage-depth.sam",
                    b"@HD\tVN:1.6\tSO:coordinate\n"
                    b"@SQ\tSN:chr1\tLN:80\n"
                    b"r1\t0\tchr1\t1\t60\t40M\t*\t0\t0\t"
                    b"ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT\t"
                    b"!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n"
                    b"r2\t0\tchr1\t1\t60\t40M\t*\t0\t0\t"
                    b"ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT\t"
                    b"!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n",
                ),
            ),
        ),
        Case(
            "coverage-ascii-histogram-uneven-tail-text",
            ("coverage", "-m", "-A", "-w", "20", "{tmp}/coverage-uneven.sam"),
            files=(
                (
                    "coverage-uneven.sam",
                    b"@HD\tVN:1.6\tSO:coordinate\n"
                    b"@SQ\tSN:chr1\tLN:83\n"
                    b"r1\t0\tchr1\t1\t60\t83M\t*\t0\t0\t"
                    b"AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA\t"
                    b"!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n",
                ),
            ),
        ),
        Case(
            "coverage-ascii-histogram-columns-text",
            ("coverage", "-m", "-A", "{tmp}/coverage-columns.sam"),
            files=(
                (
                    "coverage-columns.sam",
                    b"@HD\tVN:1.6\tSO:coordinate\n"
                    b"@SQ\tSN:chr1\tLN:80\n"
                    b"r1\t0\tchr1\t1\t60\t40M\t*\t0\t0\t"
                    b"ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT\t"
                    b"!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n"
                    b"r2\t0\tchr1\t1\t60\t40M\t*\t0\t0\t"
                    b"ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT\t"
                    b"!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n",
                ),
            ),
            env=(("COLUMNS", "70"),),
        ),
        Case(
            "depth-bam-text",
            ("depth", "-r", "ref1:1-10", str(dat / "test_input_1_a.bam")),
        ),
        Case(
            "depth-output-file",
            (
                "depth",
                "-o",
                "{side_tmp}/depth.txt",
                "-r",
                "ref1:1-10",
                str(dat / "test_input_1_a.bam"),
            ),
            output_files=("depth.txt",),
        ),
        Case(
            "head-bam-text",
            ("head", "-n", "2", str(dat / "test_input_1_a.bam")),
        ),
        Case(
            "tview-large-position-text",
            ("tview", "-d", "T", "-p", "CHROMOSOME_I:10000000000", "{tmp}/longref.sam.gz"),
            setup_external=(
                (
                    "{root}/repos/htslib-rs/repos/htslib/bgzip",
                    "-c",
                    str(large_pos / "longref.sam"),
                    ">",
                    "{tmp}/longref.sam.gz",
                ),
            ),
            setup=(("index", "-c", "{tmp}/longref.sam.gz"),),
        ),
        Case(
            "quickcheck-good-bam",
            ("quickcheck", str(dat / "test_input_1_a.bam")),
        ),
        Case(
            "quickcheck-all-fixtures-verbose",
            (
                "quickcheck",
                "-v",
                str(quickcheck / "1.quickcheck.badeof.bam"),
                str(quickcheck / "2.quickcheck.badheader.bam"),
                str(quickcheck / "3.quickcheck.ok.bam"),
                str(quickcheck / "4.quickcheck.ok.bam"),
                str(quickcheck / "5.quickcheck.scramble30.truncated.cram"),
                str(quickcheck / "6.quickcheck.cram21.ok.cram"),
                str(quickcheck / "7.quickcheck.cram30.ok.cram"),
                str(quickcheck / "8.quickcheck.cram21.truncated.cram"),
                str(quickcheck / "9.quickcheck.cram30.truncated.cram"),
                str(quickcheck / "10.quickcheck.notargets.bam"),
            ),
        ),
        Case(
            "dict-fasta-text",
            ("dict", str(dat / "mpileup.ref.fa")),
        ),
        Case(
            "dict-output-file",
            ("dict", "-o", "{side_tmp}/dict.sam", str(dat / "mpileup.ref.fa")),
            output_files=("dict.sam",),
        ),
        Case(
            "reference-embedded-cram-text",
            ("reference", "-e", "{tmp}/reference-embed.cram"),
            setup=reference_setup,
        ),
        Case(
            "reference-md-cram-text",
            ("reference", "{tmp}/reference-embed.cram"),
            setup=reference_setup,
        ),
        Case(
            "reference-md-output-file",
            ("reference", "-o", "{side_tmp}/reference.fa", "{tmp}/reference-embed.cram"),
            setup=reference_setup,
            output_files=("reference.fa",),
        ),
        Case(
            "reference-md-region-cram-text",
            ("reference", "-r", "17:1000-1500", "{tmp}/reference-embed.cram"),
            setup=reference_setup,
        ),
        Case(
            "cram-size-normal-text",
            ("cram-size", str(cram_size / "mpileup.1.cram")),
        ),
        Case(
            "cram-size-verbose-text",
            ("cram-size", "-v", str(cram_size / "mpileup.1.cram")),
        ),
        Case(
            "cram-size-encodings-text",
            ("cram-size", "-e", str(cram_size / "mpileup.1.cram")),
        ),
        Case(
            "mpileup-bam-list-baq-text",
            (
                "mpileup",
                "-b",
                "{tmp}/mpileup.bam.list",
                "-f",
                str(dat / "mpileup.ref.fa"),
                "-r17:100-150",
            ),
            files=(
                (
                    "mpileup.bam.list",
                    b"{tmp}/mpileup.1.bam\n{tmp}/mpileup.2.bam\n{tmp}/mpileup.3.bam\n",
                ),
            ),
            setup=mpileup_setup,
        ),
        Case(
            "mpileup-no-baq-filter-text",
            (
                "mpileup",
                "-B",
                "--ff",
                "0x14",
                "-f",
                str(dat / "mpileup.ref.fa"),
                "-r17:1050-1060",
                "{tmp}/mpileup.1.bam",
            ),
            setup=mpileup_setup,
        ),
        Case(
            "mpileup-output-file",
            (
                "mpileup",
                "-o",
                "{side_tmp}/mpileup.txt",
                "-B",
                "--ff",
                "0x14",
                "-f",
                str(dat / "mpileup.ref.fa"),
                "-r17:1050-1060",
                "{tmp}/mpileup.1.bam",
            ),
            setup=mpileup_setup,
            output_files=("mpileup.txt",),
        ),
        Case(
            "mpileup-overlap-removal-text",
            ("mpileup", str(test / "mpileup" / "overlap.bam")),
        ),
        Case(
            "merge-seeded-sam-text",
            (
                "merge",
                "--no-PG",
                "-s",
                "1",
                "-O",
                "sam",
                "-",
                str(dat / "test_input_1_a.sam"),
                str(dat / "test_input_1_b.sam"),
                str(dat / "test_input_1_c.sam"),
            ),
        ),
        Case(
            "merge-r-filename-rg-text",
            (
                "merge",
                "--no-PG",
                "-r",
                "-O",
                "sam",
                "-",
                str(test / "merge" / "test_no_pg_rg_co.sam"),
            ),
        ),
        Case(
            "merge-template-coordinate-text",
            (
                "merge",
                "--no-PG",
                "-O",
                "sam",
                "--template-coordinate",
                "-",
                str(test / "merge" / "test_template_coordinate.1.sam"),
                str(test / "merge" / "test_template_coordinate.2.sam"),
            ),
        ),
        Case(
            "samples-bam-text",
            ("samples", str(dat / "test_input_1_a.bam")),
        ),
        Case(
            "samples-output-file",
            ("samples", "-o", "{side_tmp}/samples.txt", str(dat / "test_input_1_a.bam")),
            output_files=("samples.txt",),
        ),
        Case(
            "flags-int-text",
            ("flags", "12"),
        ),
        Case(
            "rmdup-sam-output",
            (),
            c_args=(
                "rmdup",
                "-s",
                "{tmp}/rmdup.sam",
                "{side_tmp}/rmdup.out.sam",
            ),
            rust_args=(
                "rmdup",
                "-s",
                "--no-PG",
                "{tmp}/rmdup.sam",
                "{side_tmp}/rmdup.out.sam",
            ),
            files=(
                (
                    "rmdup.sam",
                    b"@HD\tVN:1.6\n"
                    b"@SQ\tSN:chr1\tLN:8\n"
                    b"low\t0\tchr1\t1\t10\t4M\t*\t0\t0\tACGT\t!!!!\n"
                    b"high\t0\tchr1\t1\t60\t4M\t*\t0\t0\tTGCA\t####\n",
                ),
            ),
            output_files=("rmdup.out.sam",),
        ),
        Case(
            "rmdup-paired-sam-output",
            (),
            c_args=(
                "rmdup",
                "{tmp}/rmdup-pe.sam",
                "{side_tmp}/rmdup-pe.out.sam",
            ),
            rust_args=(
                "rmdup",
                "--no-PG",
                "{tmp}/rmdup-pe.sam",
                "{side_tmp}/rmdup-pe.out.sam",
            ),
            files=(
                (
                    "rmdup-pe.sam",
                    b"@HD\tVN:1.6\tSO:coordinate\n"
                    b"@SQ\tSN:chr1\tLN:200\n"
                    b"pair_a\t99\tchr1\t1\t60\t10M\t=\t91\t100\tAAAAAAAAAA\t!!!!!!!!!!\n"
                    b"pair_a\t147\tchr1\t91\t60\t10M\t=\t1\t-100\tTTTTTTTTTT\t!!!!!!!!!!\n"
                    b"pair_b\t99\tchr1\t1\t10\t10M\t=\t91\t100\tCCCCCCCCCC\t!!!!!!!!!!\n"
                    b"pair_b\t147\tchr1\t91\t10\t10M\t=\t1\t-100\tGGGGGGGGGG\t!!!!!!!!!!\n"
                    b"pair_c\t99\tchr1\t2\t10\t10M\t=\t91\t100\tACACACACAC\t!!!!!!!!!!\n"
                    b"pair_c\t147\tchr1\t91\t10\t10M\t=\t2\t-100\tTGTGTGTGTG\t!!!!!!!!!!\n",
                ),
            ),
            output_files=("rmdup-pe.out.sam",),
        ),
        Case(
            "addreplacerg-sam-text",
            (
                "addreplacerg",
                "--no-PG",
                "-r",
                r"@RG\tID:foo\tSM:bar",
                "-O",
                "sam",
                str(addrprg / "4_fixup_norg.sam"),
            ),
        ),
        Case(
            "reset-sam-text",
            ("reset", "--no-PG", "-O", "sam", str(reset / "seq.sam")),
        ),
        Case(
            "depad-sam-text",
            (
                "depad",
                "-T",
                str(dat / "depad.001.fa"),
                "-s",
                "--no-PG",
                str(dat / "depad.001p.sam"),
            ),
        ),
        Case(
            "split-template-sam-output",
            (
                "split",
                "--no-PG",
                "-f",
                "{side_tmp}/split.%#.sam",
                "{tmp}/split-rg.bam",
            ),
            files=(
                (
                    "split-rg.sam",
                    b"@HD\tVN:1.6\n"
                    b"@SQ\tSN:chr1\tLN:20\n"
                    b"@RG\tID:g1\n"
                    b"@RG\tID:g2\n"
                    b"r1\t0\tchr1\t1\t60\t4M\t*\t0\t0\tACGT\t!!!!\tRG:Z:g1\n"
                    b"r2\t0\tchr1\t5\t60\t4M\t*\t0\t0\tTGCA\t####\tRG:Z:g2\n",
                ),
            ),
            setup=(
                (
                    "view",
                    "-b",
                    "--no-PG",
                    "-o",
                    "{tmp}/split-rg.bam",
                    "{tmp}/split-rg.sam",
                ),
            ),
            output_files=("split.0.sam", "split.1.sam"),
        ),
        Case(
            "split-missing-rg-no-unaccounted-error",
            (
                "split",
                "--no-PG",
                "--output-fmt",
                "sam",
                "-f",
                "{side_tmp}/split.%#.sam",
                "{tmp}/split-missing-rg.sam",
            ),
            files=(
                (
                    "split-missing-rg.sam",
                    b"@HD\tVN:1.6\n"
                    b"@SQ\tSN:chr1\tLN:20\n"
                    b"@RG\tID:g1\n"
                    b"r1\t0\tchr1\t1\t60\t4M\t*\t0\t0\tACGT\t!!!!\tRG:Z:g1\n"
                    b"r2\t0\tchr1\t5\t60\t4M\t*\t0\t0\tTGCA\t####\n",
                ),
            ),
            output_files=("split.0.sam",),
        ),
        Case(
            "collate-sam-text",
            (
                "collate",
                "--no-PG",
                "--output-fmt=sam",
                "-O",
                str(dat / "test_input_1_d.sam"),
            ),
        ),
        Case(
            "consensus-fasta-text",
            ("consensus", "-f", "fasta", str(consensus / "consen1.sam")),
        ),
        Case(
            "ampliconclip-sam-stdout",
            (
                "ampliconclip",
                "--no-PG",
                "-O",
                "sam",
                "-b",
                str(ampliconclip / "ac_test.bed"),
                "-o",
                "-",
                str(ampliconclip / "1_test_data.sam"),
            ),
            compare_stderr=False,
        ),
        Case(
            "ampliconstats-single-ref-output-file",
            (
                "ampliconstats",
                "-S",
                "-t",
                "50",
                "-d",
                "1,20,100",
                "-o",
                "{side_tmp}/ampliconstats.txt",
                str(ampliconclip / "ac_test.bed"),
                str(ampliconclip / "1_hard_clipped.expected.sam"),
                str(ampliconclip / "1_soft_clipped.expected.sam"),
                str(ampliconclip / "1_soft_clipped_strand.expected.sam"),
                str(ampliconclip / "2_both_clipped.expected.sam"),
            ),
            output_files=("ampliconstats.txt",),
            strip_lines_containing=("Samtools version", "Command line"),
        ),
        Case(
            "ampliconstats-mixed-output-file",
            (
                "ampliconstats",
                "-c",
                "0",
                "-o",
                "{side_tmp}/ampliconstats.txt",
                str(ampliconclip / "multi_ref.bed"),
                str(ampliconstats / "mixed_clipped.sam"),
            ),
            output_files=("ampliconstats.txt",),
            strip_lines_containing=("Samtools version", "Command line"),
        ),
        Case(
            "ampliconstats-partial-output-file",
            (
                "ampliconstats",
                "-c",
                "0",
                "-o",
                "{side_tmp}/ampliconstats.txt",
                str(ampliconclip / "ac_test.bed"),
                str(ampliconstats / "mixed_clipped.sam"),
            ),
            output_files=("ampliconstats.txt",),
            strip_lines_containing=("Samtools version", "Command line"),
        ),
        Case(
            "markdup-sam-stdout",
            (
                "markdup",
                "-O",
                "sam",
                "--no-PG",
                str(markdup / "18_primary_duplicate_count.sam"),
                "-",
            ),
        ),
        Case(
            "fixmate-sam-stdout",
            (
                "fixmate",
                "-O",
                "sam",
                "--no-PG",
                str(fixmate / "7_two_read_mapped.sam"),
                "-",
            ),
        ),
        Case(
            "calmd-default-sam-text",
            (
                "calmd",
                "--no-PG",
                str(dat / "mpileup.1.sam"),
                str(dat / "mpileup.ref.fa"),
            ),
        ),
        Case(
            "calmd-equals-sam-text",
            (
                "calmd",
                "--no-PG",
                "-e",
                str(dat / "mpileup.1.sam"),
                str(dat / "mpileup.ref.fa"),
            ),
        ),
        Case(
            "calmd-realign-equals-sam-text",
            (
                "calmd",
                "--no-PG",
                "-re",
                str(dat / "mpileup.1.sam"),
                str(dat / "mpileup.ref.fa"),
            ),
        ),
        Case(
            "calmd-drop-tags-text",
            (
                "calmd",
                "--no-PG",
                "-d",
                "{tmp}/calmd-drop.sam",
                "{tmp}/calmd-drop.fa",
            ),
            files=(
                ("calmd-drop.fa", b">chr1\nACGTACGT\n"),
                (
                    "calmd-drop.sam",
                    b"@SQ\tSN:chr1\tLN:8\n"
                    b"r1\t0\tchr1\t1\t60\t4M\t*\t0\t0\tACGT\t!!!!"
                    b"\tRG:Z:g1\tBQ:Z:!!!!\tXX:i:7\tMD:Z:4\tNM:i:0\n",
                ),
            ),
        ),
        Case(
            "calmd-max-nm-text",
            (
                "calmd",
                "--no-PG",
                "-n2",
                "{tmp}/calmd-max-nm.sam",
                "{tmp}/calmd-max-nm.fa",
            ),
            files=(
                ("calmd-max-nm.fa", b">chr1\nACGTACGT\n"),
                (
                    "calmd-max-nm.sam",
                    b"@SQ\tSN:chr1\tLN:8\n"
                    b"low\t0\tchr1\t1\t60\t8M\t*\t0\t0\tACGTTCGT\tIIIIIIII"
                    b"\tNM:i:99\tMD:Z:0A0\n"
                    b"high\t0\tchr1\t1\t60\t8M\t*\t0\t0\tACGTTCGA\tIIIIIIII"
                    b"\tNM:i:99\tMD:Z:0A0\n",
                ),
            ),
        ),
        Case(
            "calmd-quiet-max-nm-text",
            (
                "calmd",
                "--no-PG",
                "-Q",
                "-n2",
                "{tmp}/calmd-quiet-max-nm.sam",
                "{tmp}/calmd-quiet-max-nm.fa",
            ),
            files=(
                ("calmd-quiet-max-nm.fa", b">chr1\nACGTACGT\n"),
                (
                    "calmd-quiet-max-nm.sam",
                    b"@SQ\tSN:chr1\tLN:8\n"
                    b"low\t0\tchr1\t1\t60\t8M\t*\t0\t0\tACGTTCGT\tIIIIIIII"
                    b"\tNM:i:99\tMD:Z:0A0\n"
                    b"high\t0\tchr1\t1\t60\t8M\t*\t0\t0\tACGTTCGA\tIIIIIIII"
                    b"\tNM:i:99\tMD:Z:0A0\n",
                ),
            ),
        ),
        Case(
            "calmd-cap-mapq-text",
            (
                "calmd",
                "--no-PG",
                "-C40",
                "{tmp}/calmd-cap.sam",
                "{tmp}/calmd-cap.fa",
            ),
            files=(
                ("calmd-cap.fa", b">chr1\nACGTACGT\n"),
                (
                    "calmd-cap.sam",
                    b"@SQ\tSN:chr1\tLN:8\n"
                    b"perfect\t0\tchr1\t1\t60\t8M\t*\t0\t0\tACGTACGT\tIIIIIIII\n"
                    b"mismatch\t0\tchr1\t1\t60\t8M\t*\t0\t0\tACGTTCGT\tIIIIIIII\n"
                    b"softclip\t0\tchr1\t1\t60\t2S6M\t*\t0\t0\tTTACGTAC\tIIIIIIII\n",
                ),
            ),
        ),
        Case(
            "calmd-bin-quality-text",
            (
                "calmd",
                "--no-PG",
                "-q",
                "{tmp}/calmd-bin-qual.sam",
                "{tmp}/calmd-bin-qual.fa",
            ),
            files=(
                ("calmd-bin-qual.fa", b">chr1\nACGTACGT\n"),
                (
                    "calmd-bin-qual.sam",
                    b"@SQ\tSN:chr1\tLN:8\n"
                    b"r1\t0\tchr1\t1\t60\t8M\t*\t0\t0\tACGTTCGT\t!+5?IS]g"
                    b"\tNM:i:99\tMD:Z:0A0\n",
                ),
            ),
        ),
        Case(
            "calmd-no-md-nm-update-text",
            (
                "calmd",
                "--no-PG",
                "-N",
                "{tmp}/calmd-no-md-nm.sam",
                "{tmp}/calmd-no-md-nm.fa",
            ),
            files=(
                ("calmd-no-md-nm.fa", b">chr1\nACGTACGT\n"),
                (
                    "calmd-no-md-nm.sam",
                    b"@SQ\tSN:chr1\tLN:8\n"
                    b"r1\t0\tchr1\t1\t60\t8M\t*\t0\t0\tACGTTCGT\tIIIIIIII"
                    b"\tNM:i:99\tMD:Z:0A0\n",
                ),
            ),
        ),
        Case(
            "fastq-stdout",
            ("fastq", str(dat / "bam2fq.001.sam")),
        ),
        Case(
            "fastq-output-file",
            ("fastq", "-0", "{side_tmp}/reads.fq", str(dat / "bam2fq.001.sam")),
            output_files=("reads.fq",),
        ),
        Case(
            "import-paired-rg-text",
            (
                "import",
                "--no-PG",
                str(test / "bam2fq" / "1.1.fq.expected"),
                str(test / "bam2fq" / "1.2.fq.expected"),
                "-r",
                "ID:rgid",
            ),
        ),
        Case(
            "import-interleaved-aux-text",
            (
                "import",
                "--no-PG",
                str(test / "import" / "2.interleaved.fq"),
                "-T",
                "",
            ),
        ),
        Case(
            "cat-sam-input-error",
            ("cat", "--no-PG", str(dat / "view.001.sam")),
        ),
        Case(
            "reheader-sam-input-error",
            (
                "reheader",
                str(test / "reheader" / "hdr.sam"),
                str(dat / "view.001.sam"),
            ),
        ),
        Case(
            "unknown-command-error",
            ("notacommand",),
        ),
        Case(
            "view-missing-input-error",
            ("view", "{tmp}/missing.bam"),
        ),
        Case(
            "index-missing-input-error",
            ("index", "{tmp}/missing.bam"),
        ),
        Case(
            "sort-missing-input-error",
            ("sort", "{tmp}/missing.bam"),
        ),
        Case(
            "flagstat-missing-input-error",
            ("flagstat", "{tmp}/missing.bam"),
        ),
        Case(
            "idxstats-missing-input-error",
            ("idxstats", "{tmp}/missing.bam"),
        ),
        Case(
            "head-missing-input-error",
            ("head", "{tmp}/missing.bam"),
        ),
        Case(
            "quickcheck-missing-input-error",
            ("quickcheck", "{tmp}/missing.bam"),
        ),
        Case(
            "dict-missing-input-error",
            ("dict", "{tmp}/missing.fa"),
        ),
        Case(
            "faidx-missing-input-error",
            ("faidx", "{tmp}/missing.fa"),
        ),
        Case(
            "fqidx-missing-input-error",
            ("fqidx", "{tmp}/missing.fa"),
        ),
        Case(
            "checksum-missing-input-error",
            ("checksum", "{tmp}/missing.bam"),
        ),
        Case(
            "coverage-missing-input-error",
            ("coverage", "{tmp}/missing.bam"),
        ),
        Case(
            "depth-missing-input-error",
            ("depth", "{tmp}/missing.bam"),
        ),
        Case(
            "samples-missing-input-error",
            ("samples", "{tmp}/missing.bam"),
        ),
        Case(
            "addreplacerg-missing-input-error",
            ("addreplacerg", "-r", r"@RG\tID:foo\tSM:bar", "{tmp}/missing.bam"),
        ),
        Case(
            "reset-missing-input-error",
            ("reset", "{tmp}/missing.bam"),
        ),
        Case(
            "consensus-missing-input-error",
            ("consensus", "{tmp}/missing.bam"),
        ),
        Case(
            "bedcov-missing-input-error",
            ("bedcov", "{tmp}/regions.bed", "{tmp}/missing.bam"),
            files=(("regions.bed", b""),),
        ),
        Case(
            "ampliconclip-missing-input-error",
            (
                "ampliconclip",
                "-b",
                "{tmp}/amplicons.bed",
                "{tmp}/missing.bam",
                "-o",
                "{side_tmp}/out.bam",
            ),
            files=(("amplicons.bed", b"chr1\t0\t10\tamp1\t0\t+\n"),),
        ),
        Case(
            "ampliconstats-missing-input-error",
            ("ampliconstats", "{tmp}/amplicons.bed", "{tmp}/missing.bam"),
            files=(("amplicons.bed", b"chr1\t0\t10\tamp1\t0\t+\n"),),
        ),
        Case(
            "collate-missing-input-error",
            ("collate", "-O", "{tmp}/missing.bam"),
        ),
        Case(
            "cat-missing-input-error",
            ("cat", "{tmp}/missing.bam"),
        ),
        Case(
            "split-missing-input-error",
            ("split", "{tmp}/missing.bam"),
        ),
        Case(
            "stats-missing-input-error",
            ("stats", "{tmp}/missing.bam"),
        ),
        Case(
            "fixmate-missing-input-error",
            ("fixmate", "{tmp}/missing.bam", "{side_tmp}/out.bam"),
        ),
        Case(
            "markdup-missing-input-error",
            ("markdup", "{tmp}/missing.bam", "{side_tmp}/out.bam"),
        ),
        Case(
            "rmdup-missing-input-error",
            ("rmdup", "{tmp}/missing.bam", "{side_tmp}/out.bam"),
        ),
        Case(
            "calmd-missing-input-error",
            ("calmd", "{tmp}/missing.bam", "{tmp}/ref.fa"),
            files=(("ref.fa", b">ref\nACGT\n"),),
        ),
        Case(
            "fastq-missing-input-error",
            ("fastq", "{tmp}/missing.bam"),
        ),
        Case(
            "fasta-missing-input-error",
            ("fasta", "{tmp}/missing.bam"),
        ),
        Case(
            "bam2fq-missing-input-error",
            ("bam2fq", "{tmp}/missing.bam"),
        ),
        Case(
            "reheader-missing-header-error",
            ("reheader", "{tmp}/missing.hdr", "{tmp}/missing.bam"),
        ),
        Case(
            "reheader-missing-input-error",
            ("reheader", "{tmp}/header.sam", "{tmp}/missing.bam"),
            files=(("header.sam", b"@HD\tVN:1.6\n"),),
        ),
        Case(
            "import-missing-input-error",
            ("import", "{tmp}/missing.fq"),
        ),
        Case(
            "merge-missing-input-error",
            ("merge", "{side_tmp}/out.bam", "{tmp}/missing.bam"),
        ),
        Case(
            "mpileup-missing-input-error",
            ("mpileup", "{tmp}/missing.bam"),
        ),
        Case(
            "reference-missing-input-error",
            ("reference", "{tmp}/missing.bam"),
        ),
        Case(
            "cram-size-missing-input-error",
            ("cram-size", "{tmp}/missing.cram"),
        ),
        Case(
            "phase-missing-input-error",
            ("phase", "{tmp}/missing.bam"),
        ),
        Case(
            "targetcut-missing-input-error",
            ("targetcut", "-f", "{tmp}/ref.fa", "{tmp}/missing.bam"),
            files=(("ref.fa", b">ref\nACGTACGT\n"),),
        ),
        Case(
            "tview-missing-input-error",
            ("tview", "{tmp}/missing.bam"),
        ),
        Case(
            "depad-missing-input-error",
            ("depad", "{tmp}/missing.bam"),
        ),
        Case(
            "quickcheck-bad-header-error",
            ("quickcheck", str(quickcheck / "2.quickcheck.badheader.bam")),
        ),
        Case(
            "view-no-sq-header-error",
            ("view", "{tmp}/no-sq.sam"),
            files=(
                (
                    "no-sq.sam",
                    b"@HD\tVN:1.6\n"
                    b"r1\t0\tchr1\t1\t60\t1M\t*\t0\t0\tA\t!\n",
                ),
            ),
        ),
        Case(
            "quickcheck-missing-eof-bitmask",
            ("quickcheck", str(quickcheck / "1.quickcheck.badeof.bam")),
        ),
        Case(
            "quickcheck-truncated-cram-bitmask",
            ("quickcheck", str(quickcheck / "9.quickcheck.cram30.truncated.cram")),
        ),
    ]


def run(
    binary: Path,
    args: tuple[str, ...],
    cwd: Path,
    tmp: Path,
    extra_env: tuple[tuple[str, str], ...],
) -> subprocess.CompletedProcess[bytes]:
    tmp.mkdir(parents=True, exist_ok=True)
    env = os.environ.copy()
    env["TMPDIR"] = str(tmp)
    env.update(extra_env)
    return subprocess.run(
        [str(binary), *args],
        cwd=cwd,
        env=env,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        check=False,
    )


def run_external(
    args: tuple[str, ...],
    cwd: Path,
    tmp: Path,
    extra_env: tuple[tuple[str, str], ...],
) -> subprocess.CompletedProcess[bytes]:
    tmp.mkdir(parents=True, exist_ok=True)
    env = os.environ.copy()
    env["TMPDIR"] = str(tmp)
    env.update(extra_env)
    stdout = subprocess.PIPE
    if len(args) >= 2 and args[-2] == ">":
        stdout = open(args[-1], "wb")
        args = args[:-2]
    try:
        return subprocess.run(
            args,
            cwd=cwd,
            env=env,
            stdout=stdout,
            stderr=subprocess.PIPE,
            check=False,
        )
    finally:
        if hasattr(stdout, "close"):
            stdout.close()


def unified_diff(label: str, left: bytes, right: bytes) -> str:
    left_text = left.decode("utf-8", errors="replace").splitlines(keepends=True)
    right_text = right.decode("utf-8", errors="replace").splitlines(keepends=True)
    return "".join(
        difflib.unified_diff(
            left_text,
            right_text,
            fromfile=f"c/{label}",
            tofile=f"rust/{label}",
        )
    )


def strip_lines_containing(data: bytes, needles: tuple[str, ...]) -> bytes:
    if not needles:
        return data

    needle_bytes = tuple(needle.encode() for needle in needles)
    return b"".join(
        line
        for line in data.splitlines(keepends=True)
        if not any(needle in line for needle in needle_bytes)
    )


def compare_case(case: Case, c_samtools: Path, rust_samtools: Path, root: Path, keep_tmp: bool) -> bool:
    tmp_ctx = None if keep_tmp else tempfile.TemporaryDirectory(prefix=f"samtools-rs-byte-parity-{case.name}-")
    tmp = Path(
        tempfile.mkdtemp(prefix=f"samtools-rs-byte-parity-{case.name}-")
        if keep_tmp
        else tmp_ctx.name
    )
    try:
        tmp_marker = str(tmp).encode()
        for relative_path, content in case.files:
            path = tmp / relative_path
            path.parent.mkdir(parents=True, exist_ok=True)
            path.write_bytes(content.replace(b"{tmp}", tmp_marker))

        for setup_template in case.setup_external:
            setup_args = tuple(
                arg.format(root=root, tmp=tmp, side_tmp=tmp) for arg in setup_template
            )
            setup_result = run_external(setup_args, root, tmp / "setup", case.env)
            if setup_result.returncode != 0:
                print(f"{case.name}: external setup failed")
                print(unified_diff(f"{case.name}.setup.stdout", b"", setup_result.stdout or b""))
                print(unified_diff(f"{case.name}.setup.stderr", b"", setup_result.stderr))
                if keep_tmp:
                    print(f"  tmp: {tmp}")
                return False

        for setup_template in case.setup:
            setup_args = tuple(arg.format(tmp=tmp, side_tmp=tmp) for arg in setup_template)
            setup_result = run(c_samtools, setup_args, root, tmp / "setup", case.env)
            if setup_result.returncode != 0:
                print(f"{case.name}: setup failed")
                print(unified_diff(f"{case.name}.setup.stdout", b"", setup_result.stdout))
                print(unified_diff(f"{case.name}.setup.stderr", b"", setup_result.stderr))
                if keep_tmp:
                    print(f"  tmp: {tmp}")
                return False

        c_side_tmp = tmp / "c-side"
        rust_side_tmp = tmp / "rust-side"
        c_side_tmp.mkdir(parents=True, exist_ok=True)
        rust_side_tmp.mkdir(parents=True, exist_ok=True)
        for relative_path, content in case.side_files:
            for side_tmp in (c_side_tmp, rust_side_tmp):
                path = side_tmp / relative_path
                path.parent.mkdir(parents=True, exist_ok=True)
                path.write_bytes(content.replace(b"{tmp}", tmp_marker))
        c_args_template = case.c_args if case.c_args is not None else case.args
        rust_args_template = case.rust_args if case.rust_args is not None else case.args
        c_args = tuple(arg.format(tmp=tmp, side_tmp=c_side_tmp) for arg in c_args_template)
        rust_args = tuple(
            arg.format(tmp=tmp, side_tmp=rust_side_tmp) for arg in rust_args_template
        )
        c_result = run(c_samtools, c_args, root, tmp / "c", case.env)
        rust_result = run(rust_samtools, rust_args, root, tmp / "rust", case.env)
        output_matches = True
        output_diffs: list[str] = []
        for relative_path in case.output_files:
            c_path = c_side_tmp / relative_path
            rust_path = rust_side_tmp / relative_path
            if not c_path.exists() or not rust_path.exists():
                output_matches = False
                output_diffs.append(
                    f"{case.name}.{relative_path}: "
                    f"c_exists={c_path.exists()} rust_exists={rust_path.exists()}"
                )
                continue
            c_output = c_path.read_bytes()
            rust_output = rust_path.read_bytes()
            c_output = strip_lines_containing(c_output, case.strip_lines_containing)
            rust_output = strip_lines_containing(rust_output, case.strip_lines_containing)
            if c_output != rust_output:
                output_matches = False
                output_diffs.append(unified_diff(f"{case.name}.{relative_path}", c_output, rust_output))
        for relative_path, render_template in case.rendered_output_files:
            c_path = c_side_tmp / relative_path
            rust_path = rust_side_tmp / relative_path
            if not c_path.exists() or not rust_path.exists():
                output_matches = False
                output_diffs.append(
                    f"{case.name}.{relative_path}: "
                    f"c_exists={c_path.exists()} rust_exists={rust_path.exists()}"
                )
                continue
            c_render_args = tuple(
                arg.format(root=root, tmp=tmp, side_tmp=c_side_tmp, output=c_path)
                for arg in render_template
            )
            rust_render_args = tuple(
                arg.format(root=root, tmp=tmp, side_tmp=rust_side_tmp, output=rust_path)
                for arg in render_template
            )
            c_render = run(c_samtools, c_render_args, root, tmp / "c-render", case.env)
            rust_render = run(c_samtools, rust_render_args, root, tmp / "rust-render", case.env)
            c_render_stdout = strip_lines_containing(
                c_render.stdout, case.strip_lines_containing
            )
            rust_render_stdout = strip_lines_containing(
                rust_render.stdout, case.strip_lines_containing
            )
            c_render_stderr = strip_lines_containing(
                c_render.stderr, case.strip_lines_containing
            )
            rust_render_stderr = strip_lines_containing(
                rust_render.stderr, case.strip_lines_containing
            )
            if c_render.returncode != rust_render.returncode:
                output_matches = False
                output_diffs.append(
                    f"{case.name}.{relative_path}.render.exit: "
                    f"c={c_render.returncode} rust={rust_render.returncode}"
                )
            if c_render_stdout != rust_render_stdout:
                output_matches = False
                output_diffs.append(
                    unified_diff(
                        f"{case.name}.{relative_path}.render.stdout",
                        c_render_stdout,
                        rust_render_stdout,
                    )
                )
            if c_render_stderr != rust_render_stderr:
                output_matches = False
                output_diffs.append(
                    unified_diff(
                        f"{case.name}.{relative_path}.render.stderr",
                        c_render_stderr,
                        rust_render_stderr,
                    )
                )
        c_stdout = strip_lines_containing(c_result.stdout, case.strip_lines_containing)
        rust_stdout = strip_lines_containing(rust_result.stdout, case.strip_lines_containing)
        c_stderr = strip_lines_containing(c_result.stderr, case.strip_lines_containing)
        rust_stderr = strip_lines_containing(rust_result.stderr, case.strip_lines_containing)
        ok = (
            c_result.returncode == rust_result.returncode
            and c_stdout == rust_stdout
            and (not case.compare_stderr or c_stderr == rust_stderr)
            and output_matches
        )
        if ok:
            print(f"{case.name}: ok")
            return True

        print(f"{case.name}: mismatch")
        if c_result.returncode != rust_result.returncode:
            print(f"  exit: c={c_result.returncode} rust={rust_result.returncode}")
        if c_stdout != rust_stdout:
            print(unified_diff(f"{case.name}.stdout", c_stdout, rust_stdout))
        if case.compare_stderr and c_stderr != rust_stderr:
            print(unified_diff(f"{case.name}.stderr", c_stderr, rust_stderr))
        for diff in output_diffs:
            print(diff)
        if keep_tmp:
            print(f"  tmp: {tmp}")
        return False
    finally:
        if tmp_ctx is not None:
            tmp_ctx.cleanup()


def main() -> int:
    root = repo_root()
    parser = argparse.ArgumentParser()
    parser.add_argument("--c-samtools", type=Path, default=root / "samtools" / "samtools")
    parser.add_argument("--rust-samtools", type=Path, default=root / "target" / "debug" / "samtools")
    parser.add_argument("--keep-tmp", action="store_true")
    args = parser.parse_args()

    cases = default_cases(root)
    failures = 0
    for case in cases:
        if not compare_case(case, args.c_samtools, args.rust_samtools, root, args.keep_tmp):
            failures += 1

    return 1 if failures else 0


if __name__ == "__main__":
    raise SystemExit(main())
