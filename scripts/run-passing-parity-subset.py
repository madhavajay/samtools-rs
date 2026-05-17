#!/usr/bin/env python3
"""Run the upstream samtools test harness for the stable CI subset.

The upstream `samtools/test/test.pl` script has no group-selection flag. This
helper keeps the vendored harness unmodified: it writes a temporary copy beside
`test.pl`, comments out top-level `test_*($opts...)` calls that are not in the
allow-list, then executes the filtered copy from `samtools/`.
"""

from __future__ import annotations

import argparse
import os
import re
import subprocess
from pathlib import Path


DEFAULT_GROUPS = [
    "test_reference",
    "test_dict",
    "test_fqidx",
    "test_sort",
    "test_collate",
    "test_calmd",
    "test_idxstat",
    "test_quickcheck",
    "test_head",
    "test_addrprg",
    "test_markdup",
    "test_bedcov",
    "test_split",
    "test_reset",
    "test_ampliconclip",
    "test_ampliconstats",
]


TOP_LEVEL_CALL = re.compile(r"^(test_[A-Za-z0-9_]+)\(\$opts(?:,.*)?\);\s*$")
SUMMARY_START = re.compile(r'^print "\\nNumber of tests:')


def repo_root() -> Path:
    return Path(__file__).resolve().parent.parent


def filtered_harness(src: Path, dst: Path, groups: set[str]) -> None:
    in_top_level = True
    output: list[str] = []

    for line in src.read_text().splitlines(keepends=True):
        if in_top_level and SUMMARY_START.match(line):
            in_top_level = False

        match = TOP_LEVEL_CALL.match(line.rstrip("\n"))
        if in_top_level and match and match.group(1) not in groups:
            output.append(f"# samtools-rs subset skip: {line}")
        else:
            output.append(line)

    dst.write_text("".join(output))


def prepare_sort_prereqs(root: Path) -> list[Path]:
    """Stage files that upstream test_sort expects test_index to create."""

    dat = root / "samtools" / "test" / "dat"
    bam = dat / "auto_indexed.tmp.bam"
    bai = dat / "auto_indexed.tmp.bam.bai"
    existed = {path: path.exists() for path in (bam, bai)}
    subprocess.run(
        [
            str(root / "samtools" / "samtools"),
            "view",
            "--write-index",
            "-o",
            str(bam),
            str(dat / "mpileup.1.sam"),
        ],
        cwd=root / "samtools",
        check=True,
    )
    return [path for path in (bam, bai) if not existed[path]]


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "groups",
        nargs="*",
        default=DEFAULT_GROUPS,
        help="test.pl group names to run; defaults to the CI-enforced stable subset",
    )
    parser.add_argument(
        "--samtools",
        type=Path,
        default=None,
        help="Rust samtools binary to stage at samtools/samtools before running",
    )
    args = parser.parse_args()

    root = repo_root()
    samtools_dir = root / "samtools"
    test_dir = samtools_dir / "test"
    src = test_dir / "test.pl"
    dst = test_dir / ".samtools-rs-passing-subset.pl"

    if args.samtools is not None:
        target = samtools_dir / "samtools"
        target.write_bytes(args.samtools.read_bytes())
        target.chmod(0o755)

    cleanup_paths: list[Path] = []
    if "test_sort" in args.groups:
        cleanup_paths.extend(prepare_sort_prereqs(root))

    filtered_harness(src, dst, set(args.groups))
    try:
        env = os.environ.copy()
        env["SAMTOOLS_RS_PARITY_SUBSET"] = ",".join(args.groups)
        return subprocess.run(
            ["perl", f"test/{dst.name}"],
            cwd=samtools_dir,
            env=env,
            check=False,
        ).returncode
    finally:
        dst.unlink(missing_ok=True)
        for path in cleanup_paths:
            path.unlink(missing_ok=True)


if __name__ == "__main__":
    raise SystemExit(main())
