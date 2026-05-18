#!/usr/bin/env python3
"""Run stable upstream regression.sh files against the Rust samtools binary."""

from __future__ import annotations

import argparse
import subprocess
from dataclasses import dataclass
from pathlib import Path

import _parity_preflight


@dataclass(frozen=True)
class Regression:
    directory: str
    regfile: str
    cleanup_globs: tuple[str, ...] = ()


DEFAULT_REGRESSIONS = [
    Regression("test/consensus", "consensus.reg", ("*.bam.bai",)),
    Regression("test/cram_size", "cram_size.reg"),
]


def repo_root() -> Path:
    return Path(__file__).resolve().parent.parent


def parse_regression(value: str) -> Regression:
    if ":" not in value:
        raise argparse.ArgumentTypeError("expected DIR:REGFILE")
    directory, regfile = value.split(":", 1)
    if not directory or not regfile:
        raise argparse.ArgumentTypeError("expected non-empty DIR:REGFILE")
    return Regression(directory, regfile)


def stage_binary(root: Path, samtools: Path | None) -> None:
    if samtools is None:
        return

    target = root / "samtools" / "samtools"
    target.write_bytes(samtools.read_bytes())
    target.chmod(0o755)


def cleanup(workdir: Path, patterns: tuple[str, ...]) -> None:
    for pattern in patterns:
        for path in workdir.glob(pattern):
            if path.is_file() or path.is_symlink():
                path.unlink()


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "regressions",
        nargs="*",
        type=parse_regression,
        default=DEFAULT_REGRESSIONS,
        help="regression files as DIR:REGFILE; defaults to the stable CI subset",
    )
    parser.add_argument(
        "--samtools",
        type=Path,
        default=None,
        help="Rust samtools binary to stage at samtools/samtools before running",
    )
    _parity_preflight.add_preflight_arg(parser)
    args = parser.parse_args()

    _parity_preflight.enforce(args.allow_missing_bgzip)

    root = repo_root()
    samtools_dir = root / "samtools"
    stage_binary(root, args.samtools)

    for regression in args.regressions:
        workdir = samtools_dir / regression.directory
        try:
            result = subprocess.run(
                ["../regression.sh", regression.regfile],
                cwd=workdir,
                check=False,
            )
            if result.returncode != 0:
                return result.returncode
        finally:
            cleanup(workdir, regression.cleanup_globs)

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
