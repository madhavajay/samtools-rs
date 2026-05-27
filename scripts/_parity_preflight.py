#!/usr/bin/env python3
"""Shared preflight guard for the parity/regression subset runners.

The upstream `repos/samtools/test/test.pl` (and several `regression.sh` data-gen
steps) need `bgzip` and `tabix` to build their inputs. When those tools are
absent from PATH, test.pl does **not** fail the affected groups -- it aborts
their data generation on `sh: bgzip: command not found` and treats the group
as skipped. The subset runner then exits 0 even though `test_view`,
`test_cat`, `test_reheader` (and others) never executed: a *false green* that
has previously masked real CRAM/parity regressions.

This guard makes that failure mode loud: by default the runners refuse to
start unless both `bgzip` and `tabix` are on PATH. Pass `--allow-missing-bgzip`
(or set `SAMTOOLS_RS_ALLOW_MISSING_BGZIP=1`) to deliberately run a degraded,
non-authoritative subset.
"""

from __future__ import annotations

import argparse
import os
import shutil
import sys

REQUIRED_TOOLS = ("bgzip", "tabix")

_ALLOW_ENV = "SAMTOOLS_RS_ALLOW_MISSING_BGZIP"


def add_preflight_arg(parser: argparse.ArgumentParser) -> None:
    """Register the shared opt-out flag on a runner's argument parser."""

    parser.add_argument(
        "--allow-missing-bgzip",
        action="store_true",
        default=os.environ.get(_ALLOW_ENV) == "1",
        help=(
            "run even if bgzip/tabix are missing from PATH. The result is "
            "NOT authoritative: groups whose data-gen needs bgzip are "
            "silently skipped by upstream test.pl, so a 0 exit does not "
            "mean those groups passed."
        ),
    )


def missing_tools() -> list[str]:
    return [tool for tool in REQUIRED_TOOLS if shutil.which(tool) is None]


def enforce(allow_missing: bool) -> None:
    """Abort with a clear message unless bgzip+tabix are usable.

    Exits the process (code 3) when tools are missing and the caller has
    not opted out, so a misconfigured environment can never masquerade as a
    passing run.
    """

    missing = missing_tools()
    if not missing:
        return

    joined = ", ".join(missing)
    if allow_missing:
        print(
            f"WARNING: missing on PATH: {joined}. Continuing because "
            "--allow-missing-bgzip was set. Groups whose data-gen needs "
            "these tools will be SILENTLY SKIPPED by upstream test.pl; a "
            "0 exit from this run is NOT authoritative.",
            file=sys.stderr,
        )
        return

    print(
        f"ERROR: required tools missing from PATH: {joined}.\n"
        "Upstream test.pl silently skips (does not fail) groups whose "
        "input generation needs these tools, so running now would produce "
        "a FALSE GREEN.\n"
        "Fix: build and expose the vendored copies, e.g.\n"
        "  make -C repos/htslib-rs/repos/htslib tabix bgzip\n"
        '  export PATH="$PWD/repos/htslib-rs/repos/htslib:$PATH"\n'
        "then re-run. To deliberately run a degraded, non-authoritative "
        "subset anyway, pass --allow-missing-bgzip (or set "
        f"{_ALLOW_ENV}=1).",
        file=sys.stderr,
    )
    raise SystemExit(3)
