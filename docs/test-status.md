# samtools test.pl Status

This tracks the upstream `samtools/test/test.pl` groups against the Rust
`samtools-rs-cli` binary. The CI parity job currently runs the harness as a
regression watch and intentionally does not propagate its exit code. Each group
must move to `passing` before the parity gate can become required.

Status values:

- `passing`: byte-for-byte parity has been verified through the upstream Perl
  harness, or the listed subset is known to pass and the remaining cases are
  separately marked.
- `partial`: the command exists and has Rust integration coverage, but the
  upstream harness group is not fully passing.
- `not-yet-ported`: the command is still a stub or depends on unimplemented
  infrastructure.
- `blocked`: the command depends on missing `htslib-rs` support.

## Current Gate

- Rust gate: `cargo fmt --all --check`, `cargo clippy --workspace --all-targets -- -D warnings`,
  and `cargo test --workspace` are expected to pass.
- Parity gate: `.github/workflows/ci.yml` stages the Rust binary at the ignored
  `samtools/samtools` path and runs `cd samtools && perl test/test.pl || true`.
  Remove `|| true` only after the rows below are all `passing` or explicitly
  skipped with documented rationale.

## Harness Groups

| `test.pl` group | Status | Evidence / next work |
| --- | --- | --- |
| `test_reference` | not-yet-ported | `samtools reference` is a stub. Needs MD-tag/reference reconstruction work. |
| `test_bgzip` | not-yet-ported | `bgzip` is an htslib tool, not currently in the samtools-rs binary scope. Decide whether to exclude from this parity run or add an htslib-rs CLI. |
| `test_faidx` | partial | `faidx` builds `.fai` and extracts local uncompressed regions, including `-r`, `-o`, `--length`, `faidx -f` FASTQ mode, reverse-complement `-i`, and mark-strand modes. BGZI, compressed output/indexing, and exact warning text remain. |
| `test_fqidx` | partial | FASTQ index build and local uncompressed region extraction exist, including reverse-complement `-i` and mark-strand modes. BGZI, compressed output/indexing, and exact warning text remain. |
| `test_dict` | passing | TODO marks `dict.out`, `dict.alias.out`, and `dict.alt.out` byte-parity verified; covered by `crates/samtools-rs/tests/dict.rs`. |
| `test_index` | partial | BAM/CSI/CRAI/SAM index build paths exist; threads, full `view -X`, merge/index interactions, and exact binary parity still need harness verification. |
| `test_mpileup` | blocked | `mpileup` is a stub. Needs `htslib-rs` pileup iterator API. |
| `test_usage` | partial | Top-level dispatcher/help exists. Full upstream usage text for every subcommand is not yet verified. |
| `test_view` | partial | SAM/BAM basics, count/header, filters, BED, some region queries, and CRAM header-only output exist. Aux mutation, stdin, multi-input, full expressions, paired filters, and CRAM record parity remain. |
| `test_head` | partial | SAM/BAM modes and CRAM header-only modes are covered by `crates/samtools-rs/tests/head.rs`; CRAM record extraction and stdin remain. |
| `test_cat` | partial | Record-level BAM concat exists; BGZF fast path, CRAM, `-p`, `-r`, and `--no-PG` remain. |
| `test_import` | partial | Basic single FASTA/FASTQ and paired FASTQ (`-1`/`-2`, `--r1`/`--r2`, `-s` interleaved, plus two positional inputs) to SAM/BAM (`-O bam` / `--bam`) exists, including `-0` single-read alias, positional interleaved FASTQ detection from `/1`/`/2` read names, no-op `--no-PG`, CASAVA parsing (`-i`) with upstream-style reverse comments, SRA name2 (`-N`), UMI extraction with reverse comments, CASAVA barcode sequence tags, FASTQ definition aux tags (`-T`) including upstream-style float exponent spelling, explicit index reads (`--i1`/`--i2`) for `-0`, `-s`, positional interleaved, and paired `-1`/`-2` inputs with barcode sequence/quality tags and `-b`, and read-group header/tag support (`-R`/`-r`) with repeated `-r` accumulation, `-r` precedence over `-R`, and `-r` ID validation. Direct comparisons against `test/import/*.expected.sam` for the currently implemented import fixture commands pass; paired singleton/other grouping, full read-group parity, and CRAM output remain. |
| `test_bam2fq` | partial | Basic SAM/BAM FASTQ/FASTA output exists, including `-f`/`-F`/`-G` flag filters, `-0` as a single-stream output target, `-n`/`-N` read-name suffix controls, and basic flag-driven `-1`/`-2`/`-s`/`-0` split outputs; exact paired grouping, barcode/index/tag handling, and CRAM remain. |
| `test_depad` | not-yet-ported | `depad` is a stub. |
| `test_stats` | partial | Basic `SN` summary exists; histograms, insert size, GC, coverage, per-cycle, BAQ, and region support remain. |
| `test_merge` | partial | In-memory BAM merge exists; streaming merge, header reconciliation, region restriction, CRAM, and parity details remain. |
| `test_sort` | partial | In-memory coordinate/name BAM sort exists; external merge, tag/template/minimiser sorts, memory/thread caps, write-index, and CRAM remain. |
| `test_collate` | partial | In-memory BAM name grouping exists; hash-bucket/on-disk mode, random seed, record cap, and CRAM remain. |
| `test_fixmate` | partial | Basic adjacent name-sorted BAM mate fixup exists; MC/MQ tags, `-r`/`-c`/`-m`, mate rescore, and CRAM remain. |
| `test_calmd` | partial | BAQ paths for SAM input exist; MD/NM recomputation, BAM/CRAM I/O, and remaining flags remain. |
| `test_idxstat` | partial | BAM index counts exist, with streaming slow-path counts for SAM and unindexed BAM. CRAM slow path and full harness parity remain. |
| `test_quickcheck` | passing | TODO marks byte-for-byte parity against `quickcheck/all.expected`; covered by `crates/samtools-rs/tests/quickcheck.rs`. |
| `test_reheader` | partial | Basic BAM header replacement exists; BGZF fast path, CRAM in-place, command filter, and `--no-PG` remain. |
| `test_addrprg` | partial | SAM text add/replace exists; BAM/CRAM aux mutation, `orphan_first`, and mate-aware behavior remain. |
| `test_markdup` | not-yet-ported | `markdup` is a stub. |
| `test_bedcov` | partial | CIGAR-walk BAM coverage, `-H` output headers, `-c` read-count columns, and `-d` depth-threshold columns exist; exact pileup behavior remains. |
| `test_split` | partial | BAM by read group exists; `%*`, CRAM, and tag-based grouping remain. |
| `test_large_positions` | partial | Some large-position behavior routes through `view`, `index`, `merge`, and `depth`; `depth -r` and `depth -b` exist for indexed BAM. Full harness group includes unported `tview` and parity-sensitive index/query cases. |
| `test_ampliconclip` | not-yet-ported | `ampliconclip` is a stub. |
| `test_ampliconstats` | not-yet-ported | `ampliconstats` is a stub. |
| `test_reset` | partial | BAM reset exists for alignment fields/default aux strip; reverse-strand re-reversal, `--reject-PG`, `--dupflag`, SAM, and CRAM remain. |
| `test_checksum` | not-yet-ported | `checksum` is a stub. |
| `test_coverage` | partial | Basic BAM CIGAR-walk coverage and `-r` region restriction exist; base-quality mean, histogram output, and CRAM remain. |

## Rust Integration Coverage

Current Rust integration tests live under `crates/samtools-rs/tests/`:

- `dict.rs`
- `flags.rs`
- `head.rs`
- `misc.rs`
- `quickcheck.rs`
- `sort_merge.rs`
- `stats_wave_d.rs`
- `view.rs`

These are useful development checks, but they do not replace the upstream
`test.pl` parity rows above.
