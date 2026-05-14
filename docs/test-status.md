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
| `test_reference` | partial | SAM/BAM MD-tag reference reconstruction is implemented with `-o`, `-q`, basic `-r` region output, and indexed BAM region iteration when an associated BAI/CSI is present, covered by Rust integration tests. Upstream `test_reference` uses CRAM inputs and embedded-reference mode (`-e`), which remain blocked on CRAM all-record/container internals and full parity work. |
| `test_bgzip` | not-yet-ported | `bgzip` is an htslib tool, not currently in the samtools-rs binary scope. Decide whether to exclude from this parity run or add an htslib-rs CLI. |
| `test_faidx` | partial | `faidx` builds `.fai` and extracts local uncompressed regions, including `-r`, `-o`, `--length`, `faidx -f` FASTQ mode, reverse-complement `-i`, and mark-strand modes. BGZI, compressed output/indexing, and exact warning text remain. |
| `test_fqidx` | partial | FASTQ index build and local uncompressed region extraction exist, including reverse-complement `-i` and mark-strand modes. BGZI, compressed output/indexing, and exact warning text remain. |
| `test_dict` | passing | TODO marks `dict.out`, `dict.alias.out`, and `dict.alt.out` byte-parity verified; covered by `crates/samtools-rs/tests/dict.rs`. |
| `test_index` | partial | BAM/CSI/CRAI/SAM index build paths exist; threads, full `view -X`, merge/index interactions, and exact binary parity still need harness verification. |
| `test_mpileup` | blocked | `mpileup` is a stub. Needs `htslib-rs` pileup iterator API. |
| `test_usage` | partial | Top-level dispatcher/help exists. Full upstream usage text for every subcommand is not yet verified. |
| `test_view` | partial | SAM/BAM/reference-backed CRAM basics, stdin paths, count/header, filters, BED, region queries, and expression filtering are covered by `crates/samtools-rs/tests/view.rs`. Default `@PG` insertion applies to SAM-output paths (file/stdin SAM, BAM/CRAM stdin SAM, header-only) and is suppressed by `--no-PG`. `-U`/`-p` for SAM-input BAM and CRAM output is supported via a text → binary roundtrip. BAM/CRAM-output binary `@PG` insertion (requires `htslib-rs` header-injection support), BAM/CRAM-input binary aux mutation, BAM/CRAM-input `-U`/`-p`, multi-input, paired filters, and full CRAM parity remain. |
| `test_head` | partial | SAM/BAM modes, SAM/BAM/CRAM stdin, CRAM header-only modes, and reference-backed CRAM record extraction are covered by `crates/samtools-rs/tests/head.rs` plus unit tests in `commands::head`; full upstream harness parity still needs verification. |
| `test_cat` | partial | Record-level SAM and BAM concat exists with `-h` header replacement, default `@PG` insertion, `--no-PG`, and `-r region` for indexed BAM; BGZF fast path, CRAM, and `-p N/M` remain. |
| `test_import` | partial | Basic single FASTA/FASTQ and paired FASTQ (`-1`/`-2`, `--r1`/`--r2`, `-s` interleaved, plus two positional inputs) to SAM/BAM (`-O bam` / `--bam`) exists, including `-0` single-read alias, positional interleaved FASTQ detection from `/1`/`/2` read names, no-op `--no-PG`, CASAVA parsing (`-i`) with upstream-style reverse comments, SRA name2 (`-N`), UMI extraction with reverse comments, CASAVA barcode sequence tags, FASTQ definition aux tags (`-T`) including upstream-style float exponent spelling, explicit index reads (`--i1`/`--i2`) for `-0`, `-s`, positional interleaved, and paired `-1`/`-2` inputs with barcode sequence/quality tags and `-b`, and read-group header/tag support (`-R`/`-r`) with repeated `-r` accumulation, `-r` precedence over `-R`, and `-r` ID validation. Direct comparisons against `test/import/*.expected.sam` for the currently implemented import fixture commands pass; paired singleton/other grouping, full read-group parity, and CRAM output remain. |
| `test_bam2fq` | partial | Basic SAM/BAM FASTQ/FASTA output exists, including `-f`/`-F`/`-G` flag filters, `-0` as a single-stream output target, `-n`/`-N` read-name suffix controls, `-O` original-quality `OQ` tag output, and basic flag-driven `-1`/`-2`/`-s`/`-0` split outputs; exact paired grouping, barcode/index/tag handling, and CRAM remain. |
| `test_depad` | partial | SAM input with `-T` padded FASTA reference and `-s` SAM output matches the upstream `depad.001` fixture with `--no-PG`; BAM input/output, CRAM, binary output modes (`-u`/`-1`), and full upstream `test_depad` parity remain. |
| `test_stats` | partial | `SN` summary now covers basic counts plus supplementary alignments, insert size mean/standard deviation with `-i`/`--insert-size` capping and `-m`/`--most-inserts` bulk selection, inward/outward/other oriented pair counts, sequence-length lines, runtime `is sorted` for record-backed paths, record-backed `-I`/`--id` read-group/sample filtering, `-f`/`--required-flag`, `-F`/`--filtering-flag`, record-backed `-l`/`--read-length`, `-q`/`--trim-quality` BWA trim counting, percentage of properly paired reads, and target coverage percentage for `-g`/`--cov-threshold` with target-region validation. SAM/BAM and region-backed CRAM record paths emit FFQ/LFQ quality histograms, GCF/GCL GC histograms, and approximate CIGAR-walk COV coverage histograms with `-c`/`--coverage` bin ranges. SAM/indexed BAM/reference-backed CRAM positional regions and `-t` target files restrict the summary and COV positions, with overlapping BAM/CRAM regions de-duplicated. `-d` / `--remove-dups` filters duplicate-marked primary records. Exact pileup-backed COV parity, per-cycle, BAQ, and CRAM without explicit reference remain. |
| `test_depth` | partial | CIGAR-walk BAM/reference-backed CRAM depth exists, including `-r`, `-b`, `-a`/`-aa`, `-d`, `-q`, `-o`, `-H`, `-f` input lists, flag filters, `-l` minimum read length filtering, and multi-input depth columns. Exact pileup overlap/deletion behavior and CRAM without explicit reference remain. |
| `test_merge` | partial | In-memory merge exists for BAM and SAM inputs with default `@PG` insertion, `--no-PG`, `-R region` (indexed BAM only), and `-L BED` (indexed BAM only, with overlapping BED interval de-duplication); streaming merge, header reconciliation, CRAM, and parity details remain. |
| `test_sort` | partial | In-memory coordinate/name/tag sort exists for BAM and SAM inputs with default `@PG` insertion and `--no-PG`; external merge, template/minimiser sorts, memory/thread caps, and CRAM remain. |
| `test_collate` | partial | In-memory name grouping exists for BAM and SAM inputs with default `@PG` insertion and `--no-PG`; hash-bucket/on-disk mode, random seed, record cap, and CRAM remain. |
| `test_fixmate` | partial | Basic adjacent name-sorted BAM and SAM mate fixup exists with coordinate-sort rejection, TLEN recalculation, default MC/MQ tags, `-m` mate-score tags, and `-c` template-CIGAR `ct` tags; `-r` removes secondary/unmapped alignments and clears `PROPER_PAIR`/`MATE_REVERSE` on the surviving mate when its partner is unmapped. Mate rescore, sanitizer mutation, and CRAM remain. |
| `test_calmd` | partial | SAM input can recompute MD/NM tags against FASTA references, and BAQ paths for SAM input exist with default `@PG` insertion and `--no-PG`; BAM/CRAM I/O, remaining flags, and full upstream MD/BAQ parity remain. |
| `test_idxstat` | partial | BAM index counts exist, with streaming slow-path counts for SAM and unindexed BAM. CRAM slow path and full harness parity remain. |
| `test_quickcheck` | passing | TODO marks byte-for-byte parity against `quickcheck/all.expected`; covered by `crates/samtools-rs/tests/quickcheck.rs`. |
| `test_reheader` | partial | Basic BAM header replacement exists with default `@PG` insertion, `--no-PG` suppression, and BAM `-c <command>` header filtering; BGZF fast path and CRAM in-place remain. |
| `test_addrprg` | partial | SAM/BAM add/replace exists with `-O sam|bam`, default mode now matching upstream (`overwrite_all`), default `@PG` insertion, and `--no-PG`; CRAM, mate-aware behavior, and full orphan-first semantics remain. |
| `test_markdup` | partial | Single-end and paired-end markdup for SAM and BAM exists. SE records dedupe by `(tid, pos, reverse-flag)` plus optional barcode tag; PE pairs are grouped by qname and dedup by the canonical pair-of-coordinates key plus optional per-end barcode tag, with combined MAPQ as the score. Secondary/supplementary records with duplicate primary qnames inherit duplicate flags and are removed by `-r`. `-b`/`--barcode-tag`, `-c` duplicate flag/tag clearing, `-S` compatibility, `-t` duplicate-origin `do` tags, `-d` optical-distance duplicate classification with `dt:Z:SQ`/`dt:Z:LB` tags, default QCFAIL exclusion with `--include-fails` override, validated `-m t|s`/`--mode t|s` compatibility, optical-aware estimated library size in `-s` stats, `-r`, `-s` upstream-shaped summary fields, `-O`, `-o`, default `@PG` insertion, and `--no-PG` are supported. Exact upstream stats/count parity and CRAM remain. |
| `test_bedcov` | partial | CIGAR-walk BAM coverage, `-H` output headers, `-c` read-count columns, `-d` depth-threshold columns, `-g`/`-G` flag-mask controls, and `-j` deletion/refskip skipping exist; exact pileup behavior remains. |
| `test_split` | partial | BAM and SAM by read group exists with per-output `@RG` header filtering and default `@PG` insertion, plus explicit `-d TAG` string/integer aux-tag grouping with `-M`, explicit `-d RG` unknown-read-group header insertion, `-h` unaccounted SAM header override, `--no-PG`, and `--write-index` BAI generation for BAM outputs; CRAM, sorted-by-tag streaming, and deeper upstream `@PG` byte-parity for complex chains remain. |
| `test_large_positions` | partial | Some large-position behavior routes through `view`, `index`, `merge`, and `depth`; `depth -r` and `depth -b` exist for indexed BAM. Full harness group includes unported `tview` and parity-sensitive index/query cases. |
| `test_ampliconclip` | not-yet-ported | `ampliconclip` is a stub. |
| `test_ampliconstats` | not-yet-ported | `ampliconstats` is a stub. |
| `test_reset` | partial | BAM and SAM reset cover alignment-field clearing, default aux strip, reverse-strand re-reversal, `-x`/`--keep-tag`, `--dupflag`, `--no-RG`, `--reject-PG`, `--no-PG` semantics matching upstream's "skip new samtools @PG entry while preserving existing ones", and default `@PG` insertion via the shared helper. CRAM I/O remains. |
| `test_checksum` | partial | Default SAM/BAM checksum output is implemented and Rust-tested against upstream `checksum/chk1.1.expected` and `checksum/chk1.3.expected` after the harness' path-line normalization. `-T` TSV output, `-O` order-specific hashing, `-P` position columns, `-C` CIGAR columns, `-M` mate columns, `-B` bamseqchksum-compatible formatting, `-a` field-selection shorthand, wildcard/exclusion scalar/string/array aux tags with canonical integer encoding, and `-m` merging work for default/position/CIGAR/mate-column checksum reports. CRAM input, sanitizer mutation for exact `-a` parity, and full upstream group parity remain. |
| `test_coverage` | partial | BAM/reference-backed CRAM CIGAR-walk coverage exists, including `-r`, `--min-depth`, `-Q`/`--min-BQ`, read map-quality filtering, `-b`/`--bam-list`, `--ff`/`--excl-flags`, `--rf`/`--incl-flags`, `-l`/`--min-read-len`, `-d` maximum-depth capping, multi-input aggregate rows, mean depth, mean base quality, mean map quality, and a basic ASCII histogram output mode (`-m`/`--histogram` with `-A`/`--ascii` and `-D` routed through it, `-w`/`--n-bins` for column count). Byte-parity with upstream's UTF-8 + sidebar histogram remains, as does CRAM without explicit reference. |

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
