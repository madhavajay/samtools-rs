# TODO: Port samtools to Pure Rust

Goal: build a pure Rust replacement for the `samtools` C program with full subcommand parity, then port and pass the upstream `test/test.pl` suite plus add Rust-native unit/integration tests. Implementation routes through `htslib-rs` (sibling submodule); when a needed HTSlib API is not yet exposed there, extend `htslib-rs` first.

## Progress Snapshot

**Phases 0–2 complete; Wave A complete (with partials); Wave B substantially complete (with partials); Wave D partially started; Wave C + Phases 4–5 pending and blocked on htslib-rs infrastructure (aux-tag mutation, pileup APIs).**

Subcommands shipped (27 of ~40):
- ✅ byte-parity verified: `flags`, `quickcheck`, `dict`
- ✅ functional with partial-feature notes: `head`, `index`, `idxstats`, `samples` (incl. `-i`, `-f`, `-X`, stdin path lists), `flagstat`, `faidx`/`fqidx`
- 🟡 partial implementation: `view` (SAM↔SAM, SAM→BAM/CRAM, count/header, region queries, `-f/-F/-G/-q` filters, `-L` BED, `-e` for SAM, `-x/--keep-tag` aux strip), `cat`, `reheader`, `fastq`/`fasta`/`bam2fq`, `split`, `sort` (in-memory), `merge` (in-memory), `collate` (in-memory), `import`, `rmdup` (single-end), `bedcov` (CIGAR-walk), `coverage`, `depth`, `addreplacerg` (SAM text mode), `fixmate` (name-sorted BAM), `reset` (alignment field clear + default aux strip), `stats` (SN summary numbers), `calmd`/`fillmd` (BAQ paths)

Remaining subcommands and their blockers:
- **BAM aux-tag mutation is reachable** via `bam::io::Reader::read_record_buf` → `RecordBuf` (mutable) → `bam::io::Writer::write_alignment_record`. This unblocks BAM forms of `addreplacerg`, `fixmate`, `reset`, `markdup`, `calmd`, `ampliconclip`. The work is still substantial per-subcommand.
- **Need a pileup iterator in htslib-rs** (would unlock `mpileup`, `consensus`, `targetcut`, `phase`, `ampliconstats`, exact `bedcov`/`coverage`/`depth`).
- **Other complex algorithms still TODO**: `stats` (123k LOC), `checksum` (47k LOC), `reference`, `cram-size` (blocked on htslib-rs CRAM internals), `depad`.

htslib-rs extensions landed during this work:
- `AlignmentRecordSummary` accessors: `flags`, `reference_sequence_id`, `mate_reference_sequence_id`, `mapping_quality`
- `summarize_bam_records_from_path`
- BAM/SAM FASTA/FASTQ helpers: limit, flag-filter, suffix, split `-1`/`-2`/`-s`, and selected aux tag preservation paths.
- BAM/CRAM region and flag-filter writers: `write_bam_regions_from_path`, `write_bam_records_with_required_flags_from_path`, `write_cram_regions_as_bam_from_path_with_reference`, `write_cram_records_with_required_flags_as_bam_from_path_with_reference`
- FASTQ import helpers: paired FASTQ input, index FASTQ input, aux-tag allow-listing, barcode quality tags, and read group tag insertion.

Rust tests: 109 currently passing (quickcheck:6, flags:3, dict:1, view:4, head:3, sort_merge:4, misc:52, stats_wave_d:20, test_status:1, library/command unit tests:15) — `cargo test -p samtools-rs` green; `cargo fmt --check` and `cargo clippy -p samtools-rs --all-targets -- -D warnings` also clean.

## What's Next — Decision Points

### BioScript VNtyper Native API Status

The temporary VNtyper priority pass is complete and folded into the main TODO. The implemented native API surface in `crates/samtools-rs/src/native.rs` now covers the BAM path needed by BioScript VNtyper without shelling out:

P0 — required for the BioScript VNtyper BAM path:

- [x] `view_region(input_bam, region, output_bam, threads?, reference_fasta?)`: indexed BAM + BAI region slicing to BAM output for `samtools view -P -b input.bam chr:start-end -o sliced.bam`.
- [x] `view_bed(input_bam, bed_file, output_bam, threads?, reference_fasta?)`: BED `-L` slicing to BAM output for `samtools view -P -b input.bam -L regions.bed -o sliced.bam`.
- [x] `index(input_bam, output_bai?)`: BAM index wrapper with implicit and explicit `.bai` output paths.
- [x] `bam_to_fastq_pair(...)` / `fastq_native(...)`: optional integrated name-sort plus paired FASTQ output to `-1`, `-2`, `-0`, and `-s` paths, including `.fastq.gz` output support.
- [x] `depth(input_bam, region, include_zero, threads?)`: structured per-base depth values for a region.
- [x] `depth_summary(input_bam, region, include_zero, threads?)`: mean, median, min, max, and uncovered count for VNTR coverage QC.

P1 — needed for faithful upstream behavior:

- [x] `merge(output_bam, input_bams, force, threads?)`: in-memory BAM merge wrapper with `-f` overwrite semantics for sliced + unmapped read workflows.
- [x] `sort(input_bam, output_bam, by_name, threads?)`: in-memory coordinate sort and name sort wrappers. Coordinate-sorted output remains indexable.
- [x] `quickcheck(input_alignment, verbose)`: BAM/CRAM validation wrapper returning structured errors.

P2 — CRAM and edge compatibility:

- [x] CRAM `view_region` with required `reference_fasta`, associated CRAI/CSI index use, and BAM output compatible with `index`, `sort`, `merge`, and `fastq`.
- [x] `extract_unmapped_pairs(input_alignment, output_bam, flag = 12, threads?, reference_fasta?)`: BAM and reference-backed CRAM flag-filter extraction without shell pipes.
- [x] Tests cover BAM/CRAM region slicing, CRAM reference-backed `flag = 12` unmapped-pair extraction, FASTQ conversion, depth summaries, sort/merge/index, and quickcheck.

Follow-up polish after the MVP:

- [ ] Propagate `threads` into BGZF/noodles worker pools instead of accepting it as an API-compatible no-op in several wrappers.
- [ ] Replace in-memory `sort`/`merge` implementations with streaming/external algorithms for large BAMs.
- [ ] Deepen `@PG` parity for native-generated BAM headers.
- [ ] Add broader real-world CRAM fixtures when available from VNtyper/BioScript workflows.

Three roughly orthogonal directions, each substantial:

1. **Unblock the pileup-dependent subcommands.** Add a `bam_plp_*`-shaped iterator API to `htslib-rs::alignment_compat`. Unlocks: `mpileup`, `consensus`, `targetcut`, `phase`, `ampliconstats`, and the exact (byte-parity) versions of `depth`/`coverage`/`bedcov`. Estimate: days of careful design in `htslib-rs`, then each subcommand is its own piece on top.
2. **Deepen the existing 27 partials toward byte-parity.** Pick a subcommand (probably `view` as the anchor) and drive every flag combination from `test.pl` to byte-for-byte. This is where `@PG` insertion, full filter expression support, and CRAM I/O need to land. Estimate: each `test_<name>` group is hours-to-a-day.
3. **Wire the upstream `test.pl` as a gating CI run.** Currently parity-gate fires `|| true`. Flip it as subcommands land, using `docs/test-status.md` to track per-test status (passing / skipped / not-yet-ported / cosmetic-diff). Forces parity attention.

Pick a direction and the per-subcommand TODOs below get rearranged accordingly.

## Current Inputs

- `samtools/`: upstream C samtools source and test suite. 54 C files (~42k LOC), ~40 subcommands dispatched from `bamtk.c:227`, and a 4159-line Perl test harness (`samtools/test/test.pl`) with expected-output fixtures under `samtools/test/<subcommand>/`.
- `htslib-rs/`: sibling pure-Rust HTSlib compatibility workspace. Re-exports `noodles` and ships HTSlib-shaped adapters under `crates/htslib-rs/src/*_compat.rs`. Currently has 41 passing HTSlib test groups.

## Pinned Scope Decisions

The following are decided up front and shape every phase below:

- **Subcommands**: target full parity with all upstream subcommands except those explicitly deferred (see *Out of Scope*).
- **Layout**: workspace mirroring `htslib-rs`:
  - `crates/samtools-rs` — library, one module per subcommand
  - `crates/samtools-rs-cli` — the `samtools` binary (dispatch + main)
- **HTSlib API gaps**: when samtools-rs needs an HTSlib-shaped API that `htslib-rs` does not yet expose, add it to `htslib-rs` first and route through it. Do not bypass `htslib-rs` for HTSlib-shaped APIs. (Direct `noodles` use from samtools-rs is acceptable only for code that has no HTSlib analogue.)
- **Tests — two gates**:
  1. **Parity gate**: upstream `samtools/test/test.pl` run against the Rust binary. Expected outputs are the checked-in files under `samtools/test/`. Used as a regression gate in CI.
  2. **Rust unit/integration tests**: per-subcommand `tests/` under each crate using `cargo test`. Used for fine-grained development feedback and Rust-native edge cases.
- **Parity level**:
  - **Strict (byte-for-byte)**: BAM/CRAM/VCF/BCF binary outputs, FASTQ/FASTA outputs, BED/depth/coverage/stats text outputs, idxstats, flagstat, sort order, index file bytes, exit codes.
  - **Semantic**: `@PG` header lines (use `ignore_pg_header` where needed, otherwise emit a `@PG` that matches the upstream `ID:samtools VN:<version> CL:<...>` shape — see *@PG strategy* below), stderr error messages (same key information, wording may differ), usage/help text.
- **C oracle**: local dev only. Devs MAY build upstream `samtools` (in `samtools/`) and use `test.pl --redo-outputs` to refresh expected fixtures. CI does NOT build or run C samtools — it only diffs against the checked-in fixtures.
- **Binary name**: `samtools` (the upstream test harness invokes `samtools` by default; we pass `-e samtools=<path>` to point it at our Rust build).

## Out of Scope (deferred)

- `tview` subcommand and the curses/HTML viewers (`bam_tview*.c`).
- Remote I/O backends: `https://`, `s3://`, `ftp://`, `gs://`. Local-file paths only. (`htslib-rs` also defers these.)
- `misc/` programs (`wgsim`, `md5fa`, `ace2sam`, `maq2sam-*`) and `misc/` Perl/Python scripts. These are auxiliary, not part of the samtools binary.
- `lz4/` vendored library — only needed by `bam_sort.c`'s temp-file compression path. Use a Rust crate (`lz4_flex` or similar) when sort temp-file compatibility is needed.
- `test/maintainer/` checks — these are tied to the C build/release process.
- C ABI exposure. samtools-rs is a Rust binary, not a library callable from C.
- Examples under `samtools/examples/`.

## Porting Principles

- Stay pure Rust. No `bindgen`, no `cc` crate, no linking to HTSlib C.
- Default to `htslib-rs` for HTSlib-shaped helpers (header manipulation, aux tags, pileup, region parsing, format detection, BGZF, index I/O). When `htslib-rs` lacks the API, file a task in `htslib-rs/TODO.md` and add it there before consuming from samtools-rs.
- Preserve observable behavior under the parity rules above. Treat each `test.pl` test case as an acceptance test; do not mark a subcommand complete until both its `test.pl` cases and its Rust integration tests pass.
- Each subcommand is one module under `crates/samtools-rs/src/commands/<name>.rs`, exposing `pub fn main(args: &[OsString]) -> ExitCode` (or similar). The CLI crate dispatches on `argv[1]` exactly like `bamtk.c:246`.
- Use `clap` for arg parsing but configure it to accept upstream's flag forms (short flags, long flags, value layout). Aliases and synonyms (`stat`/`stats`, `flag`/`flags`, `fastq`/`fasta`/`bam2fq`, `idxstat`/`idxstats`, `pad2unpad`/`depad`, `bamshuf`/`collate`) must be preserved.
- Errors: prefer `Result<T, E>` internally with a samtools-rs error type; surface via `print_error` / `print_error_errno` equivalents that match upstream's "[subcommand] message" stderr format.

## @PG Strategy

Upstream samtools writes a `@PG ID:samtools PN:samtools VN:<version> CL:<command-line>` line into output SAM/BAM/CRAM headers (`sam_view.c:1412` and equivalent in other subcommands). To stay close to byte parity:

- Emit `@PG ID:samtools PN:samtools VN:<samtools-upstream-version> CL:<reconstructed command-line>` where the VN matches the upstream samtools version we are tracking (pin this in `version.rs`). PP chaining must match upstream.
- The CL string must be reconstructed from `argv` using upstream's quoting rules (see `sam_hdr_add_pg` in HTSlib `header.c`).
- Where exact VN match is impossible or tests assert on the literal string, set `ignore_pg_header => 1` in the test harness invocation (already supported in `test.pl`).

## Phase 0: Workspace Skeleton

- [x] Create root `Cargo.toml` workspace mirroring `htslib-rs/Cargo.toml`:
  - members: `crates/samtools-rs`, `crates/samtools-rs-cli`
  - workspace deps include `htslib-rs = { path = "htslib-rs/crates/htslib-rs" }`
  - rust-version + edition matched to htslib-rs (`1.89.0`, `2024`)
- [x] Create the two crate skeletons with empty `lib.rs` / `main.rs` and a placeholder dispatcher.
- [x] Wire up `cargo fmt`, `cargo clippy --workspace --all-targets -- -D warnings`, `cargo test --workspace` as the Rust gate.
- [x] Add a top-level GitHub Actions workflow with two jobs:
  - Rust gate (fmt + clippy + test).
  - Parity gate: build release binary, then run `cd samtools/test && perl test.pl -e samtools=$WORKSPACE/target/release/samtools`.
- [x] Document the project goal, scope, and CI gates in `README.md`.

## Phase 1: Shared Infrastructure

These are used by nearly every subcommand and must exist before subcommands can land.

- [x] **CLI dispatcher** (`samtools-rs/src/dispatch.rs`): ports `bamtk.c:227-319`. Subcommand table with aliases. `samtools --version`, `samtools --version-only`, `samtools help [cmd]`. Unknown-subcommand error matches upstream wording.
- [x] **Version + feature string** (`samtools-rs/src/version.rs`): exports `SAMTOOLS_VERSION` constant (currently pinned to `1.23.1`, the upstream tag tracked by the `samtools/` submodule). `long_version` prints both samtools and htslib-rs versions.
- [ ] **Global args** (`samtools-rs/src/sam_global.rs`): port `sam_opts.h` / `sam_opts.c`. The shared `--input-fmt`, `--input-fmt-option`, `--output-fmt`, `--output-fmt-option`, `--reference`, `--threads`, `--write-index`, `--verbosity` long options. Re-usable clap fragment that subcommands inject. **(Currently each subcommand parses args locally; promote to a shared helper as Wave A grows.)**
- [ ] **Open/close helpers** (`samtools-rs/src/io.rs`):
  - `sam_open_format` / `sam_open_mode` equivalents over `htslib-rs::format` (extension- and content-based format detection).
  - `auto_index` — write an associated index alongside a writer (BAI/CSI/CRAI/TBI).
  - `autoflush_if_stdout` + `release_autoflush` semantics, so error paths flush stdout before logging.
  - `check_sam_close` — propagate close errors into exit code.
- [x] **`print_error` / `print_error_errno`** (`samtools-rs/src/diagnostics.rs`): subcommand-prefixed stderr printers matching upstream's `samtools <subcommand>: message` format and errno appending.
- [x] **BAM flag bits + parse/format** (`samtools-rs/src/bam_flag.rs`): the `BAM_F*` constants plus `bam_str2flag` / `bam_flag2str` equivalents. Used by `flags` and (later) `view` / `flagstat`.
- [x] **Raw header text extractor** (`samtools-rs/src/header_text.rs`): preserves original `@`-line order for SAM, BAM, and CRAM inputs. Required by `head`, `view -h`, `samples`, etc. — noodles' canonical writer reorders headers, which breaks byte parity.
- [ ] **BAM sanitizer** (`samtools-rs/src/sanitize.rs`): port `samtools.h`'s `FIX_*` flags and `bam_sanitize` / `bam_sanitize_options`. Used by sort/merge/view/markdup.
- [ ] **@PG add helper** (`samtools-rs/src/pg.rs`): builds `@PG ID:samtools VN:... CL:...` lines with PP chaining, mirroring `sam_hdr_add_pg`. Command-line reconstruction must replicate HTSlib's quoting.
- [x] **Aux-tag list parser** (`samtools-rs/src/aux_list.rs`): port `parse_aux_list` from `sam_utils.c`. Used by `view`, `reset`, `fastq`, and future aux-aware commands.
- [ ] **BED index** (`samtools-rs/src/bedidx.rs`): port `bedidx.c` (interval tree over BED regions). Used by `view -L`, `bedcov`, `ampliconclip`, `mpileup`.
- [ ] **Reference helpers** (`samtools-rs/src/reference.rs`): mmap FASTA / `--reference` handling shared by `calmd`, `consensus`, `mpileup`, `phase`, `import`.
- [ ] **Temp file helper** (`samtools-rs/src/tmp_file.rs`): port `tmp_file.c` (BAM record temp spooling with optional lz4 compression). Used by `sort`, `collate`.
- [ ] **Logging passthrough**: bridge to `htslib-rs::log_compat` so `--verbosity` flows correctly.

## Phase 2: Subcommand Surface Mapping

Mapping document exists at [`docs/subcommand-coverage.md`](docs/subcommand-coverage.md). It lists every subcommand, the HTSlib APIs it depends on, the `htslib-rs` coverage status, and a rolled-up list of extensions needed in `htslib-rs`.

- [x] Per-subcommand HTSlib API surface enumerated.
- [x] `htslib-rs` coverage status per API (already exposed / needs adding / out of scope).
- [x] Gap list rolled up at the end of `docs/subcommand-coverage.md`.

## Phase 3: Subcommand Implementation Waves

Each subcommand below maps to: (a) one Rust module under `crates/samtools-rs/src/commands/`, (b) `test_<name>` in `samtools/test/test.pl` passing against the Rust binary, (c) at least one Rust integration test under `crates/samtools-rs/tests/<name>.rs`.

The waves are ordered to land foundational machinery first (read/write/index) and unblock the rest.

### Wave A — Read/Write/Index Foundation

- [~] `view` (`sam_view.c`, 68k) — partial: SAM↔SAM passthrough, SAM→BAM/CRAM, header-only / count modes (including CRAM `-H`), `-h` `-H` `-c` `-b` `-C` `-T` `-o` `--no-PG`, filter flags `-f`/`-F`/`-G`/`-q`, region queries (`<chr:start-end>`), `-L FILE` BED restrict, `-e EXPR` filter expression (SAM input + count mode), `-O FORMAT` output-fmt option. **Pending:** aux-tag manipulation, unmapped output, multi-file inputs, stdin, full `-e EXPR` for BAM/CRAM, paired-aware filters, CRAM record decoding/parity.
- [x] `head` (`sam_view.c` shared) — SAM and BAM input; CRAM header-only modes; `-h N`, `-n N`, all-default. CRAM record extraction and stdin still TODO.
- [x] `quickcheck` (`bam_quickcheck.c`) — passes byte-for-byte against `quickcheck/all.expected`.
- [x] `index` (`bam_index.c`) — BAI/CSI/CRAI build, `-c` CSI mode, `--min-shift`, `-M`, `-o`, legacy `<in> <out.idx>` synopsis. **Pending:** `-@` threads not yet propagated to noodles workers.
- [x] `idxstats` (`bam_stat.c`) — index-based per-reference counts for BAM, with streaming slow-path counts for SAM and unindexed BAM. **Pending:** CRAM slow-path fallback.
- [~] `faidx` / `fqidx` (`faidx.c`) — index-build mode works (`samtools faidx file.fa` produces `file.fa.fai`); local uncompressed region extraction works for positional regions, `-r` region files, `-o`, `--length`, FASTQ mode via `fqidx` and `faidx -f`, reverse-complement `-i` with mark-strand modes, and `--continue`-style missing-region tolerance. **Pending:** BGZI support, compressed output, output indexing, and full warning text parity.
- [x] `dict` (`dict.c`) — sequence dictionary builder. Passes byte-for-byte against `dict.out`, `dict.alias.out`, `dict.alt.out` (run via test.pl-style stdin/file invocations).
- [x] `flagstat` / `flagstats` (`bam_stat.c`) — SAM and BAM input. Default + `-O json` + `-O tsv` output modes. Required extending `htslib-rs::alignment_compat::AlignmentRecordSummary` with `flags_u16` / `reference_sequence_id` / `mate_reference_sequence_id` / `mapping_quality` accessors, plus a new `summarize_bam_records_from_path`. **Pending:** CRAM input.

### Wave B — File Ops

- [~] `sort` (`bam_sort.c`, 138k — the largest single file) — basic in-memory coordinate sort and `-n` name sort for BAM. Sets `@HD SO` correctly. **Pending:** on-disk external merge for large inputs, tag sort (`-t`), template-coordinate sort (`-M`), minimiser sort (`-N`), thread/memory caps, write-index, CRAM.
- [~] `merge` (`bam_sort.c` shared) — basic in-memory multi-input BAM merge with coordinate or name sort. `-f` overwrite, `-o`/`--output-fmt`. **Pending:** k-way streaming merge, header merging for inputs with differing `@SQ`, region restriction, CRAM.
- [~] `collate` / `bamshuf` (`bamshuf.c`) — basic in-memory name-sort grouping for BAM. `-o PREFIX`, `-O` stdout, `--output-fmt sam|bam`. **Pending:** on-disk hash-bucket implementation for inputs larger than memory, `-r` random seed, `-n` in-memory record cap, CRAM.
- [~] `cat` (`bam_cat.c`) — basic BAM concatenation works (record-level decompress + re-encode). Supports `-o`, `-h` (header replacement). **Pending:** BGZF block-level fast path, CRAM, `-p N/M`, `-r`, `--no-PG` (currently silently ignored).
- [~] `split` (`bam_split.c`) — basic BAM-by-`@RG` splitting with `-f` template (`%!`, `%#`, `%.`), `-u` unaccounted, `--output-fmt sam|bam`, `-p N` padding. **Pending:** `%*` (RG SM lookup), CRAM, `-d`/`-D` tag-based grouping.
- [~] `reheader` (`bam_reheader.c`) — basic BAM header replacement (record-level rewrite). **Pending:** BGZF block-level fast path, CRAM `--in-place`, `-c <command>` external filter, `--no-PG` (currently silently ignored).
- [~] `addreplacerg` (`bam_addrprg.c`) — SAM→SAM text-level add/replace for `@RG` header line and `RG:Z` aux tag. `-r SPEC`, `-R ID`, `-m overwrite_all|orphan_only`, `-O sam`, `-o FILE`. **Pending:** BAM/CRAM input/output (blocked on aux mutation), `orphan_first` semantics, mate-aware updates.
- [~] `fastq` / `fasta` / `bam2fq` (`bam_fastq.c`) — basic single-stream output works for SAM and BAM (records written to stdout, `-o FILE`, or `-0 FILE`), with `-f`/`-F`/`-G` flag filters, read-name suffix controls (`-n`/`-N`), and basic flag-driven paired split outputs (`-1`/`-2`/`-s`/`-0`). **Pending:** barcode/index/tag handling and exact name-grouped paired/singleton/other semantics.
- [~] `import` (`bam_import.c`) — basic single FASTA/FASTQ and paired FASTQ (`-1`/`-2`, `--r1`/`--r2`, `-s` interleaved, plus two positional inputs) → SAM/BAM (`-O bam` / `--bam`), including positional single input plus `-0` single-read alias, positional interleaved FASTQ detection from `/1`/`/2` read names, no-op `--no-PG`, CASAVA parsing (`-i`) with upstream-style reverse comments, SRA name2 (`-N`), UMI extraction (`-U`/`--UMI-tag`) with reverse comments, CASAVA barcode sequence tags (`--barcode-tag`), FASTQ definition aux tags (`-T`) including upstream-style float exponent spelling, explicit index reads (`--i1`/`--i2`) for `-0`, `-s`, positional interleaved, and paired `-1`/`-2` inputs with barcode sequence/quality tags (`--barcode-tag`/`--quality-tag`) and `-b`, and read-group header/tag support (`-R`/`-r`) with repeated `-r` accumulation, `-r` precedence over `-R`, and `-r` ID validation. Direct comparisons against `test/import/*.expected.sam` for the currently implemented import fixture commands pass. **Pending:** paired singleton/other grouping, full read-group parity, CRAM output.

### Wave C — Editing / Mate-aware

- [~] `fixmate` (`bam_mate.c`) — basic mate flag/pos fixup for adjacent paired records in name-sorted BAM (`FMUNMAP`, `FMATE_REVERSE_COMPLEMENTED`, `mate_reference_sequence_id`, `mate_alignment_start`). **Pending:** MC/MQ aux tags, `-r`/`-c`/`-m` modes, CRAM, mate-rescore.
- [ ] `markdup` (`bam_markdup.c`, 89k) — duplicate marking (single + paired), barcode-aware, `-r`/`-l`/`-s` modes, opt distance.
- [~] `rmdup` (`bam_rmdup.c` + `bam_rmdupse.c`) — single-end duplicate removal (by `(tid, pos, reverse-flag)`, keeping highest MAPQ). **Pending:** paired-end mode, CRAM/SAM.
- [~] `calmd` / `fillmd` (`bam_md.c`) — BAQ paths (`-r`, `-r -e`, `-E`) wired through `htslib_rs::alignment_compat::recalculate_baq_*` and `apply_existing_baq_from_sam_path`. SAM input only. **Pending:** MD/NM tag recomputation against the reference (per-base diff), BAM/CRAM I/O, `-A`/`-d`/`-C cap`.
- [ ] `targetcut` (`cut_target.c`) — fosmid pool target cutting.
- [~] `reset` (`reset.c`) — strip alignment fields (`tid`/`pos`/`cigar`/`mapq`/`mate_*`/`template_length`), drop a default set of aligner aux tags (NM, MD, AS, XS, SA, MC, MQ, NH, HI, ms), clear `PROPER_PAIR`/`SECONDARY`/`SUPPLEMENTARY`/`DUP`/`MATE_UNMAPPED`/`REVERSE`/`MATE_REVERSE` flag bits, set `UNMAPPED`. `-x`/`--keep-tag` honored. **Pending:** reverse-strand seq/qual re-reversal, `--reject-PG`, `--dupflag`, SAM/CRAM I/O.
- [ ] `ampliconclip` (`bam_ampliconclip.c`, 40k) — primer clipping with BED amplicon spec.

### Wave D — Stats / Pileup

- [~] `depth` (`bam2depth.c`) — per-position depth via CIGAR walks for BAM. `-a`/`-aa`/`-d`/`-q`/`-o`, `-r` region restriction, and `-b` BED restriction are supported. **Pending:** pileup-based exact handling of overlaps/deletions, multi-input columnar output, CRAM.
- [~] `coverage` (`coverage.c`) — per-reference/`-r` region `numreads`, `covbases`, `coverage`, `meandepth`, `meanmapq` via CIGAR walks for BAM. **Pending:** `meanbaseq` (needs per-base qualities), histogram/ASCII-plot output modes, CRAM.
- [~] `bedcov` (`bedcov.c`) — total aligned-base coverage per BED region, walking each record's CIGAR. `-Q` mapq filter, `-H` output headers, `-c` read-count columns, and `-d` depth-threshold columns are supported. **Pending:** pileup-based exact coverage for byte parity.
- [~] `stats` (`stats.c` + `stats_isize.c`, 123k + 8k) — basic `SN` (Summary Numbers) section: raw total / filtered / 1st & last fragments / mapped / paired / properly paired / unmapped / duplicated / MQ0 / QC-failed / non-primary / singletons / diffchr pairs. **Pending:** FFQ/LFQ, IS (insert sizes), GCF, COV histogram, per-cycle, BAQ, region restriction.
- [ ] `mpileup` (`bam_plcmd.c`, 49k) — multi-way pileup with `htslib-rs::alignment_compat` pileup support. Output formats including VCF, depth-aware.
- [ ] `consensus` (`bam_consensus.c`, 126k, + `consensus_pileup.c`) — consensus FASTA/FASTQ/pileup builder.
- [ ] `phase` (`phase.c`) — heterozygote phasing.
- [ ] `depad` / `pad2unpad` (`padding.c`) — padded → unpadded BAM.
- [ ] `ampliconstats` (`amplicon_stats.c`, 65k) — per-amplicon stats over the same BED model as `ampliconclip`.
- [ ] `cram-size` (`cram_size.c`) — CRAM Content-ID and Data-Series byte breakdown. Depends on `htslib-rs` CRAM internals (currently out of scope in htslib-rs — see *htslib-rs Extensions Needed*).
- [ ] `checksum` (`bam_checksum.c`, 47k) — order-agnostic sequence-content checksums.
- [x] `samples` (`bam_samples.c`) — list `@RG SM:` samples across inputs. Header-driven dedup, `-T TAG`, `-o`, `-h`, `-i` index-presence column, `-f`/`-F` FASTA dictionary matching, stdin path lists, `-X` custom index pairs, and CRAM headers are implemented.
- [ ] `reference` (`reference.c`) — generate a reference from aligned data + MD tags.
- [x] `flags` (`bam_flags.c`) — explain a numeric BAM flag. Byte-for-byte parity with upstream.

## Phase 4: Test Harness Integration

- [ ] **Parity gate setup**: confirm `samtools/test/test.pl` can be driven via `-e samtools=<rust-binary-path>` without modifying the harness. Verify `regression.sh` still works.
- [ ] **`@PG` handling**: where upstream expected outputs include `@PG` lines with a specific VN that we cannot reproduce, set `ignore_pg_header => 1` in those tests. Avoid touching the actual expected output files.
- [x] **Status ledger**: `docs/test-status.md` tracks the upstream `test.pl` groups as `passing` / `partial` / `not-yet-ported` / `blocked`, including why CI still runs the parity harness with `|| true`.
- [ ] **Run progressively**: as each subcommand lands in Phase 3, enable its `test_<name>` in CI. Disabled tests should be tracked in `docs/test-status.md` as `not-yet-ported` (NOT just commented out).
- [ ] **Rust integration tests per subcommand**: under `crates/samtools-rs/tests/<name>.rs`, write at least: happy path, error path, region/`-L`/format-flag variants, threaded variant where applicable. These run on every PR independently of the Perl gate.
- [ ] **Compile-side test binaries**: `samtools/test/merge/test_bam_translate`, `test_rtrans_build`, `test_trans_tbl_init`, `samtools/test/split/test_*`, `samtools/test/vcf-miniview.c` — port to Rust integration tests under the relevant subcommand crate.

## Phase 5: Parity Polishing

- [ ] **Diff every `test_<name>` output byte-for-byte** against the C samtools outputs on a known fixture corpus (locally, dev-only). For each diff: classify (real bug / acceptable cosmetic / `@PG` only) and either fix or document.
- [ ] **Threads**: verify `-@ N` propagates to `htslib-rs`/noodles worker pools and matches upstream's parallelism behavior. Verify the `--threads` global flag short-circuits to BGZF workers where appropriate.
- [ ] **Exit codes**: confirm exit code matches upstream for invalid inputs, missing files, truncated BGZF, malformed CIGAR, etc.
- [ ] **Performance triage**: measure each subcommand on a representative dataset vs C samtools. Goal is "within 2x" initially; performance fixes come after parity.
- [ ] **Bench harness**: criterion or custom timing harness under `benches/` for `view`, `sort`, `markdup`, `stats`, `mpileup`.

## htslib-rs Extensions Needed (rolling list)

This list is filled in during Phase 2 as the subcommand surface mapping uncovers gaps. Each entry creates a tracked task in `htslib-rs/TODO.md`.

- [ ] **`sam_hdr_add_pg`** — programmatic `@PG` chain insertion with PP linkage (currently `htslib-rs::alignment_compat` exposes header read/write but the `@PG` chain helper is the samtools workhorse).
- [ ] **`bam_aux_update_*`** — string/int/array aux updates with re-sizing semantics.
- [ ] **`sam_pileup` / `bam_plp_*` API surface** — multi-input pileup iterator. `htslib-rs::alignment_compat` has fixture-level pileup but the iterator API surface needs auditing before mpileup/consensus land.
- [ ] **`hts_set_threads`** — wire-up to BGZF worker count for samtools `-@`.
- [ ] **`auto_index` / index save during write** — write BAI/CSI/CRAI alongside writer.
- [ ] **CRAM internals for `cram-size`** — currently marked out of scope in `htslib-rs/README.md`. Either expose a minimal block/container/codec inventory API in `htslib-rs`, or drop `cram-size` from samtools-rs scope.
- [ ] **`htslib-rs::region`** — confirm coverage of HTSlib's full region-string grammar including `*` (unmapped) and `.` (everything else).
- [ ] **`probaln_glocal` and BAQ recalculation** — exposed by `htslib-rs::probaln`; verify wiring for `calmd` and `mpileup`.

## Submodule Pinning

- [x] Pin `samtools/` to a specific upstream release tag once Phase 0 lands (record tag + commit in `README.md` and `version.rs`). Current pin: upstream tag `1.23.1`, commit `6efb9b6da35224cf804921dedecf9fb8f411365d`.
- [x] Pin `htslib-rs/` to a known-green commit when Phase 0 lands. Current pin: `88bd29f5f0d5e87d3f5d28da1f106a4b518e3926`.

## Repository Map (target end state)

- `crates/samtools-rs/` — library with one module per subcommand plus shared infra.
- `crates/samtools-rs-cli/` — the `samtools` binary.
- `samtools/` — upstream C source + tests, used as fixture and reference only.
- `htslib-rs/` — sibling Rust HTSlib compatibility workspace consumed via path dep.
- `docs/subcommand-coverage.md` — per-subcommand HTSlib API surface and `htslib-rs` coverage status.
- `docs/test-status.md` — per-test pass/skip/not-yet-ported status.
- `TODO.md` — this file.
- `README.md` — project overview, scope decisions, build/test instructions.

## Development Workflow

```sh
# clone with submodules
git clone --recurse-submodules <repo>

# rust gate
cargo fmt --all
cargo fmt --all --check
cargo clippy --workspace --all-targets -- -D warnings
cargo test --workspace

# parity gate (against checked-in expected outputs)
cargo build --release
cd samtools/test && perl test.pl -e samtools=$PWD/../../target/release/samtools

# optional: refresh expected outputs from C samtools (local dev only)
cd samtools && autoreconf -i && ./configure && make
cd test && perl test.pl --redo-outputs
```
