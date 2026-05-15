# TODO: Port samtools to Pure Rust

Goal: build a pure Rust replacement for the `samtools` C program with full subcommand parity, then port and pass the upstream `test/test.pl` suite plus add Rust-native unit/integration tests. Implementation routes through `htslib-rs` (sibling submodule); long-term, when a needed HTSlib API is not yet exposed there, extend `htslib-rs` first.

## Active Goal

For the current pass, keep working only in `samtools-rs`. Do not change the
`htslib-rs` or noodles submodules. If a TODO item requires an underlying-library
change, move or keep it at the end of this file under **htslib-rs Extensions
Needed** or **noodles Extensions Needed**, skip it for now, and continue with the
next samtools-only item. Stop once the only remaining actionable work requires
`htslib-rs` or noodles changes, and report those blockers instead of modifying
the underlying libraries.

## Current Handoff — 2026-05-15

Merged baseline:
- PR #7, https://github.com/madhavajay/samtools-rs/pull/7, is merged into `main`.
- Merge commit: `ae15e4ad603892912fe0c5175491e1d0e3f210eb`.
- The next work should start from `main` (or one of the in-flight PR branches below).

In-flight PRs (samtools-rs-only, ready for review/merge):
- PR #8, https://github.com/madhavajay/samtools-rs/pull/8, branch `fastq-index-files` — `fastq --i1 FILE` / `--i2 FILE` per-record index FASTQ extraction with `--index-format` (default `i*i*`), `--quality-tag` (default `QT`), and `--barcode-tag`. Emits one index record per primary non-READ2 record; exact upstream name-grouped one-record-per-qname-pair emission remains pending.
- PR #9, https://github.com/madhavajay/samtools-rs/pull/9, branch `fastq-index-paired-grouping` — five stacked commits:
  1. Upstream-style name-grouped fastq split routing (paired R1+R2 → `-1`/`-2`; R1-only / R2-only singleton → `-s` with fallback to `-1`/`-2`; READ_OTHER → `-0` with fallback to `-s`).
  2. Accumulating `-t` / `-T` aux-tag selections (union rather than override).
  3. Per-record interleaved output when `-1` and `-2` paths alias to the same file.
  4. Repeated `-d` / `-D` value-union for the same tag with mismatched-tag rejection.
  5. Route all FASTA paths through the local renderer so reverse-strand records are reverse-complemented (the htslib-rs FASTA fast paths skipped that step).
  Brings `bam2fq/{1,2,3,4,6,7,9,11,13,15,16,17,18,19,20}.{1,2,s}.fq.expected` and `bam2fq/11.fa.expected` upstream fixtures to byte parity for SAM input.
- PR #10, https://github.com/madhavajay/samtools-rs/pull/10, branch `view-n-qname-filter` — `view -N FILE` / `--qname-file FILE` allow/deny qname list with `^FILE` negation; wired into the shared `line_passes` filter pass.
- PR #11, https://github.com/madhavajay/samtools-rs/pull/11, branch `view-r-rg-filter` (stacked on PR #10) — `view -r STR` / `-R FILE` accumulating read-group ID filter plus `-n` exclude-no-RG.
- PR #13, https://github.com/madhavajay/samtools-rs/pull/13, branch `view-d-aux-tag-filter` (stacked on PR #11) — `view -d TAG[:VAL]` / `-D TAG:FILE` aux-tag presence/value filtering with shared-tag validation.

Latest known validation (against the tip of PR #13 / latest fastq-index-paired-grouping):
- Rust tests: 408 passing (407 on PR #13, plus the FASTA revcomp test on PR #9 tip).
- Full gate: `cargo fmt --all --check`, `cargo clippy -p samtools-rs --all-targets -- -D warnings`, `cargo test -p samtools-rs`, and `cargo test -p samtools-rs -- --list | rg ': test$' | wc -l` all green.
- New focused tests added across PRs: `fastq_index_files_extract_from_barcode_tag`, `fastq_routes_r1_only_singletons_to_singleton_output`, `fastq_dash_t_and_dash_cap_t_combine_aux_tags`, `fastq_interleaves_read1_read2_when_paths_alias`, `fastq_repeated_dash_d_unions_same_tag_values`, `fasta_reverse_strand_record_reverse_complemented_in_output`, `view_qname_file_filters_records_by_name`, `view_r_and_dash_cap_r_filter_by_read_group`, and `view_d_and_dash_cap_d_filter_by_aux_tag`.

Estimated whole-project completion:
- Roughly 60–65% complete toward the full `samtools` replacement goal once the five in-flight PRs land.
- Rationale: core workspace/subcommand layout, common I/O, many read/write/index/file-operation/statistics/editing commands, the upstream-style fastq split routing, an extended view filter suite (qname/RG/aux-tag), and 407 Rust tests are in place. The remaining risk is concentrated in byte-for-byte upstream parity for the higher-complexity subcommands (view binary aux mutation, sort external merge, markdup/stats per-cycle, full checksum/reference/depad), pileup-dependent commands, full CRAM streaming, and large external algorithms.

What to do next:
1. Land PRs #8–#13 in order (#8, #10, #9, #11, #13; #12 is documentation-only and merges any time).
2. Pick the next bounded slice from **Remaining tractable samtools-rs-only items** below.
3. Run the full gate, update `TODO.md`, `docs/subcommand-coverage.md`, and `docs/test-status.md`, then commit, push, and open the next PR.

Remaining tractable samtools-rs-only items (no htslib-rs / noodles changes required):
- **`fastq` index extraction × name-grouping interaction.** With both PR #8 and PR #9 merged, fold the per-record index emission into the name-grouped flush so each qname-group emits at most one index record per `--i1` / `--i2`. This matches upstream's `flush_rec` → `output_index` and is required for `bam2fq/{5,8,10,12}` parity.
- **`view -d TAG[:VAL]` / `-D TAG:FILE` aux-tag filter.** Same pattern as PR #11 but matches an arbitrary aux tag value rather than `RG`. Already implemented in `fastq`; can be ported into `view::line_passes` with a small helper.
- **`view --library` (`-l`)** library filter via `@RG LB:` aux lookup. Builds on PR #11's read-group infrastructure.
- **`view -X` legacy custom-index synopsis** for BAM/CRAM; the second positional becomes the index path. Already supported as a no-op in `idxstats`.
- **`merge -s SEED`** random-seed acceptance for `-n` mode (currently the option is parsed but the seed is unused).
- **`reheader -c COMMAND`** alternate header source paths (`samtools reheader -c "sed s/foo/bar/" in.bam`) — already partial; rounding out edge cases.
- **`samples` BAM index path verification** for the `-i` index-presence column when index files are at non-default locations.
- **`addreplacerg --output-fmt=cram`** with a `-T` reference — needs reference-backed CRAM writer path (already used in `view`).
- **`stats -d` / `--remove-dups` edge cases**: ensure histogram contributions are excluded for primary duplicates across CRAM record paths.

Items blocked on htslib-rs / noodles extensions (see the rolling list at the end of this file):
- All pileup-dependent commands (`mpileup`, `consensus`, `targetcut`, `phase`, `ampliconstats`, exact pileup-based `bedcov`/`coverage`/`depth`).
- Full `stats` and `checksum` for CRAM input without a region (needs CRAM all-record iterator).
- `cram-size` (needs CRAM container/block API exposure).
- `view --no-PG` for BAM/CRAM output (needs `sam_hdr_add_pg` equivalent that writes through the binary header).
- `flagstat` / `idxstats` for CRAM input without an explicit reference (needs CRAM index meta accessor in `htslib-rs::index_compat`).
- CSI query robustness for very large references (noodles `index out of bounds` panic on the upstream `test_index` `large_chrom.bam ref2` query).

## Progress Snapshot

**Phases 0–2 complete; Wave A complete (with partials); Wave B substantially complete (with partials); Wave D substantially complete (with partials); Wave C in progress (`fixmate`, `rmdup`, `calmd`, `reset`, `markdup`, `depad` partials landed); Phases 4–5 pending. Remaining work largely blocks on htslib-rs infrastructure (pileup APIs, CRAM all-record iterator, custom-index paths) or on substantial per-subcommand efforts (full `mpileup`/`consensus`/`phase`/`ampliconstats`/`targetcut`, full `depad`, full `checksum`, full `reference`, `markdup` full stats parity, `fastq` exact paired/singleton/other semantics, `sort` external merge / template-coordinate / minimiser sorts).**

Subcommands shipped (30 of ~40):
- ✅ byte-parity verified: `flags`, `quickcheck`, `dict`
- ✅ functional with partial-feature notes: `head`, `index`, `idxstats`, `samples` (incl. `-i`, `-f`, `-X`, stdin path lists), `flagstat`, `faidx`/`fqidx`
- 🟡 partial implementation: `view` (SAM↔SAM, SAM→BAM/CRAM, count/header, region queries, `-f/-F/-G/-q` filters, `-L` BED, `-e` filter expression, `-x/--keep-tag` aux strip, `-z` sanitizer mutation, `-p`/`-U` for SAM-input binary output, SAM-output `@PG`/`--no-PG`), `cat` (SAM/BAM record-level concat, `-h`, `-b FILE` input lists, `-r region` for indexed BAM, `@PG`/`--no-PG`), `reheader` (SAM/BAM with `-c` filter, `@PG`/`--no-PG`), `fastq`/`fasta`/`bam2fq` (including `-O` original-quality tags, `-v INT` missing-quality defaults, `-U`/`--UMI-tag` UMI read-name suffixes, and `-i`/`--barcode-tag` CASAVA barcode fields), `split` (with `--no-PG`, `--write-index`), `sort` (in-memory coordinate/name/tag for SAM/BAM/reference-backed CRAM + `@PG`/`--no-PG`), `merge` (in-memory coordinate/name/tag + differing `@SQ` union/remap + `-R region`/`-L BED` + `@PG`/`--no-PG`), `collate` (in-memory name grouping plus `-f` fast primary-pair mode, `-n INT` temp-count compatibility, and legacy positional output prefixes for SAM/BAM/reference-backed CRAM + `@PG`/`--no-PG`), `import`, `rmdup` (single-end + paired-end + `@PG`/`--no-PG`), `markdup` (single-end + paired-end + barcode key + optical-distance `dt` tags + QCFAIL inclusion control + `--mode` compatibility + secondary/supplementary qname propagation + `-r`/`-s`/`-O`/`-o`/`@PG`/`--no-PG`), `bedcov` (CIGAR-walk), `coverage` (CIGAR-walk + ASCII histogram), `depth`, `addreplacerg` (SAM/BAM `-O sam|bam`, `overwrite_all` default, `@PG`/`--no-PG`), `fixmate` (name-sorted BAM/SAM, coordinate-sort rejection, mate TLEN recalculation, MC/MQ, `-m` mate-score tags, `-c` template-CIGAR `ct` tags, default sanitizer mutation, `-r` mode, `@PG`/`--no-PG`), `reset` (alignment field clear, default aux strip, `--reject-PG`/`--no-RG`/`--no-PG` matching upstream `noPGentry` semantics, `@PG` insertion), `depad`/`pad2unpad` (SAM `-T` padded reference to `-s` SAM output), `stats` (extensive SN coverage plus `-f`/`-F` flag filters, `-i` insert-size cap, `-m` insert-size bulk selection, `-l` read-length filtering, `-q` BWA trim counting, FFQ/LFQ quality histograms, GCF/GCL GC histograms, and approximate COV coverage histogram), `calmd`/`fillmd` (SAM/BAM/reference-backed CRAM text MD/NM + SAM BAQ paths + `-d` + `@PG`/`--no-PG`), `reference` (SAM/BAM MD-tag reconstruction + indexed BAM `-r` + `-o`/`-q`)

Remaining subcommands and their blockers:
- **BAM aux-tag mutation is reachable** via `bam::io::Reader::read_record_buf` → `RecordBuf` (mutable) → `bam::io::Writer::write_alignment_record`. This has landed for `addreplacerg` SAM/BAM paths, `fixmate` MC/MQ tags, `markdup` barcode-key grouping, and `reset` aux-keep semantics, and unblocks further BAM aux-rewriting paths in `calmd` BAM MD/NM recompute and `ampliconclip`. The work is still substantial per-subcommand.
- **Need a pileup iterator in htslib-rs** (would unlock `mpileup`, `consensus`, `targetcut`, `phase`, `ampliconstats`, and exact pileup-based `bedcov`/`coverage`/`depth` byte parity).
- **Need a CRAM all-record iterator in htslib-rs** to compute `stats` sequence-length / quality / NM lines on CRAM input without a region (the current `summarize_*` path discards those fields).
- **Other complex algorithms still TODO** (samtools-rs-only but each is a multi-hundred-LOC implementation): full `depad` (623 LOC C; SAM `-T -s` path is partial), full `checksum` (1324 LOC C; default SAM/BAM checksum path is partial), full `reference` (598 LOC C; SAM/BAM MD path is partial), `cram-size` (blocked on htslib-rs CRAM internals), `markdup` full stats parity, `fastq` exact name-grouped paired/singleton/other routing, full `import` read-group/CRAM parity, and `sort` external/template-coordinate/minimiser sorts.

htslib-rs extensions landed during this work:
- `AlignmentRecordSummary` accessors: `flags`, `reference_sequence_id`, `mate_reference_sequence_id`, `mapping_quality`
- `summarize_bam_records_from_path`
- SAM/BAM/reference-backed CRAM filter-expression helpers for `view -c -e`, SAM-output `view -e`, BAM-output `view -b -e`, and CRAM-output `view -C -e` (including indexed BAM/CRAM regions).
- SAM/BAM/reference-backed CRAM stdin reader helpers for `view` count/text/BAM/CRAM paths, including filter-expression support for stdin SAM/BAM/CRAM.
- BAM/SAM FASTA/FASTQ helpers: limit, flag-filter, suffix, split `-1`/`-2`/`-s`, and selected aux tag preservation paths.
- BAM/CRAM region and flag-filter writers: `write_bam_regions_from_path`, `write_bam_regions_as_cram_from_path_with_reference`, `write_bam_records_with_required_flags_from_path`, `write_cram_regions_as_bam_from_path_with_reference`, `write_cram_regions_from_path_with_reference`, `write_cram_records_with_required_flags_as_bam_from_path_with_reference`
- BAM to CRAM writer: `write_cram_from_bam_path_with_reference`
- FASTQ import helpers: paired FASTQ input, index FASTQ input, aux-tag allow-listing, barcode quality tags, and read group tag insertion.

Rust tests: 404 currently passing. `cargo fmt --all --check`, `cargo clippy -p samtools-rs --all-targets -- -D warnings`, and `cargo test -p samtools-rs` are green after the most recent additions: shared `@PG` insertion via `pg::add_samtools_pg_to_header` integrated into `sort`, `merge`, `collate`, `addreplacerg`, `reset`, `fixmate`, `rmdup`, `calmd`, and `view` SAM-output paths; `view -p`/`-U` for SAM-input BAM and CRAM output via text-roundtrip; `view -z` sanitizer mutation for SAM text paths and sanitizer-triggered text roundtrips into BAM/CRAM output; reference-backed CRAM input for in-memory `sort` and `collate`; `merge` differing `@SQ` union/remap, compatible `@SQ` metadata union with conflict rejection, compatible `@HD` metadata union with conflict rejection, compatible `@RG` and `@PG` union, `@CO` comment preservation, `-t TAG` tag ordering with coordinate/name secondary keys, `-s` compatibility, stdout `-` output, and `--output-fmt=FORMAT` parsing, and `-b FILE` input lists; `collate -f` fast primary-pair mode with `-r` working-read cap, `-n INT` temp-count compatibility, legacy positional output prefixes, `-o`/`-O` conflict validation, and `--output-fmt=FORMAT` parsing; `merge -R region` and `-L BED` indexed-BAM restrictions; `cat` gained SAM record-level concatenation, `-b FILE` input lists, and `-r region` indexed-BAM restriction; `depth -H` header output, `-f` input file lists, flag filters (`-g`, `-G`/`--excl-flags`, `--incl-flags`, `--require-flags`), and `-l` minimum read length filtering; SAM/BAM `reheader` with command-filtered headers; `coverage -m`/`-A`/`-w` ASCII histogram output, `-b`/`--bam-list`, `--ff`/`--excl-flags`, `--rf`/`--incl-flags`, `-l` minimum read length filtering, and `-d` maximum-depth capping; `bedcov -g`/`-G` filter-mask controls and `-j` deletion/refskip skipping; `fastq -O` / `bam2fq -O` original-quality tag output, `fastq -v INT` missing-quality defaults, `fastq -U`/`--UMI-tag` UMI read-name suffixes, and `fastq -i`/`--barcode-tag` CASAVA barcode fields for SAM/BAM paths; SAM `depad -T -s` padded-reference conversion now matches the upstream `depad.001` fixture; `import -0` singleton FASTQ input now works alongside paired `-1`/`-2` inputs; `fixmate` now applies default sanitizer mutation against the upstream `fixmate/sanitize.sam.expected` fixture in addition to `-r` mode, coordinate-sort rejection, mate TLEN recalculation, default MC/MQ mate tags, `-m` mate-score tags, and `-c` template-CIGAR `ct` tags; `rmdup` gained paired-end duplicate removal; SE + PE `markdup` (qname-paired groups, combined MAPQ score, barcode-key grouping with `-b`/`--barcode-tag`, `-c` duplicate flag/tag clearing, `-S` compatibility, duplicate-origin `do` tags with `-t`, optical-distance `dt` duplicate-type tags with `-d`, default QCFAIL exclusion with `--include-fails` override, validated `-m`/`--mode` compatibility, optical-aware estimated library size in `-s` stats, secondary/supplementary qname propagation, upstream-shaped `-s` summary fields, `-r`/`-O`/`-o`/`@PG`/`--no-PG`); `calmd` gained SAM/BAM/reference-backed CRAM text MD/NM recomputation against FASTA references and `-d` BQ-tag removal; `stats` extended with `-I`/`--id` read-group/sample filtering, `-f`/`--required-flag`, `-F`/`--filtering-flag`, `-i`/`--insert-size`, `-m`/`--most-inserts`, `-l`/`--read-length`, and `-q`/`--trim-quality`, runtime `is sorted`, supplementary, insert size mean/stddev, inward/outward/other oriented pair counts, total/average/maximum sequence length (per fragment-1/fragment-2), bases mapped (incl. cigar), mismatches, error rate, bases duplicated, bases trimmed, average quality, FFQ/LFQ quality histograms, GCF/GCL GC histograms, approximate COV coverage histograms with `-c`/`--coverage` bin ranges, `-g`/`--cov-threshold` target-percentage SN lines with target-region validation, and percentage properly paired; `checksum` gained SAM/BAM default output plus `-P`/`-C`/`-M` columns, `-B`, wildcard scalar/string/array aux tags, `-a` field-selection shorthand, `-z` sanitizer mutation, and report merging; `reference` gained SAM/BAM MD-tag reconstruction with indexed BAM `-r` and `-o`/`-q`; `sort` gained in-memory `-t TAG` ordering with coordinate/name secondary keys and upstream-style `SS`; `reset --no-PG` semantics fixed to match upstream (preserve existing `@PG`, skip new entry); `addreplacerg` gained SAM/BAM `-O sam|bam` record-path RG header/tag rewriting, and its default mode changed to upstream's `overwrite_all`.

Progress snapshot PR chain:
- noodles: https://github.com/madhavajay/noodles/pull/1
- htslib-rs: https://github.com/madhavajay/htslib-rs/pull/5
- samtools-rs: https://github.com/madhavajay/samtools-rs/pull/7

## What's Next — Decision Points

### Current Goal Constraint

See **Active Goal** above: this pass is samtools-rs-only. Underlying-library
blockers belong at the end of this file, and the pass stops when those are the
only actionable tasks left.

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
2. **Deepen the existing 28 partials toward byte-parity.** Pick a subcommand (probably `view` as the anchor) and drive every flag combination from `test.pl` to byte-for-byte. This is where the remaining byte-parity work, full filter expression support, and CRAM I/O need to land. Estimate: each `test_<name>` group is hours-to-a-day.
3. **Wire the upstream `test.pl` as a gating CI run.** The parity job stages the Rust binary at the harness' ignored `samtools/samtools` path and still fires `|| true`. Flip it as subcommands land, using `docs/test-status.md` to track per-test status (passing / skipped / not-yet-ported / cosmetic-diff). Forces parity attention.

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
- **HTSlib API gaps**: long-term, when samtools-rs needs an HTSlib-shaped API that `htslib-rs` does not yet expose, add it to `htslib-rs` first and route through it. During the current samtools-only pass, defer these gaps to the end of this file and keep working on items that do not require underlying-library changes. Do not bypass `htslib-rs` for HTSlib-shaped APIs. (Direct `noodles` use from samtools-rs is acceptable only for code that has no HTSlib analogue.)
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
- Default to `htslib-rs` for HTSlib-shaped helpers (header manipulation, aux tags, pileup, region parsing, format detection, BGZF, index I/O). For the current samtools-only pass, when `htslib-rs` lacks the API, add the blocker to the end of this `TODO.md` and continue with other samtools-rs work instead of editing the underlying library.
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
- [~] **Global args** (`samtools-rs/src/sam_global.rs`): partial port of `sam_opts.h` / `sam_opts.c`. Top-level parsing now strips and records shared `--input-fmt`, `--input-fmt-option`, `--output-fmt`, `--output-fmt-option`, `--reference`, `--threads`, `--write-index`, and `--verbosity` long options before dispatch, with `--verbosity` applied through `htslib-rs::log_compat`; parsed globals are now stored for command I/O, and `view` / `head` / `idxstats` / `flagstat` / `stats` / `depth` / `coverage` / `bedcov` consume top-level `--reference` for supported CRAM decoding paths. **Pending:** broader thread/reference/format/write-index injection into command I/O paths, and subcommands that accept these options locally still parse them themselves.
- [~] **Open/close helpers** (`samtools-rs/src/io.rs`): shared helpers now wrap `htslib-rs::format::detect_path` as `sam_open_format`, resolve a partial output `sam_open_mode` from explicit format / output extension / default, open text outputs to file/stdout, expose stdout autoflush, and provide `check_sam_close` / `write_all_and_close` flush-error propagation. `calmd`, `view` text modes, `dict`, `faidx`/`fqidx` retrieval, `fastq`/`fasta` single-output text mode, `import`, `addreplacerg`, `depth`, `coverage`, `stats`, and `samples` use the shared text-output close path. `view`, `head`, `calmd`, `addreplacerg`, `fastq`, `depth`, `coverage`, `bedcov`, `stats`, `sort`, `merge`, `collate`, `cat`, `split`, `reheader`, `index`, `fixmate`, `rmdup`, `reset`, `idxstats`, `flagstat`, raw header extraction, and the native API wrappers use the shared format detector; `view` uses the shared output-mode resolver. **Pending:** full `sam_open_mode` writer state/options, auto-index writers (BAI/CSI/CRAI/TBI), complete stdout autoflush/release semantics, and broader binary/subcommand integration.
- [x] **`print_error` / `print_error_errno`** (`samtools-rs/src/diagnostics.rs`): subcommand-prefixed stderr printers matching upstream's `samtools <subcommand>: message` format and errno appending.
- [x] **BAM flag bits + parse/format** (`samtools-rs/src/bam_flag.rs`): the `BAM_F*` constants plus `bam_str2flag` / `bam_flag2str` equivalents. Used by `flags` and (later) `view` / `flagstat`.
- [x] **Raw header text extractor** (`samtools-rs/src/header_text.rs`): preserves original `@`-line order for SAM, BAM, and CRAM inputs. Required by `head`, `view -h`, `samples`, etc. — noodles' canonical writer reorders headers, which breaks byte parity.
- [~] **BAM sanitizer** (`samtools-rs/src/sanitize.rs`): partial port of `samtools.h`'s `FIX_*` flags and `bam_sanitize_options`, including upstream reset semantics for `all` / `none` / `off`, `cigarx` implying `cigdup`, and the parser's special `on` behavior. Record-level mutation now covers position/unmap correction, mapping-quality reset, unmapped CIGAR clearing, NM/MD/CG/SM aux stripping, overhang CIGAR trimming, CIGAR `=`/`X` conversion, and duplicate CIGAR op merging; `fixmate`, `checksum`, and `view -z` consume the shared parser/mutator, with `fixmate` matching the upstream `fixmate/sanitize.sam.expected` fixture under `--no-PG`. **Pending:** broaden edge-case parity against upstream sanitizer fixtures and direct binary-record sanitizer paths where text roundtrips are currently used.
- [~] **@PG add helper** (`samtools-rs/src/pg.rs`): shared helper now builds raw-header `@PG` lines with HTSlib-style argv stringification, generated unique IDs, `PN`, `VN`, `CL`, and `PP` links for terminal program chains. `cat`, `split`, `reheader`, `sort`, `merge`, `collate`, `addreplacerg`, `reset`, `fixmate`, `rmdup`, and `view`'s SAM-output paths (file-input header-only, SAM output with `-h`, plus BAM/CRAM stdin SAM/header-only) use it for default output headers and honor `--no-PG`. **Pending:** integrate `pg::add_samtools_pg_to_header` into `view`'s binary BAM/CRAM output paths (currently emitted by `htslib-rs` internal writers — see *htslib-rs Extensions Needed*) and verify byte-parity against upstream `sam_hdr_add_pg` for complex merge/split/reheader cases.
- [x] **Aux-tag list parser** (`samtools-rs/src/aux_list.rs`): port `parse_aux_list` from `sam_utils.c`. Used by `view`, `reset`, `fastq`, and future aux-aware commands.
- [~] **BED index** (`samtools-rs/src/bedidx.rs`): shared BED parser/index now stores 0-based half-open intervals by reference, skips comments/UCSC metadata, emits HTSlib-style 1-based inclusive region strings, supports overlap queries, and is used by `view -L`, `depth -b`, `bedcov`, and native `view_bed`. **Pending:** interval-tree acceleration/parity with `bedidx.c`, stricter upstream diagnostics where needed, and integration into `ampliconclip` and future `mpileup`.
- [~] **Reference helpers** (`samtools-rs/src/reference.rs`): shared FASTA helper now derives associated `.fai` paths, builds missing FASTA indexes through `htslib-rs::faidx_compat`, loads `(SN, LN)` dictionaries, and matches candidate FASTA references against BAM/CRAM `@SQ` dictionaries for `samples -f/-F`. **Pending:** mmap/FASTA sequence cache, common `--reference` option plumbing, CRAM reference resolution, and integration into `calmd`, `consensus`, `mpileup`, `phase`, and `import`.
- [~] **Temp file helper** (`samtools-rs/src/tmp_file.rs`): shared temp path helper now creates collision-resistant temp files, owns best-effort cleanup on drop, supports explicit persist/close, and is used by native name-sort FASTQ conversion instead of ad hoc temp names. **Pending:** BAM record temp spooling, compression support, and integration into external `sort` / `collate` algorithms.
- [x] **Logging passthrough**: bridge to `htslib-rs::log_compat` so top-level `--verbosity` flows correctly.

## Phase 2: Subcommand Surface Mapping

Mapping document exists at [`docs/subcommand-coverage.md`](docs/subcommand-coverage.md). It lists every subcommand, the HTSlib APIs it depends on, the `htslib-rs` coverage status, and a rolled-up list of extensions needed in `htslib-rs`.

- [x] Per-subcommand HTSlib API surface enumerated.
- [x] `htslib-rs` coverage status per API (already exposed / needs adding / out of scope).
- [x] Gap list rolled up at the end of `docs/subcommand-coverage.md`.

## Phase 3: Subcommand Implementation Waves

Each subcommand below maps to: (a) one Rust module under `crates/samtools-rs/src/commands/`, (b) `test_<name>` in `samtools/test/test.pl` passing against the Rust binary, (c) at least one Rust integration test under `crates/samtools-rs/tests/<name>.rs`.

The waves are ordered to land foundational machinery first (read/write/index) and unblock the rest.

### Wave A — Read/Write/Index Foundation

- [~] `view` (`sam_view.c`, 68k) — partial: SAM↔SAM passthrough, SAM→BAM/CRAM, SAM/BAM/reference-backed CRAM stdin count/text/BAM/CRAM paths, reference-backed CRAM→SAM text/count paths including flag/MAPQ filtered count mode, reference-backed CRAM→BAM full-file and region output, reference-backed BAM→CRAM and CRAM→CRAM full-file and region output, header-only / count modes (including CRAM `-H`), `-h` `-H` `-c` `-b` `-C` `-T` `-o` `--no-PG`, filter flags `-f`/`-F`/`-G`/`-q` for SAM output/count modes, SAM-input BAM/CRAM output, BAM/CRAM-input binary output, and SAM/BAM/reference-backed CRAM stdin binary output, region queries (`<chr:start-end>`), `-L FILE` BED restrict, `-U FILE` unselected SAM-output splitting for flag/MAPQ and expression filters, `-p/--unmap` SAM-output marking for records failing flag/MAPQ and expression filters (sets UNMAP, MAPQ=0, CIGAR=`*`, TLEN=0), `-U FILE` and `-p/--unmap` for SAM-input BAM output (text → BAM roundtrip via `build_split_sam_text` and `write_bam_from_sam_reader`), `-e EXPR` filter expression count/SAM/BAM/CRAM output modes for SAM/BAM/reference-backed CRAM (including indexed BAM/CRAM regions and SAM/BAM/reference-backed CRAM stdin), `-x/--remove-tag` and `--keep-tag` aux stripping for SAM output and SAM-input BAM/CRAM output, `-O FORMAT` output-fmt option, default `@PG` insertion on SAM-output paths (header-only/SAM/BAM-stdin/CRAM-stdin), and `--no-PG`. **Pending:** BAM/CRAM-input binary aux-tag manipulation, BAM/CRAM-input `-U`/`-p` binary output (needs aux mutation), BAM/CRAM-output binary `@PG` insertion (`htslib-rs` writer extension), multi-file inputs, paired-aware filters, full CRAM parity.
- [x] `head` (`sam_view.c` shared) — SAM and BAM input; SAM/BAM/CRAM stdin header/record output; CRAM header-only modes; reference-backed CRAM record extraction for `-n N`; `-h N`, `-n N`, all-default.
- [x] `quickcheck` (`bam_quickcheck.c`) — passes byte-for-byte against `quickcheck/all.expected`.
- [x] `index` (`bam_index.c`) — BAI/CSI/CRAI build, `-c` CSI mode, `--min-shift`, `-M`, `-o`, legacy `<in> <out.idx>` synopsis. **Pending:** `-@` threads not yet propagated to noodles workers.
- [x] `idxstats` (`bam_stat.c`) — index-based per-reference counts for BAM, with streaming slow-path counts for SAM, reference-backed CRAM, and unindexed BAM; tests cover both successful reference-backed CRAM and clean missing-reference failure. **Pending:** index-derived CRAM counting path for CRAM inputs without a reference.
- [~] `faidx` / `fqidx` (`faidx.c`) — index-build mode works (`samtools faidx file.fa` produces `file.fa.fai`); BGZF FASTA/FASTQ input now writes `.gzi` and can be indexed/retrieved; local region extraction works for positional regions, `-r` region files, `-o`, `.gz`/`.bgz`/`.bgzf` BGZF output, `--length`, `--write-index` for file outputs, FASTQ mode via `fqidx` and `faidx -f`, reverse-complement `-i` with mark-strand modes, `--continue`-style missing-region tolerance, and upstream-style zero/truncated region warning keywords. The upstream `test_faidx`/`test_fqidx` section now progresses through its checked commands in the local parity harness. **Pending:** exact warning text parity, compression-level/thread option effects, and broader BGZI edge cases.
- [x] `dict` (`dict.c`) — sequence dictionary builder. Passes byte-for-byte against `dict.out`, `dict.alias.out`, `dict.alt.out` (run via test.pl-style stdin/file invocations).
- [x] `flagstat` / `flagstats` (`bam_stat.c`) — SAM, BAM, and reference-backed CRAM input. Default + `-O json` + `-O tsv` output modes. Tests cover both successful reference-backed CRAM and clean missing-reference failure. Required extending `htslib-rs::alignment_compat::AlignmentRecordSummary` with `flags_u16` / `reference_sequence_id` / `mate_reference_sequence_id` / `mapping_quality` accessors, plus BAM and reference-backed CRAM summary paths. **Pending:** CRAM input without an explicit reference remains unsupported.

### Wave B — File Ops

- [~] `sort` (`bam_sort.c`, 138k — the largest single file) — basic in-memory coordinate sort, `-n` name sort, and `-t TAG` tag sort for BAM, SAM, and reference-backed CRAM inputs. Supports `-o`, `-O sam|bam`, `--output-fmt sam|bam`, coordinate-sort BAM `--write-index`, default `@PG` insertion via `pg::add_samtools_pg_to_header`, `--no-PG`, and sets `@HD SO`/tag-sort `SS`. **Pending:** on-disk external merge for large inputs, template-coordinate sort (`-M`), minimiser sort (`-N`), thread/memory caps, and CRAM output.
- [~] `merge` (`bam_sort.c` shared) — basic in-memory multi-input merge for BAM and SAM inputs with coordinate, name, or `-t TAG` tag ordering. It unions differing `@SQ` dictionaries, remaps reference and mate-reference ids into the merged header, unions compatible same-name `@SQ` metadata fields, preserves and unions compatible `@HD` metadata fields while leaving output `SO`/`SS` controlled by the selected merge order, unions compatible `@RG` and `@PG` definitions, appends `@CO` comments, and rejects conflicting sequence lengths, sequence metadata fields, header metadata fields, read-group definitions, or program definitions. `-f` overwrite, `-o`, stdout `-`, `-O sam|bam`, `--output-fmt sam|bam` / `--output-fmt=sam|bam`, accepted `-s INT`, `-b FILE` input lists, coordinate-sort BAM `--write-index`, default `@PG` insertion via `pg::add_samtools_pg_to_header`, `--no-PG`, `-R region` (indexed BAM only), and `-L bed` (indexed BAM only, de-duplicates records from overlapping BED intervals) are supported. **Pending:** k-way streaming merge, broader header reconciliation beyond `@HD`/`@SQ`/`@RG`/`@PG`/`@CO`, CRAM.
- [~] `collate` / `bamshuf` (`bamshuf.c`) — basic in-memory name-sort grouping for BAM, SAM, and reference-backed CRAM inputs. `-f` fast mode outputs strict primary READ1/READ2 pairs early and omits secondary/supplementary records, with `-r` controlling the working-read cap before unmatched records are deferred to grouped output. `-o FILE`, `-O` stdout, `-n INT` temp-count compatibility, legacy positional output prefixes, `-o`/`-O` conflict validation, validated `--output-fmt sam|bam` / `--output-fmt=sam|bam`, default `@PG` insertion via `pg::add_samtools_pg_to_header`, upstream-style `@HD SO:unsorted GO:query`, and `--no-PG` are supported. **Pending:** on-disk hash-bucket implementation for inputs larger than memory and CRAM output.
- [~] `cat` (`bam_cat.c`) — basic SAM and BAM concatenation works (record-level decompress + re-encode). Supports `-o`, `-h` (header replacement), `-b FILE` input lists (expanded before positional inputs), default `@PG` insertion, `--no-PG`, and `-r region` (indexed BAM only — restricts each input to records overlapping the region via `query_bam_records_from_path`). **Pending:** BGZF block-level fast path, CRAM, `-p N/M`.
- [~] `split` (`bam_split.c`) — basic BAM/SAM-by-`@RG` splitting with per-output `@RG` header filtering and default `@PG` insertion; explicit `-d TAG` string/integer aux-tag splitting with on-demand outputs; explicit `-d RG` unknown-read-group header insertion; `-M`/`--max-split`, `-f` template (`%*`, `%!`, `%#`, `%.`), `-u` unaccounted, `-h` unaccounted SAM header override, `--output-fmt sam|bam`, `--no-PG`, `--write-index` BAI generation for BAM outputs, and `-p N` padding. **Pending:** CRAM, sorted-by-tag streaming mode, and deeper upstream `@PG` byte-parity for complex chains.
- [~] `reheader` (`bam_reheader.c`) — basic SAM/BAM header replacement (record-level rewrite) with default `@PG` insertion, `--no-PG` suppression, and `-c <command>` external header filtering. **Pending:** BGZF block-level BAM fast path and CRAM `--in-place`.
- [~] `addreplacerg` (`bam_addrprg.c`) — SAM/BAM add/replace for the `@RG` header line and `RG:Z` aux tag. `-r SPEC`, `-R ID`, `-m overwrite_all|orphan_only` (default `overwrite_all`, matching upstream), `-O sam|bam`, `-o FILE`, default `@PG` insertion, and `--no-PG`. The SAM→SAM path remains streaming text rewrite; other SAM/BAM combinations use mutable `RecordBuf` records. **Pending:** CRAM input/output, mate-aware updates, full orphan-first semantics.
- [~] `fastq` / `fasta` / `bam2fq` (`bam_fastq.c`) — basic single-stream output works for SAM and BAM (records written to stdout, `-o FILE`, or `-0 FILE`), with `-f`/`--require-flags`, `--rf`/`--include-flags`, `-F`/`--exclude-flags`, `-G`, the upstream default `0x900` secondary/supplementary exclusion, read-name suffix controls (`-n`/`-N`), `-O` original-quality `OQ` tag output, `-v INT` missing-quality defaults for FASTQ, `-U`/`--UMI-tag` UMI read-name suffixes, `-i`/`--barcode-tag` CASAVA barcode fields, basic flag-driven paired split outputs (`-1`/`-2`/`-s`/`-0`), SAM/BAM selected aux comments via `-T` in single and split output modes, all-tag SAM/BAM comments via `-T ''` / `-T '*'` in single and split-output FASTQ mode, SAM/BAM `B` array aux comment formatting, SAM/BAM single and split-output FASTQ tag filtering via `-d`/`--tag TAG[:VALUE]` and `-D`/`--tag-file TAG:FILE`, and `-t` as the upstream shortcut for `RG,BC,QT`. The upstream SAM-input all-tags fixture `bam2fq/15.fq.expected` now matches for `-T ''` and `-t -T '*'`. **Pending:** index FASTQ file extraction and exact name-grouped paired/singleton/other semantics.
- [~] `import` (`bam_import.c`) — basic single FASTA/FASTQ and paired FASTQ (`-1`/`-2`, `--r1`/`--r2`, `-s` interleaved, plus two positional inputs) → SAM/BAM (`-O bam` / `--bam`), including positional single input plus `-0` single-read alias, `-0` singleton input alongside paired `-1`/`-2`, positional interleaved FASTQ detection from `/1`/`/2` read names, no-op `--no-PG`, CASAVA parsing (`-i`) with upstream-style reverse comments, SRA name2 (`-N`), UMI extraction (`-U`/`--UMI-tag`) with reverse comments, CASAVA barcode sequence tags (`--barcode-tag`), FASTQ definition aux tags (`-T`) including upstream-style float exponent spelling, explicit index reads (`--i1`/`--i2`) for `-0`, `-s`, positional interleaved, and paired `-1`/`-2` inputs with barcode sequence/quality tags (`--barcode-tag`/`--quality-tag`) and `-b`, and read-group header/tag support (`-R`/`-r`) with repeated `-r` accumulation, `-r` precedence over `-R`, and `-r` ID validation. Direct comparisons against `test/import/*.expected.sam` for the currently implemented import fixture commands pass. **Pending:** full paired singleton/other grouping parity, full read-group parity, CRAM output.

### Wave C — Editing / Mate-aware

- [~] `fixmate` (`bam_mate.c`) — basic mate flag/pos fixup for adjacent paired records in name-sorted BAM and SAM inputs (`FMUNMAP`, `FMATE_REVERSE_COMPLEMENTED`, `mate_reference_sequence_id`, `mate_alignment_start`) and rejects `@HD SO:coordinate` input like upstream. TLEN is recalculated from mate 5-prime positions, including large coordinate inputs where the resulting template length still fits. Default MC/MQ mate aux tags are added for mapped mates and cleared when the mate is unmapped. `-m` adds `ms:i` mate-score tags for markdup. `-c` adds lowercase template-CIGAR `ct:Z` tags to the earlier mapped mate and clears stale `ct` tags from both mates. `-z`/`--sanitize` now parses and validates through the shared sanitizer option parser, and the default sanitizer mutates records to match the upstream `sanitize.sam` fixture. Default `@PG` insertion via `pg::add_samtools_pg_to_header` and `--no-PG` are supported. `-r` removes secondary and unmapped alignments and clears `PROPER_PAIR`/`MATE_REVERSE` on the surviving mate when its partner is unmapped, matching upstream's `remove_reads` semantics. **Pending:** CRAM, mate-rescore, base-modification `-M` parity.
- [~] `markdup` (`bam_markdup.c`, 89k) — single-end and paired-end duplicate marking for SAM and BAM inputs. SE records are keyed by `(tid, pos, reverse-flag)` plus optional barcode tag; PE records pair by qname and are keyed by the canonical pair of `(tid, pos, strand)` triples plus optional per-end barcode tag. Within a group the entry with highest (combined) MAPQ stays primary and the rest receive `BAM_FDUP`; secondary and supplementary alignments with duplicate primary qnames inherit the duplicate flag and are removed by `-r`. Supports `-b TAG`/`--barcode-tag TAG`, `-c` (clear existing duplicate flags and duplicate metadata tags), `-S` (accepted; propagation is always on), `-t` duplicate-origin `do` tags, `-d DISTANCE` optical-distance duplicate classification with `dt:Z:SQ`/`dt:Z:LB` tags, default QCFAIL exclusion with `--include-fails` override, validated `-m t|s`/`--mode t|s` compatibility, optical-aware estimated library size in `-s` stats, `-r` (remove duplicates from output), `-s` (upstream-shaped summary counts to stderr), `-O sam|bam`, `-o FILE`, default `@PG` insertion, and `--no-PG`. **Pending:** exact upstream stats output/count parity, CRAM.
- [~] `rmdup` (`bam_rmdup.c` + `bam_rmdupse.c`) — single-end and paired-end duplicate removal for BAM and SAM inputs. SE records are keyed by `(tid, pos, reverse-flag)`; PE records pair by qname and are keyed by the canonical pair of `(tid, pos, strand)` triples, retaining the highest MAPQ/combined MAPQ record or pair. `-s`/`-S` force single-end treatment. Default `@PG` insertion via `pg::add_samtools_pg_to_header` and `--no-PG` are supported. **Pending:** CRAM, full upstream deprecated-command parity.
- [~] `calmd` / `fillmd` (`bam_md.c`) — SAM, BAM, and reference-backed CRAM input can emit SAM text with MD/NM tags recomputed against a FASTA reference via CIGAR/reference walking. BAQ paths (`-r`, `-r -e`, `-E`) are wired through `htslib_rs::alignment_compat::recalculate_baq_*` and `apply_existing_baq_from_sam_path` for SAM input, and `-d` drops existing `BQ` tags from the SAM-text output. Default `@PG` insertion via `pg::add_samtools_pg` (text-level) and `--no-PG` are supported. **Pending:** BAM/CRAM output, BAM/CRAM BAQ paths, `-A`/`-C cap`, full upstream MD/BAQ parity.
- [ ] `targetcut` (`cut_target.c`) — fosmid pool target cutting.
- [~] `reset` (`reset.c`) — strip alignment fields (`tid`/`pos`/`cigar`/`mate_*`/`template_length`) for BAM and SAM inputs, set MAPQ to `0`, drop a default set of aligner aux tags (NM, MD, AS, XS, SA, MC, MQ, NH, HI, ms), clear `PROPER_PAIR`/`SECONDARY`/`SUPPLEMENTARY`/`REVERSE`/`MATE_REVERSE`, set `UNMAPPED`, set `MATE_UNMAPPED` for paired reads, reverse-restore reverse-strand sequence/quality, preserve duplicate flags with `--dupflag`, remove read-group headers/tags with `--no-RG`, remove program header chains with `--reject-PG`, add a new samtools `@PG` chain entry by default (via the shared `pg::add_samtools_pg_to_header` helper), suppress the new `@PG` with `--no-PG` while preserving existing entries (matching upstream's `noPGentry` semantics), accept SAM/BAM input from stdin/no positional input/`-`, and tolerate legacy SAM `@HD VN:1` headers. `-x`/`--keep-tag` honored, with `--no-RG` taking precedence over keeping `RG`. **Pending:** CRAM I/O.
- [ ] `ampliconclip` (`bam_ampliconclip.c`, 40k) — primer clipping with BED amplicon spec.

### Wave D — Stats / Pileup

- [~] `depth` (`bam2depth.c`) — per-position depth via CIGAR walks for SAM, BAM, and reference-backed indexed CRAM. `-a`/`-aa`/`-d`/`-q`/`-o`, `-H` header output, `-f` input file lists, flag filters (`-g`, `-G`/`--excl-flags`, `--incl-flags`, `--require-flags`), `-l` minimum read length filtering, `-r` region restriction, `-b` BED restriction, and multi-input columnar output are supported. Tests cover SAM region depth, successful CRAM region depth with top-level `--reference`, multi-input BAM columns/header/list files, flag-filter controls, minimum read length filtering, and clean missing-reference failure. **Pending:** pileup-based exact handling of overlaps/deletions, CRAM input without an explicit reference.
- [~] `coverage` (`coverage.c`) — per-reference/`-r` region `numreads`, `covbases`, `coverage`, `meandepth`, `meanbaseq`, and `meanmapq` via CIGAR walks for SAM, BAM, and reference-backed indexed CRAM. `--min-depth` thresholds covered-base counts, `-Q`/`--min-BQ` filters low-quality bases, `-q`/`--min-MQ` filters reads, `-b`/`--bam-list` expands input filename lists, `--ff`/`--excl-flags` replaces the default filter-out flags, `--rf`/`--incl-flags` requires at least one selected flag, `-l`/`--min-read-len` filters short alignments, `-d` caps per-position depth for reported coverage/depth metrics, and multiple inputs aggregate into one row per reference/region. `-m`/`--histogram` (with `-A`/`--ascii` and `-D`/`--plot-depth` routed through the same path) emits a 10-row ASCII histogram with `-w`/`--n-bins` controlling the column count. Tests cover SAM region coverage, successful CRAM region coverage with top-level `--reference`, non-zero BAM/CRAM mean base quality, base/depth filtering, input list expansion, flag/read-length filtering, maximum-depth capping, multi-input BAM aggregation, ASCII histogram output, and clean missing-reference failure. **Pending:** byte-parity match against upstream's UTF-8 box-drawing histogram + sidebar text, `-D` true depth-plot semantics, CRAM input without an explicit reference.
- [~] `bedcov` (`bedcov.c`) — total aligned-base coverage per BED region, walking each record's CIGAR for SAM, BAM, and reference-backed indexed CRAM. `-Q` mapq filter, `-g`/`-G` filter-mask controls, `-j` deletion/refskip skipping, `-H` output headers, `-c` read-count columns, and `-d` depth-threshold columns are supported. Tests cover SAM depth/count columns, flag-mask controls, deletion/refskip handling, successful CRAM BED coverage with top-level `--reference`, and clean missing-reference failure. **Pending:** pileup-based exact coverage for byte parity.
- [~] `stats` (`stats.c` + `stats_isize.c`, 123k + 8k) — `SN` (Summary Numbers) section plus FFQ/LFQ first/last fragment quality histograms, GCF/GCL first/last fragment GC histograms, and approximate CIGAR-walk COV coverage histograms for SAM, BAM, and reference-backed CRAM region paths, including record-backed `-I`/`--id` read-group/sample filtering, `-f`/`--required-flag`, `-F`/`--filtering-flag`, `-i`/`--insert-size` insert-size cap, `-m`/`--most-inserts` insert-size bulk selection, `-l`/`--read-length`, `-q`/`--trim-quality` BWA trim counting, `-c`/`--coverage MIN,MAX,STEP` COV binning, and `-g`/`--cov-threshold` target percentage lines with target-region validation. SAM and BAM iterate records directly to populate sequence-length, quality, GC, CIGAR, NM, COV, and runtime coordinate-order accumulators (CRAM non-region remains on the `summarize_*` summary path and therefore reports zeros for those fields, with `is sorted` still falling back to header `SO`). The emitted lines now cover: raw total / filtered / sequences / runtime is sorted / 1st & last fragments / mapped / mapped+paired / unmapped / properly paired / paired / duplicated / MQ0 / QC-failed / non-primary / supplementary / total length / total first & last fragment length / bases mapped / bases mapped (cigar) / mismatches (NM aux) / error rate / average length & per-fragment / maximum length & per-fragment / bases trimmed / average quality / singletons / insert size mean & stddev / inward, outward, and other oriented pair counts / pairs on different chromosomes / percentage of properly paired reads / target bases / target genome coverage above threshold. SAM, indexed BAM, and reference-backed CRAM positional region arguments and `-t` target files restrict the summary and COV positions, with overlapping BAM/CRAM regions de-duplicated. `-d` / `--remove-dups` filters duplicate-marked primary records and their quality/GC/COV histogram contributions. Missing CRAM references fail cleanly. **Pending:** exact pileup-backed COV byte parity, per-cycle, BAQ, CRAM non-region read-group/read-length/trim-quality filtering, and CRAM non-region sequence-length/quality/GC/COV stats (requires a CRAM all-record iterator in `htslib-rs`).
- [ ] `mpileup` (`bam_plcmd.c`, 49k) — multi-way pileup with `htslib-rs::alignment_compat` pileup support. Output formats including VCF, depth-aware.
- [ ] `consensus` (`bam_consensus.c`, 126k, + `consensus_pileup.c`) — consensus FASTA/FASTQ/pileup builder.
- [ ] `phase` (`phase.c`) — heterozygote phasing.
- [~] `depad` / `pad2unpad` (`padding.c`) — SAM input with `-T` padded FASTA reference and `-s` SAM output converts padded reference columns to unpadded coordinates/CIGAR (`I`/`P`) and matches the upstream `depad.001` fixture with `--no-PG`. **Pending:** BAM input/output, CRAM, binary output modes (`-u`/`-1`), and full upstream `test_depad` parity.
- [ ] `ampliconstats` (`amplicon_stats.c`, 65k) — per-amplicon stats over the same BED model as `ampliconclip`.
- [ ] `cram-size` (`cram_size.c`) — CRAM Content-ID and Data-Series byte breakdown. Depends on `htslib-rs` CRAM internals (currently out of scope in htslib-rs — see *htslib-rs Extensions Needed*).
- [~] `checksum` (`bam_checksum.c`, 47k) — default order-agnostic checksum output for SAM/BAM input is implemented, including `-o`, `-f`/`-F`/`-b`, `-c`, `-N`, `-q`, `-v`, `-T`, `-O`, `-P`, `-C`, `-M`, `-B`, `-a` field-selection shorthand with upstream-style sanitizer defaults, `-z`/`--sanitize` record mutation, `-m` checksum-output merging for default/position/CIGAR/mate-column reports, selected and wildcard/exclusion aux-tag hashing for scalar/string/array tags with canonical integer encoding, read-group grouping, and Rust tests matching upstream `chk1.1.expected` and `chk1.3.expected` after the harness' checksum path-line normalization. **Pending:** CRAM input (needs htslib-rs CRAM all-record iterator) and full upstream `test_checksum`.
- [x] `samples` (`bam_samples.c`) — list `@RG SM:` samples across inputs. Header-driven dedup, `-T TAG`, `-o`, `-h`, `-i` index-presence column, `-f`/`-F` FASTA dictionary matching, stdin path lists, `-X` custom index pairs, and CRAM headers are implemented.
- [~] `reference` (`reference.c`) — SAM/BAM MD-tag reconstruction to FASTA is implemented, including `-o`, `-q`, basic `-r region` output, and indexed BAM region iteration when an associated BAI/CSI is present. **Pending:** embedded-reference CRAM mode (`-e`, blocked on CRAM container/block internals), CRAM input MD path (needs CRAM all-record iteration/reference handling), and full upstream `test_reference` parity.
- [x] `flags` (`bam_flags.c`) — explain a numeric BAM flag. Byte-for-byte parity with upstream.

## Phase 4: Test Harness Integration

- [~] **Parity gate setup**: confirmed the pinned upstream harness does not honor `-e samtools=<rust-binary-path>` for most commands because it constructs commands from `$$opts{bin}/samtools` after option parsing. CI now stages the Rust binary at the ignored `samtools/samtools` path and runs `cd samtools && perl test/test.pl || true` without modifying the harness. **Pending:** verify `regression.sh` and flip the parity job from advisory to required once the tracked groups pass.
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
- [ ] **Indexing BAMs without `@HD SO:coordinate` metadata** — upstream `samtools index` indexes fixtures such as `test/dat/test_input_1_a.bam` and `test_input_1_b.bam` even when the header has no coordinate sort-order tag. Current `htslib-rs` index creation rejects these with `invalid sort order: expected coordinate, got None`; defer this here rather than changing `htslib-rs` during the current samtools-only pass.
- [ ] **CRAM internals for `cram-size`** — currently marked out of scope in `htslib-rs/README.md`. Either expose a minimal block/container/codec inventory API in `htslib-rs`, or drop `cram-size` from samtools-rs scope.
- [ ] **`htslib-rs::region`** — confirm coverage of HTSlib's full region-string grammar including `*` (unmapped) and `.` (everything else).
- [ ] **`probaln_glocal` and BAQ recalculation** — exposed by `htslib-rs::probaln`; verify wiring for `calmd` and `mpileup`.
- [ ] **CRAM all-record iterator** — `htslib-rs::alignment_compat` exposes `iter_cram_records_from_path_with_reference` for *indexed region* iteration but no non-region streaming iterator. Needed for `stats` to populate sequence-length and quality SN lines on CRAM input that does not request a region, and for `checksum` to process whole CRAM inputs such as the upstream `test_checksum` fixtures.

## noodles Extensions Needed (rolling list)

Keep these at the end during the current samtools-only pass. Do not modify the
noodles submodule for these blockers until explicitly switching back to
underlying-library work.

- [ ] **CSI query robustness for very large references/regions** — the local parity harness now reaches `test_index`, where `samtools view large_chrom.bam ref2` panics inside `noodles-csi/src/binning_index/index/reference_sequence.rs` with `index out of bounds`, and `ref2:1-541556283` reports `invalid end bound`. Defer to noodles/htslib-rs region/index handling rather than patching noodles from this samtools-rs pass.

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
cp target/release/samtools samtools/samtools
cd samtools && perl test/test.pl

# optional: refresh expected outputs from C samtools (local dev only)
cd samtools && autoreconf -i && ./configure && make
cd test && perl test.pl --redo-outputs
```
