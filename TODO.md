# TODO: Port samtools to Pure Rust

Goal: build a pure Rust replacement for the `samtools` C program with full subcommand parity, then port and pass the upstream `test/test.pl` suite plus add Rust-native unit/integration tests. Implementation routes through `htslib-rs` (sibling submodule); long-term, when a needed HTSlib API is not yet exposed there, extend `htslib-rs` first.

## Active Goal

**`TODO-NEXT.md` is COMPLETE (all 12 library/infra items done,
byte/fixture-verified, committed).** The earlier samtools-only
constraint no longer applies: minimal patches to the **owned vendored
noodles fork** (`madhavajay/noodles`, an `htslib-rs` submodule) are
sanctioned by the Ground rule's "(and carry minimal patches)" clause,
and every htslib-rs / noodles gap has been closed that way (pins
bumped). The remaining project work is **only the two fixtureless
subcommands `phase` + `targetcut`** (no upstream `test/*` expected
outputs → faithful ports verifiable by unit tests, not the byte-exact
harness) and **Phase 4/5 polish**. Per-subcommand integration tests
and the byte-exact upstream harness already cover every other
subcommand.

## Current Handoff — 2026-05-17

> **`TODO-NEXT.md` COMPLETE (12/12).** All library/infra blockers
> resolved via the owned vendored noodles fork; every upstream-fixtured
> subcommand byte-exact vs its full harness. Only `phase` +
> `targetcut` (no upstream fixtures) and Phase 4/5 polish remain.
> The dated narrative below this banner is historical.

Merged baseline:
- PR #7, https://github.com/madhavajay/samtools-rs/pull/7, is merged into `main`.
- Merge commit: `ae15e4ad603892912fe0c5175491e1d0e3f210eb`.
- PRs #8–#15 (the work below) were consolidated into `integration-all-prs` and merged back to `main`.
- The next work should start from `main` on a new short-lived branch.

Merged PRs (samtools-rs-only, consolidated via PR #16 `integration-all-prs` → `main` as `b312c99`, CI green):
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
- PR #14, https://github.com/madhavajay/samtools-rs/pull/14, branch `addreplacerg-default-rg` (from `main`, two commits) — `addreplacerg` defaults to the first header `@RG` ID when neither `-r`/`-R` is given, upstream `@RG` header reconciliation (`-r` + `overwrite_all` strips other `@RG` lines; `-w` overwrites a same-ID line), and `-R ID` rejection when the ID is absent from the header. Brings the whole upstream `test_addrprg` group (`addrprg/{1,2,3,4,5}`, #3 = expected failure) to parity modulo `@PG`.
- PR #15, https://github.com/madhavajay/samtools-rs/pull/15, branch `reheader-parity` (from `main`) — reorders the shared `pg::push_pg_line` output to upstream's `@PG` field order `ID, PN, PP, VN, CL`. Benefits every command that inserts a samtools `@PG`; the upstream harness strips `\tVN:.*`, so `PP` must precede `VN`. Brings the `reheader/{1,4}` header section to parity after harness reordering.

Active working branch (not yet merged): `work-sam-float-renderer` (off `main`).
Landed slices on this branch:
- Shared `sam_render` module with htslib-style aux float formatting.
- `view` BAM/CRAM→SAM aux float spelling via `sam_render`.
- `split` SAM output routed through `sam_render::write_record`.
- **SAM aux float formatting — remaining commands (this slice):** `reheader`
  SAM→SAM, `sort`, `merge`, `collate`, `addreplacerg`, `reset`, `fixmate`,
  `rmdup`, `markdup`, and `cat` SAM-output sinks now wrap a plain
  `File`/`Stdout` and render through `sam_render::write_record` /
  `write_header`. Full gate green after this slice (samtools-rs: 422
  passing, 0 failing; `cargo test --workspace`: 2593 passing; fmt + clippy
  `-D warnings` clean). New test:
  `sort_sam_output_uses_htslib_float_aux_spelling`.
- **`view -X` legacy custom-index synopsis:** `view -X` /
  `--customized-index` accepts `in.bam in.bam.bai [region…]` (index
  positional accepted as a no-op). New test:
  `view_dash_cap_x_accepts_legacy_custom_index_synopsis`.
- **`view --library` / `-l`:** resolves `@RG LB:STR` → RG-ID set from
  the header (path + SAM/BAM/CRAM stdin) and filters records by
  `RG:Z:` membership. New test:
  `view_dash_l_filters_by_read_group_library`. (`merge -s SEED` was
  examined and skipped: it's already accepted/consumed, and its only
  upstream effect — random RG/PG-ID collision suffixing — would require
  reworking merge header reconciliation, out of scope for a bounded
  slice.)
- **`samples -i` custom-index resolution:** `-X` index path now
  resolves an exact file, a directory holding the index, or a prefix
  (matching `sam_index_load3`). New test:
  `samples_custom_index_directory_reports_index_presence`.
- **`addreplacerg` CRAM output:** `-O cram` / `--output-fmt[=]cram`
  with `-T`/`--reference` writes reference-backed CRAM (SAM/BAM input)
  via a temp-BAM → shared CRAM-writer path. New test:
  `addreplacerg_writes_cram_output_with_reference`.
- **`stats -d` CRAM region path verified:** the CRAM region path
  shares the SAM/BAM `update` chokepoint, so `--remove-dups` already
  excludes `BAM_FDUP` records from histograms there. New test:
  `stats_remove_dups_excludes_duplicates_on_cram_region_path`. The
  no-region CRAM summarize path remains blocked on the htslib-rs CRAM
  all-record iterator.
- **`fastq` index × name-grouping (complete):** one index record per
  qname-group, htslib-exact CASAVA barcode normalization
  (`ac-gt` → `AC+GT`), the CASAVA comment on `-i` index records, and
  cross-mate barcode propagation (R2/other inherits the R1 mate's
  `BC`). **`bam2fq/{5,8,10,12}` now byte parity on every output.** New
  tests: `fastq_index_emits_one_record_per_qname_group_with_casava_comment`,
  `fastq_casava_barcode_propagates_from_r1_to_r2_mate`.

Batch status: **all "Remaining tractable samtools-rs-only items" are now
done or explicitly deferred** (only `merge -s SEED`, whose sole upstream
effect needs a merge header-reconciliation rework — a substantial
samtools-rs effort, not a bounded slice). Per the **Active Goal**, the
remaining actionable work requires htslib-rs / noodles changes (see
*Items blocked …* below); this samtools-only pass stops here. Final
validation on `work-sam-float-renderer`: `cargo fmt --all --check`,
`cargo clippy --workspace --all-targets -- -D warnings`, and
`cargo test --workspace` (**2601 passing, 0 failing**) all green; the
upstream `bam2fq/{5,8,10,12}` fixtures now match byte-for-byte. Next:
open one PR for this branch, get CI green, merge to `main`, then a new
working branch (the next batch is htslib-rs/noodles-blocked work).

Latest known validation (on `main` at `b312c99`, post-merge):
- Rust tests: 416 `samtools-rs` passing, 0 failing (`cargo test --workspace`: 2587 passing).
- Full gate green in CI on PR #16: `cargo fmt --all --check`, `cargo clippy --workspace --all-targets -- -D warnings`, `cargo test --workspace`, and the advisory parity gate.
- New focused tests added across PRs: `fastq_index_files_extract_from_barcode_tag`, `fastq_routes_r1_only_singletons_to_singleton_output`, `fastq_dash_t_and_dash_cap_t_combine_aux_tags`, `fastq_interleaves_read1_read2_when_paths_alias`, `fastq_repeated_dash_d_unions_same_tag_values`, `fasta_reverse_strand_record_reverse_complemented_in_output`, `view_qname_file_filters_records_by_name`, `view_r_and_dash_cap_r_filter_by_read_group`, `view_d_and_dash_cap_d_filter_by_aux_tag`, `addreplacerg_defaults_to_first_header_rg_and_preserves_lines`, `addreplacerg_r_overwrite_all_removes_other_header_rg_lines`, and `addreplacerg_dash_cap_r_unknown_id_is_rejected`.

Estimated whole-project completion:
- **Roughly 95%+** toward the full `samtools` replacement goal. All 12
  `TODO-NEXT.md` library/infra items are done and committed, and every
  upstream-fixtured subcommand is byte-exact vs its entire harness:
  `consensus` (77/77 `consensus.reg`), `sort` (all `test_sort` incl.
  minimiser/template-coordinate), `cram-size` (3/3), `reference` (all
  `test_reference` incl. embed_ref read+write), `calmd`, `stats`,
  `markdup`, `ampliconclip`, `ampliconstats`, `merge`, `fixmate`,
  `addreplacerg`, `reset`, `split`, `view`, `idxstats`/`flagstat`
  (incl. CRAM-no-ref), plus the `coverage`/`bedcov`/`depth`/`mpileup`
  tabular suites.
- Remaining risk is confined to: `phase`/`targetcut` (no upstream
  fixtures — faithful ports only), the optional non-parity niceties
  (CRAM-NM recompute for exact `stats` error-rate; UTF-8 `coverage`
  histogram; mpileup BAQ quals/VCF; external-merge perf), and Phase
  4/5 (exit-code/thread/perf triage; the per-subcommand integration
  tests are already in place).

### Workflow rule: one large PR branch at a time

Do **not** open many small per-slice PRs. Use a single long-lived working
branch off `main` (e.g. `work-<topic>`), land multiple related bounded
slices onto it as separate commits, keep the full gate green after each
commit, and open **one** large PR for that batch. Only start the next
working branch after that PR has merged to `main`. (PRs #8–#15 were the
fragmented anti-pattern; #16 consolidated them — keep it to one large PR
branch going forward.)

What to do next:
1. Create one working branch off `main` for the next batch of bounded slices.
2. Land each slice from **Remaining tractable samtools-rs-only items** below as its own commit on that branch, running the full gate (`cargo fmt --all --check`, `cargo clippy --workspace --all-targets -- -D warnings`, `cargo test --workspace`) green after every commit, and updating `TODO.md`, `docs/subcommand-coverage.md`, and `docs/test-status.md` as you go.
3. Open a single PR for the whole batch, get CI green, merge to `main`, then start the next working branch.

Remaining tractable samtools-rs-only items (no htslib-rs / noodles changes required):
- ~~**SAM aux float formatting — remaining commands.**~~ **Done.** The shared `samtools_rs::sam_render` module (`format_aux_float`, `format_htslib_exponent`, `fix_sam_aux_floats`, `fix_sam_text`, `write_record`, `write_header`) now backs every noodles-`sam::io::Writer` SAM-output path: `view`, `split`, plus `reheader` SAM→SAM, `sort`, `merge`, `collate`, `addreplacerg`, `reset`, `fixmate`, `rmdup`, `markdup`, and `cat` (their `SamFile`/`SamStdout`/`Sam*Sink` sinks now wrap a plain `File`/`Stdout` and render through `sam_render`). So every SAM-text output path emits htslib `%g`-style float aux spelling. Regression covered by `sort_sam_output_uses_htslib_float_aux_spelling` plus the existing `view`/`split` fixtures.
- ~~**`fastq` index extraction × name-grouping interaction.**~~ **Done.** `emit_index_files` dedupes to **one index record per adjacent qname-group** (matching upstream `flush_rec` → `output_index`). The CASAVA barcode field is normalized exactly like htslib `fastq_format1` (`casava_barcode_field`: absent/non-sequence-first → `0`; otherwise non-alpha → `+`, lowercase → upper, e.g. `ac-gt` → `AC+GT`), `-i` index records carry the ` <rnum>:<filt>:0:<barcode>` CASAVA comment, and `GroupedSplitWriter` now **propagates the group barcode across mates** so an R2 (or other) record lacking its own `BC` gets the R1 mate's barcode in its CASAVA comment (`fill_casava_barcode`; upstream `bam_fastq.c:952`). **`bam2fq/{5,8,10,12}` now match byte-for-byte on every output (`1.fq`/`2.fq`/`s.fq` and the index/`bc` files).** New tests: `fastq_index_emits_one_record_per_qname_group_with_casava_comment`, `fastq_casava_barcode_propagates_from_r1_to_r2_mate` (plus the 45 existing fastq tests still green).
- ~~**`view --library` (`-l`)** library filter via `@RG LB:` aux lookup.~~ **Done.** `view -l STR` / `--library STR` resolves the requested library to the set of `@RG` IDs whose `LB:` equals STR (scanned from the input header for path, SAM/BAM/CRAM stdin), then a record passes iff its `RG:Z:` value is in that set (no-RG / non-matching RG excluded, matching upstream `bam_get_library`). Regression: `view_dash_l_filters_by_read_group_library`.
- ~~**`view -X` legacy custom-index synopsis**~~ **Done.** `view -X` / `--customized-index` accepts the legacy synopsis where the second positional is the explicit index path (`view -X in.bam in.bam.bai [region…]`); accepted as a no-op (our region queries build/find the index themselves), matching `idxstats -X`. Regression: `view_dash_cap_x_accepts_legacy_custom_index_synopsis`.
- **`merge -s SEED`** — *deferred (not a bounded slice).* The option is already parsed and its value consumed. Upstream's only use of the seed is `hts_srand48` feeding `lrand48()` for random `@RG`/`@PG`-ID collision suffixes during header merge (`bam_sort.c:408`). Our merge reconciles headers by *rejecting* ID conflicts rather than random-suffixing, so the seed has no observable effect until that suffixing path is implemented — a header-reconciliation rework, larger than a bounded slice.
- ~~**`samples` BAM index path verification**~~ **Done.** `samples -i` with a custom `-X` index path now mirrors `sam_index_load3`: an exact index file, a *directory* holding the index (`<dir>/<data-name>.bai`), or a suffix-less prefix all resolve via the shared `locate_associated_index` resolver, so index files at non-default locations register `Y`. Regression: `samples_custom_index_directory_reports_index_presence` (and the existing exact-file/pair test still passes).
- ~~**`addreplacerg --output-fmt=cram`** with a `-T` reference~~ **Done.** `addreplacerg` accepts `-O cram` / `--output-fmt cram` / `--output-fmt=cram` and `-T`/`--reference[=]FILE`; SAM/BAM input → CRAM output spools rewritten records to a temp BAM and converts via the shared `write_cram_from_bam_path_with_reference` (the `.fai` is built if missing). CRAM output without `-T` errors. Regression: `addreplacerg_writes_cram_output_with_reference`.
- ~~**`stats -d` / `--remove-dups` edge cases**~~ **Done (tractable part).** The CRAM *region* path iterates real records through the same `update_record_with_targets` → `update` chokepoint as SAM/BAM, which already gates all histogram/seq/quality accumulation on `self.total` increasing (and `--remove-dups` filters `BAM_FDUP` before `total` is bumped). Verified end-to-end by `stats_remove_dups_excludes_duplicates_on_cram_region_path` (SAM→CRAM→indexed→region stats with/without `-d`). The CRAM *no-region* path uses the `summarize_cram_records_from_path_with_reference` summary path, which discards per-record seq/quality — that remains **blocked on the htslib-rs CRAM all-record iterator** (already tracked in the blocked list).

Items previously blocked on htslib-rs / noodles extensions — **ALL
RESOLVED** (via `TODO-NEXT.md` #1–#12; the "htslib-rs Extensions
Needed" list at the end of this file is fully checked off):
- ~~pileup-dependent commands~~ — pileup iterator done (#1); `mpileup`/
  `consensus`/`coverage`/`bedcov`/`depth` byte-exact (`consensus`
  77/77); `ampliconstats` done; only `phase`/`targetcut` remain
  (no fixtures).
- ~~`stats`/`checksum` CRAM no-region~~ — CRAM all-record iterator
  done (#2); wired.
- ~~`cram-size`~~ — done (#3), all 3 `cram_size.reg` byte-exact.
- ~~binary `@PG`~~ — done (#4), `view` SAM→BAM/CRAM + BAM→BAM.
- ~~`flagstat`/`idxstats` CRAM-no-ref~~ — done (#5), byte-exact.
- ~~large-reference CSI~~ — done (#12).
- ~~SAM aux float formatting~~ — resolved via `sam_render`.

## Progress Snapshot

**Phases 0–2 complete; Waves A/B/C/D substantially complete and byte-exact for the upstream-fixtured subcommands. `TODO-NEXT.md` is COMPLETE — all 12 numbered library/infra items done, byte/fixture-verified and committed (pileup iterator, CRAM all-record iterator, CRAM container/codec inventory → `cram-size`, binary `@PG`, aux mutation, threads, write-index, SO-less index, region grammar, BAQ, large-ref CSI, embed_ref read+write).** Subcommands now byte-exact vs their entire upstream harness include: `consensus` (all 77 `consensus.reg`), `sort` (all `test_sort`), `cram-size` (all 3), `reference` (all `test_reference`), `calmd` (`test_calmd`), `stats`, `markdup`, `ampliconclip`, `ampliconstats`, `merge`, `fixmate`, `addreplacerg`, `reset`, plus `coverage`/`bedcov`/`depth`/`mpileup` tabular suites. **Remaining: only `phase` (`phase.c`, 843 LOC) and `targetcut` (`cut_target.c`, 257 LOC) — neither has upstream test fixtures (dense numerical HMM / errmod ports, verifiable only by faithful port + unit tests, not the byte-exact harness) — plus Phase 4/5 polish (per-subcommand integration tests largely already in place; thread/exit-code/perf triage).**

Subcommands shipped (30 of ~40):
- ✅ byte-parity verified: `flags`, `quickcheck`, `dict`
- ✅ functional with partial-feature notes: `head`, `index`, `idxstats`, `samples` (incl. `-i`, `-f`, `-X`, stdin path lists), `flagstat`, `faidx`/`fqidx`
- 🟡 partial implementation: `view` (SAM↔SAM, SAM→BAM/CRAM, count/header, region queries, `-f/-F/-G/-q` filters, `-L` BED, `-e` filter expression, `-x/--keep-tag` aux strip, `-z` sanitizer mutation, `-p`/`-U` for SAM-input binary output, SAM-output `@PG`/`--no-PG`), `cat` (SAM/BAM record-level concat, `-h`, `-b FILE` input lists, `-r region` for indexed BAM, `@PG`/`--no-PG`), `reheader` (SAM/BAM with `-c` filter, `@PG`/`--no-PG`), `fastq`/`fasta`/`bam2fq` (including `-O` original-quality tags, `-v INT` missing-quality defaults, `-U`/`--UMI-tag` UMI read-name suffixes, and `-i`/`--barcode-tag` CASAVA barcode fields), `split` (with `--no-PG`, `--write-index`), `sort` (in-memory coordinate/name/tag for SAM/BAM/reference-backed CRAM + `@PG`/`--no-PG`), `merge` (in-memory coordinate/name/tag + differing `@SQ` union/remap + `-R region`/`-L BED` + `@PG`/`--no-PG`), `collate` (in-memory name grouping plus `-f` fast primary-pair mode, `-n INT` temp-count compatibility, and legacy positional output prefixes for SAM/BAM/reference-backed CRAM + `@PG`/`--no-PG`), `import`, `rmdup` (single-end + paired-end + `@PG`/`--no-PG`), `markdup` (single-end + paired-end + barcode key + optical-distance `dt` tags + QCFAIL inclusion control + `--mode` compatibility + secondary/supplementary qname propagation + `-r`/`-s`/`-O`/`-o`/`@PG`/`--no-PG`), `bedcov` (CIGAR-walk), `coverage` (CIGAR-walk + ASCII histogram), `depth`, `addreplacerg` (SAM/BAM `-O sam|bam`, `overwrite_all` default, `@PG`/`--no-PG`), `fixmate` (name-sorted BAM/SAM, coordinate-sort rejection, mate TLEN recalculation, MC/MQ, `-m` mate-score tags, `-c` template-CIGAR `ct` tags, default sanitizer mutation, `-r` mode, `@PG`/`--no-PG`), `reset` (alignment field clear, default aux strip, `--reject-PG`/`--no-RG`/`--no-PG` matching upstream `noPGentry` semantics, `@PG` insertion), `depad`/`pad2unpad` (SAM `-T` padded reference to `-s` SAM output), `stats` (extensive SN coverage plus `-f`/`-F` flag filters, `-i` insert-size cap, `-m` insert-size bulk selection, `-l` read-length filtering, `-q` BWA trim counting, FFQ/LFQ quality histograms, GCF/GCL GC histograms, and approximate COV coverage histogram), `calmd`/`fillmd` (SAM/BAM/reference-backed CRAM text MD/NM + SAM BAQ paths + `-d` + `@PG`/`--no-PG`), `reference` (SAM/BAM MD-tag reconstruction + indexed BAM `-r` + `-o`/`-q`)

Remaining subcommands and their blockers — **resolved** (kept for
history; current truth in the banner + Progress Snapshot above):
- ~~BAM aux-tag mutation~~ — done; all parity consumers byte-exact.
- ~~pileup iterator in htslib-rs~~ — done (#1); `mpileup`/`consensus`
  (77/77) / `ampliconstats` / pileup-based `bedcov`/`coverage`/`depth`
  byte-exact. Only `phase`/`targetcut` remain (no upstream fixtures).
- ~~CRAM all-record iterator~~ — done (#2); `stats`/`checksum`
  no-region CRAM + `reference` MD path wired.
- ~~Other complex algorithms~~ — `cram-size` ✅ (3/3), `reference` ✅
  (full `test_reference`, embed_ref read+write), `checksum`/`markdup`/
  `sort` (incl. minimiser + template-coordinate) byte-exact. The only
  not-fully-ported subcommands are `phase` + `targetcut` (no
  fixtures); `depad` BAM/CRAM and a few non-parity niceties remain
  cosmetic.

htslib-rs extensions landed during this work:
- `AlignmentRecordSummary` accessors: `flags`, `reference_sequence_id`, `mate_reference_sequence_id`, `mapping_quality`
- `summarize_bam_records_from_path`
- SAM/BAM/reference-backed CRAM filter-expression helpers for `view -c -e`, SAM-output `view -e`, BAM-output `view -b -e`, and CRAM-output `view -C -e` (including indexed BAM/CRAM regions).
- SAM/BAM/reference-backed CRAM stdin reader helpers for `view` count/text/BAM/CRAM paths, including filter-expression support for stdin SAM/BAM/CRAM.
- BAM/SAM FASTA/FASTQ helpers: limit, flag-filter, suffix, split `-1`/`-2`/`-s`, and selected aux tag preservation paths.
- BAM/CRAM region and flag-filter writers: `write_bam_regions_from_path`, `write_bam_regions_as_cram_from_path_with_reference`, `write_bam_records_with_required_flags_from_path`, `write_cram_regions_as_bam_from_path_with_reference`, `write_cram_regions_from_path_with_reference`, `write_cram_records_with_required_flags_as_bam_from_path_with_reference`
- BAM to CRAM writer: `write_cram_from_bam_path_with_reference`
- FASTQ import helpers: paired FASTQ input, index FASTQ input, aux-tag allow-listing, barcode quality tags, and read group tag insertion.

Rust tests: 405 currently passing. `cargo fmt --all --check`, `cargo clippy -p samtools-rs --all-targets -- -D warnings`, and `cargo test -p samtools-rs` are green after the most recent additions: shared `@PG` insertion via `pg::add_samtools_pg_to_header` integrated into `sort`, `merge`, `collate`, `addreplacerg`, `reset`, `fixmate`, `rmdup`, `calmd`, and `view` SAM-output paths; `view -p`/`-U` for SAM-input BAM and CRAM output via text-roundtrip; `view -z` sanitizer mutation for SAM text paths and sanitizer-triggered text roundtrips into BAM/CRAM output; reference-backed CRAM input for in-memory `sort` and `collate`; `merge` differing `@SQ` union/remap, compatible `@SQ` metadata union with conflict rejection, compatible `@HD` metadata union with conflict rejection, compatible `@RG` and `@PG` union, `@CO` comment preservation, `-t TAG` tag ordering with coordinate/name secondary keys, `-s` compatibility, stdout `-` output, and `--output-fmt=FORMAT` parsing, and `-b FILE` input lists; `collate -f` fast primary-pair mode with `-r` working-read cap, `-n INT` temp-count compatibility, legacy positional output prefixes, `-o`/`-O` conflict validation, and `--output-fmt=FORMAT` parsing; `merge -R region` and `-L BED` indexed-BAM restrictions; `cat` gained SAM record-level concatenation, `-b FILE` input lists, and `-r region` indexed-BAM restriction; `depth -H` header output, `-f` input file lists, flag filters (`-g`, `-G`/`--excl-flags`, `--incl-flags`, `--require-flags`), and `-l` minimum read length filtering; SAM/BAM `reheader` with command-filtered headers; `coverage -m`/`-A`/`-w` ASCII histogram output, `-b`/`--bam-list`, `--ff`/`--excl-flags`, `--rf`/`--incl-flags`, `-l` minimum read length filtering, and `-d` maximum-depth capping; `bedcov -g`/`-G` filter-mask controls and `-j` deletion/refskip skipping; `fastq -O` / `bam2fq -O` original-quality tag output, `fastq -v INT` missing-quality defaults, `fastq -U`/`--UMI-tag` UMI read-name suffixes, `fastq -i`/`--barcode-tag` CASAVA barcode fields for SAM/BAM paths, and `fastq --i1`/`--i2` per-record index FASTQ extraction with `--index-format`/`--quality-tag`; SAM `depad -T -s` padded-reference conversion now matches the upstream `depad.001` fixture; `import -0` singleton FASTQ input now works alongside paired `-1`/`-2` inputs; `fixmate` now applies default sanitizer mutation against the upstream `fixmate/sanitize.sam.expected` fixture in addition to `-r` mode, coordinate-sort rejection, mate TLEN recalculation, default MC/MQ mate tags, `-m` mate-score tags, and `-c` template-CIGAR `ct` tags; `rmdup` gained paired-end duplicate removal; SE + PE `markdup` (qname-paired groups, combined MAPQ score, barcode-key grouping with `-b`/`--barcode-tag`, `-c` duplicate flag/tag clearing, `-S` compatibility, duplicate-origin `do` tags with `-t`, optical-distance `dt` duplicate-type tags with `-d`, default QCFAIL exclusion with `--include-fails` override, validated `-m`/`--mode` compatibility, optical-aware estimated library size in `-s` stats, secondary/supplementary qname propagation, upstream-shaped `-s` summary fields, `-r`/`-O`/`-o`/`@PG`/`--no-PG`); `calmd` gained SAM/BAM/reference-backed CRAM text MD/NM recomputation against FASTA references and `-d` BQ-tag removal; `stats` extended with `-I`/`--id` read-group/sample filtering, `-f`/`--required-flag`, `-F`/`--filtering-flag`, `-i`/`--insert-size`, `-m`/`--most-inserts`, `-l`/`--read-length`, and `-q`/`--trim-quality`, runtime `is sorted`, supplementary, insert size mean/stddev, inward/outward/other oriented pair counts, total/average/maximum sequence length (per fragment-1/fragment-2), bases mapped (incl. cigar), mismatches, error rate, bases duplicated, bases trimmed, average quality, FFQ/LFQ quality histograms, GCF/GCL GC histograms, approximate COV coverage histograms with `-c`/`--coverage` bin ranges, `-g`/`--cov-threshold` target-percentage SN lines with target-region validation, and percentage properly paired; `checksum` gained SAM/BAM default output plus `-P`/`-C`/`-M` columns, `-B`, wildcard scalar/string/array aux tags, `-a` field-selection shorthand, `-z` sanitizer mutation, and report merging; `reference` gained SAM/BAM MD-tag reconstruction with indexed BAM `-r` and `-o`/`-q`; `sort` gained in-memory `-t TAG` ordering with coordinate/name secondary keys and upstream-style `SS`; `reset --no-PG` semantics fixed to match upstream (preserve existing `@PG`, skip new entry); `addreplacerg` gained SAM/BAM `-O sam|bam` record-path RG header/tag rewriting, and its default mode changed to upstream's `overwrite_all`.

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
- [~] **@PG add helper** (`samtools-rs/src/pg.rs`): shared helper now builds raw-header `@PG` lines with HTSlib-style argv stringification, generated unique IDs, and upstream `sam_hdr_add_pg` field order `ID, PN, PP, VN, CL` (PP precedes VN/CL so the upstream harness' `s/\tVN:.*//` normalization keeps `PP`), with `PP` links for terminal program chains. `cat`, `split`, `reheader`, `sort`, `merge`, `collate`, `addreplacerg`, `reset`, `fixmate`, `rmdup`, and `view`'s SAM-output paths (file-input header-only, SAM output with `-h`, plus BAM/CRAM stdin SAM/header-only) use it for default output headers and honor `--no-PG`. The upstream `reheader/{1,4}` header section now matches after harness reordering. **Pending:** integrate `pg::add_samtools_pg_to_header` into `view`'s binary BAM/CRAM output paths (currently emitted by `htslib-rs` internal writers — see *htslib-rs Extensions Needed*) and verify byte-parity against upstream `sam_hdr_add_pg` for complex merge/split/reheader cases.
- [x] **Aux-tag list parser** (`samtools-rs/src/aux_list.rs`): port `parse_aux_list` from `sam_utils.c`. Used by `view`, `reset`, `fastq`, and future aux-aware commands.
- [~] **BED index** (`samtools-rs/src/bedidx.rs`): shared BED parser/index now stores 0-based half-open intervals by reference, skips comments/UCSC metadata, emits HTSlib-style 1-based inclusive region strings, supports overlap queries, and is used by `view -L`, `depth -b`, `bedcov`, and native `view_bed`. **Pending:** interval-tree acceleration/parity with `bedidx.c`, stricter upstream diagnostics where needed, and integration into `ampliconclip` and future `mpileup`.
- [~] **Reference helpers** (`samtools-rs/src/reference.rs`): shared FASTA helper now derives associated `.fai` paths, builds missing FASTA indexes through `htslib-rs::faidx_compat`, loads `(SN, LN)` dictionaries, and matches candidate FASTA references against BAM/CRAM `@SQ` dictionaries for `samples -f/-F`. **Pending:** mmap/FASTA sequence cache, common `--reference` option plumbing, CRAM reference resolution, and integration into `calmd`, `consensus`, `mpileup`, `phase`, and `import`.
- [~] **Temp file helper** (`samtools-rs/src/tmp_file.rs`): shared temp path helper now creates collision-resistant temp files, owns best-effort cleanup on drop, supports explicit persist/close, and is used by native name-sort FASTQ conversion instead of ad hoc temp names. **Pending:** BAM record temp spooling, compression support, and integration into external `sort` / `collate` algorithms.
- [x] **Logging passthrough**: bridge to `htslib-rs::log_compat` so top-level `--verbosity` flows correctly.
- [x] **SAM render helper** (`samtools-rs/src/sam_render.rs`): shared htslib-style aux float formatting (`format_aux_float`/`format_htslib_exponent`), SAM-line/SAM-text float fixers (`fix_sam_aux_floats`/`fix_sam_text`), and noodles `sam::io::Writer` drop-ins (`write_record`/`write_header`). Used by `fastq` (float helper) and now every SAM-text output path: `view`, `split`, `reheader` SAM→SAM, `sort`, `merge`, `collate`, `addreplacerg`, `reset`, `fixmate`, `rmdup`, `markdup`, and `cat` route through `write_record`/`write_header`, so binary→SAM and SAM→SAM both get htslib `%g` float spelling.

## Phase 2: Subcommand Surface Mapping

Mapping document exists at [`docs/subcommand-coverage.md`](docs/subcommand-coverage.md). It lists every subcommand, the HTSlib APIs it depends on, the `htslib-rs` coverage status, and a rolled-up list of extensions needed in `htslib-rs`.

- [x] Per-subcommand HTSlib API surface enumerated.
- [x] `htslib-rs` coverage status per API (already exposed / needs adding / out of scope).
- [x] Gap list rolled up at the end of `docs/subcommand-coverage.md`.

## Phase 3: Subcommand Implementation Waves

Each subcommand below maps to: (a) one Rust module under `crates/samtools-rs/src/commands/`, (b) `test_<name>` in `samtools/test/test.pl` passing against the Rust binary, (c) at least one Rust integration test under `crates/samtools-rs/tests/<name>.rs`.

The waves are ordered to land foundational machinery first (read/write/index) and unblock the rest.

### Wave A — Read/Write/Index Foundation

- [~] `view` (`sam_view.c`, 68k) — partial: SAM↔SAM passthrough, SAM→BAM/CRAM, SAM/BAM/reference-backed CRAM stdin count/text/BAM/CRAM paths, reference-backed CRAM→SAM text/count paths including flag/MAPQ filtered count mode, reference-backed CRAM→BAM full-file and region output, reference-backed BAM→CRAM and CRAM→CRAM full-file and region output, header-only / count modes (including CRAM `-H`), `-h` `-H` `-c` `-b` `-C` `-T` `-o` `--no-PG`, filter flags `-f`/`-F`/`-G`/`-q` for SAM output/count modes, SAM-input BAM/CRAM output, BAM/CRAM-input binary output, and SAM/BAM/reference-backed CRAM stdin binary output, region queries (`<chr:start-end>`), `-L FILE` BED restrict, `-U FILE` unselected SAM-output splitting for flag/MAPQ and expression filters, `-p/--unmap` SAM-output marking for records failing flag/MAPQ and expression filters (sets UNMAP, MAPQ=0, CIGAR=`*`, TLEN=0), `-U FILE` and `-p/--unmap` for SAM-input BAM output (text → BAM roundtrip via `build_split_sam_text` and `write_bam_from_sam_reader`), `-e EXPR` filter expression count/SAM/BAM/CRAM output modes for SAM/BAM/reference-backed CRAM (including indexed BAM/CRAM regions and SAM/BAM/reference-backed CRAM stdin), `-x/--remove-tag` and `--keep-tag` aux stripping for SAM output and SAM-input BAM/CRAM output, `-O FORMAT` output-fmt option, default `@PG` insertion on SAM-output paths (header-only/SAM/BAM-stdin/CRAM-stdin), `--no-PG`, `-N`/`--qname-file FILE` read-name allow/deny filtering with `^FILE` negation, accumulating `-r STR` / `-R FILE` read-group filtering, `-n` exclude-no-read-group filtering, and `-d TAG[:VAL]` / `-D TAG:FILE` aux-tag presence/value filtering with shared-tag validation on every SAM-line-based path. **Pending:** BAM/CRAM-input binary aux-tag manipulation, BAM/CRAM-input `-U`/`-p` binary output (needs aux mutation), BAM/CRAM-output binary `@PG` insertion (`htslib-rs` writer extension), multi-file inputs, paired-aware filters, full CRAM parity.
- [x] `head` (`sam_view.c` shared) — SAM and BAM input; SAM/BAM/CRAM stdin header/record output; CRAM header-only modes; reference-backed CRAM record extraction for `-n N`; `-h N`, `-n N`, all-default.
- [x] `quickcheck` (`bam_quickcheck.c`) — passes byte-for-byte against `quickcheck/all.expected`.
- [x] `index` (`bam_index.c`) — BAI/CSI/CRAI build, `-c` CSI mode, `--min-shift`, `-M`, `-o`, legacy `<in> <out.idx>` synopsis. **Pending:** `-@` threads not yet propagated to noodles workers.
- [x] `idxstats` (`bam_stat.c`) — index-based per-reference counts for BAM, with streaming slow-path counts for SAM, reference-backed CRAM, and unindexed BAM; tests cover both successful reference-backed CRAM and clean missing-reference failure. **Pending:** index-derived CRAM counting path for CRAM inputs without a reference.
- [~] `faidx` / `fqidx` (`faidx.c`) — index-build mode works (`samtools faidx file.fa` produces `file.fa.fai`); BGZF FASTA/FASTQ input now writes `.gzi` and can be indexed/retrieved; local region extraction works for positional regions, `-r` region files, `-o`, `.gz`/`.bgz`/`.bgzf` BGZF output, `--length`, `--write-index` for file outputs, FASTQ mode via `fqidx` and `faidx -f`, reverse-complement `-i` with mark-strand modes, `--continue`-style missing-region tolerance, and upstream-style zero/truncated region warning keywords. The upstream `test_faidx`/`test_fqidx` section now progresses through its checked commands in the local parity harness. **Pending:** exact warning text parity, compression-level/thread option effects, and broader BGZI edge cases.
- [x] `dict` (`dict.c`) — sequence dictionary builder. Passes byte-for-byte against `dict.out`, `dict.alias.out`, `dict.alt.out` (run via test.pl-style stdin/file invocations).
- [x] `flagstat` / `flagstats` (`bam_stat.c`) — SAM, BAM, and reference-backed CRAM input. Default + `-O json` + `-O tsv` output modes. Tests cover both successful reference-backed CRAM and clean missing-reference failure. Required extending `htslib-rs::alignment_compat::AlignmentRecordSummary` with `flags_u16` / `reference_sequence_id` / `mate_reference_sequence_id` / `mapping_quality` accessors, plus BAM and reference-backed CRAM summary paths. **Pending:** CRAM input without an explicit reference remains unsupported.

### Wave B — File Ops

- [~] `sort` (`bam_sort.c`, 138k — in-memory coordinate/`-n` natural-name/`-t TAG` sort for BAM/SAM/reference-backed CRAM. `-o`/`-o -`(stdout)/`-O sam|bam`/`--write-index`/`--no-PG`. Emits the **raw input header** (preserving @SQ/@RG field order & @CO) with @HD `SO`/`SS` applied + text @PG; name sort uses the exact `bam_sort.c` comparator (`strnum_cmp` natural order + the `flag&0xc0/0x100/0x800` READ1/READ2/supp/sec tiebreak). Tolerates SAM `c/C/s/S/I` scalar aux integer synonyms (htslib-compatible). **Byte-exact vs upstream `sort/{pos,name,name3,tag.rg,tag.rg.n,tag.as}.sort.expected.sam`** (modulo @PG); test `sort_matches_upstream_test_sort_fixtures`. **Pending:** on-disk external merge, template-coordinate (`-M`), minimiser (`-N`/`-K`), CRAM output.
- [x] `merge` (`bam_sort.c` shared) — in-memory multi-input merge (BAM/SAM) with coordinate/`-n` natural-name/`-t TAG` order. **Upstream `-s SEED` @RG/@PG reconciliation implemented**: seeded `gen_unique_id` (`crate::rand48` glibc LCG) suffixes colliding IDs in header-line order per file; raw merged header (input[0] @HD verbatim, no forced SO for coordinate; SO/SS for name/tag; @SQ unioned by SN; @RG/@PG ID/PG:/PP: remapped; @CO appended); records' `RG:Z:`/`PG:Z:` remapped, dropped when unresolved (`bam_translate` 'tag lost'). `-r` attaches a filename-stem @RG to every record (order-preserving RG:Z delete-then-append); `-c`/`-p` combine identical @RG/@PG IDs (grouped short opts `-cp`/`-rp`). `-f`/`-o`/`-o -`/`-O sam|bam`/`-b`/`-R`/`-L`/`--write-index`/`--no-PG`. **Byte-exact vs ALL upstream `merge/{2,4,5,6,7}.merge.expected.sam`** (modulo @PG); test `merge_reconciles_rg_pg_byte_exact_vs_upstream`. **Pending (not fixture-covered):** k-way streaming merge, CRAM output, `--template-coordinate`.
- [x] `collate` / `bamshuf` (`bamshuf.c`) — in-memory grouping for BAM, SAM (tolerant `c/C/s/S/I` aux reader), and reference-backed CRAM. **Non-fast order is the exact upstream `bamshuf` order**: bucket by `hash_X31_Wang(qname) % 64`, then sort each bucket by `(hash, qname, flag>>6&3)` (ported `hash_Wang` bit-mix). Fast `-f` mode mirrors the ring buffer: evict-the-oldest-after-insert so a read whose mate is further than `-r` away is deferred. SAM output emits the **raw input header** with `@HD SO:unsorted GO:query` applied (preserving input `@RG`/`@SQ` field order) + records via `sam_render`. Output format inferred from the `-o` filename extension when `--output-fmt` is absent. `-o`/`-O`/`-n`/positional prefix/`--no-PG`. **Byte-exact vs the ENTIRE upstream `test_collate` harness (6/6)**; tests `collate_*` + `collate_matches_upstream_test_collate_fixtures` in `sort_merge.rs`. **Pending:** on-disk hash-bucket for inputs larger than memory, CRAM output.
- [~] `cat` (`bam_cat.c`) — basic SAM and BAM concatenation works (record-level decompress + re-encode). Supports `-o`, `-h` (header replacement), `-b FILE` input lists (expanded before positional inputs), default `@PG` insertion, `--no-PG`, and `-r region` (indexed BAM only — restricts each input to records overlapping the region via `query_bam_records_from_path`). **Pending:** BGZF block-level fast path, CRAM, `-p N/M`.
- [~] `split` (`bam_split.c`) — basic BAM/SAM-by-`@RG` splitting with per-output `@RG` header filtering and default `@PG` insertion; explicit `-d TAG` string/integer aux-tag splitting with on-demand outputs; explicit `-d RG` unknown-read-group header insertion; `-M`/`--max-split`, `-f` template (`%*`, `%!`, `%#`, `%.`), `-u` unaccounted, `-h` unaccounted SAM header override, `--output-fmt sam|bam`, `--no-PG`, `--write-index` BAI generation for BAM outputs, and `-p N` padding. **Pending:** CRAM, sorted-by-tag streaming mode, and deeper upstream `@PG` byte-parity for complex chains.
- [~] `reheader` (`bam_reheader.c`) — basic SAM/BAM header replacement (record-level rewrite) with default `@PG` insertion, `--no-PG` suppression, and `-c <command>` external header filtering. **Pending:** BGZF block-level BAM fast path and CRAM `--in-place`.
- [~] `addreplacerg` (`bam_addrprg.c`) — SAM/BAM add/replace `@RG` + `RG:Z`. `-r` now unescapes `\t`/`\n` so a full `@RG\tID:..\tCN:..` spec works; incremental `-r KEY:VAL`, `-R ID` (rejected if absent), default-first-`@RG`-ID, `-m overwrite_all|orphan_only`, `-w` edit, `-O sam|bam`, `-o`, `@PG`/`--no-PG`. **Byte-exact vs the whole upstream `test_addrprg` group** (`addrprg/{1,2,4,5}` + `-R` overwrite, modulo `@PG`; `addrprg/3` = expected `-R` failure); integration test `addreplacerg_matches_upstream_group`. **Pending:** CRAM input/output, mate-aware updates, full orphan-first semantics.
- [~] `fastq` / `fasta` / `bam2fq` (`bam_fastq.c`) — basic single-stream output works for SAM and BAM (records written to stdout, `-o FILE`, or `-0 FILE`), with `-f`/`--require-flags`, `--rf`/`--include-flags`, `-F`/`--exclude-flags`, `-G`, the upstream default `0x900` secondary/supplementary exclusion, read-name suffix controls (`-n`/`-N`), `-O` original-quality `OQ` tag output, `-v INT` missing-quality defaults for FASTQ, `-U`/`--UMI-tag` UMI read-name suffixes, `-i`/`--barcode-tag` CASAVA barcode fields, upstream-style name-grouped paired split outputs (`-1`/`-2`/`-s`/`-0`) that pick the best per-readpart record per qname-group and route R1+R2 to `-1`/`-2`, R1-only or R2-only singletons to `-s` (falling back to `-1`/`-2` when `-s` is absent), and READ_OTHER to `-0` (falling back to `-s` when `-0` is absent), per-record interleaved output when `-1` and `-2` paths alias to the same file, SAM/BAM selected aux comments via `-T` in single and split output modes, all-tag SAM/BAM comments via `-T ''` / `-T '*'` in single and split-output FASTQ mode, SAM/BAM `B` array aux comment formatting, SAM/BAM single and split-output FASTQ tag filtering via `-d`/`--tag TAG[:VALUE]` and `-D`/`--tag-file TAG:FILE`, accumulating `-t` (`RG,BC,QT` upstream shortcut) and `-T TAG,...` selections that union rather than override, repeated `-d TAG[:VAL]` / `-D TAG:FILE` invocations that union value sets for the same tag and reject mismatched tags, FASTA reverse-complement of reverse-strand records, per-record `--i1`/`--i2` index FASTQ extraction with `--index-format` (default `i*i*`), `--quality-tag` (default `QT`), and `--barcode-tag`, and the upstream SAM-input all-tags fixture `bam2fq/15.fq.expected` matches for `-T ''` and `-t -T '*'`. `bam2fq/{1,2,3,4,6,7,9,11,13,15,16,17,18,19,20}.{1,2,s}.fq.expected` and `bam2fq/11.fa.expected` fixtures pass against the current Rust binary. **Pending:** exact upstream name-grouped one-record-per-qname index emission, index emission for stdin input, CASAVA paired-end barcode propagation, exact upstream behavior for `-i`/`--index-format` interaction with split mode, CRAM.
- [~] `import` (`bam_import.c`) — basic single FASTA/FASTQ and paired FASTQ (`-1`/`-2`, `--r1`/`--r2`, `-s` interleaved, plus two positional inputs) → SAM/BAM (`-O bam` / `--bam`), including positional single input plus `-0` single-read alias, `-0` singleton input alongside paired `-1`/`-2`, positional interleaved FASTQ detection from `/1`/`/2` read names, no-op `--no-PG`, CASAVA parsing (`-i`) with upstream-style reverse comments, SRA name2 (`-N`), UMI extraction (`-U`/`--UMI-tag`) with reverse comments, CASAVA barcode sequence tags (`--barcode-tag`), FASTQ definition aux tags (`-T`) including upstream-style float exponent spelling, explicit index reads (`--i1`/`--i2`) for `-0`, `-s`, positional interleaved, and paired `-1`/`-2` inputs with barcode sequence/quality tags (`--barcode-tag`/`--quality-tag`) and `-b`, and read-group header/tag support (`-R`/`-r`) with repeated `-r` accumulation, `-r` precedence over `-R`, and `-r` ID validation. Direct comparisons against `test/import/*.expected.sam` for the currently implemented import fixture commands pass. **Pending:** full paired singleton/other grouping parity, full read-group parity, CRAM output.

### Wave C — Editing / Mate-aware

- [~] `fixmate` (`bam_mate.c`) — basic mate flag/pos fixup for adjacent paired records in name-sorted BAM and SAM inputs (`FMUNMAP`, `FMATE_REVERSE_COMPLEMENTED`, `mate_reference_sequence_id`, `mate_alignment_start`) and rejects `@HD SO:coordinate` input like upstream. TLEN is recalculated from mate 5-prime positions, including large coordinate inputs where the resulting template length still fits. Default MC/MQ mate aux tags are added for mapped mates and cleared when the mate is unmapped. `-m` adds `ms:i` mate-score tags for markdup. `-c` adds lowercase template-CIGAR `ct:Z` tags to the earlier mapped mate and clears stale `ct` tags from both mates. `-z`/`--sanitize` parses/validates through the shared sanitizer parser. A `-` output operand now means stdout. Aux updates use an order-preserving `aux_del`/`aux_set_append` (mirroring HTSlib `bam_aux_del`+`bam_aux_append`: drop-then-append-to-tail), MQ before MC, and `MC:Z:*` is added when *either* read is mapped (`bam_mate.c:197`); the `-c` `ct` removal is order-preserving too. Default `@PG`/`--no-PG` supported; `-r` matches upstream `remove_reads`. **Byte-exact vs the entire upstream `test_fixmate` group** (`fixmate/{2,3,4,5,6,7,8}*` + `sanitize`, modulo `@PG`); integration test `fixmate_matches_upstream_group`. **Pending:** CRAM, mate-rescore, base-modification `-M`.
- [x] `markdup` (`bam_markdup.c`, 89k) — **faithful upstream key/score port**. SE/PE duplicate marking for SAM/BAM. PE reads build the upstream `make_pair_key` (template default + `--mode s` sequence; unclipped coords from CIGAR & `MC` tag; `R_LE`/`R_RI` left/right discriminator so a template's two mates get distinct keys and only corresponding mates of duplicate templates collide) plus a shared `make_single_key`; the kept read of a colliding key is the one with the higher `calc_score` = Σ(base qual ≥ 15) + `ms` mate-score tag, with the QCFAIL-asymmetry override and qname `strcmp` tie-break. `-S` seeds a qname `dup_hash` from marked-duplicate reads carrying `SA`/`XA` or an unmapped mate and flags matching supplementary/secondary/unmapped records (gated on `-S`, as upstream). `-b`/`--barcode-tag`, `-c`, `-t` `do` tags, `-d` `dt:Z:SQ|LB` with the full `find_duplicate_chains` optical re-tagging (per-read `original`/`duplicate` chain links + `check_chain_against_original` + `check_duplicate_chain`), `get_coordinates_colons` optical-name parse, `--use-read-groups` (rg-keyed), `--duplicate-count` (`dc:i`), `--include-fails`, `-m`/`--mode t|s`, `-r`, `-s`, `-O`, `-o`, regex `--read-coords`/`--coords-order`/`--barcode-rgx`/`--barcode-name` (via the `regex` crate; capture-span-bounded coord parse), raw-header SAM output (preserves input `@RG`/`@SQ` order), `@PG`/`--no-PG`. **Byte-exact vs the ENTIRE upstream `test_markdup` SAM harness — `markdup/{5..18}.expected.sam` (all 14 fixtures)**; test `markdup_matches_upstream_test_markdup_fixtures`. **Pending:** exact `-s` stats counts, CRAM, the `1..4` expect-fail error-message cases.
- [~] `rmdup` (`bam_rmdup.c` + `bam_rmdupse.c`) — single-end and paired-end duplicate removal for BAM and SAM inputs. SE records are keyed by `(tid, pos, reverse-flag)`; PE records pair by qname and are keyed by the canonical pair of `(tid, pos, strand)` triples, retaining the highest MAPQ/combined MAPQ record or pair. `-s`/`-S` force single-end treatment. Default `@PG` insertion via `pg::add_samtools_pg_to_header` and `--no-PG` are supported. **Pending:** CRAM, full upstream deprecated-command parity.
- [~] `calmd` / `fillmd` (`bam_md.c`) — SAM, BAM, and reference-backed CRAM input can emit SAM text with MD/NM tags recomputed against a FASTA reference via CIGAR/reference walking. BAQ paths (`-r`, `-r -e`, `-E`) are wired through `htslib_rs::alignment_compat::recalculate_baq_*` and `apply_existing_baq_from_sam_path` for SAM input, and `-d` drops existing `BQ` tags from the SAM-text output. Default `@PG` insertion via `pg::add_samtools_pg` (text-level) and `--no-PG` are supported. **Pending:** BAM/CRAM output, BAM/CRAM BAQ paths, `-A`/`-C cap`, full upstream MD/BAQ parity.
- [ ] `targetcut` (`cut_target.c`) — fosmid pool target cutting.
- [x] `reset` (`reset.c`) — strip alignment fields (`tid`/`pos`/`cigar`/`mate_*`/`template_length`) for BAM and SAM inputs, set MAPQ to `0`, drop a default set of aligner aux tags (NM, MD, AS, XS, SA, MC, MQ, NH, HI, ms), clear `PROPER_PAIR`/`SECONDARY`/`SUPPLEMENTARY`/`REVERSE`/`MATE_REVERSE`, set `UNMAPPED`, set `MATE_UNMAPPED` for paired reads, reverse-restore reverse-strand sequence/quality, preserve duplicate flags with `--dupflag`, remove read-group headers/tags with `--no-RG`, remove program header chains with `--reject-PG`, add a new samtools `@PG` chain entry by default (via the shared `pg::add_samtools_pg_to_header` helper), suppress the new `@PG` with `--no-PG` while preserving existing entries (matching upstream's `noPGentry` semantics), accept SAM/BAM input from stdin/no positional input/`-`, and tolerate legacy SAM `@HD VN:1` headers. `-x`/`--keep-tag` honored, with `--no-RG` taking precedence over keeping `RG`. **Order-preserving aux drop** (input field order, HTSlib `bam_aux_del` semantics) + **raw-header SAM output** (keep `@HD`, drop `@SQ`/`@CO`, `@RG` verbatim, `--reject-PG` removes the named `@PG` + its `PP`-chain descendants) + format inferred from the `-o` extension (`sam_open_mode`: SAM unless `.bam`). `--reject-PG` uses the upstream positional rule (`reset.c:223`: keep `@PG` until the first matching `ID`, drop it and all subsequent `@PG`). **Byte-exact vs every tested upstream `reset` fixture: `basic.1.mp.1` (reset\|view, stdin, file), `basic.output.mp.1` (`-o` SAM in), `basic.bam.input` (`-o` BAM in), `output.nRG.1` (`--no-RG`), `reject.1`, `reject.2`** (harness `hskip=1` + `ignore_pg_header`); test `reset_matches_upstream_test_reset_fixtures`. **Pending:** CRAM I/O.
- [~] `ampliconclip` (`bam_ampliconclip.c`, 40k) — **faithful port**. Per-reference BED primer sites (sorted by `right`), `matching_clip_site` (binary-search + `--tolerance`/`--strand` overlap pick), `bam_trim_left`/`bam_trim_right` soft/hard clip (CIGAR/POS/SEQ/QUAL rewrite, hardclip merge, full-consume→empty), `active_query_len`-gated `--filter-len`/`--fail-len`/`--unmap-len`, `--both-ends`, `--original` (`OA` tag), `--keep-tag` (default deletes `NM`/`MD`, order-preserving), `--clipped`, `--no-excluded`, `--rejects-file`, `--primer-counts` TSV, `-f` stats, `-o`/`-O sam|bam`, `-b`, raw-header `@HD SO:coordinate→unknown`, default `@PG`/`--no-PG`. **Byte-exact vs the entire upstream `test_ampliconclip` harness** (10 SAM fixtures + 3 primer-counts TSVs); test `ampliconclip_matches_upstream_test_ampliconclip_fixtures`. **Pending:** CRAM, BGZF block fast path, the unused `3_multi_ref_both_clip` edge.

### Wave D — Stats / Pileup

- [~] `depth` (`bam2depth.c`) — per-position depth via CIGAR walks for SAM, BAM, and reference-backed indexed CRAM. `-a`/`-aa`/`-d`/`-q`/`-o`, `-H` header output, `-f` input file lists, flag filters (`-g`, `-G`/`--excl-flags`, `--incl-flags`, `--require-flags`), `-l` minimum read length filtering, `-r` region restriction, `-b` BED restriction, and multi-input columnar output are supported. Sparse per-position storage (no OOM on huge references). **Byte-exact** vs the upstream `large_pos` fixtures `depth.expected.out` and `depth_bed.expected.out` (deletions excluded like pileup; bedidx now whitespace-tolerant for space-delimited BEDs). Tests cover SAM/CRAM region depth, multi-input columns/header/list files, flag filters, read-length filtering, large-reference depth + BED, and clean missing-reference failure. **Pending:** pileup-based overlap handling for paired double-counting, CRAM input without an explicit reference.
- [~] `coverage` (`coverage.c`) — per-reference/`-r` region `numreads`, `covbases`, `coverage`, `meandepth`, `meanbaseq`, and `meanmapq` via CIGAR walks for SAM, BAM, and reference-backed indexed CRAM. `--min-depth` thresholds covered-base counts, `-Q`/`--min-BQ` filters low-quality bases, `-q`/`--min-MQ` filters reads, `-b`/`--bam-list` expands input filename lists, `--ff`/`--excl-flags` replaces the default filter-out flags, `--rf`/`--incl-flags` requires at least one selected flag, `-l`/`--min-read-len` filters short alignments, `-d` caps per-position depth for reported coverage/depth metrics, and multiple inputs aggregate into one row per reference/region. `-m`/`--histogram` (with `-A`/`--ascii` and `-D`/`--plot-depth` routed through the same path) emits a 10-row ASCII histogram with `-w`/`--n-bins` controlling the column count. **Byte-exact** vs the entire upstream `test_coverage` tabular suite (`coverage/{1..5}.expected`, incl. multi-input and `-Q`/`-q`): C `printf %g`/`%.3g` formatting (`c_printf_g`, `coverage.c:211`), `min_depth`-gated `meandepth`/`meanbaseq` accumulators (per-position baseq vecs), and pileup-arrival reference row ordering. Test `coverage_matches_upstream_tabular_fixtures`. **Pending:** byte-parity for the UTF-8 box-drawing histogram + sidebar, `-D` true depth-plot, CRAM without explicit reference.
- [~] `bedcov` (`bedcov.c`) — total aligned-base coverage per BED region, walking each record's CIGAR for SAM, BAM, and reference-backed indexed CRAM. `-Q` mapq filter, `-g`/`-G` filter-mask controls, `-j` deletion/refskip skipping, `-H` output headers, `-c` read-count columns, and `-d` depth-threshold columns are supported. **Byte-exact** vs all four upstream `test_bedcov` fixtures (`bedcov/bedcov{,_j,_gG,_c}.expected`), including attached `-g512 -G2048`. **Pending:** none for the tabular suite; CRAM without explicit reference.
- [~] `stats` (`stats.c` + `stats_isize.c`, 123k + 8k) — `SN` (Summary Numbers) section plus FFQ/LFQ first/last fragment quality histograms, GCF/GCL first/last fragment GC histograms, and approximate CIGAR-walk COV coverage histograms for SAM, BAM, and reference-backed CRAM region paths, including record-backed `-I`/`--id` read-group/sample filtering, `-f`/`--required-flag`, `-F`/`--filtering-flag`, `-i`/`--insert-size` insert-size cap, `-m`/`--most-inserts` insert-size bulk selection, `-l`/`--read-length`, `-q`/`--trim-quality` BWA trim counting, `-c`/`--coverage MIN,MAX,STEP` COV binning, and `-g`/`--cov-threshold` target percentage lines with target-region validation. SAM and BAM iterate records directly to populate sequence-length, quality, GC, CIGAR, NM, COV, and runtime coordinate-order accumulators (CRAM non-region remains on the `summarize_*` summary path and therefore reports zeros for those fields, with `is sorted` still falling back to header `SO`). The emitted lines now cover: raw total / filtered / sequences / runtime is sorted / 1st & last fragments / mapped / mapped+paired / unmapped / properly paired / paired / duplicated / MQ0 / QC-failed / non-primary / supplementary / total length / total first & last fragment length / bases mapped / bases mapped (cigar) / mismatches (NM aux) / error rate / average length & per-fragment / maximum length & per-fragment / bases trimmed / average quality / singletons / insert size mean & stddev / inward, outward, and other oriented pair counts / pairs on different chromosomes / percentage of properly paired reads / target bases / target genome coverage above threshold. SAM, indexed BAM, and reference-backed CRAM positional region arguments and `-t` target files restrict the summary and COV positions, with overlapping BAM/CRAM regions de-duplicated. `-d` / `--remove-dups` filters duplicate-marked primary records and their quality/GC/COV histogram contributions. Missing CRAM references fail cleanly. **Pending:** exact pileup-backed COV byte parity, per-cycle, BAQ, CRAM non-region read-group/read-length/trim-quality filtering, and CRAM non-region sequence-length/quality/GC/COV stats (requires a CRAM all-record iterator in `htslib-rs`).
- [~] `mpileup` (`bam_plcmd.c`, 49k) — default **text** pileup implemented on the new `htslib-rs::alignment_compat` pileup iterator (`pileup_from_alignment_paths_with[_reference]_and_options` + `PileupColumn`/`PileupRead`). Supports multi-input + `-b` list (incl. `file://`), `-f` reference (plain or bgzipped, for the ref-base column and CRAM decode), `-r region` (incl. attached `-r17:1-2` form), `-Q`/`--min-BQ` (default 13), `-q`/`--min-MQ`, `--ff`/`--excl-flags` (mask **replace**, default `0x704`), `--rf`/`--incl-flags`, `-A`/`--count-orphans`, `-x`/`--ignore-overlaps`, `-o`. Faithful `pileup_seq` byte encoding (`.`/`,`, upper/lower mismatch, `^`+mapq head, `$` tail, `*` deletion, `<`/`>` ref-skip, `+`/`-` indels), HTSlib `MPLP_SMART_OVERLAPS` overlap removal and `MPLP_NO_ORPHAN` orphan filter (both in the htslib-rs pileup engine), and the `[mpileup] N samples in M input files` stderr line. **Byte parity:** upstream `mpileup.out.3` (`-B --ff 0x14`) and `mpileup.out.5` (overlap) match exactly; `mpileup.out.1` matches on depth + read bases (a few quality chars differ only where HTSlib applies BAQ). Integration tests `mpileup_minus_b_ff_matches_upstream_out3`, `mpileup_overlap_removal_matches_upstream_out5`. **Pending:** BAQ-adjusted qualities (TODO-NEXT #11), `@RG`-`SM` sample grouping (currently one sample per file), VCF/BCF (`-g`/`-v`) output, `-a`/`-aa` all-positions, base-modification columns, per-position mods/qpos/qname extra columns, CRAM-without-index-via-region efficiency.
- [x] `consensus` (`bam_consensus.c`, 126k, + `consensus_pileup.c`) — **byte-exact vs ALL 77 `test/consensus/consensus.reg` cases**: both `--mode simple` (freq/score) and the default `bayesian`/`recall` Gap5 model (`calculate_consensus_gap5` + the `consensus_init` probability tables / `fast_exp`/`fast_log2`/`ph_log` math, fed by the htslib-rs pileup `nm_init` precompute `PileupRead::bayes_poly`/`bayes_nm_local`), fasta/fastq/pileup, `-a`/`-aa`, `-r`, `-T`/`--ref-qual`, `--min-MQ`/`--min-BQ`, show-del/ins, glued short options. In-process harness test `consensus_matches_upstream_consensus_reg` (77/77).
- [ ] `phase` (`phase.c`) — heterozygote phasing.
- [~] `depad` / `pad2unpad` (`padding.c`) — SAM input with `-T` padded FASTA reference and `-s` SAM output converts padded reference columns to unpadded coordinates/CIGAR (`I`/`P`) and matches the upstream `depad.001` fixture with `--no-PG`. **Pending:** BAM input/output, CRAM, binary output modes (`-u`/`-1`), and full upstream `test_depad` parity.
- [~] `ampliconstats` (`amplicon_stats.c`, 1776 LOC) — **faithful port**. Per-ref BED (file order, strand), `count_amplicon`/`bed2amplicon`, ±`pos-margin` position→amplicon lookup, `accumulate_stats` (flag filter, qname read-pair overlap removal, `depth_all`/`depth_valid`, `nreads`/`nbases`/`coverage`, `amp_dist` via TLEN±`tlen-adjust`, `tcoord` freq map), `append_lstats` (sum+sum² s.d.), full `dump_stats` (`SS`/`AMPLICON`/`F·`/`C·` incl. `depth_bin` RLE, `TCOORD` ≥ `tcoord-min-count`, `FAMP`, `COMBINED` MEAN/STDDEV). `-S`/`-c`/`-t`/`-d d1,d2,d3`/`-m`/`-D`/`-b`/`-a`/`-l`/`-f`/`-F`/`-o`. **Byte-exact vs the entire upstream `test_ampliconstats` harness** (`stats`, `stats_mixed`, `stats_partial`, modulo the harness-stripped version/command-line lines); test `ampliconstats_matches_upstream_test_ampliconstats_fixtures`. **Pending:** `--tcoord-bin` aggregation, CRAM, `--use-sample-name`.
- [x] `cram-size` (`cram_size.c`) — **byte-exact vs the entire upstream `test/cram_size/cram_size.reg`** (`normal.out`, `verbose.out`, `encodings.out`). Faithful `cram_expand_method`/`comp_method2expanded` method decoder + verbatim tables, `Container::blocks()` walk, `cram_cid2ds` content_id→data-series map, normal (aggregate-by-cid) / `-v` (by cid+method) reports with the exact `BLOCK …%6.2f%% %-Ns` formatting + ratio/`>999%`/summary, and `-e` `cram_describe_encodings`/`cram_codec_describe` with htslib's exact DS + `tag_encoding_map` ordering. Built on the vendored-noodles CRAM inventory surface (`CompressionHeader` encodings public, `Container::blocks()`, ordered `TagEncodings`). Test `cram_size_matches_upstream_cram_size_reg` (all 3).
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

## htslib-rs Extensions Needed (rolling list) — ✅ ALL RESOLVED

**Every entry below is DONE.** Tracked and completed as `TODO-NEXT.md`
items #1–#12 (all 12 byte- or fixture-verified, committed, with the
htslib-rs / vendored-noodles-fork pins bumped). The "CRAM internals"
and "CRAM all-record iterator" gaps were closed via minimal patches to
the **owned vendored noodles fork** (`madhavajay/noodles`, htslib-rs
submodule) per the Ground rule's "(and carry minimal patches)" clause.
Kept for historical traceability:

- [x] **`sam_hdr_add_pg`** — binary `@PG` for `view` SAM→BAM/CRAM + BAM→BAM (TODO-NEXT #4).
- [x] **`bam_aux_update_*`** — aux mutation via `RecordBuf`; all parity consumers byte-exact (TODO-NEXT #7).
- [x] **`sam_pileup` / `bam_plp_*` API surface** — pileup iterator (TODO-NEXT #1); unlocked mpileup/consensus/coverage/bedcov/depth.
- [x] **`hts_set_threads`** — `-@` accepted, byte-identical output (TODO-NEXT #8).
- [x] **`auto_index` / index save during write** — `sort`/`view`/`merge --write-index` BAI == post-pass, byte-verified (TODO-NEXT #9).
- [x] **Indexing BAMs without `@HD SO:coordinate`** — DONE (TODO-NEXT #6, htslib-rs `530b27c`).
- [x] **CRAM internals for `cram-size`** — full inventory surface in the vendored noodles fork; `cram-size` byte-exact on all 3 fixtures (TODO-NEXT #3).
- [x] **`htslib-rs::region`** — `*`/`.` grammar done (TODO-NEXT #10).
- [x] **`probaln_glocal` and BAQ recalculation** — verified + `calmd` BAM output (TODO-NEXT #11).
- [x] **CRAM all-record iterator** — `query_cram_records_all_from_path[_with_reference]`; `stats`/`checksum` no-region CRAM, `reference` MD path, embed_ref read+write (TODO-NEXT #2).

## noodles Extensions Needed (rolling list)

Keep these at the end during the current samtools-only pass. Do not modify the
noodles submodule for these blockers until explicitly switching back to
underlying-library work.

- [x] **CSI query robustness for very large references/regions** — **DONE** (TODO-NEXT #12, htslib-rs `8372873`). Root cause was `build_bam_csi_with_min_shift` using a fixed CSI depth of 5; it now auto-sizes depth from the largest reference via `alignment_csi_depth_for_header` (as the SAM-CSI builder already did). `samtools view large_chrom.bam ref2` and `ref2:1-541556283` are byte-exact vs `dat/large_chrom.out`, no panic / no `invalid end bound`. Fixed entirely in `htslib-rs` — **noodles unpatched**.
- [~] **SAM aux float formatting (`f:` scalars and `B:f` arrays)** — RESOLVED via samtools-rs option (b): `samtools_rs::sam_render` reuses the htslib-style `format_aux_float` / `format_htslib_exponent` and adds `fix_sam_aux_floats` / `fix_sam_text` (post-process noodles SAM text) plus `write_record` / `write_header` (drop-in for noodles `sam::io::Writer`). Wired into `view` (all binary→SAM text paths + `record_to_sam_line`) and `split` SAM output, bringing `reheader/1` (via `… | view -h`) and `split.expected.grp{1,2}.sam` to byte parity. Remaining commands (`reheader` SAM→SAM, `sort`, `merge`, `collate`, `addreplacerg`, `reset`, `fixmate`, `rmdup`, `markdup`, `cat`) still write through noodles directly — tracked as a bounded follow-up in the *Remaining tractable items*. The noodles-side option (a) is no longer required.

## Submodule Pinning

- [x] Pin `samtools/` to a specific upstream release tag once Phase 0 lands (record tag + commit in `README.md` and `version.rs`). Current pin: upstream tag `1.23.1`, commit `6efb9b6da35224cf804921dedecf9fb8f411365d`.
- [x] Pin `htslib-rs/` to a known-green commit when Phase 0 lands. Current pin: `83728732f88b53e0db817583e0fc7157cb6795d6` (TODO-NEXT #1 pileup engine + #2 whole-CRAM iterator + #6 build_bai no-SO + #12 header-aware BAM-CSI depth; prior: `530b27c`, `ca812dd`, `9cf30b3`, `e25f3929`, `5b25622`, `da4d3319`, `6bd6fb0`, `88bd29f`).

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
