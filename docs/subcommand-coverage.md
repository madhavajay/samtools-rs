# Subcommand Coverage Map

For each upstream `samtools` subcommand, this document lists:

- the C source file(s) that implement it (under `samtools/`)
- the HTSlib APIs it depends on
- the corresponding `htslib-rs` API coverage status
- the samtools-rs implementation status

`htslib-rs` API status legend:
- ✅ already exposed in `htslib-rs` (specify module)
- ⚠️ partial — covers some but not all behavior we need
- ❌ missing — must be added to `htslib-rs` first
- 🟦 not an HTSlib-shaped API (would call `noodles` directly)

samtools-rs status legend:
- ✅ implemented + passes targeted tests
- 🟡 stubbed; partial implementation lands here
- ⬜ not yet implemented

## Wave A — Read/Write/Index Foundation

### view

- **C source:** `sam_view.c` (~68k LOC)
- **HTSlib APIs used:** `sam_open_format`, `sam_index_load2`, `sam_hdr_read`, `sam_hdr_add_pg`, `sam_read1`, `sam_write1`, `hts_set_threads`, `hts_set_opt`, `hts_itr_querys`, `bam_str2flag`, `bam_aux_*`, `hts_expr_*`, `fai_path`
- **htslib-rs coverage:** ⚠️
  - `htslib_rs::format::detect_path` — ✅
  - `htslib_rs::alignment_compat::write_bam_from_path` / `write_bam_from_sam_path` — ✅
  - `htslib_rs::alignment_compat::write_cram_from_sam_path_with_reference` — ✅
  - Region queries via `htslib_rs::region::parse_region` + `query_*_records_from_path` — ✅
  - `sam_hdr_add_pg` equivalent — ❌ (needed for non-`--no-PG` paths)
  - `hts_expr_*` filter evaluation — ✅ via `htslib_rs::expr` (audit needed)
  - Aux-tag filter / remove (`-x`, `--keep-tag`) — needs `bam_aux_*` helpers — ⚠️
- **samtools-rs status:** 🟡 (basic SAM↔SAM and SAM→BAM paths work; many flags and region/BED queries TODO)

### head

- **C source:** `sam_view.c::main_head` (lines 1772+)
- **HTSlib APIs used:** `sam_open_format`, `sam_hdr_read`, `sam_hdr_str`, `sam_read1`, `sam_format1`
- **htslib-rs coverage:** ⚠️
  - Raw header byte access — handled locally via `samtools_rs::header_text` (BAM read-text helper); could be promoted to `htslib-rs`
- **samtools-rs status:** ✅ (SAM and BAM input; CRAM and stdin TODO)

### quickcheck

- **C source:** `bam_quickcheck.c`
- **HTSlib APIs used:** `hts_open`, `hts_get_format`, `sam_hdr_read`, `sam_hdr_nref`, `hts_check_EOF`, `hts_close`
- **htslib-rs coverage:** ⚠️
  - `htslib_rs::format::detect_path` — ✅
  - `htslib_rs::alignment_compat::read_*_header_from_path` — ✅
  - BGZF/CRAM EOF marker check — implemented locally; should be promoted to `htslib-rs::format` or `htslib-rs::bgzf_compat`
- **samtools-rs status:** ✅ (parity output matches `quickcheck/all.expected`)

### index

- **C source:** `bam_index.c`
- **HTSlib APIs used:** `sam_index_build3`, `sam_index_build`, `sam_index_save`, `hts_set_threads`, `bgzf_check_EOF`
- **htslib-rs coverage:** ✅ via `htslib_rs::index_compat::build_bai` / `build_bam_csi` / `build_cram_crai`
- **samtools-rs status:** ⬜

### idxstats

- **C source:** `bam_index.c::bam_idxstats`
- **HTSlib APIs used:** `sam_open_format`, `sam_hdr_read`, `sam_index_load2`, `hts_idx_get_stat`, `hts_idx_get_n_no_coor`, `sam_hdr_tid2name`, `sam_hdr_tid2len`
- **htslib-rs coverage:** ⚠️
  - `read_associated_bam_index` returns `Box<dyn csi::BinningIndex>` — ✅
  - Per-reference mapped/unmapped counts from index meta — needs accessor — ❌ (must extend `htslib-rs::index_compat`)
- **samtools-rs status:** ⬜

### faidx / fqidx

- **C source:** `faidx.c`
- **HTSlib APIs used:** `fai_build`, `fai_load`, `fai_fetch`, `faidx_fetch_seq64`, `fai_destroy`
- **htslib-rs coverage:** ✅ via `htslib_rs::faidx_compat::build_index`, `read_index`, `write_index`, `fetch_sequence`, `fetch_region_sequence`
- **samtools-rs status:** ⬜

### dict

- **C source:** `dict.c`
- **HTSlib APIs used:** `gzopen`, `kseq_*`, `hts_md5_*` — only HTSlib's md5 wrapper is used; everything else is its own FASTA reader
- **htslib-rs coverage:** 🟦 (use `noodles-fasta` for FASTA iteration + `md-5` crate for MD5)
- **samtools-rs status:** ⬜

### flagstat

- **C source:** `bam_stat.c`
- **HTSlib APIs used:** `sam_open_format`, `sam_hdr_read`, `sam_read1`, `hts_set_threads`, `hts_set_opt` (CRAM_OPT_REQUIRED_FIELDS)
- **htslib-rs coverage:** ⚠️
  - Full BAM record iteration with access to `flag`, `tid`, `mtid`, `qual` — not yet exposed via a public accessor on `AlignmentRecordSummary` (private fields). Must extend `htslib-rs::alignment_compat` to expose accessors or add a streaming record iterator. — ❌
- **samtools-rs status:** ⬜

## Wave B — File Operations

### sort

- **C source:** `bam_sort.c` (~138k LOC — largest file in samtools)
- **HTSlib APIs used:** `sam_open_format`, `sam_index_*`, `sam_read1`, `sam_write1`, `sam_hdr_*`, `hts_set_threads`, BGZF temp file I/O, lz4 compression for temps
- **htslib-rs coverage:** ⚠️
  - Streaming write — ✅
  - Custom per-record sort key extraction — needs accessor APIs — ❌
  - Multi-way merge — ❌
- **samtools-rs status:** ⬜

### merge

- **C source:** `bam_sort.c::bam_merge_core`
- **HTSlib APIs used:** as above
- **samtools-rs status:** ⬜

### collate / bamshuf

- **C source:** `bamshuf.c`
- **HTSlib APIs used:** `sam_open_format`, `sam_read1`, `sam_write1`, BGZF temp file I/O
- **samtools-rs status:** ⬜

### cat

- **C source:** `bam_cat.c`
- **HTSlib APIs used:** BGZF passthrough — fast concatenation without re-decompression
- **htslib-rs coverage:** ⚠️ — BGZF block-level concatenation is not yet a first-class API in `htslib-rs::bgzf_compat`. Must be added. — ❌
- **samtools-rs status:** ⬜

### split

- **C source:** `bam_split.c`
- **HTSlib APIs used:** `sam_open_format`, `sam_read1`, `sam_write1`, `sam_hdr_*`, RG-tag inspection
- **samtools-rs status:** ⬜

### reheader

- **C source:** `bam_reheader.c`
- **HTSlib APIs used:** BAM/CRAM header rewriting, in-place mode for CRAM
- **htslib-rs coverage:** ❌ — in-place CRAM header replacement is not exposed
- **samtools-rs status:** ⬜

### addreplacerg

- **C source:** `bam_addrprg.c`
- **HTSlib APIs used:** `sam_hdr_add_line`, `bam_aux_update_str`, record iteration
- **htslib-rs coverage:** ❌ — `bam_aux_update_str` not yet exposed
- **samtools-rs status:** ⬜

### fastq / fasta / bam2fq

- **C source:** `bam_fastq.c` (~48k LOC)
- **HTSlib APIs used:** record iteration, aux tag access for barcodes/QT/RX/QX
- **htslib-rs coverage:** ⚠️ — `view_sam_as_fastq_text_from_path_with_limit` exists; full feature set (paired/single, barcode-aware) needs more
- **samtools-rs status:** ⬜

### import

- **C source:** `bam_import.c`
- **HTSlib APIs used:** FASTA/FASTQ reading + SAM/BAM/CRAM writing
- **samtools-rs status:** ⬜

## Wave C — Editing / Mate-aware

### fixmate

- **C source:** `bam_mate.c` (~43k LOC)
- **HTSlib APIs used:** record iteration, mate-flag/pos rewriting, MC/MQ tag updates
- **htslib-rs coverage:** ⚠️ — needs `bam_aux_*` and mutable record APIs
- **samtools-rs status:** ⬜

### markdup

- **C source:** `bam_markdup.c` (~89k LOC)
- **HTSlib APIs used:** name-sort grouping, position-sort grouping, barcode parsing, duplicate marking via flag updates
- **htslib-rs coverage:** ⚠️ — needs mutable record APIs and indexed I/O
- **samtools-rs status:** ⬜

### rmdup

- **C source:** `bam_rmdup.c` + `bam_rmdupse.c`
- **samtools-rs status:** ⬜

### calmd / fillmd

- **C source:** `bam_md.c`
- **HTSlib APIs used:** record iteration, BAQ via `probaln_glocal`, MD/NM recomputation
- **htslib-rs coverage:** ✅ partial via `htslib_rs::probaln` and `htslib_rs::alignment_compat::recalculate_baq_*`
- **samtools-rs status:** ⬜

### targetcut

- **C source:** `cut_target.c`
- **samtools-rs status:** ⬜

### reset

- **C source:** `reset.c`
- **HTSlib APIs used:** record iteration, aux-tag stripping, flag/CIGAR/pos resets
- **samtools-rs status:** ⬜

### ampliconclip

- **C source:** `bam_ampliconclip.c` (~40k LOC)
- **HTSlib APIs used:** record iteration, CIGAR rewriting, BED parsing
- **samtools-rs status:** ⬜

## Wave D — Stats / Pileup

### depth

- **C source:** `bam2depth.c`
- **HTSlib APIs used:** pileup iteration (`bam_plp_*`)
- **htslib-rs coverage:** ❌ — `bam_plp_*` API surface needs to be exposed
- **samtools-rs status:** ⬜

### coverage

- **C source:** `coverage.c` (~30k LOC)
- **samtools-rs status:** ⬜

### bedcov

- **C source:** `bedcov.c`
- **samtools-rs status:** ⬜

### stats

- **C source:** `stats.c` (~123k LOC) + `stats_isize.c`
- **HTSlib APIs used:** record iteration, CIGAR analysis, base/qual histograms, GC bias, insert-size distribution
- **samtools-rs status:** ⬜

### mpileup

- **C source:** `bam_plcmd.c` (~49k LOC)
- **HTSlib APIs used:** multi-input pileup, VCF output, regions
- **htslib-rs coverage:** ❌ — multi-input pileup iterator not yet exposed
- **samtools-rs status:** ⬜

### consensus

- **C source:** `bam_consensus.c` (~126k LOC) + `consensus_pileup.c`
- **samtools-rs status:** ⬜

### phase

- **C source:** `phase.c`
- **samtools-rs status:** ⬜

### depad / pad2unpad

- **C source:** `padding.c`
- **samtools-rs status:** ⬜

### ampliconstats

- **C source:** `amplicon_stats.c` (~65k LOC)
- **samtools-rs status:** ⬜

### cram-size

- **C source:** `cram_size.c`
- **HTSlib APIs used:** CRAM internal block/container/codec inspection
- **htslib-rs coverage:** ❌ — explicitly out-of-scope in `htslib-rs`. May need to drop this subcommand or expose minimal CRAM internals.
- **samtools-rs status:** ⬜

### checksum

- **C source:** `bam_checksum.c` (~47k LOC)
- **samtools-rs status:** ⬜

### samples

- **C source:** `bam_samples.c`
- **HTSlib APIs used:** `sam_hdr_*` to list `@RG` SM values
- **htslib-rs coverage:** ✅
- **samtools-rs status:** ⬜

### reference

- **C source:** `reference.c`
- **samtools-rs status:** ⬜

### flags

- **C source:** `bam_flags.c`
- **HTSlib APIs used:** `bam_str2flag`, `bam_flag2str`
- **htslib-rs coverage:** 🟦 — implemented locally in `samtools_rs::bam_flag` (could be promoted to `htslib-rs`)
- **samtools-rs status:** ✅

## htslib-rs Extensions Required (rolled up from above)

The following items are blockers for one or more samtools-rs subcommands. They should be added to `htslib-rs/TODO.md`:

1. **`sam_hdr_add_pg`** — programmatic `@PG` chain insertion with PP linkage. Required by every subcommand that writes output without `--no-PG`.
2. **`bam_aux_update_*`** — string/int/array aux updates. Required by `addreplacerg`, `fixmate`, `reset`.
3. **`AlignmentRecordSummary` accessors** — public getters for `flags`, `reference_sequence_id`, `mate_reference_sequence_id`, `mapping_quality`, `alignment_start`. Required by `flagstat`, `stats`, and many others. (Alternatively, expose a streaming record iterator that returns full `bam::Record`/`sam::RecordBuf`.)
4. **`hts_idx_get_stat` equivalent** — per-reference mapped/unmapped counts from `Box<dyn csi::BinningIndex>`. Required by `idxstats`.
5. **`bam_plp_*` pileup iterator** — multi-input pileup with overlap detection. Required by `mpileup`, `depth`, `consensus`, `coverage`.
6. **BGZF block-level concatenation** — fast `bam_cat` requires concatenating BGZF blocks without re-decompression.
7. **In-place CRAM header rewrite** — required by `samtools reheader --in-place` for CRAM files.
8. **`hts_set_threads` wiring** — propagate `-@ N` to BGZF worker count consistently.
9. **CRAM EOF marker check in `htslib-rs::format`** — currently duplicated in `samtools_rs::commands::quickcheck`.
10. **Raw BAM header text** — `read_bam_header_text_from_path` (duplicated in `samtools_rs::header_text`).
11. **CRAM internals for `cram-size`** — either expose minimal block/codec inspection or drop the subcommand.
12. **Region grammar coverage audit** — confirm `htslib_rs::region::parse_region` handles `*` (unmapped), `.` (everything), and the `chr:from-to` / `chr:from` forms.
