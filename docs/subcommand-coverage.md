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
  - Aux-tag filter / remove (`-x`, `--keep-tag`) — implemented for SAM output and SAM-input BAM/CRAM output; BAM/CRAM-input binary aux mutation still needs deeper mutable-record support — ⚠️
- **samtools-rs status:** 🟡 (SAM/BAM/reference-backed CRAM count/text/BAM/CRAM paths, stdin paths, region/BED queries, simple filters, expression filters, SAM-output `-U`/`-p`, and SAM-input aux stripping work; BAM/CRAM-input binary aux mutation, binary `-U`/`-p`, multi-file inputs, paired filters, and full CRAM parity remain)

### head

- **C source:** `sam_view.c::main_head` (lines 1772+)
- **HTSlib APIs used:** `sam_open_format`, `sam_hdr_read`, `sam_hdr_str`, `sam_read1`, `sam_format1`
- **htslib-rs coverage:** ⚠️
  - Raw header byte access — handled locally via `samtools_rs::header_text` and command-local stdin helpers; could be promoted to `htslib-rs`
- **samtools-rs status:** ✅ (SAM/BAM/CRAM file and stdin header/record output; reference-backed CRAM record extraction)

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
- **samtools-rs status:** 🟡 — `-Q`, `-g`/`-G`, `-H`, `-c`, and `-d` implemented with SAM/BAM/reference-backed CRAM input; exact pileup parity remains.

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
  - Custom per-record sort key extraction — ✅ partial for coordinate, query-name, and aux-tag keys
  - Multi-way merge — ❌
- **samtools-rs status:** 🟡 — in-memory coordinate, query-name, and aux-tag sort works for BAM and SAM inputs; external merge, template/minimiser sorts, write-index, thread/memory caps, and CRAM remain.

### merge

- **C source:** `bam_sort.c::bam_merge_core`
- **HTSlib APIs used:** as above
- **samtools-rs status:** 🟡 — in-memory coordinate and query-name merge works for BAM and SAM inputs, with `-R region` and `-L BED` restriction for indexed BAM; streaming k-way merge, header reconciliation for differing `@SQ`, and CRAM remain.

### collate / bamshuf

- **C source:** `bamshuf.c`
- **HTSlib APIs used:** `sam_open_format`, `sam_read1`, `sam_write1`, BGZF temp file I/O
- **samtools-rs status:** 🟡 — in-memory name grouping works for BAM and SAM inputs; on-disk hash bucketing, random seed, record cap, and CRAM remain.

### cat

- **C source:** `bam_cat.c`
- **HTSlib APIs used:** BGZF passthrough — fast concatenation without re-decompression
- **htslib-rs coverage:** ⚠️ — BGZF block-level concatenation is not yet a first-class API in `htslib-rs::bgzf_compat`. Must be added. — ❌
- **samtools-rs status:** 🟡 — record-level SAM and BAM concatenation works with `-o`, `-h`, default `@PG` insertion, `--no-PG`, and `-r region` for indexed BAM; BGZF block-level fast path, CRAM, and `-p` remain.

### split

- **C source:** `bam_split.c`
- **HTSlib APIs used:** `sam_open_format`, `sam_read1`, `sam_write1`, `sam_hdr_*`, RG-tag inspection
- **samtools-rs status:** 🟡 — basic `@RG` splitting with per-output `@RG` header filtering and default `@PG` insertion, plus explicit `-d TAG` string/integer aux-tag splitting, work for BAM and SAM inputs with `-f`, `-u`, `-h`, `--output-fmt sam|bam`, `--no-PG`, `--write-index` BAI generation for BAM outputs, `-M`, and `-p`; explicit `-d RG` can add missing output `@RG` headers. CRAM, sorted-by-tag streaming, and deeper upstream `@PG` byte-parity for complex chains remain.

### reheader

- **C source:** `bam_reheader.c`
- **HTSlib APIs used:** BAM/CRAM header rewriting, in-place mode for CRAM
- **htslib-rs coverage:** ❌ — in-place CRAM header replacement is not exposed
- **samtools-rs status:** 🟡 — record-level SAM/BAM header replacement works with default `@PG` insertion, `--no-PG`, and `-c <command>` external header filtering; BAM BGZF block-level fast path and CRAM `--in-place` remain.

### addreplacerg

- **C source:** `bam_addrprg.c`
- **HTSlib APIs used:** `sam_hdr_add_line`, `bam_aux_update_str`, record iteration
- **htslib-rs coverage:** ⚠️ — mutable SAM/BAM `RecordBuf` paths cover current RG string replacement; direct `bam_aux_update_str` parity and CRAM remain unavailable
- **samtools-rs status:** 🟡 — SAM/BAM add/replace exists with `-O sam|bam`, default `@PG` insertion, and `--no-PG`; CRAM, mate-aware updates, and full orphan-first semantics remain.

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
- **htslib-rs coverage:** ⚠️ — basic mutable record rewriting works through `RecordBuf`; direct `bam_aux_*` parity is still useful for deeper aux-tag behavior.
- **samtools-rs status:** 🟡 — basic adjacent name-sorted mate flag/reference/position fixup works for BAM and SAM inputs, including coordinate-sort rejection, TLEN recalculation, default MC/MQ mate tags, `-m` mate-score tags, `-c` template-CIGAR `ct` tags, and `-r`; mate rescore, sanitizer mutation, and CRAM remain.

### markdup

- **C source:** `bam_markdup.c` (~89k LOC)
- **HTSlib APIs used:** name-sort grouping, position-sort grouping, barcode parsing, duplicate marking via flag updates
- **htslib-rs coverage:** ⚠️ — mutable SAM/BAM `RecordBuf` paths cover current flag and aux-tag reads; indexed/streaming parity and CRAM remain.
- **samtools-rs status:** 🟡 — single-end and paired-end duplicate marking for SAM/BAM exists with optional barcode-key grouping (`-b`/`--barcode-tag`), duplicate flag/tag clearing (`-c`), `-S` compatibility for supplementary propagation, `-t` duplicate-origin `do` tags, `-d` optical-distance duplicate classification with `dt:Z:SQ`/`dt:Z:LB` tags, default QCFAIL exclusion with `--include-fails` override, validated `-m t|s`/`--mode t|s` compatibility, optical-aware estimated library size in `-s` stats, secondary/supplementary qname propagation, `-r`, upstream-shaped `-s` summary fields, `-O`, `-o`, default `@PG`, and `--no-PG`; exact stats/count parity and CRAM remain.

### rmdup

- **C source:** `bam_rmdup.c` + `bam_rmdupse.c`
- **samtools-rs status:** 🟡 — single-end and paired-end duplicate removal for BAM and SAM inputs exists, with `-s`/`-S`, default `@PG`, and `--no-PG`; CRAM and full deprecated-command parity remain.

### calmd / fillmd

- **C source:** `bam_md.c`
- **HTSlib APIs used:** record iteration, BAQ via `probaln_glocal`, MD/NM recomputation
- **htslib-rs coverage:** ✅ partial via `htslib_rs::probaln` and `htslib_rs::alignment_compat::recalculate_baq_*`
- **samtools-rs status:** 🟡 — SAM, BAM, and reference-backed CRAM input can emit SAM text with recomputed MD/NM tags against FASTA references; SAM input can also run BAQ paths, and `-d` drops existing `BQ` tags, with default `@PG`/`--no-PG`. BAM/CRAM output, BAM/CRAM BAQ paths, remaining flags, and full upstream MD/BAQ parity remain.

### targetcut

- **C source:** `cut_target.c`
- **samtools-rs status:** ⬜

### reset

- **C source:** `reset.c`
- **HTSlib APIs used:** record iteration, aux-tag stripping, flag/CIGAR/pos resets
- **samtools-rs status:** 🟡 — BAM and SAM reset paths clear alignment fields, default aux tags, and alignment-dependent flags; reverse-strand sequence/quality re-reversal, `-x`/`--keep-tag`, `--no-RG`, `--reject-PG`, `--dupflag`, default `@PG`, and `--no-PG` are supported. CRAM remains.

### ampliconclip

- **C source:** `bam_ampliconclip.c` (~40k LOC)
- **HTSlib APIs used:** record iteration, CIGAR rewriting, BED parsing
- **samtools-rs status:** ⬜

## Wave D — Stats / Pileup

### depth

- **C source:** `bam2depth.c`
- **HTSlib APIs used:** pileup iteration (`bam_plp_*`)
- **htslib-rs coverage:** ⚠️ — exact pileup parity still needs `bam_plp_*`; current implementation uses alignment record iteration/CIGAR walks for BAM and reference-backed CRAM.
- **samtools-rs status:** 🟡 — `-a`/`-aa`/`-d`/`-q`/`-o`, `-H`, `-f` input lists, flag filters, `-l` minimum read length filtering, `-r`, `-b`, and multi-input depth columns are implemented; overlap/deletion parity and CRAM without explicit reference remain.

### coverage

- **C source:** `coverage.c` (~30k LOC)
- **HTSlib APIs used:** pileup iteration (`bam_plp_*`), CIGAR/quality access
- **htslib-rs coverage:** ⚠️ — exact pileup parity still needs `bam_plp_*`; current implementation uses alignment record iteration/CIGAR walks for BAM and reference-backed CRAM.
- **samtools-rs status:** 🟡 — `numreads`, `covbases`, percent coverage, mean depth, mean base quality, mean map quality, `-r` regions, `--min-depth`, `-Q`/`--min-BQ`, read map-quality filtering, `-b`/`--bam-list`, `--ff`/`--excl-flags`, `--rf`/`--incl-flags`, `-l`/`--min-read-len`, `-d` maximum-depth capping, multi-input aggregate rows, and ASCII histogram mode are implemented; byte-parity histogram/depth-plot semantics and CRAM without explicit reference remain.

### bedcov

- **C source:** `bedcov.c`
- **HTSlib APIs used:** pileup iteration (`bam_plp_*`), BED parsing, record filters
- **htslib-rs coverage:** ⚠️ — exact pileup parity still needs `bam_plp_*`; current implementation uses alignment record iteration/CIGAR walks for SAM, BAM, and reference-backed CRAM.
- **samtools-rs status:** 🟡 — CIGAR-walk total coverage is implemented with `-Q`, `-g`/`-G` flag-mask controls, `-j` deletion/refskip skipping, `-H` output headers, `-c` read-count columns, and `-d` depth-threshold columns; exact pileup behavior remains.

### stats

- **C source:** `stats.c` (~123k LOC) + `stats_isize.c`
- **HTSlib APIs used:** record iteration, CIGAR analysis, base/qual histograms, GC bias, insert-size distribution
- **samtools-rs status:** 🟡 — `SN` summary numbers, runtime `is sorted` for record-backed paths, record-backed `-I`/`--id` read-group/sample filtering, `-f`/`--required-flag`, `-F`/`--filtering-flag`, `-i`/`--insert-size` insert-size capping, `-m`/`--most-inserts` insert-size bulk selection, record-backed `-l`/`--read-length`, `-q`/`--trim-quality` BWA trim counting, FFQ/LFQ quality histograms, GCF/GCL GC histograms, and approximate CIGAR-walk COV coverage histograms with `-c`/`--coverage` bin ranges and `-g`/`--cov-threshold` target percentage lines with target-region validation are implemented for SAM, BAM, reference-backed CRAM, SAM/indexed BAM/reference-backed CRAM positional regions, and SAM/indexed BAM/reference-backed CRAM `-t` target files, with overlapping BAM/CRAM regions de-duplicated. `-d` / `--remove-dups` filters duplicate-marked primary records from the summary and record-level histograms; exact pileup-backed COV parity, per-cycle metrics, and CRAM without explicit reference remain.

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
- **HTSlib APIs used:** `sam_open_format`, `sam_hdr_read`, `sam_read1`, `bam_aux_*`, `bam_sanitize`, `hts_crc32`, `hts_set_threads`
- **htslib-rs coverage:** ⚠️ — SAM/BAM record iteration is available; whole-CRAM record iteration and lower-level raw aux/CIGAR byte access still need coverage for full parity.
- **samtools-rs status:** 🟡 — default SAM/BAM checksum output works with read-group grouping, flag filters/masks, reverse-complement handling, selected and wildcard/exclusion scalar/string/array aux tags with canonical integer encoding, `-N`, `-o`, `-q`, `-v`, `-T`, `-O`, `-P` position columns, `-C` CIGAR columns, `-M` mate columns, `-B` bamseqchksum-compatible formatting, and `-m` merging for default/position/CIGAR/mate-column checksum reports; CRAM, all-field mode, and sanitizer mutation remain.

### samples

- **C source:** `bam_samples.c`
- **HTSlib APIs used:** `sam_hdr_*` to list `@RG` SM values
- **htslib-rs coverage:** ✅
- **samtools-rs status:** ⬜

### reference

- **C source:** `reference.c`
- **samtools-rs status:** 🟡 — SAM/BAM MD-tag reconstruction to FASTA works with `-o`, `-q`, and basic `-r` region output; CRAM input, embedded-reference extraction (`-e`), indexed region iteration, and full upstream parity remain.

### flags

- **C source:** `bam_flags.c`
- **HTSlib APIs used:** `bam_str2flag`, `bam_flag2str`
- **htslib-rs coverage:** 🟦 — implemented locally in `samtools_rs::bam_flag` (could be promoted to `htslib-rs`)
- **samtools-rs status:** ✅

## htslib-rs Extensions Required (rolled up from above)

The following items are blockers for one or more samtools-rs subcommands. They should be added to `htslib-rs/TODO.md`:

1. **`sam_hdr_add_pg`** — programmatic `@PG` chain insertion with PP linkage. Required by every subcommand that writes output without `--no-PG`.
2. **`bam_aux_update_*`** — string/int/array aux updates. Still useful for closer HTSlib parity; remaining blockers include deeper `fixmate`, `markdup`, and BAM/CRAM-output `calmd` aux-tag behavior.
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
