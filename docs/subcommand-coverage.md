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
  - `sam_hdr_add_pg` equivalent — ✅ via shared header transforms / `pg` helpers for the currently covered paths
  - `hts_expr_*` filter evaluation — ✅ via `htslib_rs::expr` (audit needed)
  - Aux-tag filter / remove (`-x`, `--keep-tag`) — implemented for SAM/BAM/CRAM input to SAM/BAM/reference-backed CRAM output via record rendering / binary roundtrips — ✅
  - Shared sanitizer mutation (`-z`/`--sanitize`) — implemented through text-record rewrites / text roundtrips for supported view output paths — ⚠️
- **samtools-rs status:** ✅ — the full upstream `test_view` group passes (445 total: 427 passed, 18 expected failures). Covered surface includes SAM/BAM/reference-backed CRAM count/text/BAM/CRAM paths, stdin paths, region/BED queries including sequential SAM/BGZF-SAM region fallback for large coordinates, `chr:pos` open-ended regions, two-column BED point semantics, `-L` intersection with positional regions, simple filters, expression filters, `-U`/`-p` split/unmap output for SAM/BAM/CRAM input to SAM/BAM/reference-backed CRAM output, `-x`/`--keep-tag` aux stripping for SAM/BAM/CRAM input to SAM/BAM/reference-backed CRAM output, `-z` sanitizer mutation, `--save-counts FILE` JSON counters, count-only no-reference CRAM `--save-counts` for simple summary-backed filters, reference-independent MAPQ/flag expressions, and read-group/library/aux-tag filters, `-m INT`, `-N`/`--qname-file FILE`, read-group/library filters, aux-tag filters, `-t`/`-T`, `-B`, `-s`, `--remove-flags`, `--exclude-flags`, upstream `--fetch-pairs` fixture behavior, `-X` legacy index synopsis, htslib-style aux float spelling for binary→SAM / CRAM stdin SAM / filtered SAM-input output via `sam_render`, and BAM/CRAM→CRAM plus CRAM→BAM file output. Direct C-vs-Rust smoke locks SAM text output, BAM `-h` text output, count-only output, MAPQ/required-flag/filtering-flag count filters, indexed BAM region text output, MAPQ, `flag.proper_pair`, CIGAR-derived (`endpos`/`rlen`/`sclen`/`hclen`), and mate/reference (`mpos`/`pnext`/`tlen`/`rnext`/`mrname`/`refid`/`mrefid`/`ncigar`) expression filters, qname allow/deny filters, read-group filters including upstream's no-`RG` pass-through unless `-n`, library filters, aux-tag value/file filters, rendered BAM/CRAM binary aux-tag rewrite output, and `view -c --save-counts` JSON sidecar output, including no-reference CRAM MAPQ/flag expression and read-group/library/aux-tag filter counts. Non-fixture follow-up: no-reference CRAM save-counts for text/binary output and reference-dependent expression counts, multi-file inputs, paired filters, and deeper CRAM performance/streaming parity.

### head

- **C source:** `sam_view.c::main_head` (lines 1772+)
- **HTSlib APIs used:** `sam_open_format`, `sam_hdr_read`, `sam_hdr_str`, `sam_read1`, `sam_format1`
- **htslib-rs coverage:** ⚠️
  - Raw header byte access — handled locally via `samtools_rs::header_text` and command-local stdin helpers; could be promoted to `htslib-rs`
- **samtools-rs status:** ✅ (SAM/BAM/CRAM file and stdin header/record output; reference-backed CRAM record extraction)

### tview

- **C source:** `bam_tview.c`, `bam_tview.h`, `bam_tview_html.c`, `bam_lpileup.c`
- **HTSlib APIs used:** `sam_open_format`, `sam_index_load`, pileup buffering, optional FASTA reference lookup, curses/html/text display callbacks
- **htslib-rs coverage:** ⚠️ — SAM/BGZF-SAM text iteration is available; full indexed pileup-backed interactive display remains broader work.
- **samtools-rs status:** 🟡 — the noninteractive text path used by the upstream large-position harness is implemented: `tview -d T -p REGION` over SAM/BGZF-SAM emits the expected 80-column padded text view and passes `test_large_positions`. Direct C-vs-Rust smoke now locks the indexed BGZF-SAM text output. `tests/tview.rs` covers top-level text-mode dispatch, width/region parsing, help, and unsupported/malformed invocation errors. Interactive curses, HTML output, sample/read-group selection, reference FASTA display, and general BAM/CRAM indexed pileup rendering remain.

### quickcheck

- **C source:** `bam_quickcheck.c`
- **HTSlib APIs used:** `hts_open`, `hts_get_format`, `sam_hdr_read`, `sam_hdr_nref`, `hts_check_EOF`, `hts_close`
- **htslib-rs coverage:** ⚠️
  - `htslib_rs::format::detect_path` — ✅
  - `htslib_rs::alignment_compat::read_*_header_from_path` — ✅
  - BGZF/CRAM EOF marker check — implemented locally; should be promoted to `htslib-rs::format` or `htslib-rs::bgzf_compat`
- **samtools-rs status:** ✅ — parity output matches `quickcheck/all.expected`; CLI coverage locks the full `quickcheck -v` fixture set, stderr diagnostics, and exit bitmask, while library integration coverage locks representative per-file status codes.

### index

- **C source:** `bam_index.c`
- **HTSlib APIs used:** `sam_index_build3`, `sam_index_build`, `sam_index_save`, `hts_set_threads`, `bgzf_check_EOF`
- **htslib-rs coverage:** ✅ via `htslib_rs::index_compat::build_bai` / `build_bam_csi` / `build_sam_csi` / `build_cram_crai`, plus explicit-index lookup via `##idx##`.
- **samtools-rs status:** ✅ — the full upstream `test_index` group passes (26/26): BAI/CSI/CRAI/BGZF-SAM CSI build, `-c` CSI mode, `--min-shift`, `-M`, `-o`, legacy `<in> <out.idx>` synopsis, `view -X`, `merge -X -R`, and `view --write-index` auto-indexing. Direct C-vs-Rust smoke now locks explicit `index -o` BAI bytes for a stable BAM fixture. Local `-@ N`, attached `-@N`, `--threads N`, `--threads=N`, and top-level `--threads` now route BAM/BGZF-SAM BAI/CSI construction through noodles multithreaded BGZF readers when nonzero. Non-fixture follow-up: CRAI thread parity and broader throughput comparison against C samtools.

### idxstats

- **C source:** `bam_index.c::bam_idxstats`
- **HTSlib APIs used:** `sam_open_format`, `sam_hdr_read`, `sam_index_load2`, `hts_idx_get_stat`, `hts_idx_get_n_no_coor`, `sam_hdr_tid2name`, `sam_hdr_tid2len`
- **htslib-rs coverage:** ⚠️
  - `read_associated_bam_index` returns `Box<dyn csi::BinningIndex>` — ✅
  - Per-reference mapped/unmapped counts from index meta — needs accessor — ❌ (must extend `htslib-rs::index_compat`)
- **samtools-rs status:** 🟡 — BAM index counts exist, with streaming slow-path counts for SAM, BAM, and CRAM. CRAM works **with or without** an explicit reference (no-reference path uses the synthesizing-reference summary, since idxstats only needs the reference-independent reference id + flags); `samtools idxstats dat/test_input_1_a.cram` is byte-exact vs `idxstats/test_input_1_a.bam.expected` (test `idxstats_cram_without_reference_succeeds`). Full harness parity (index-meta fast path) remains.

### faidx / fqidx

- **C source:** `faidx.c`
- **HTSlib APIs used:** `fai_build`, `fai_load`, `fai_fetch`, `faidx_fetch_seq64`, `fai_destroy`
- **htslib-rs coverage:** ✅ via `htslib_rs::faidx_compat::build_index`, `read_index`, `write_index`, `fetch_sequence`, `fetch_region_sequence`
- **samtools-rs status:** 🟡 — FASTA/FASTQ index build and local region extraction work, including default C-style retrieval headers, BGZF input with `.gzi`, region files, `-o`, BGZF output, `--length` line wrapping, `--write-index`, `faidx -f`, reverse-complement `-i`, mark-strand modes, and missing/truncated-region handling. Local `-@` / `--threads` and top-level `--threads` now route BGZF read/write through worker-count BGZF paths; direct C-vs-Rust smoke locks `.fai` index bytes, default retrieval headers, zero-length/truncated-region output, and `--continue` missing-region warning/output text for both `faidx` and `fqidx`. Compression-level effects and broader BGZI edge cases remain.

### dict

- **C source:** `dict.c`
- **HTSlib APIs used:** `gzopen`, `kseq_*`, `hts_md5_*` — only HTSlib's md5 wrapper is used; everything else is its own FASTA reader
- **htslib-rs coverage:** 🟦 (use `noodles-fasta` for FASTA iteration + `md-5` crate for MD5)
- **samtools-rs status:** ✅ — sequence dictionary output matches upstream `dict.out`, `dict.alias.out`, and `dict.alt.out` fixtures byte-for-byte. Direct C-vs-Rust smoke locks both stdout and `-o FILE` output.

### flagstat

- **C source:** `bam_stat.c`
- **HTSlib APIs used:** `sam_open_format`, `sam_hdr_read`, `sam_read1`, `hts_set_threads`, `hts_set_opt` (CRAM_OPT_REQUIRED_FIELDS)
- **htslib-rs coverage:** ✅ — `AlignmentRecordSummary` exposes flag/reference/mate/mapq accessors, with SAM/BAM/reference-backed CRAM summary helpers.
- **samtools-rs status:** 🟡 — SAM, BAM, and CRAM input are implemented with default text, `-O json`, and `-O tsv` output. CRAM works **with or without** an explicit reference (no-reference path uses the synthesizing-reference summary, since flagstat only needs reference-independent flags); CRAM `flagstat` output equals the BAM equivalent. Direct C-vs-Rust smoke locks BAM text, JSON, TSV, and CRAM no-reference text output. Test `flagstat_cram_without_reference_succeeds`.

## Wave B — File Operations

### sort

- **C source:** `bam_sort.c` (~138k LOC — largest file in samtools)
- **HTSlib APIs used:** `sam_open_format`, `sam_index_*`, `sam_read1`, `sam_write1`, `sam_hdr_*`, `hts_set_threads`, BGZF temp file I/O, lz4 compression for temps
- **htslib-rs coverage:** ⚠️
  - Streaming write — ✅
  - Custom per-record sort key extraction — ✅ partial for coordinate, query-name, and aux-tag keys
  - Multi-way merge — ❌
- **samtools-rs status:** 🟡 — in-memory coordinate, query-name (`-n` natural / `-N` lexicographical), aux-tag (`-t`), **minimiser (`-M`/`-K`/`-H`/`-R`/`-I`)**, and **`--template-coordinate`** sort work for BAM, SAM, and reference-backed CRAM inputs; SAM/BAM/CRAM output is supported, with CRAM output requiring top-level `--reference` and using a temporary BAM conversion. **Every upstream `test_sort` fixture is byte-exact** (pos/name/name2/name3/tag.rg/tag.rg.n/tag.as/tag.fi/minimiser-{basic,indexed,indexed-poly}/template-coordinate; tests `sort_matches_upstream_test_sort_fixtures` + `sort_minimiser_all_variants_match_upstream`). Remaining: external/temp-file merge for very large inputs (perf) and thread/memory caps.

### merge

- **C source:** `bam_sort.c::bam_merge_core`
- **HTSlib APIs used:** as above
- **samtools-rs status:** ✅ — the full upstream `test_merge` group passes
  (28/28). In-memory coordinate, natural query-name, `-t TAG`, and
  `--template-coordinate` merge work for BAM and SAM inputs, including
  seeded `@RG`/`@PG` reconciliation, `-r`, `-c`/`-p`, raw-header SAM output,
  unresolved aux-tag warnings, `-b FILE` input lists, stdout `-`,
  `--output-fmt=FORMAT`, `--no-PG`, `-R region` / `-L BED` restriction for
  indexed BAM, and CRAM output with top-level `--reference` via temporary
  BAM conversion. Direct C-vs-Rust smoke locks `--no-PG` SAM output for
  seeded three-input merge, `-r` filename-derived read groups, and
  template-coordinate merge. Remaining non-fixture follow-ups: streaming k-way merge
  and broader header reconciliation beyond `@HD`/`@SQ`/`@RG`/`@PG`/`@CO`.

### collate / bamshuf

- **C source:** `bamshuf.c`
- **HTSlib APIs used:** `sam_open_format`, `sam_read1`, `sam_write1`, BGZF temp file I/O
- **samtools-rs status:** 🟡 — in-memory name grouping works for BAM, SAM, and reference-backed CRAM inputs, including `-f` fast primary-pair mode, `-r` working-read cap, accepted `-n INT` temp-count compatibility, legacy positional output prefixes, `-o`/`-O` conflict validation, `--output-fmt=cram` / `.cram` output with top-level `--reference`, and upstream-style `@HD SO:unsorted GO:query`; on-disk hash bucketing remains.

### cat

- **C source:** `bam_cat.c`
- **HTSlib APIs used:** BGZF passthrough — fast concatenation without re-decompression
- **htslib-rs coverage:** ⚠️ — BGZF block-level concatenation is not yet a first-class API in `htslib-rs::bgzf_compat`; CRAM fixture parity currently uses a SAM-visible stream rewrite rather than true CRAM block concatenation.
- **samtools-rs status:** ✅ — the full upstream `test_cat` group passes (26/26): BAM, recompressed BAM, CRAM visible concatenation, CRAM region paths, stdout redirection, `-p 1/2` + `-p 2/2`, and `-h` replacement headers. Direct C-vs-Rust smoke locks upstream-shaped SAM-input rejection. Non-fixture follow-up: true CRAM-preserving concatenation and BGZF block-level BAM fast path.

### split

- **C source:** `bam_split.c`
- **HTSlib APIs used:** `sam_open_format`, `sam_read1`, `sam_write1`, `sam_hdr_*`, RG-tag inspection
- **samtools-rs status:** 🟡 — basic `@RG` splitting with per-output `@RG` header filtering and default `@PG` insertion, plus explicit `-d TAG` string/integer aux-tag splitting, work for BAM, SAM, and whole-file CRAM inputs with `-f`, `-u`, `-h`, extension-based output-format inference, `--output-fmt sam|bam|cram`, `--no-PG`, `--write-index` BAI generation for BAM outputs, `-M`, and `-p`; explicit `-d RG` can add missing output `@RG` headers. The full upstream `test_split` group passes (18/18), CRAM input covers the `test_checksum` split/merge path, direct C-vs-Rust smoke locks `.sam` template inference on a stable RG-only BAM and the missing-RG/no-`-u` error path, and `split_sam_input_by_rg_to_cram_outputs_with_reference` covers CRAM output with top-level `--reference`. Sorted-by-tag streaming and deeper upstream `@PG` byte-parity for complex chains remain.

### reheader

- **C source:** `bam_reheader.c`
- **HTSlib APIs used:** BAM/CRAM header rewriting, in-place mode for CRAM
- **htslib-rs coverage:** ⚠️ — in-place CRAM header replacement is not exposed; the current Rust path uses a harness-visible SAM rewrite for CRAM fixtures
- **samtools-rs status:** ✅ — the full upstream `test_reheader` group passes (7/7): BAM header replacement, CRAM v2.1/v3.0 visible reheader, in-place harness paths, and `-c <command>` external header filtering. Direct C-vs-Rust smoke locks upstream-shaped SAM-input rejection. Non-fixture follow-up: true CRAM-preserving in-place/binary rewrite and BAM BGZF block-level fast path.

### addreplacerg

- **C source:** `bam_addrprg.c`
- **HTSlib APIs used:** `sam_hdr_add_line`, `bam_aux_update_str`, record iteration
- **htslib-rs coverage:** ⚠️ — mutable `RecordBuf` paths cover current RG string replacement for SAM/BAM/reference-backed CRAM; direct `bam_aux_update_str` parity remains unavailable
- **samtools-rs status:** 🟡 — SAM/BAM/reference-backed CRAM add/replace exists with `-O sam|bam|cram` (`cram` requires `-T`/`--reference`), default `@PG` insertion, and `--no-PG`; tests cover CRAM input and CRAM output. Mate-aware updates and full orphan-first semantics remain.

### fastq / fasta / bam2fq

- **C source:** `bam_fastq.c` (~48k LOC)
- **HTSlib APIs used:** record iteration, aux tag access for barcodes/QT/RX/QX
- **htslib-rs coverage:** ⚠️ — `view_sam_as_fastq_text_from_path_with_limit` exists; full feature set (paired/single, barcode-aware) needs more
- **samtools-rs status:** ✅ — the full upstream `test_bam2fq` group passes (84/84), including threaded duplicates. SAM/BAM/CRAM FASTQ/FASTA conversion supports single-output and upstream-style name-grouped split-output paths (paired R1+R2 to `-1`/`-2`, R1-only or R2-only singletons to `-s` with fallback to `-1`/`-2`, READ_OTHER to `-0` with fallback to `-s`), per-record interleaved output when `-1`/`-2` paths alias, flag filters, read-name suffix controls, selected/all aux comments, accumulating `-t`/`-T` aux-tag selections including compact `-Tfoo`, `-d`/`-D` value-union aux-tag filtering, FASTQ `-O` original-quality `OQ` tags, `-v INT` missing-quality defaults, `-U`/`--UMI-tag` UMI read-name suffixes, `-i`/`--barcode-tag` CASAVA barcode fields, FASTA/FASTQ reverse-complement of reverse-strand records, headerless `.sam` index extraction, `--no-sc` soft-clip trimming with optional backup aux fields, per-record `--i1`/`--i2` index FASTQ extraction with `--index-format` (default `i*i*`) and `--quality-tag` (default `QT`), and CRAM input with explicit or discoverable references. Direct C-vs-Rust smoke locks default FASTQ stdout bytes, the `[M::bam2fq_mainloop]` processed/discarded informational stderr summary, and `-0 FILE` routing where paired reads remain on suffixed stdout while READ_OTHER records go to the file. Non-fixture follow-up: broader CRAM reference-discovery edge cases and worker-thread propagation beyond accepted `-@`.

### import

- **C source:** `bam_import.c`
- **HTSlib APIs used:** FASTA/FASTQ reading + SAM/BAM/CRAM writing
- **samtools-rs status:** ✅ — the full upstream `test_import` group passes (21/21): single FASTA/FASTQ, paired FASTQ, `-0` singleton input alongside paired `-1`/`-2`, positional interleaved FASTQ, index reads, CASAVA/SRA name parsing, UMI/barcode/comment aux tags, read-group header/tag support, and SAM/BAM/CRAM output. CRAM output is supported via `--cram`, `-O cram`, `--output-fmt=cram`, or a `.cram` output path by encoding the unmapped import stream against an empty temporary FASTA. Direct C-vs-Rust smoke locks paired read-group SAM output and interleaved `-T ""` aux-output SAM bytes.

## Wave C — Editing / Mate-aware

### fixmate

- **C source:** `bam_mate.c` (~43k LOC)
- **HTSlib APIs used:** record iteration, mate-flag/pos rewriting, MC/MQ tag updates, base-modification tag maintenance
- **htslib-rs coverage:** ✅ — mutable `RecordBuf` paths now cover the upstream fixture surface, including order-preserving aux rewrites and MM/ML/MN base-modification trimming/validation.
- **samtools-rs status:** ✅ — the full upstream `test_fixmate` group passes (42 total: 40 passed, 2 expected failures), including threaded duplicates. Name-grouped BAM/SAM/reference-backed CRAM mate fixup supports coordinate-sort rejection, TLEN recalculation, default MC/MQ mate tags, `-m` mate-score tags, `-c` template-CIGAR `ct` tags, default sanitizer mutation, `-r`, raw-header SAM output, CRAM output with top-level `--reference` (`-O cram` / `--output-fmt=cram`), and `-M` base-modification parity for draft tags, hard-clipped secondary/supplementary records, invalid ML/MN cases, and missing-sequence cases. Rust regression `fixmate_accepts_reference_backed_cram_input` covers CRAM input and CRAM output.

### markdup

- **C source:** `bam_markdup.c` (~89k LOC)
- **HTSlib APIs used:** name-sort grouping, position-sort grouping, barcode parsing, duplicate marking via flag updates
- **htslib-rs coverage:** ⚠️ — mutable SAM/BAM `RecordBuf` paths cover current flag and aux-tag reads; CRAM is handled through reference-backed decode plus temporary-BAM encode bridges; indexed/streaming parity remains.
- **samtools-rs status:** ✅ — single-end and paired-end duplicate marking for SAM/BAM/reference-backed CRAM input exists with optional barcode-key grouping (`-b`/`--barcode-tag`), duplicate flag/tag clearing (`-c`), `-S` compatibility for supplementary propagation, `-t` duplicate-origin `do` tags, `-d` optical-distance duplicate classification with `dt:Z:SQ`/`dt:Z:LB` tags, default QCFAIL exclusion with `--include-fails` override, validated `-m t|s`/`--mode t|s` compatibility, optical-aware estimated library size in `-s` stats, secondary/supplementary qname propagation, upstream-shaped expect-fail exits for queryname sort, bad coordinate order, missing `MC`, and missing `ms`, `-r`, upstream-shaped `-s` summary fields with fixture-verified counts for the promoted `5..18` matrix, `-O sam|bam|cram`, `-o`, local `-T`/`--reference` for CRAM input/output, default `@PG`, and `--no-PG`.

### rmdup

- **C source:** `bam_rmdup.c` + `bam_rmdupse.c`
- **samtools-rs status:** 🟡 — single-end and paired-end duplicate removal for BAM, SAM, and reference-backed CRAM inputs exists, with `-s`/`-S`, `-O sam|bam|cram`, `.cram` output inference, local `-T`/`--reference` or top-level `--reference` for CRAM input/output, default `@PG`, and `--no-PG`; Rust regression `rmdup_accepts_reference_backed_cram_input_and_output` covers CRAM input to CRAM output through the shared temporary-BAM conversion path. The dev byte-parity smoke now includes direct C-vs-Rust output-file and stderr diagnostic cases for single-end SAM output and paired SAM output. Broader deprecated-command parity for binary output, CRAM, and non-smoke diagnostic/stat output remains.

### calmd / fillmd

- **C source:** `bam_md.c`
- **HTSlib APIs used:** record iteration, BAQ via `probaln_glocal`, MD/NM recomputation
- **htslib-rs coverage:** ✅ partial via `htslib_rs::probaln` and `htslib_rs::alignment_compat::recalculate_baq_*`
- **samtools-rs status:** 🟡 — SAM, BAM, and reference-backed CRAM input can emit SAM text with recomputed MD/NM tags against FASTA references; SAM input can also run BAQ paths (`-r`/`-E`/`-A`) directly, while BAM/CRAM input runs the same BAQ helpers through a temporary SAM bridge. `-e` changes matching bases to `=`, including on unmapped records that still have reference/CIGAR fields, `-d` drops all aux tags except `RG`, `-q` bins base qualities, `-N` suppresses MD/NM aux updates and diagnostics, `-C cap` caps MAPQ through the upstream `sam_cap_mapq` port when `cap > 10`, `-n max_nm` masks matching bases and zeroes their qualities when the recomputed NM reaches the threshold, `-b`/`-u` emit BGZF BAM output, `.cram` / `-O cram` / `--output-fmt=cram` emit reference-backed CRAM output, and getopt-style glued short clusters (`-uAr`, `-C40`, `-n2`) are split, with default `@PG`/`--no-PG`. The upstream `test_calmd` invocation (`calmd -uAr mpileup.1.sam mpileup.ref.fa` → BGZF) passes (integration test `calmd_dash_u_a_r_emits_bgzf_bam_like_upstream`), direct C-vs-Rust smoke locks default `calmd --no-PG mpileup.1.sam mpileup.ref.fa` SAM stdout/stderr/exit parity plus `-e` matching-base conversion, `-r -e` BAQ plus matching-base conversion, `-d` tag dropping, `-n` masking diagnostics, `-Q` quiet suppression, `-C` MAPQ capping, `-q` quality binning, and `-N` no-MD/NM-update behavior, `calmd_writes_cram_output_with_reference` covers CRAM output, `calmd_cap_mapping_quality_uses_sam_cap_mapq` covers `-C`, `calmd_max_nm_masks_matching_bases_and_qualities` covers `-n`, and `calmd_baq_accepts_bam_and_cram_input` covers BAM/CRAM BAQ input. Remaining: broader upstream MD/BAQ option byte parity beyond the promoted harness and direct smoke cases.

### targetcut

- **C source:** `cut_target.c`
- **HTSlib APIs used:** pileup, revised MAQ error model (`errmod_cal`), optional reference-backed BAQ realignment
- **samtools-rs status:** ✅ — faithful port of the upstream pileup consensus, revised MAQ error-model scoring, and two-state target HMM. Supports `-Q`, `-i`, `-0`, `-1`, `-2`, `-f`/`--reference`, and `-o`; SAM/BAM input works via the htslib-rs pileup engine, and CRAM requires `-f`. Upstream ships no `test_targetcut` fixtures, so coverage is focused Rust unit tests for long supported interval emission, min-baseQ filtering, and attached scoring-option parsing, plus `tests/targetcut.rs` CLI coverage for dispatch, `-o`, quality filtering, and error exits. Optional exact `sam_prob_realn` BAQ side effects for `-f` remain a non-fixture polish item shared with the broader BAQ parity work.

### reset

- **C source:** `reset.c`
- **HTSlib APIs used:** record iteration, aux-tag stripping, flag/CIGAR/pos resets
- **samtools-rs status:** 🟡 — BAM, SAM, and CRAM reset paths clear alignment fields, default aux tags, and alignment-dependent flags; reverse-strand sequence/quality re-reversal, `-x`/`--keep-tag`, `--no-RG`, `--reject-PG`, `--dupflag`, default `@PG`, and `--no-PG` are supported. The output header is rebuilt faithfully per `reset.c:307-324` — a fresh `@HD VN:1.6` + verbatim `@RG`/`@PG` (no `@SQ`/`@CO`), for SAM, BAM, and CRAM output. CRAM output uses a temporary BAM plus empty temporary FASTA because all reset records are emitted unmapped; regression `reset_writes_cram_output` covers `.cram` inference and `--output-fmt=cram`. Broader CRAM reference-discovery parity remains.

### ampliconclip

- **C source:** `bam_ampliconclip.c` (~40k LOC)
- **HTSlib APIs used:** record iteration, CIGAR rewriting, BED parsing
- **samtools-rs status:** 🟡 — faithful SAM/BAM/reference-backed CRAM port with soft/hard primer clipping, `--both-ends`, `--original` `OA` tags, `--keep-tag` / default `NM`+`MD` deletion, length-gated filtering/failing/unmapping, strand/tolerance controls, rejects output, primer-count TSVs, stats output, `-O sam|bam|cram`, default `@PG`, and `--no-PG`. **Byte-exact vs the full upstream `test_ampliconclip` harness** for SAM fixtures, plus the dormant current-upstream `3_multi_ref_both_clip` `--both-ends` multi-reference edge; `ampliconclip_accepts_reference_backed_cram_input_and_output` covers CRAM input and output through the shared temporary-BAM conversion path. BGZF block fast path remains.

## Wave D — Stats / Pileup

### depth

- **C source:** `bam2depth.c`
- **HTSlib APIs used:** pileup iteration (`bam_plp_*`)
- **htslib-rs coverage:** ⚠️ — exact pileup parity still needs `bam_plp_*`; current implementation uses alignment record iteration/CIGAR walks for BAM and CRAM, including a synthetic-reference CRAM region-query path for reference-independent metrics.
- **samtools-rs status:** 🟡 — `-a`/`-aa`/`-d`/`-q`/`-o`, `-H`, `-f` input lists, flag filters, `-l` minimum read length filtering, `-r`, `-b`, and multi-input depth columns are implemented; direct C-vs-Rust smoke locks both stdout and `-o FILE` output for a stable region. CRAM region depth works with an explicit reference or the no-reference synthetic-query path for this CIGAR-only metric. Broader pileup-edge parity remains.

### coverage

- **C source:** `coverage.c` (~30k LOC)
- **HTSlib APIs used:** pileup iteration (`bam_plp_*`), CIGAR/quality access
- **htslib-rs coverage:** ⚠️ — exact pileup parity still needs `bam_plp_*`; current implementation uses alignment record iteration/CIGAR walks for BAM and CRAM, including a synthetic-reference CRAM region-query path for reference-independent CIGAR/quality metrics.
- **samtools-rs status:** 🟡 — `numreads`, `covbases`, percent coverage, mean depth, mean base quality, mean map quality, `-r` regions, `--min-depth`, `-Q`/`--min-BQ`, read map-quality filtering, `-b`/`--bam-list`, `--ff`/`--excl-flags`, `--rf`/`--incl-flags`, `-l`/`--min-read-len`, `-d` maximum-depth capping, multi-input aggregate rows, and upstream-shaped histogram/depth-plot output are implemented; missing quality scores contribute HTSlib's `0xff` base-quality sentinel, which is locked by direct C-vs-Rust smoke over `test_input_1_a.bam`; direct smoke also locks `-o FILE`, ASCII `-m` breadth histograms, `-D` mean-depth plots, an uneven-bin-tail histogram, and `$COLUMNS`-driven default-bin selection, including C-style x-axis tick placement. `-m`/`-A` and `-D` use UTF-8/ASCII glyphs, sidebars, x-axis labels, and `-w` bin control. CRAM region coverage works with an explicit reference or the no-reference synthetic-query path for current CIGAR/quality metrics. Broader interactive TTY-width histogram byte-parity edge cases remain.

### bedcov

- **C source:** `bedcov.c`
- **HTSlib APIs used:** pileup iteration (`bam_plp_*`), BED parsing, record filters
- **htslib-rs coverage:** ⚠️ — exact pileup parity still needs `bam_plp_*`; current implementation uses alignment record iteration/CIGAR walks for SAM, BAM, and CRAM, including a synthetic-reference CRAM region-query path for reference-independent CIGAR metrics.
- **samtools-rs status:** ✅ — the full upstream `test_bedcov` group passes (8/8): BAM coverage, `-j`, attached `-g512 -G2048`, `-c`, and all `-H` header cases including custom headers, empty source header fields, and BED12-derived placeholder columns. Direct C-vs-Rust smoke locks default output, `-j`, `-g512 -G2048`, and `-c -H` text output. CRAM bedcov works with an explicit reference or the no-reference synthetic-query path for the current CIGAR-only metric.

### stats

- **C source:** `stats.c` (~123k LOC) + `stats_isize.c`
- **HTSlib APIs used:** record iteration, CIGAR analysis, base/qual histograms, GC bias, insert-size distribution
- **samtools-rs status:** 🟡 — `SN` summary numbers, runtime `is sorted` for record-backed paths, record-backed `-I`/`--id` read-group/sample filtering, `-f`/`--required-flag`, `-F`/`--filtering-flag`, `-i`/`--insert-size` insert-size capping, `-m`/`--most-inserts` insert-size bulk selection, record-backed `-l`/`--read-length`, `-q`/`--trim-quality` BWA trim counting, FFQ/LFQ quality histograms, GCF/GCL GC histograms, no-reference and reference-backed GCD bins, and approximate CIGAR-walk COV coverage histograms with `-c`/`--coverage` bin ranges and `-g`/`--cov-threshold` target percentage lines with target-region validation are implemented for SAM, BAM, reference-backed CRAM, SAM/indexed BAM/reference-backed CRAM positional regions, and SAM/indexed BAM/reference-backed CRAM `-t` target files, with overlapping BAM/CRAM regions de-duplicated. READ_OTHER records match upstream for average-quality accumulation, the properly-paired percentage denominator, no-reference GCD zero-GC contribution behavior, and reference-backed MPC/GCD handling of ambiguous reference bases. `-d` / `--remove-dups` filters duplicate-marked primary records from the summary and record-level histograms; direct C-vs-Rust smoke now locks raw default BAM output plus representative `-d`, `-f`, `-F`, `-i`, `-m`, `-l`, `-c`, `-g`, `-t`, `-q`, reference-backed, `--ref-stats`, reference-backed `--ref-stats`, and positional-region paths, including the upstream-compatible stats banner and command-line header. Broader pileup COV/GCD parity, per-cycle metrics, and CRAM without explicit reference remain.

### mpileup

- **C source:** `bam_plcmd.c` (~49k LOC)
- **HTSlib APIs used:** multi-input pileup, VCF output, regions
- **htslib-rs coverage:** 🟡 — multi-input text pileup iterator exposed with overlap/orphan filtering, BAQ-adjusted mpileup qualities, and BGZF FASTA reference repositories for CRAM decode; advanced VCF/BCF and extra-column modes remain.
- **samtools-rs status:** ✅ — full upstream `test_mpileup` group passes (7/7): BAM/CRAM `-b` lists, `file://` lists, default BAQ qualities against BGZF FASTA, `-B --ff`, and overlap fixture parity. `test/mpileup:depth.reg` mpileup-simulation lines also pass. Direct C-vs-Rust smoke locks successful BAQ multi-BAM-list output, `-B --ff 0x14` stdout and `-o FILE`, overlap-removal output, stderr sample-count lines, and the missing-input error path. Remaining non-fixture work: sample grouping, VCF/BCF, and extra output columns.

### consensus

- **C source:** `bam_consensus.c` (~126k LOC) + `consensus_pileup.c`
- **samtools-rs status:** ✅ — byte-exact vs all 77 upstream
  `test/consensus/consensus.reg` cases: simple + bayesian/recall
  (Gap5) modes, fasta/fastq/pileup, `-a`/`-aa`, `-r`, `-T`/`--ref-qual`,
  `--min-MQ`/`--min-BQ`, show-del/ins, glued short options. Locked by
  `consensus_matches_upstream_consensus_reg`.

### phase

- **C source:** `phase.c`
- **HTSlib APIs used:** pileup, revised MAQ error model (`errmod_cal`), SAM/BAM/CRAM read/write, optional reference-backed CRAM decoding
- **samtools-rs status:** ✅ — faithful port of heterozygote discovery, local haplotype dynamic programming, fragment phasing, ambiguity masks, optional chimera fixing, site-list controls, and split-BAM output (`-b` prefix). Supports `-Q`/`--min-BQ`, `-q`, `-k`, `-D`, `-F`, `-A`, `-l`, `-e`, `--no-PG`, and reference-backed CRAM via `-f`/`--reference`. Upstream ships no `test_phase` fixtures, so coverage is focused Rust unit tests for phase-set marker/evidence output, min-baseQ filtering, and split-BAM creation, plus `tests/phase.rs` CLI coverage for dispatch, split-BAM output, and error exits.

### depad / pad2unpad

- **C source:** `padding.c`
- **samtools-rs status:** ✅ — SAM, BAM, and reference-backed CRAM input with `-T` padded FASTA reference convert padded reference columns to unpadded coordinates/CIGAR. SAM output (`-s`), BAM output (default/`-u`/`-1`/`.bam`), and CRAM output (`-O cram`, `--output-fmt=cram`, or `.cram`) are supported; CRAM output derives an unpadded FASTA for encoding. The full upstream `test_depad` group passes (9/9) against the `depad.001` fixture with `--no-PG`; direct C-vs-Rust smoke locks `depad -s` SAM output, including upstream's preserved `@SQ M5` field while `LN` is rewritten. Rust regression `depad_cram_input_and_output_roundtrip` covers CRAM input/output.

### ampliconstats

- **C source:** `amplicon_stats.c` (~65k LOC)
- **HTSlib APIs used:** SAM/BAM record iteration, header `@SQ` and `@RG SM:` lookup
- **samtools-rs status:** ✅ — faithful SAM/BAM/reference-backed CRAM port with byte-exact upstream `test_ampliconstats` harness coverage (`stats`, `stats_mixed`, `stats_partial`, modulo version/command-line filtering). Covers BED primer pairing, position-to-amplicon lookup, overlap-aware accumulation, depth/coverage/read/template-coordinate sections, `--tcoord-bin` template-coordinate aggregation, `COMBINED` mean/stddev output, and `-s`/`--use-sample-name` sample labels from the first `@RG SM:` header line. Direct C-vs-Rust smoke now locks all three output-file shapes with metadata lines stripped. Tests: `ampliconstats_matches_upstream_test_ampliconstats_fixtures`, `ampliconstats_use_sample_name_uses_first_read_group_sample`, `ampliconstats_accepts_bam_and_reference_backed_cram_input`, `aggregate_tcoord_merges_nearby_same_status_into_most_frequent_site`.

### cram-size

- **C source:** `cram_size.c`
- **HTSlib APIs used:** CRAM internal block/container/codec inspection
- **htslib-rs coverage:** ✅ — the vendored noodles fork exposes `Container::blocks()` (raw per-block content_id/type/method/sizes) + the public `CompressionHeader` encodings/preservation-map inventory.
- **samtools-rs status:** ✅ — default, `-v` (verbose), **and `-e` (encodings)** are all **byte-exact** vs the entire upstream `test/cram_size/cram_size.reg` (`normal.out`/`verbose.out`/`encodings.out`): faithful `cram_expand_method`/`comp_method2expanded` method decoder, block walk, `cram_cid2ds` map, aggregation, summary, and `cram_describe_encodings`/`cram_codec_describe` with htslib's exact DS + `tag_encoding_map` ordering. Test `cram_size_matches_upstream_cram_size_reg` (all 3); direct C-vs-Rust smoke now locks normal, verbose, and encodings text output.

### checksum

- **C source:** `bam_checksum.c` (~47k LOC)
- **HTSlib APIs used:** `sam_open_format`, `sam_hdr_read`, `sam_read1`, `bam_aux_*`, `bam_sanitize`, `hts_crc32`, `hts_set_threads`
- **htslib-rs coverage:** ⚠️ — SAM/BAM record iteration and whole-CRAM record iteration are available; lower-level raw aux/CIGAR byte access remains partial for broader non-fixture parity.
- **samtools-rs status:** ✅ — the full upstream `test_checksum` group passes (14/14), including threaded duplicates. Default SAM/BAM/CRAM checksum output works with read-group grouping, flag filters/masks, reverse-complement handling, selected and wildcard/exclusion scalar/string/array aux tags with canonical integer encoding, `-N`, `-o`, `-q`, `-v`, `-T`, `-O`, `-P` position columns, `-C` CIGAR columns, `-M` mate columns, `-B` bamseqchksum-compatible formatting, `-a` all-field shorthand with upstream-style sanitizer defaults, `-z`/`--sanitize` record mutation, and `-m` merging for default/position/CIGAR/mate-column checksum reports. Direct C-vs-Rust smoke locks both stdout and `-o FILE` checksum reports. CRAM input uses explicit-reference and no-reference all-record iterator paths; broader reference-compressed CRAM coverage remains a non-fixture follow-up.

### samples

- **C source:** `bam_samples.c`
- **HTSlib APIs used:** `sam_hdr_*` to list `@RG` SM values
- **htslib-rs coverage:** ✅
- **samtools-rs status:** ✅ — lists `@RG SM:` samples across inputs with header-driven de-duplication, `-T`, `-o`, `-h`, `-i`, `-f`/`-F`, stdin path lists, `-X` custom index pairs (exact file, directory, or prefix — `sam_index_load3`-style resolution), and CRAM headers. Direct C-vs-Rust smoke locks both stdout and `-o FILE` output.

### reference

- **C source:** `reference.c`
- **samtools-rs status:** ✅ — **entire upstream `test_reference` byte-exact**. SAM/BAM/CRAM MD-tag reconstruction (`-o`, `-q`, `-r`, indexed BAM iteration); CRAM MD mode discovers `@SQ UR:` references for C-generated embedded-reference CRAMs; embed_ref CRAM **read + write** in the vendored noodles fork; `view -O cram,embed_ref=1`; `-e`/`--embedded` faithful `cram2ref`. All four upstream invocations (`reference` MD no-`-T`, `-e`, and both `-r 17:1000-1500` variants) match `reference/mpileup.{MD,embed}.fa{,.reg}.expected` (tests `reference_embed_ref_full_test_reference_byte_exact`, `reference_cram_md_path_with_reference_matches_upstream`). Direct C-vs-Rust smoke locks MD full/region output, MD `-o FILE`, and `reference -e` over a C-generated embedded-reference CRAM.

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
8. **`hts_set_threads` wiring** — `faidx` / `fqidx` BGZF paths now propagate `-@ N` and global `--threads`; `index` propagates local/global thread counts to multithreaded BGZF readers for BAM/BGZF-SAM BAI/CSI construction; native BAM-output wrappers for region slicing and required-flag extraction use multithreaded BGZF writers when `threads` is nonzero. Propagate the same worker-count plumbing through the remaining alignment readers/writers and subcommands.
9. **CRAM EOF marker check in `htslib-rs::format`** — currently duplicated in `samtools_rs::commands::quickcheck`.
10. **Raw BAM header text** — `read_bam_header_text_from_path` (duplicated in `samtools_rs::header_text`).
11. **CRAM internals for `cram-size`** — either expose minimal block/codec inspection or drop the subcommand.
12. **Region grammar coverage audit** — confirm `htslib_rs::region::parse_region` handles `*` (unmapped), `.` (everything), and the `chr:from-to` / `chr:from` forms.
