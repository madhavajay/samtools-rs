# TODO-NEXT: Library-Blocked Work (htslib-rs first, then samtools-rs)

This file tracks the remaining `samtools-rs` parity work that is **blocked on
underlying-library extensions**. It is the successor to the "Items blocked on
htslib-rs / noodles extensions" and "Extensions Needed" sections of `TODO.md`,
re-scoped for the next pass where the samtools-only constraint is lifted.

## Batch progress (working branch `work-htslib-blocked-batch`)

htslib-rs pinned forward through several known-green commits (current pin
in `TODO.md` "Submodule Pinning"); every commit keeps both gates green.

- **#6** ✅ done — `build_bai` no longer needs `@HD SO:coordinate`.
- **#12** ✅ done — header-aware BAM-CSI depth fixes the `large_chrom`
  `> 2^29` reference panic; byte-exact vs `dat/large_chrom.out`.
- **#1** ✅ core + most wiring — public pileup iterator engine
  (`PileupColumn`/`PileupRead` + `PileupOptions` flag/mapq/overlap/orphan
  + smart overlap removal). **Byte-exact vs upstream fixtures:**
  `mpileup` (`mpileup.out.3`/`out.5`/stderr; `out.1` depth+bases, only
  BAQ quals differ → #11), `consensus --mode simple` (every
  `test/consensus/expected` in `consensus.reg`), `coverage`
  (`coverage/{1..5}`), `bedcov` (all four `test_bedcov`), `depth`
  (`large_pos/depth{,_bed}`), and **consensus `recall`/Bayesian modes**
  (all 77 `consensus.reg`). *Remaining:* `targetcut`, `phase`,
  `ampliconstats`.
- **#2** ✅ core + wiring — whole-CRAM all-record iterator; `stats` and
  `checksum` no-region CRAM byte-identical to BAM (bar NM-derived
  lines). `reference` CRAM MD path now wired
  (`query_cram_records_all_from_path[_with_reference]`, `fefc4ff`):
  with `-T/--reference`, `samtools reference` (whole-file + `-r`) on
  the upstream embed_ref test CRAM is **byte-exact** vs
  `reference/mpileup.MD.fa{,.reg}.expected`. ⛔ *Blocked on noodles:*
  the upstream no-reference / `-e` invocations need embedded-reference
  decoding, which **noodles-cram 0.93.0 does not implement** (it
  treats embed_ref slices as requiring an external reference and
  `expect()`-panics on an empty repository at
  `io/reader/container/slice.rs:446`) — same noodles-internals family
  as #3; raised for a decision (no noodles patch). Optional CRAM NM
  recompute for exact `stats` mismatch/error-rate also still open.
- **#8** correctness ✅ — `-@`/`--threads` accepted everywhere, output
  byte-identical regardless of count (perf worker-pool wiring deferred).
- **#9** partial — `sort`/`view --write-index` BAI == post-pass BAI.
  *Remaining:* `merge --write-index`, CSI/CRAI auto-write.
- **#10** ✅ core — `.` (everything) and `*` (unplaced) both done for
  the upstream-tested SAM/count paths. *Remaining:* binary `*` output.
- **#11** ✅ done — htslib-rs probaln/BAQ verified (realn fixtures) and
  `calmd` BAM output wired: `calmd -uAr in.sam ref.fa` emits BGZF
  (upstream `test_calmd` acceptance), getopt-style `-uAr` cluster split,
  `-A` applies recalculated BAQ to QUAL. Integration test
  `calmd_dash_u_a_r_emits_bgzf_bam_like_upstream`.
- **#5** ✅ done — `idxstats`/`flagstat` on CRAM without an explicit
  reference, via the synthesizing-reference summary (`4aea535` /
  `5f871a5`); byte-exact vs the BAM equivalents.
- **#4** 🟢 mostly done — binary `@PG` for `view` SAM→BAM, SAM→CRAM,
  and BAM→BAM (`5f9c643`/`741e368`/`0604d24`, htslib-rs `96e9e46`);
  remaining: CRAM-input + filtered/region binary-copy sub-paths.
- **#3** ⛔ blocked on noodles-cram `pub(crate)` internals (decision
  required). **#7** assessed (#7 aux
  mutation is functionally unblocked via `RecordBuf` — and an
  order-preserving `aux_del`/`aux_set_append` now exists in `fixmate`,
  the exact `bam_aux_del`+`bam_aux_append` semantics #7 wanted).

**§13 byte-parity also achieved this batch** (beyond the pileup family):
`fixmate` (entire `test_fixmate` group), `addreplacerg` (entire
`test_addrprg` group), and `sort` (`pos/name/name3/tag.rg/tag.rg.n/tag.as`) are
byte-exact modulo `@PG` (sort via raw-header preservation + the exact
`strnum_cmp` natural comparator + `-o -`→stdout; the `-`→stdout bug also
fixed in `fixmate`/`markdup`/`rmdup`). So upstream byte-exact now:
`flags`/`quickcheck`/`dict`/`idxstats`(BAM/SAM)/`coverage`/`bedcov`/
`depth`/`mpileup`/`consensus`/`fixmate`/`addreplacerg`/`sort`(basic);
`split` matches `test_split` under the harness' `reorder_header`.

**Precisely-characterized remaining blockers** (each a large port, not a
quick fix — verified by probing):
- `merge` → ✅ **DONE (fully)**: @RG/@PG `-s SEED` reconciliation
  (crate::rand48 + raw-header), plus `-r` (filename-stem @RG attached to
  every record, order-preserving RG:Z delete-then-append) and `-c`/`-p`
  (combine identical @RG/@PG IDs; grouped short opts `-cp`/`-rp`).
  **Byte-exact vs ALL upstream test_merge fixtures merge/{2,4,5,6,7}**
  (modulo @PG); integration test `merge_reconciles_rg_pg_byte_exact_vs_upstream`
  covers all five. Remaining (not fixture-covered): k-way streaming
  merge, CRAM output, `--template-coordinate`.
- `markdup` → ✅ **DONE — all 14 `test_markdup` SAM fixtures
  byte-exact** (`markdup/{5..18}`): faithful `make_pair_key`
  (template + sequence), `make_single_key`, `calc_score`+`ms`,
  QCFAIL/qname tie-break, `-S` `dup_hash` propagation,
  `get_coordinates_colons` + regex `get_coordinates` optical-name
  parse, the full `find_duplicate_chains` optical-chain re-tagging,
  `--use-read-groups`, `--duplicate-count`, `--read-coords`/
  `--coords-order`/`--barcode-rgx`/`--barcode-name` (regex crate;
  capture-span-bounded coord parse), raw-header SAM output. Remaining:
  exact `-s` stats counts, CRAM, the `1..4` expect-fail cases.
- `sort -M` minimiser → ✅ **DONE — all 3 upstream fixtures
  byte-exact** (`e9e1717`/`1a2c9d4` + the indexed slice): faithful
  `worker_minhash` + `bam1_cmp_by_minhash` + `build_minhash_index` +
  `minhash_with_idx[_squash]` port (fwd/rev strand, `-H` squash, `-R`,
  `-K` kmer, `-I FILE` indexed reference), plus the `reset.c:307-324`
  fresh-header rebuild (`@HD VN:1.6` + `@RG`, no `@SQ`, for SAM **and**
  BAM output). The upstream `reset --dupflag | sort -m 10M …` pipeline
  is **byte-identical** to `sort/minimiser-{basic,indexed,indexed-poly}
  .sam` under the harness' `ignore_pg_header` — full `@HD`/`@RG`
  header **and** all 569 records — locked by
  `sort_minimiser_all_variants_match_upstream`.
- `sort` → ✅ **all upstream `test_sort` fixtures byte-exact**: `-N`
  lexicographical name sort (`name2`), `-t FI` (`tag.fi`), and
  `--template-coordinate` (full `template_coordinate_key` +
  `bam1_cmp_template_coordinate` + `unclipped_*` + `lookup_libraries`
  port, `@HD GO:query`) all land byte-exact. Only external/temp-file
  merge (large-input perf, not fixture-blocking) remains.
- `stats` → ✅ **DONE — all 20 fixture groups byte-exact** (stat/1–19
  plus the four stat/12 `-t`/`-p` variants), incl. the `-p`/
  `--remove-overlaps` paired-overlap chunk subtraction and the f32
  error-rate cast. Out of scope (no fixtures / library-blocked):
  multi-bin (>20 kbp) GC-depth, exact pileup-backed COV, CRAM without
  explicit reference (the htslib-rs CRAM all-record iterator). Original
  detail retained below.
  Byte-exact end to end vs
  upstream `stat/{1..11,13,14,15,16,17,18,19}` (CHK checksum, all SN
  lines + comments, FFQ/LFQ with `max_qual` width + reverse-strand
  quality orientation, full MPC reference-mismatch engine, GCF/GCL,
  GCC/GCT/FBC/FTC/LBC/LTC, region-clipped `bases mapped (cigar)`,
  barcode BCC/QTQ, IS, RL/FRL/LRL, MAPQ, ID/IC indel engine +
  `nindels`=300 cap, COV, single-bin GCD; supplementary handling;
  `-S RG`/`-P` per-tag `.bamstat`; `-I` read-group filter; streaming
  `-t` BAM path; upstream unpaired-read `order`; `--ref-stats` RFS
  with/without reference incl. region-merged targets; insert-size
  avg/sd from the integer-halved per-size arrays so a lone in-region
  read whose mate was filtered drops out; stat/12 `3reads.overlap` +
  `2reads.overlap` byte-exact). **Remaining — only the two stat/12
  `-p`/`--remove-overlaps` `nooverlap` variants:** the paired-overlap
  chunk subtraction. Port `remove_overlaps` (stats.c:1057-1170) and
  `cleanup_overlaps` (stats.c:1023): a `qname -> pair_t` map keyed by
  template, storing the first mate's mapped-ref `[pmin,pmax]` chunks
  by `order`; when the other mate arrives, subtract the overlap from
  `nbases_mapped_cigar` and the coverage round-buffer
  (`nbases_mapped_cigar -= (pmax-pmin)` etc. at stats.c:1140-1170),
  gated by the `-p` flag (currently a no-op in our arg parser).
  Touches `bases mapped (cigar)`, `error rate`, and `COV` (both are
  byte-exact across the other 19 fixture groups, so the change must be
  `-p`-only). Plus (non-fixture-blocking): multi-bin (>20 kbp)
  GC-depth, exact pileup-backed COV, CRAM without explicit reference
  (blocked on the htslib-rs CRAM all-record iterator).
- `ampliconclip` → ✅ **DONE**: full port, byte-exact vs the entire
  upstream `test_ampliconclip` harness (10 SAM + 3 primer-count TSVs).
- `ampliconstats` (`amplicon_stats.c`, 1776 LOC) → ✅ **DONE**:
  full faithful port — byte-exact vs the entire upstream
  `test_ampliconstats` harness (`stats`, `stats_mixed`,
  `stats_partial`, modulo the harness-stripped version/command-line
  lines). Remaining: `--tcoord-bin` aggregation, CRAM,
  `--use-sample-name`.
- `phase`/`targetcut` (no upstream fixtures + dense numerical HMM),
  `reference` CRAM (needs MD recompute / embedded-ref internals).
- **consensus Bayesian/default mode** — fixture-backed port, in
  progress. Default = `MODE_RECALL`, single-pass
  `calculate_consensus_gap5(cp_r)` (mixed/precise are opt-in).
  Sub-steps, in order, each a verifiable commit:
  1. ✅ **DONE** (`c95c8cd`): `consensus_init` probability tables —
     `consensus.rs::bayes::{ConsProbs, cons_probs_init, Qcal::flat,
     default_recall}`, unit-tested. Defaults: P_HET=1e-3,
     P_INDEL=2e-4, het_scale=1.0, poly_mul=0.01, flat qcal,
     MODE_RECALL.
  2. ✅ **DONE** (`024d489`): math accel helpers —
     `consensus.rs::bayes::{q2p_table, mqual_pow_1m_table, ETab
     (fast_exp), fast_log2, ph_log, TENLOG2OVERLOG10}`, unit-tested
     vs the `bam_consensus_tab.h` literals and libm.
  3. ✅ **DONE** (`2a6d97d`): the S[15] accumulation + call —
     `consensus.rs::bayes::{Gap5Obs, Gap5Opts, Gap5Cons, gap5_call,
     L, MAP_SING, MAP_HET}`, the faithful default-build port (no K2/
     DO_FRACT/DO_HDW/DO_POLY_DIST/DISCREP) incl. the `nm_adjust` /
     (nm_local+1) + `td` depth fudge, `qual2` poly_mul, the 6-case
     switch, +lprior15, pure/het argmax, shift-normalise, phred,
     het_logodd. Unit-tested (pure A/G, empty→N, A/C het).
  4. ✅ **resolved**: `calculate_consensus_gap5m` (:1797) — all
     fixture cases are `MODE_RECALL` (`bayesian`/`bayesian_r`/
     default; verified at bam_consensus.c:3124-3142). For non-MIXED
     it is exactly `gap5_call` with `cp_recall`, so step 3 already
     covers every fixture; the experimental MIXED combination
     (`-m bayesian_m`, no fixtures) is deferred.
  5. **Library-blocked (htslib-rs first).** The numerical engine
     (steps 1–4) is done and unit-tested, but wiring it byte-exact is
     blocked: `htslib_rs::alignment_compat::PileupRead` only carries
     the base at the column (`base`, `qpos`, `mapping_quality`, …) —
     not the aligned read's sequence/CIGAR/`NM`. Upstream defaults
     enable `use_mqual` **and** `nm_adjust`, and `poly_len` is applied
     unconditionally, so byte-exact output needs, per pileup read,
     **(a)** `poly_len()` (bam_consensus.c:989 — homopolymer run
     length around `qpos` in the read sequence) and **(b)**
     `nm_local()` (bam_consensus.c:976 — local mismatch count within
     `nm_halo` of the position, from the read + its `NM`/MD). Neither
     was derivable from the old `PileupRead`. **htslib-rs extension
     ✅ DONE** (htslib-rs `0d81dec`/`abd52e9`/`5385da8`, pinned via
     samtools-rs `f588f17`/`b54b021`/`d484beb`): (a) `TestPileupRecord`
     now captures per-read CIGAR + `MD`/`NM`; (b) `compute_local_nm`
     faithfully ports `nm_init`'s default path (adj_qual deficit,
     homopolymer high-8-bits, soft-clip cost, MD-walk mismatch halo —
     `homopoly_fix` opt-in path deferred) + `local_nm_poly`/
     `local_nm_score`; (c) `PileupRead` now exposes `bayes_poly` and
     `bayes_nm_local` (precomputed once per record, indexed at
     `qpos+1`). All unit-tested; htslib-rs 133/0, workspace green.
  5b. **Next — the samtools-rs wiring (now unblocked):** in
     `consensus.rs`, dispatch `--mode bayesian`/default (currently
     rejected at consensus.rs:150) to a `consensus_bayes(reads,cfg)`
     paralleling `consensus_simple`: per column build `Vec<bayes::
     Gap5Obs>` from `&[PileupRead]` — filter `!is_refskip` and
     `qual >= min_qual`; `base4` = ASCII→SAM-4bit (A1 C2 G4 T8 N15,
     deletion `*`→16); `qual` = `quality.or(qpos_quality)` with the
     `255 || (0 && raw==255)`→`default_qual` rule; `mqual` =
     `mapping_quality`; `nm_local` = `bayes_nm_local`; `poly` =
     `bayes_poly`. Call `bayes::gap5_call` (td = column depth), then
     the `consensus_base` thresholds (bam_consensus.c:2137-2167):
     `cons.depth < min_depth && call!=4`→`N`/0; else if
     `het_logodd>0 && ambig` → the 25-char "AMRWa MCSYc RSGKg WYKTt
     acgt*"[het_call], cq=het_logodd; else `"ACGT*"[call]`, cq=phred;
     then `cq<cons_cutoff && cb!='*' && het_call%5!=4 && het_call/5
     !=4` → `N`/0. Reuse the simple path's ref-column + insertion
     sub-column loop and FASTA/FASTQ/pileup/line-len writer. Defaults:
     min_depth=1, call_fract=0.75, het_fract=0.5, cons_cutoff=10,
     default_qual=10, use_mqual=1, nm_adjust=1, scale_mqual=1,
     low/high_mqual=1/60, line_len=70, show_ins=1, show_del=0
     (bam_consensus.c:2985+).
     Fixtures: `samtools/test/consensus/consensus.reg`. ✅ **DONE —
     all 77 cases byte-exact** (`4dfef81`; engine wired end-to-end
     from `4ff63fa`/`a0d4d53`). Delivered, by feature:
     • **`-a`/`-aa` all-bases** — faithful `basic_fasta`/`empty_pileup2`
       gap rule: per-contig `last_pos` init, advance even on skipped
       `*`, fill `[last_pos+1,pos)` only when
       `pos>last_pos && (last_pos>0 || all_bases)`; `-a` pads the
       covered span, `-aa` the whole @SQ contig; `-aa` pileup re-emits
       every header @SQ in header order (uncovered → full-length
       empty-row block).
     • **`-r region`** — `parse_region_spec` + column-loop region clip.
     • **`--min-MQ <n>`** — `min_mqual` into `PileupOptions`.
     • **`--min-BQ`/pileup-format nuances** — `empty_pileup2` row format
       + non-refskip depth.
     • **30/31/32 bayesian `--show-del yes --show-ins no -C0`** — glued
       short-opt parse + show-del in the bayesian path.
     • **`--ref-qual`/`-T` reference** — `load_ref_seqs` + `gap_base`
       (reference base + `--ref-qual` for uncovered positions).
     Locked by integration test
     `consensus_matches_upstream_consensus_reg` (in-process INIT+P
     harness, 77/77); clippy + full workspace test suite green.
- TODO-NEXT #3/#4/#5 (CRAM internals / binary-`@PG` / CRAM index meta).

**Honest remaining scope:** the library-blocked *foundations* are all
shipped, tested and pinned; the rest is ordinary (if large) samtools-rs
porting + Phases 4–5 — a multi-week/multi-engineer effort, not
single-session work.

## Ground rules for this pass

- **All library changes go in `htslib-rs`. Do not patch `noodles`.**
  `noodles` is a third-party library we do not own; we only consume it (and
  carry minimal patches) as the backing implementation for parts of
  `htslib-rs`. Every gap below — *including the CSI large-reference issue that
  surfaces in `noodles-csi`* — must be addressed inside `htslib-rs` (a guard,
  clamp, fallback, or alternate code path in `htslib-rs`'s region/index layer),
  **not** by editing the `noodles` submodule. If a fix genuinely cannot be done
  without a `noodles` change, stop and raise it for a decision rather than
  modifying `noodles`.

- **Test-first, in `htslib-rs`.** For each item: add/extend the API in
  `htslib-rs`, and write the htslib-rs-level tests there (unit/integration
  under `htslib-rs/crates/htslib-rs`) proving the new surface works against
  fixtures. Keep the `htslib-rs` gate green
  (`cargo fmt --all --check`, `cargo clippy --workspace --all-targets -- -D
  warnings`, `cargo test --workspace`).

- **Then wire `samtools-rs` and test here.** Once the `htslib-rs` API exists
  and is tested, make the final consuming changes in `samtools-rs`, and add the
  `samtools-rs` integration tests (per-subcommand `tests/`) plus the relevant
  upstream `test.pl` fixtures. Keep the `samtools-rs` gate green.

- **Workflow unchanged:** one long-lived working branch per batch, multiple
  commits, one PR; gate green after every commit; update `TODO.md` /
  `docs/test-status.md` / `docs/subcommand-coverage.md` as work lands. If an
  `htslib-rs` change is needed, land it (with its tests) and pin the new
  `htslib-rs` commit before the dependent `samtools-rs` commit.

- **Pinning:** bump the `htslib-rs/` submodule pin in `TODO.md`
  ("Submodule Pinning") to the known-green commit once each htslib-rs change
  merges, so the dependent samtools-rs work builds reproducibly.

---

## Items (each: htslib-rs change + htslib-rs tests → samtools-rs wiring + samtools-rs tests)

> **Status (in progress):** htslib-rs side **DONE + pinned** (`9cf30b3`):
> public `pileup_from_alignment_paths[_with_reference][_and_options]` +
> `PileupColumn`/`PileupRead` (`bam_pileup1_t`-shaped, incl. `qpos_quality`,
> indel, head/tail, refskip) + `PileupOptions` (flag/mapq/`detect_overlaps`/
> `discard_orphans`) + ported `tweak_overlap_quality`/`overlap_push` smart
> overlap removal + `MPLP_NO_ORPHAN` filter; unit tests + CRAM==BAM equality.
> samtools-rs side: **`mpileup`** default text output byte-exact vs
> upstream `mpileup.out.3`/`out.5`/stderr (`out.1` depth+bases match;
> quals differ only under BAQ → #11). **`consensus --mode simple`**
> implemented, byte-exact vs every `test/consensus/expected` fixture in
> `consensus.reg` (FASTA/FASTQ/pileup, show-del/ins, call/het-fract).
> **`coverage`/`bedcov`/`depth` byte-exact** vs their full upstream
> tabular suites: coverage (`%g`/`%.3g`, min_depth-gated means, row
> ordering — `coverage/{1..5}`), bedcov (all four `test_bedcov` incl.
> attached `-g512`), depth (`large_pos/depth{,_bed}.expected.out` — sparse
> storage fixes the `LN:10001009800` OOM; bedidx now whitespace-tolerant).
> **consensus `recall`/Bayesian modes DONE** (all 77 `consensus.reg`
> byte-exact). Remaining pileup-dependent ports: `targetcut`, `phase`,
> `ampliconstats`.

### 1. Pileup iterator — highest leverage

- **htslib-rs:** expose a multi-input pileup iterator API surface
  (`bam_plp_*` / `sam_pileup` shaped). `htslib-rs::alignment_compat` has
  fixture-level pileup but no public iterator; audit and expose it.
- **htslib-rs tests:** per-column pileup over known BAM/CRAM fixtures
  (depth, base/qual, indel start/len, ref-skip) including multi-input merge.
- **samtools-rs wiring:** `mpileup`, `consensus`, `targetcut`, `phase`,
  `ampliconstats`, and the byte-exact pileup-based paths of `bedcov` /
  `coverage` / `depth`.
- **samtools-rs tests:** per-subcommand integration tests + the upstream
  `test_mpileup` / `test_consensus` / `test_ampliconstats` / etc. fixtures.
- **Unblocks:** ~5 whole subcommands plus exact `bedcov`/`coverage`/`depth`.
  Do this first.

### 2. CRAM all-record (non-region) streaming iterator — 🟡 mostly DONE

> **htslib-rs DONE + pinned** (`ca812dd`):
> `query_cram_records_all_from_path_with_reference` /
> `iter_cram_records_all_from_path_with_reference` (full RecordBuf;
> count/name/flags/start/seq/qual == range.bam equivalent).
> **samtools-rs `stats` DONE:** no-region CRAM now uses the full-record
> iterator (`collect_cram_full_stats`) + wires `stats -r/--ref-seq`;
> SN lines match BAM except NM-derived `mismatches`/`error rate` (CRAM
> stores no NM and noodles does not synthesize it — recompute-from-ref is
> a separate follow-up). Test `stats_cram_without_region_matches_bam_*`.
> **samtools-rs `checksum` DONE:** whole-CRAM via the iterator,
> byte-identical to the BAM checksum (test
> `checksum_cram_matches_bam_via_all_record_iterator`). **Remaining:**
> `reference` CRAM-input MD path on the same iterator; optional CRAM NM
> recompute for exact `stats` mismatch/error-rate parity.

- **htslib-rs:** add a non-region streaming record iterator for CRAM
  (today only `iter_cram_records_from_path_with_reference` for indexed
  *regions* exists; the non-region `summarize_*` path discards per-record
  sequence/quality/NM).
- **htslib-rs tests:** iterate a whole CRAM fixture and assert record
  count + per-record sequence/quality/NM/flags match the BAM equivalent.
- **samtools-rs wiring:** `stats` seq-length/quality/NM SN lines for CRAM
  without a region; `stats -d` no-region CRAM dedup; `checksum` whole-CRAM
  inputs; `reference` CRAM-input MD path.
- **samtools-rs tests:** `stats`/`checksum` CRAM-without-region integration
  tests + upstream `test_stats` / `test_checksum` CRAM fixtures.
- **Unblocks:** several CRAM-without-region gaps at once. Do this second.

### 3. CRAM container / block / codec inventory API — ⛔ BLOCKED (decision required)

> **Assessed (precise blocker).** `cram-size` needs, per container:
> the block table (Content-ID, uncompressed/compressed size,
> compression method, backing data-series) and the "Container
> encodings" map (DataSeries → codec: `EXTERNAL(id)`,
> `HUFFMAN(codes,lengths)`, `BYTE_ARRAY_STOP(stop,id)`,
> `BYTE_ARRAY_LEN(len_codec,val_codec)`), plus tag encodings. In the
> vendored **noodles-cram 0.93.0** these are all `pub(crate)` and not
> re-exported, so they are unreachable from `htslib-rs` (which only
> *consumes* noodles) **and** from `samtools-rs`:
> - `CompressionHeader::{preservation_map,data_series_encodings,
>   tag_encodings}` — `pub(crate)`.
> - `DataSeriesEncodings`, `TagEncodings`, `PreservationMap` — type is
>   `pub(crate)`.
> - `io::reader::container::block::Block` (has public
>   `compression_method`/`content_type`/`content_id`/
>   `uncompressed_size`) — **not exported**; `Slice::decode_blocks()`
>   only returns `(ContentId, Cow<[u8]>)` with no size/method metadata.
>
> Per the **Ground rules** ("Do not patch `noodles` … if a fix
> genuinely cannot be done without a `noodles` change, stop and raise
> it for a decision"), `cram-size` is **blocked on a scope decision**:
> either (a) accept a minimal upstreamed/forked noodles-cram public
> inventory surface, or (b) **drop `cram-size` from samtools-rs scope**
> (the documented fallback in this item). Not half-implemented; no
> noodles patch made. The dependent `reference -e` embedded-reference
> mode shares this blocker.

- **htslib-rs:** (blocked) would expose a CRAM container/block/codec
  inventory — requires the noodles-cram accessors above to be public.
- **samtools-rs wiring:** `cram-size` (entirely); `reference -e`.
- **samtools-rs tests:** upstream `test/cram_size/cram_size.reg`
  (`normal.out`/`verbose.out`/`encodings.out`).

### 4. `@PG` through the binary header — 🟢 mostly DONE

> **Done** (`5f9c643`/`741e368`/`0604d24` + htslib-rs `96e9e46`):
> `view` injects the samtools `@PG` into the **binary** header for
> **SAM-input → BAM**, **SAM-input → CRAM**, and **BAM-input → BAM**
> (the common, no-filter/no-region path) — `--no-PG` (or no captured
> argv) stays a no-op / fast path. Done at the samtools-rs level:
> `sam_bytes_with_pg` injects into the SAM-text header before
> SAM→binary conversion across every SAM-input sub-path; BAM→BAM uses
> the new htslib-rs `write_bam_from_path_transforming_header` (header
> serialized → SAM text → `apply_pg_to_header` → parsed back; records
> streamed unchanged, `PP`-chained correctly). Tests
> `view_b_embeds_pg_in_binary_bam_header` (SAM→BAM + BAM→BAM),
> `view_c_embeds_pg_in_binary_cram_header` (SAM→CRAM), and htslib-rs
> `write_bam_from_path_transforming_header_rewrites_header_keeps_records`.
> Note: upstream `test_view` deliberately uses `--no-PG` for binary
> comparisons, so this is correctness/completeness, locked by the
> samtools-rs integration tests rather than an upstream fixture.
>
> **Remaining (smaller):** CRAM-*input* → binary output, and the
> filtered/region BAM/CRAM binary-copy sub-paths
> (`write_bam_matching_filter_from_path`,
> `write_cram_*_from_*_path_with_reference`) — each would need an
> analogous header-transform variant; lower-frequency, no fixture.

### 5. CRAM `idxstats`/`flagstat` without an explicit reference — ✅ DONE

> Solved without a CRAM-index-meta accessor: `idxstats`/`flagstat`
> only use per-record **flags + reference id** (CRAM core/external
> data, reference-independent); only the *sequence* is reconstructed
> against the reference. htslib-rs `4aea535` adds
> `summarize_cram_records_from_path_synthesizing_reference`, which
> builds an all-`N` `fasta::Repository` sized from the CRAM header
> `@SQ` lines so the noodles decoder runs without erroring while the
> reference-independent fields stay byte-identical (the plain no-ref
> path errors — noodles eagerly resolves). samtools-rs `5f871a5`
> wires `idxstats`/`flagstat` to it when no `--reference` is given.
> `samtools idxstats dat/test_input_1_a.cram` is byte-exact vs
> `idxstats/test_input_1_a.bam.expected`; `flagstat` CRAM == BAM.
> Tests: htslib-rs `cram_summaries_without_reference_match_bam_flags_and_tids`,
> samtools-rs `idxstats_cram_without_reference_succeeds` /
> `flagstat_cram_without_reference_succeeds`.

- ~~**htslib-rs:** accessor to read per-reference counts from the CRAM
  index meta~~ — superseded by the synthesizing-reference summary
  (simpler, exact, no `.crai`-meta parsing needed).

### 6. Index BAMs lacking `@HD SO:coordinate` — ✅ DONE

> Completed: htslib-rs `530b27c` (`build_bai` no longer requires the SO
> header tag; `range.bam.bai` byte regression preserved; new no-SO test),
> samtools-rs `index` works on `test_input_1_{a,b}.bam` with integration
> test `index_bam_without_so_coordinate_header`. No noodles patch.

- **htslib-rs:** allow BAI/CSI creation for coordinate-ordered data whose
  header has no `SO:coordinate` tag (currently rejected with
  `invalid sort order: expected coordinate, got None`); upstream
  `samtools index` indexes such fixtures anyway.
- **htslib-rs tests:** build an index for `test_input_1_a.bam` /
  `test_input_1_b.bam`-shaped fixtures and query it.
- **samtools-rs wiring:** `samtools index` on those fixtures.
- **samtools-rs tests:** `index` integration test + the `test_index`
  fixtures that currently can't be indexed.

### 7. `bam_aux_update_*` (string / int / array, with resize) — 🟡 unblocked

> Assessed: htslib-rs already exposes `sam_aux_get` / `sam_aux_insert` /
> `sam_aux_remove` on `RecordBuf`, and `addreplacerg` already rewrites
> aux via mutable `RecordBuf::data_mut()` (its upstream fixture group
> passes). So aux mutation is **functionally unblocked**; the remaining
> work is per-subcommand wiring (`calmd` BAM MD/NM recompute,
> `ampliconclip`) and, optionally, true in-place binary-resize primitives
> for performance — not a blocking library gap.

- **htslib-rs:** binary aux update primitives with re-sizing semantics
  (the proper path; today partially worked around via mutable `RecordBuf`).
- **htslib-rs tests:** update scalar/string/array aux on a `RecordBuf`,
  re-serialize, and assert byte layout matches expectation.
- **samtools-rs wiring:** byte-exact `addreplacerg` (`bam_aux_update_str`),
  broader BAM aux rewrite in `calmd` BAM MD/NM recompute, `ampliconclip`.
- **samtools-rs tests:** `addreplacerg`/`calmd` BAM aux integration tests +
  the relevant `test.pl` groups.

### 8. `hts_set_threads` wiring to BGZF worker pools — 🟡 correctness DONE

> Correctness deliverable met: `-@`/`--threads` (incl. attached `-@4`,
> `-m768M`, `--threads=4`, …) accepted by `view`/`sort` (`index` already
> did); `-@ N` output is byte-identical to `-@ 1` — test
> `threads_flag_is_byte_identical_for_view_and_sort` (`--no-PG` isolates
> the @PG CL which legitimately embeds the arg). **Remaining (perf only,
> explicitly out of correctness scope):** actually wire the thread count
> into noodles BGZF worker pools so `-@` speeds up I/O.

- **htslib-rs:** wire a thread-count option into the BGZF/noodles worker
  pools so `-@` is honored (currently an API-compatible no-op everywhere).
- **htslib-rs tests:** functional test that multi-threaded and
  single-threaded writes produce identical bytes.
- **samtools-rs wiring:** propagate `-@` / `--threads` into `index`,
  `sort`, `view`, native API, etc.
- **samtools-rs tests:** assert `-@ N` output is byte-identical to `-@ 1`
  (correctness, not perf).

### 9. `auto_index` / index-save-during-write — 🟡 PARTIAL

> `sort --write-index` (BAI) and now `view --write-index` (BAM file
> output) both work via a post-write `index_compat::build_bai` pass that
> is **byte-identical** to a separate `samtools index` run (test
> `view_write_index_matches_post_pass_index`). The TODO deliverable
> ("writer-produced index == post-pass byte-for-byte") is met for these.
> **Remaining:** `merge --write-index`, CSI/CRAI auto-write, and true
> inline (during-write) index emission instead of a post-write pass
> (perf/streaming — not a correctness gap).

- **htslib-rs:** write BAI/CSI/CRAI alongside the writer (some samtools-rs
  paths currently do a separate post-write index pass).
- **htslib-rs tests:** writer-produced index matches a post-pass-built
  index byte-for-byte.
- **samtools-rs wiring:** `--write-index` for all writer paths.
- **samtools-rs tests:** `--write-index` integration tests across
  `view`/`sort`/`merge`.

### 10. Region-string grammar coverage (`htslib-rs::region`) — ✅ core DONE

> `.` (everything) and `*` (unplaced/no-coordinate, RNAME `*`) both
> implemented in `view`: normalized out of the region list, `*` adds an
> `only_unplaced` filter via `line_passes`/`has_filters`. Tests
> `view_dot_region_means_whole_file`,
> `view_star_region_selects_unplaced_reads` (matches upstream
> `reg_unmapped1`). **Minor remaining:** binary BAM/CRAM-output of `*`
> isn't filtered (shares view's other binary-output filter limitation);
> upstream `reg_unmapped1` compares SAM, which passes.

- **htslib-rs:** confirm/extend coverage of HTSlib's full region grammar,
  notably `*` (unmapped) and `.` (everything else).
- **htslib-rs tests:** parse + query tests for `*`, `.`, and edge spans.
- **samtools-rs wiring:** region-query edge cases across `view`/`stats`/etc.
- **samtools-rs tests:** region-grammar integration tests + `test.pl`
  region cases.

### 11. `probaln_glocal` / BAQ wiring verification — ✅ DONE

> htslib-rs side **verified**: `htslib_rs::probaln::probaln_glocal` is
> implemented with unit tests, and the BAQ surface
> (`recalculate_baq_from_sam_path`, `apply_existing_baq_from_sam_path`,
> `revert_existing_baq_from_sam_path`, `force_recalculate_baq_from_sam_path`,
> extended BAQ) is wired and passes the upstream realn fixtures
> (`ports_test_realn_*` in `alignment_io.rs`). **samtools-rs DONE**
> (`5298222`): `calmd -uAr mpileup.1.sam mpileup.ref.fa` emits a BGZF
> (BAM) stream — the exact upstream `test_calmd` acceptance check —
> via `expand_short_clusters` (getopt-style `-uAr` split), `-A`
> (apply recalculated BAQ to QUAL), and `-b`/`-u` BAM output through
> `write_bam_from_sam_reader[_with_compression_level]`. Integration
> test `calmd_dash_u_a_r_emits_bgzf_bam_like_upstream` (BGZF magic +
> 569/569 record round-trip). Remaining (not fixture-blocking):
> `-C cap`, `-n max_nm`, CRAM output, BAQ over BAM/CRAM input,
> and feeding `mpileup` out.1 BAQ-adjusted qualities.

- **htslib-rs:** verify `htslib-rs::probaln` (`probaln_glocal`) is wired and
  correct for BAQ recalculation (likely verification, not new API).
- **htslib-rs tests:** BAQ output on a known fixture matches expected.
- **samtools-rs wiring:** `calmd` BAQ paths (and later `mpileup`).
- **samtools-rs tests:** `calmd -b`/BAQ integration test + `test_calmd`
  BAQ fixtures.

### 12. CSI robustness for very large references/regions — ✅ DONE

> Completed: htslib-rs `8372873`. `build_bam_csi_with_min_shift` was
> using a fixed CSI depth of 5; it now auto-sizes depth from the largest
> reference (`alignment_csi_depth_for_header`, as the SAM-CSI builder
> already did). `view large_chrom.bam ref2` and `ref2:1-541556283` are
> byte-exact vs `dat/large_chrom.out`, no panic / no `invalid end bound`.
> **Fixed entirely in htslib-rs; noodles unpatched.** Tests:
> `builds_csi_for_very_large_reference_and_queries_it`,
> `view_large_chrom_csi_region_matches_upstream`.

- **htslib-rs:** make large-reference CSI queries robust **inside
  `htslib-rs`'s region/index handling layer** — e.g. validate/clamp the
  region end bound and guard the bin computation before it reaches
  `noodles-csi`. **Do not patch `noodles`.** Symptoms today:
  `samtools view large_chrom.bam ref2` panics `index out of bounds` in
  `noodles-csi/.../reference_sequence.rs`, and `ref2:1-541556283` reports
  `invalid end bound`.
- **htslib-rs tests:** query `large_chrom.bam ref2` and a near-INT_MAX
  span without panicking; assert correct (possibly empty) results.
- **samtools-rs wiring:** none beyond consuming the fixed API.
- **samtools-rs tests:** the upstream `test_index` `large_chrom.bam ref2`
  case passes.
- **Note:** if it provably cannot be fixed in `htslib-rs` without a
  `noodles` change, stop and raise it for a decision.

---

## Suggested sequencing

1. **#1 Pileup iterator** — by far the biggest unlock (≈5 subcommands).
2. **#2 CRAM all-record iterator** — closes several CRAM-without-region gaps.
3. (#4 🟢 mostly done — binary `@PG` for SAM→BAM/CRAM + BAM→BAM;
   #5 ✅ done, #6 ✅ done.)
4. **#7–#12** — incremental hardening / correctness.

## Out of scope here (NOT library-blocked — large samtools-rs ports)

For context only — these need no library change, just substantial
samtools-rs work, and are tracked in `TODO.md` proper: full `depad`
(~623 LOC C), full `checksum` (~1324 LOC C), full `reference` (~598 LOC C),
`markdup` full stats parity, `sort` external-merge / template-coordinate /
minimiser sorts, full `import` read-group/CRAM parity, and the deferred
`merge -s SEED` header-reconciliation rework.
