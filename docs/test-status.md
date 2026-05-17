# samtools test.pl Status

This tracks the upstream `samtools/test/test.pl` groups against the Rust
`samtools-rs-cli` binary. The CI parity job now fails on a stable upstream
subset and also runs the full harness as an advisory regression watch without
propagating its exit code. Each remaining group must move to `passing` or be
explicitly excluded before the full parity gate can become required.

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
- Required parity subset: `.github/workflows/ci.yml` runs
  `scripts/run-passing-parity-subset.py`, which generates a temporary filtered
  copy of upstream `test.pl`. This is a stable CI subset rather than every
  row marked `passing`; the enforced groups are currently
  `test_reference`, `test_dict`, `test_faidx`, `test_fqidx`, `test_import`, `test_sort`, `test_collate`, `test_calmd`,
  `test_idxstat`, `test_quickcheck`, `test_head`, `test_addrprg`,
  `test_markdup`, `test_bedcov`, `test_split`, `test_coverage`, `test_stats`, `test_depad`, `test_reset`,
  `test_ampliconclip`, and `test_ampliconstats`.
- Required regression subset: CI also runs
  `scripts/run-passing-regression-subset.py`, which executes stable upstream
  `regression.sh` files. The enforced files are currently
  `test/consensus/consensus.reg` and `test/cram_size/cram_size.reg`.
- Full parity watch: CI still stages the Rust binary at the ignored
  `samtools/samtools` path and runs `cd samtools && perl test/test.pl || true`.
  Remove `|| true` only after all rows below are `passing` or explicitly
  skipped with documented rationale.

## Harness Groups

| `test.pl` group | Status | Evidence / next work |
| --- | --- | --- |
| `test_reference` | passing | **Entire upstream `test_reference` byte-exact.** SAM/BAM/CRAM MD-tag reconstruction; CRAM MD path with `-T`; embed_ref CRAM read+write in the vendored noodles fork + `view -O cram,embed_ref=1`; `-e`/`--embedded` faithful `cram2ref` extraction. Building the embed_ref CRAM via `view -e EXPR -O cram,embed_ref=1 -T ref`, all four invocations — `reference` (MD, no `-T`), `reference -e`, and both `-r 17:1000-1500` variants — match `reference/mpileup.{MD,embed}.fa{,.reg}.expected` (tests `reference_embed_ref_full_test_reference_byte_exact`, `reference_cram_md_path_with_reference_matches_upstream`). |
| `test_bgzip` | not-yet-ported | `bgzip` is an htslib tool, not currently in the samtools-rs binary scope. Decide whether to exclude from this parity run or add an htslib-rs CLI. |
| `test_faidx` | passing | The full upstream `test_faidx` group passes (8/8): FASTA index creation plus uncompressed and gzip-compressed region retrieval for the checked local ranges. Non-fixture follow-up: broader BGZI edge cases, compressed output/indexing options, thread/compression-level effects, and exact warning text beyond the harness cases. |
| `test_fqidx` | passing | The full upstream `test_fqidx` group passes (16/16): FASTQ index creation plus uncompressed and gzip-compressed region retrieval for the checked local ranges. Non-fixture follow-up: broader BGZI edge cases, compressed output/indexing options, thread/compression-level effects, and exact warning text beyond the harness cases. |
| `test_dict` | passing | The full upstream `test_dict` group passes (5/5): local FASTA, bgzip-compressed FASTA, stdin FASTA, alias header output, and alt-location header output. Covered by `crates/samtools-rs/tests/dict.rs`. |
| `test_index` | partial | BAM/CSI/CRAI/SAM index build paths exist; threads, full `view -X`, merge/index interactions, and exact binary parity still need harness verification. |
| `test_mpileup` | partial | Default text pileup is implemented on the `htslib-rs` pileup iterator (no longer a stub): multi-input + `-b`, `-f` reference, `-r` region, `-Q`/`-q`/`--ff`/`--rf`/`-A`/`-x`/`-o`, faithful `pileup_seq` encoding, HTSlib smart-overlap removal + orphan filter. Byte-exact vs upstream `mpileup.out.3` (`-B --ff 0x14`) and `mpileup.out.5` (overlap); `mpileup.out.1` matches on depth + read bases (quality chars differ only where HTSlib applies BAQ). Tests `mpileup_minus_b_ff_matches_upstream_out3`, `mpileup_overlap_removal_matches_upstream_out5`. BAQ-adjusted qualities from the completed BAQ/probaln surface, `@RG`-`SM` sample grouping, VCF/BCF (`-g`/`-v`), `-a`/`-aa` remain. |
| `test_usage` | partial | Top-level dispatcher/help exists. Full upstream usage text for every subcommand is not yet verified. |
| `test_view` | partial | SAM/BAM/reference-backed CRAM basics, stdin paths, count/header, filters, BED, region queries, expression filtering, `-z`/`--sanitize` record mutation, `-N`/`--qname-file FILE` qname allow/deny lists (with `^FILE` negation), `-r STR` / `-R FILE` read-group filtering, `-n` exclude-no-read-group filtering, and `-d TAG[:VAL]` / `-D TAG:FILE` aux-tag presence/value filtering are covered by `crates/samtools-rs/tests/view.rs`. Default `@PG` insertion applies to SAM-output paths and now also to **binary output**: SAM→BAM, SAM→CRAM, and BAM→BAM (no-filter/no-region) carry the samtools `@PG` in the binary header (via `sam_bytes_with_pg` / htslib-rs `write_bam_from_path_transforming_header`), suppressed by `--no-PG`; tests `view_b_embeds_pg_in_binary_bam_header`, `view_c_embeds_pg_in_binary_cram_header`. CRAM-producing paths support stdout redirection, e.g. `view -C -T ref in.sam > out.cram`, as well as `-o FILE`. `-U`/`-p` for SAM-input BAM and CRAM output is supported via a text → binary roundtrip. `-X`/`--customized-index` accepts the legacy `in.bam in.bam.bai [region…]` synopsis. `-l`/`--library STR` filters by `@RG LB:`. Remaining: binary `@PG` for CRAM-input and the filtered/region binary-copy sub-paths, BAM/CRAM-input binary aux mutation, BAM/CRAM-input `-U`/`-p`, multi-input, paired filters, and full CRAM parity. |
| `test_head` | passing | The full upstream `test_head` group passes (10/10): default output, `-h` 0/1/5/29/30, `-n` 0/1/5, and stdin with `-h 5 -n 5`. Rust coverage remains in `crates/samtools-rs/tests/head.rs` plus `commands::head` unit tests. |
| `test_cat` | partial | Record-level SAM and BAM concat exists with `-h` header replacement, `-b FILE` input lists, default `@PG` insertion, `--no-PG`, and `-r region` for indexed BAM; SAM output routes through `sam_render` for htslib `%g` float aux spelling; BGZF fast path, CRAM, and `-p N/M` remain. |
| `test_import` | passing | The full upstream `test_import` group passes (21/21): single `-0` and interleaved `-s` FASTQ roundtrips through `fastq`, paired positional and `-1`/`-2` imports, CASAVA `-i` and alternate barcode tags, selected/all FASTQ definition aux tags, `-R`/`-r` read-group forms, positional interleaved FASTQ detection, explicit `--i1`/`--i2` index reads with barcode/quality tags, and UMI extraction. Non-fixture follow-up: CRAM output. |
| `test_bam2fq` | partial | Basic SAM/BAM FASTQ/FASTA output exists, including `-f`/`-F`/`-G` flag filters, `-0` as a single-stream output target, `-n`/`-N` read-name suffix controls, `-O` original-quality `OQ` tag output, `-v INT` missing-quality defaults, `-U`/`--UMI-tag` UMI read-name suffixes, `-i`/`--barcode-tag` CASAVA barcode fields, upstream-style name-grouped `-1`/`-2`/`-s`/`-0` split outputs (paired R1+R2 → `-1`/`-2`; R1-only or R2-only singletons → `-s` with fallback to `-1`/`-2`; READ_OTHER → `-0` with fallback to `-s`), per-record interleaved output when `-1` and `-2` paths alias to the same file, accumulating `-t` plus `-T TAG,…` aux-tag selections, repeated `-d TAG:VAL` / `-D TAG:FILE` invocations that union value sets for the same tag (with mismatched-tag rejection), and `--i1`/`--i2` index FASTQ extraction with `--index-format` (default `i*i*`) and `--quality-tag` (default `QT`), now emitting one index record per adjacent qname-group with htslib-exact CASAVA barcode normalization (`ac-gt` → `AC+GT`) and, under `-i`, the CASAVA comment on index records. `bam2fq/{1,2,3,4,6,7,9,11,13,15,16,17,18,19,20}.{1,2,s}.fq.expected`, `bam2fq/11.fa.expected`, and `bam2fq/{5,8,10,12}` (every output, incl. `2.fq` via cross-mate barcode propagation) pass against the current Rust binary. CRAM remains. |
| `test_depad` | passing | The full upstream `test_depad` group passes (9/9): SAM and BAM input, default BAM output, `-s` SAM output, and `-u`/`-1` BAM output modes against the padded `depad.001` fixture with `-T` and `--no-PG`. The implementation reuses the locked SAM depadding transform and converts the depadded SAM stream to BAM for binary output. Non-fixture follow-up: CRAM input/output. |
| `test_stats` | passing | The full upstream `test_stats` group passes (42 total: 38 passed, 4 expected failures): SN, FFQ/LFQ, GCF/GCL, target and positional regions, overlap handling, barcode histograms, read-group/sample filtering, big-deletion and ref-stats RFS fixtures. The final blocker was the `--ref-stats -t targets ref1` harness variant; RFS reporting now prefers the target file's intervals while positional regions still drive the no-`-t` case. Test `stats_matches_upstream_stat_fixtures` covers the locked fixture set. Non-fixture follow-up: broader CRAM no-region per-cycle/quality/COV parity and exact pileup-backed COV edge cases. |
| `test_depth` | partial | CIGAR-walk BAM/reference-backed CRAM depth exists, including `-r`, `-b`, `-a`/`-aa`, `-d`, `-q`, `-o`, `-H`, `-f` input lists, flag filters, `-l` minimum read length filtering, and multi-input depth columns. Exact pileup overlap/deletion behavior and CRAM without explicit reference remain. |
| `test_merge` | partial | Core merge fixtures are byte-exact (`merge/{2,4,5,6,7}.merge.expected.sam`, modulo harness-stripped `@PG`): glibc-LCG `-s SEED` `@RG`/`@PG` reconciliation, raw-header preservation, `-r` (filename `@RG`), `-c`/`-p` (combine identical IDs). `-L BED` now accepts regions named by `@SQ AN` aliases after the SAM-to-BAM alias writer fix (`merge_l_bed_resolves_reference_aliases`). Current full upstream group still has diffs in `merge/3` and PG-tag sort cases, plus unsupported `--template-coordinate`, so it is not in the required CI subset yet. Test `merge_reconciles_rg_pg_byte_exact_vs_upstream`. Streaming/k-way merge and CRAM remain. |
| `test_sort` | passing | **Every upstream `test_sort` fixture is byte-exact (modulo the harness' `@PG`/`VN` strip)**: coordinate (`pos`), `-n` natural name (`name`,`name3`), `-N` lexicographical name (`name2`), `-t` aux-tag (`tag.rg`,`tag.rg.n`,`tag.as`,`tag.fi`), `-M` minimiser non-indexed/`-I` indexed/`-MH` squash (`minimiser-{basic,indexed,indexed-poly}`, incl. the `reset --dupflag` fresh-header rebuild), and `--template-coordinate` (full `template_coordinate_key`/`bam1_cmp_template_coordinate`/`unclipped_*`/`lookup_libraries` port, `@HD GO:query`). Tests `sort_matches_upstream_test_sort_fixtures` + `sort_minimiser_all_variants_match_upstream`. Remaining (not fixture-blocking): external/temp-file merge for very large inputs, memory/thread caps, CRAM output. |
| `test_collate` | passing | Byte-exact vs the ENTIRE upstream test_collate harness (6/6): tolerant SAM reader, the exact bamshuf order (bucket by hash_X31_Wang(qname)%64, sort by hash,qname,flag>>6&3), raw-header SAM with @HD SO:unsorted GO:query, -o-extension format inference, and fast-mode ring eviction (evict-after-insert deferral). Tests collate_* + collate_matches_upstream_test_collate_fixtures in sort_merge.rs. On-disk hash-bucket for >memory inputs and CRAM output remain. |
| `test_fixmate` | partial | Basic adjacent name-sorted BAM and SAM mate fixup exists with coordinate-sort rejection, TLEN recalculation, default MC/MQ tags, `-m` mate-score tags, `-c` template-CIGAR `ct` tags, and default sanitizer mutation matching the upstream `sanitize.sam` fixture; `-r` removes secondary/unmapped alignments and clears `PROPER_PAIR`/`MATE_REVERSE` on the surviving mate when its partner is unmapped. SAM output routes through `sam_render` for htslib `%g` float aux spelling. Mate rescore, base-modification `-M` parity, and CRAM remain. |
| `test_calmd` | passing | The upstream `test_calmd` invocation (`calmd -uAr mpileup.1.sam mpileup.ref.fa` → BGZF) passes: getopt-style `-uAr` cluster split, `-A` applies recalculated BAQ to QUAL, `-b`/`-u` emit BGZF BAM output; SAM-input MD/NM recompute and BAQ paths with `@PG`/`--no-PG`. Test `calmd_dash_u_a_r_emits_bgzf_bam_like_upstream` (BGZF magic + 569/569 round-trip). Remaining (not fixture-blocking): `-C`/`-n`, CRAM I/O, BAQ over BAM/CRAM input, full upstream MD/BAQ byte parity. |
| `test_idxstat` | passing | Byte-exact (stdout + stderr) vs `idxstats/test_input_1_a.bam.expected` for `test_input_1_a.bam`, `.sam`, and `.cram`. CRAM without an explicit reference is supported via the htslib-rs synthesizing-reference summary path. |
| `test_quickcheck` | passing | TODO marks byte-for-byte parity against `quickcheck/all.expected`; covered by `crates/samtools-rs/tests/quickcheck.rs`. |
| `test_reheader` | partial | Basic BAM header replacement exists with default `@PG` insertion, `--no-PG` suppression, and BAM `-c <command>` header filtering. The `reheader/1_view1.sam.expected` pipeline (`reheader hdr.sam in.bam \| view -h --no-PG`) is byte-for-byte after the harness' VN-strip + header reorder (`@PG` field order and aux float spelling now match upstream). BGZF fast path and CRAM in-place remain. |
| `test_addrprg` | passing | The full upstream `test_addrprg` group passes (14/14), including threaded duplicates of every case: `overwrite_all`, `orphan_only`, unknown `-R` rejection stderr, full `@RG` `-r`, `-R` existing ID reuse, incremental `-r ID`/`-r CN`, and `-w` header edit. SAM output routes through `sam_render` for htslib `%g` float aux spelling. `-O cram`/`--output-fmt[=]cram` with `-T`/`--reference` writes reference-backed CRAM (SAM/BAM input, via temp-BAM → shared CRAM writer). Non-fixture follow-up: CRAM input and mate-aware behavior. |
| `test_markdup` | passing | **Byte-exact vs the entire upstream `test_markdup` SAM harness — all 14 fixtures `markdup/{5..18}`**. Faithful `bam_markdup.c` port: `make_pair_key` (template + `--mode s` sequence, unclipped coords via CIGAR/`MC`, `R_LE`/`R_RI` mate discriminator), `make_single_key`, `calc_score`+`ms` with QCFAIL/qname tie-break, `-S` `dup_hash` propagation, `get_coordinates_colons` + regex `get_coordinates`, the full `find_duplicate_chains` optical-chain re-tagging, `--use-read-groups`, `--duplicate-count`, `--read-coords`/`--coords-order`/`--barcode-rgx`/`--barcode-name`, raw-header SAM output. Test `markdup_matches_upstream_test_markdup_fixtures`. Exact `-s` stats counts, the `1..4` expect-fail cases, and CRAM remain. |
| `test_bedcov` | passing | The full upstream `test_bedcov` group passes (8/8): BAM coverage, `-j`, attached `-g512 -G2048`, `-c`, and all `-H` header cases including custom headers, empty source header fields, and BED12-derived placeholder columns. Non-fixture follow-up: CRAM without explicit reference. |
| `test_split` | passing | The full upstream `test_split` group passes (18/18), including threaded duplicates of every case: `%#`/`%!` template expansion, padded output indexes, `-d RG`, string/integer aux-tag splitting, `-M` overflow routing, sorted-by-tag stdin pipelines (`sort -t nn | split -`), unaccounted output, and harness header reordering / `@PG` normalization. SAM output routes through `sam_render::write_record`. Non-fixture follow-up: CRAM and sorted-by-tag streaming beyond the in-memory harness path. |
| `test_large_positions` | partial | Some large-position behavior routes through `view`, `index`, `merge`, and `depth`; `depth -r` and `depth -b` exist for indexed BAM. Full harness group includes unported `tview` and parity-sensitive index/query cases. |
| `test_ampliconclip` | passing | Full `bam_ampliconclip.c` port (from a stub). **Byte-exact vs the entire upstream `test_ampliconclip` harness** (10 SAM fixtures + 3 primer-count TSVs): per-ref BED, `matching_clip_site`, soft/hard `bam_trim_left`/`right`, `--both-ends`, `--original` `OA`, `--keep-tag`/`NM`-`MD` drop, `--filter-len`/`--fail-len`/`--unmap-len`, `--strand`, `--primer-counts`, raw-header `@HD SO:coordinate→unknown`. Test `ampliconclip_matches_upstream_test_ampliconclip_fixtures`. CRAM/BGZF fast path remain. |
| `test_ampliconstats` | passing | Full `amplicon_stats.c` port (~1776 lines, from a stub). **Byte-exact (modulo harness-stripped version/command-line) vs the entire upstream `test_ampliconstats` harness** (`stats`, `stats_mixed`, `stats_partial`): `count_amplicon`/`bed2amplicon`, ±`pos-margin` lookup, `accumulate_stats` (qname overlap removal, depth, amplicon classification, tcoord freq), `append_lstats`, the full multi-section `dump_stats` with `COMBINED` MEAN/STDDEV. Test `ampliconstats_matches_upstream_test_ampliconstats_fixtures`. `--tcoord-bin` aggregation, CRAM, `--use-sample-name` remain. |
| `test_reset` | passing | The full upstream `test_reset` group passes (18/18): basic.1.mp.1 (reset\|view, stdin, file), basic.output.mp.1 (-o SAM from stdin), basic.bam.input (-o SAM from BAM), basic.cram.input (-o SAM from CRAM with adjacent FASTA discovery when no explicit reference is supplied), output.nRG.* (`--no-RG` plus keep-tag precedence), output.keep.* (`--keep-tag`, `--remove-tag`, `--remove-tag ^...` unioned with `--keep-tag`), output.flg.* (flag update + reverse flip), and reject.1/reject.2 (the positional `--reject-PG` "onwards" removal). `reset_matches_upstream_test_reset_fixtures` locks the fixture set. Non-fixture follow-up: broader CRAM reference-discovery parity and binary CRAM output. |
| `test_checksum` | partial | Default SAM/BAM checksum output is implemented and Rust-tested against upstream `checksum/chk1.1.expected` and `checksum/chk1.3.expected` after the harness' path-line normalization. `-T` TSV output, `-O` order-specific hashing, `-P` position columns, `-C` CIGAR columns, `-M` mate columns, `-B` bamseqchksum-compatible formatting, `-a` field-selection shorthand with upstream-style sanitizer defaults, `-z`/`--sanitize` record mutation, wildcard/exclusion scalar/string/array aux tags with canonical integer encoding, and `-m` merging work for default/position/CIGAR/mate-column checksum reports. CRAM input and full upstream group parity remain. |
| `test_coverage` | passing | The full upstream `test_coverage` group passes (6/6): default sample coverage, `--min-depth`, `-Q`/`-q`, and multi-input aggregate rows. Non-fixture follow-up: exact UTF-8/sidebar histogram parity and CRAM without explicit reference. |

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

The required CI subset is maintained in
`scripts/run-passing-parity-subset.py` and
`scripts/run-passing-regression-subset.py`; update their default lists whenever
a row or regression file is stable enough to gate.
