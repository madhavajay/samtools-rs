# TODO: Port samtools to Pure Rust

Goal: build a pure Rust replacement for the `samtools` C program with full subcommand parity, then port and pass the upstream `test/test.pl` suite plus add Rust-native unit/integration tests. Implementation routes through `htslib-rs` (sibling submodule); long-term, when a needed HTSlib API is not yet exposed there, extend `htslib-rs` first.

## Status: COMPLETE ✅ (2026-05-18)

**The samtools→Rust port is functionally complete and parity-green.**

- All ~40 in-scope subcommands implemented (per the kickoff scope:
  excludes interactive `tview` curses/HTML, remote I/O, `misc/`,
  `lz4/`, C ABI — these were never in scope).
- The full upstream `samtools/test/test.pl` passes, bgzip-honest:
  **998 total / 966 passed / 0 failed / 32 expected-failure /
  0 unexpected-pass.** Every harness group is `passing`. The 32
  expected-failures are upstream `test.pl`'s own `expect_fail`
  negative/edge tests (e.g. `view -d`/`-D` invalid-tag-syntax cases,
  `stats` barcode-fail fixtures); samtools-rs fails them identically
  to C samtools — 0 unexpected-pass confirms exact parity.
- The prioritized backlog (Tasks 0–3) is **fully done and merged**:
  trust reset + authoritative count reconciliation, full
  `seqs_per_slice`/`slices_per_slice` (all write paths), exit-code
  breadth, paired/mate count-mode expr, thread byte-identity, and the
  `cat`/`reheader` BAM BGZF block-level fast paths. CI badges on both
  the samtools-rs and htslib-rs READMEs.
- Required gates green: `cargo fmt`/`clippy`/`test --workspace`, the
  stable parity + regression subsets, and the full `test.pl`.

**Remaining is optional, no failing fixtures** — perf and C
edge-parity only: CRAM container-level `cat`/`reheader` fast path (BAM
done), deeper internal thread *parallelism* (correctness already
byte-identical), `mpileup` VCF/BCF (`-g`/`-v`) and `@RG`→`SM`
grouping, `sort` tag-streaming beyond the in-memory path, `tview`
indexed-SAM seeking, `rmdup`/pileup COV/GCD edges, `faidx`/`fqidx`
BGZI/compression-level edges. These are tracked under
`What to do next` as scope-with-user-first, not open backlog.

<details><summary>Original "Active Goal" framing (pre-completion, kept
for history)</summary>

The dependency/library blocker batch was completed and folded into
this file: all 12 htslib-rs / vendored-noodles work items
byte/fixture-verified, committed, pinned, summarized below under
**Completed Library / Infra Batch**. All major subcommands had an
implementation path or focused tests where upstream ships no fixture
group (`phase`, `targetcut`); the full upstream `test.pl` harness is a
required green gate.

</details>

## Current Handoff — 2026-05-18 (session 2: bgzip false-green + roll-up merged)

> **Merged baseline:** samtools-rs PR #43 merged to `main`
> (`39d22d5`). htslib-rs submodule pinned at `815428b` (htslib-rs PR #8
> + PR #9 merged to its `main`). Vendored noodles at `f998c0f`
> (noodles PR #8 merged to `madhava/bioscript`). All three repos'
> CIs were green at merge.

**Critical discovery — false-green parity gate.** The local parity
gate had been **silently passing while skipping `test_view`,
`test_cat`, `test_reheader`** (and any other `bgzip`-dependent group)
because `bgzip` was not on `PATH`; upstream `test.pl` treats the
missing tool as a skip, not a failure. This invalidated prior
"promoted/passing" claims in this file and hid **7 real bugs**, all
fixed and merged in PR #43:

1. `view -X` ignored the explicit index path (only "worked" via a
   stale `dat/*.bai`); now stages a temp dir under the default name.
2. CRAM region **count** over-counted ~2× — `query_cram_records_*`
   now applies `record_intersects_region` (noodles CRAM `query` is
   slice-granular). [htslib-rs #9]
3. CRAM `-m`/min-qlen used sequence length, not CIGAR query length —
   added `AlignmentRecordSummary::cigar_query_len()`. [htslib-rs #9]
4. `view -t`/`-T` rejected no-`@SQ` SAM instead of injecting `@SQ`.
5. `cat` wrote SAM text into `.cram` outputs (broke re-`cat`); now
   re-encodes a real CRAM.
6. `reheader` wrote SAM text into CRAM (broke `-c`/`--in-place`); now
   re-encodes a real CRAM.
7. Two rustfmt-version-sensitive lines (CI stable rustfmt ≠ local
   1.9.0); rewritten to format identically across versions.

After the fixes the **full** parity + regression subset and the full
upstream `test.pl` pass in CI (and locally with `bgzip`/`tabix`).
See `What to do next` for the resulting remaining tasks (trust reset,
TODO reconciliation, `seqs_per_slice`, hardening/perf).

## Current Handoff — 2026-05-18

> **Library/infra batch COMPLETE (12/12) and rolled into `TODO.md`.**
> All blockers resolved via htslib-rs plus the owned vendored noodles
> fork. This does **not** mean samtools-rs is finished: full upstream
> harness parity is now required, but non-fixture hardening, explicit
> follow-up triage, and performance work remain.

Merged baseline:
- PR #7, https://github.com/madhavajay/samtools-rs/pull/7, is merged into `main`.
- Merge commit: `ae15e4ad603892912fe0c5175491e1d0e3f210eb`.
- PRs #8–#15 (the work below) were consolidated into `integration-all-prs` and merged back to `main`.
- The next work should start from `main` on a new short-lived branch.

Dependency changes for the current Phase 4/5 branch:
- htslib-rs PR #8, https://github.com/madhavajay/htslib-rs/pull/8, branch `reheader-cram-unmapped-md-nm` — merged; skips MD/NM recomputation for unmapped CRAM records with no reference id, unblocking CRAM→SAM rendering in the promoted `test_reheader` path. Local validation: `cargo test -p htslib-rs` (137 unit tests + 1 doctest).
- Current branch carries an additional vendored noodles CRAM writer/reader fix in
  `htslib-rs`: CRAM writing now preserves explicit mate fields instead of
  recomputing TLEN through downstream mate links, and CRAM reading restores
  unmapped MAPQ as `0`. This is required for the upstream `test_bam2fq`
  generated-CRAM setup roundtrip.
- Current branch also carries the vendored noodles CRAM aux parity fix:
  the CRAM writer no longer stores `RG` twice (as both the CRAM read-group
  data series and an aux tag), and explicit `MD`/`NM` tags are serialized at
  the CRAM-visible aux tail in upstream order.

Current Phase 4/5 branch work:
- Promoted full upstream groups into the required parity subset
  (per-group counts below are **pre-bgzip, illustrative only**;
  superseded by the authoritative bgzip-honest whole-suite figure —
  **998 total / 966 passed / 0 failed / 32 xfail** on `main`
  `81b4d87`, see "Authoritative whole-suite parity". The per-group
  splits are kept for historical context, not as verified counts;
  `test_view`/`test_cat`/`test_reheader` in particular were silently
  skipped when first promoted and only genuinely passed after the PR
  #43 fixes):
  `test_merge` (28/28), `test_fixmate` (42 total: 40 passed,
  2 expected failures), `test_reheader` (7/7), `test_cat` (26/26),
  `test_index` (26/26), `test_checksum` (14/14),
  `test_large_positions` (9/9), `test_mpileup` (7/7), and
  `test_view` (445 total: 427 passed, 18 expected failures).
- `view` filtered SAM-input output now uses the shared SAM renderer for
  htslib-style aux float spelling, including `B:f` arrays.
- `view -r` / `-R` SAM output prunes `@RG` header lines to the selected
  read groups; `view -l` / `--library` still filters records by
  `@RG LB:` but deliberately keeps all `@RG` header lines, matching
  upstream.
- `view -C` BAM/CRAM output, CRAM→BAM file output, CRAM stdin text/BAM
  output, and CRAM count/save-count paths now discover references from
  explicit `-T`, header `@SQ UR:`, or the M5/`REF_PATH` lookup path.
  No-reference CRAM count/save-count paths now cover summary-backed
  MAPQ/flag expressions plus read-group, library, and aux-tag filters.
  No-reference CRAM non-region text and BAM output now also support
  reference-independent expression filtering with `--save-counts`; verified
  byte-for-byte vs C samtools for `view --no-PG --save-counts ... -e
  'mapq>=20'` to SAM output and BAM output rendered back to SAM.
  BAM/CRAM input binary output paths now support `-U`/`-p` split/unmap
  routing and `-x`/`--keep-tag` aux rewriting for BAM and reference-backed
  CRAM output, including CRAM `RG`/`MD`/`NM` roundtrip parity with C
  samtools. The full upstream `test_view` group is now promoted; remaining
  `view` follow-up is non-fixture hardening around reference-dependent
  expression counts, multi-file inputs, paired-aware filters, and deeper
  CRAM performance/streaming parity.
- `checksum` is promoted into the required parity subset. It now handles
  embedded-reference CRAM without explicit `--reference` for the upstream
  checksum fixtures, and `split` can split whole-file CRAM input to SAM/BAM
  outputs for the checksum split/merge harness path.
- `bam2fq` / `fastq` is now promoted: the full upstream `test_bam2fq`
  group passes (84/84), including threaded duplicates. This branch fixes
  SAM/BAM default single-output reverse-complementing, paired stdout routing
  when only discard files are specified, headerless `.sam` index extraction
  (`bam2fq.014.sam`), compact `-Tfoo` parsing, `--no-sc` /
  `--no-sc-bkp` / `--sc-aux` soft-clip trimming with backup aux fields,
  CRAM input to `fastq` / `fastq -n`, and the generated `view -C` CRAM
  setup roundtrip via the vendored noodles CRAM mate/MAPQ preservation fix.
- `test_usage` is now verified under a real TTY: 40/40 passing. The
  dispatcher prints stdout usage for no-argument subcommands invoked from
  a terminal, and `tview` is registered so every advertised top-level
  command is dispatchable.
  This row is not in the default non-TTY parity subset because upstream
  `test.pl` skips usage tests when no TTY/PTY is available.

Additional local validation after the no-reference CRAM / quickcheck
updates: `cargo fmt --all --check`;
`cargo test -p samtools-rs --test quickcheck`;
`cargo test -p samtools-rs-cli --test quickcheck`;
`cargo test -p htslib-rs alignment_compat::tests::cram_summaries_without_reference_match_bam_flags_and_tids`;
`cargo test -p samtools-rs --test view view_count_save_counts_no_reference_cram_supports_summary_expr_filters`;
`cargo clippy -p samtools-rs --all-targets -- -D warnings`;
`cargo clippy -p htslib-rs --lib -- -D warnings`;
direct C-vs-Rust probes for no-reference CRAM `view --save-counts` SAM/BAM
output; and the full `scripts/run-byte-parity-smoke.py` suite using the
rebuilt Rust binary under `/dev/shm`.

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

Historical working branch: `work-sam-float-renderer` (merged via the
samtools-rs-only batch).
Landed slices on that branch:
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
  `RG:Z:` membership while keeping all `@RG` header lines. New test:
  `view_dash_l_filters_by_read_group_library`.
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
  no-region CRAM path was later unblocked by the completed htslib-rs
  all-record iterator.
- **`fastq` index × name-grouping (complete):** one index record per
  qname-group, htslib-exact CASAVA barcode normalization
  (`ac-gt` → `AC+GT`), the CASAVA comment on `-i` index records, and
  cross-mate barcode propagation (R2/other inherits the R1 mate's
  `BC`). **`bam2fq/{5,8,10,12}` now byte parity on every output.** New
  tests: `fastq_index_emits_one_record_per_qname_group_with_casava_comment`,
  `fastq_casava_barcode_propagates_from_r1_to_r2_mate`.

Batch status: **all "Remaining tractable samtools-rs-only items" are now
done.** The one item deferred during this historical branch,
`merge -s SEED`, was later completed with upstream-style seeded
`@RG`/`@PG` reconciliation and is now byte-exact vs the upstream merge
fixtures. Per the **Active Goal**, the remaining actionable work moved
to the now-completed library/infra batch summarized below. Final
validation on `work-sam-float-renderer`:
`cargo fmt --all --check`,
`cargo clippy --workspace --all-targets -- -D warnings`, and
`cargo test --workspace` (**2601 passing, 0 failing**) all green; the
upstream `bam2fq/{5,8,10,12}` fixtures now match byte-for-byte.

**Authoritative whole-suite parity (2026-05-18, `main` at `81b4d87`,
post-PR-#44, release binary, `bgzip`+`tabix` on `PATH` — bgzip-honest):**
full upstream `perl test/test.pl`:
- total **998**, passed **966**, failed **0**, expected failure **32**,
  unexpected pass **0**, exit 0.
- This supersedes every pre-bgzip per-group tally below as the
  trustworthy project-level signal. `test.pl` emits no per-group
  breakdown, so individual "NNN total: NNN passed" lines are *not*
  separately re-measured — treat the whole-suite figure as
  authoritative and the per-group splits as historical/illustrative.
- Reproduced via `scripts/run-passing-parity-subset.py` /
  `run-passing-regression-subset.py`, which now hard-fail without
  `bgzip`/`tabix` (Task 0 preflight guard), so this number cannot
  silently regress to a false green again.

Latest known validation (on `main` at `39d22d5`, post-PR-#43-merge):
- Full samtools-rs CI green on PR #43: rust-gate (`cargo fmt --all
  --check`, `cargo clippy --workspace --all-targets -- -D warnings`,
  `cargo test --workspace`) **and** the parity-gate
  (`run-passing-parity-subset.py`, `run-passing-regression-subset.py`,
  and the full upstream `perl test/test.pl`) — the parity-gate is now a
  required, honest signal (CI installs `tabix`).
- Reproduced locally with `bgzip`/`tabix` on `PATH`: full parity +
  regression subset rc=0; `cargo test --workspace` green for both
  samtools-rs and htslib-rs.
- Older note (stale, kept for history): "on `main` at `b312c99`: 416
  `samtools-rs` passing, `cargo test --workspace` 2587 passing; PR #16
  parity gate advisory at the time." Pre-dates the bgzip trust reset —
  do not rely on its counts.
- New focused tests added across PRs: `fastq_index_files_extract_from_barcode_tag`, `fastq_routes_r1_only_singletons_to_singleton_output`, `fastq_dash_t_and_dash_cap_t_combine_aux_tags`, `fastq_interleaves_read1_read2_when_paths_alias`, `fastq_repeated_dash_d_unions_same_tag_values`, `fasta_reverse_strand_record_reverse_complemented_in_output`, `view_qname_file_filters_records_by_name`, `view_r_and_dash_cap_r_filter_by_read_group`, `view_d_and_dash_cap_d_filter_by_aux_tag`, `addreplacerg_defaults_to_first_header_rg_and_preserves_lines`, `addreplacerg_r_overwrite_all_removes_other_header_rg_lines`, and `addreplacerg_dash_cap_r_unknown_id_is_rejected`.

Estimated whole-project completion: **~88–92% (good confidence)**.
- Breadth ~95% (all ~40 subcommands implemented). Parity is now a
  trustworthy signal: the authoritative bgzip-honest whole-suite run
  above is **998/966 passed, 0 failed, 32 expected failure** on current
  `main`, and the preflight guard (Task 0) makes a future false green
  impossible. The trust reset that previously capped confidence is
  resolved.
- Confidence raised from *moderate* to *good*: the count uncertainty is
  retired at the project level (Task 1 done). Residual gap (~8–12%) is
  genuine non-fixture hardening + perf + a few latent encoder gaps
  (`seqs_per_slice`, Task 2), not unknown parity risk.
- **Not done.** The library/infra blockers and Phase 3 command
  implementation are complete, and the full upstream `test/test.pl` harness
  is now required, but `TODO.md` still tracks non-fixture hardening and
  exit-code/thread/perf triage before the port can be considered complete.
- Current required gates cover the stable parity subset in
  `scripts/run-passing-parity-subset.py`, the stable regression files in
  `scripts/run-passing-regression-subset.py`, and the full upstream harness.
  Remaining project work is Phase 4/5: close or explicitly defer the
  non-fixture pending items and finish performance triage.

### Workflow rule: one large PR branch at a time

Do **not** open many small per-slice PRs. Use a single long-lived working
branch off `main` (e.g. `work-<topic>`), land multiple related bounded
slices onto it as separate commits, keep the full gate green after each
commit, and open **one** large PR for that batch. Only start the next
working branch after that PR has merged to `main`. (PRs #8–#15 were the
fragmented anti-pattern; #16 consolidated them — keep it to one large PR
branch going forward.)

What to do next (remaining tasks, priority order):

> **Status 2026-05-18:** Tasks 0, 1, 2, **and 3 are complete.**
> Task 0 + 1 merged in samtools-rs PR #45; Task 2 via noodles #9 +
> htslib-rs #10 + samtools-rs #46. Task 3's tracked follow-ups all
> landed: filter/region `seqs_per_slice` (htslib-rs #11, samtools-rs
> #48), exit-code breadth (#49), paired/mate count-mode expr (#50),
> `stats -@` thread byte-identity (#51), and the cat/reheader BGZF
> block-level fast path (noodles #10 → htslib-rs #12 → samtools-rs
> #53 (cat) + #54 (reheader)). The prioritized Tasks 0–3 backlog is
> **exhausted**; what remains below is open-ended deeper hardening
> with no fixture failures, to be scoped with the user. Tasks 0–3
> are kept as a closed record, not open work.

0. **[DONE 2026-05-18] Trust reset — local gate honesty.** The local parity
   gate silently false-greens: `scripts/run-passing-parity-subset.py`
   and `run-passing-regression-subset.py` drive upstream `test.pl`,
   which **skips** any group whose data-gen needs `bgzip` when `bgzip`
   is missing from `PATH`, then returns 0. This hid 7 real bugs (see
   the 2026-05-18 session-2 handoff). CI is honest (installs `tabix`);
   only local runs lie.
   - [x] Add a preflight check to both subset scripts that hard-errors
     (non-zero, loud message) if `bgzip`/`tabix` are not on `PATH`,
     instead of letting test.pl skip silently. Done: shared
     `scripts/_parity_preflight.py`, wired into both runners; exits 3
     with a remediation message, `--allow-missing-bgzip` /
     `SAMTOOLS_RS_ALLOW_MISSING_BGZIP=1` opt-out for degraded runs.
   - [x] Document the requirement (`README`/dev setup): prebuilt at
     `htslib-rs/htslib/{bgzip,tabix}` via
     `make -C htslib-rs/htslib tabix bgzip`; export that dir on `PATH`.
     Done: README "Parity testing" section.

1. **[DONE 2026-05-18] Reconcile TODO.md / `docs/test-status.md` with
   reality.** Resolved at the project level: the authoritative
   bgzip-honest whole-suite run (see "Authoritative whole-suite parity"
   above) is **998 total / 966 passed / 0 failed / 32 expected
   failure** on `main` `81b4d87`. `test.pl` emits no per-group
   breakdown, so the individual pre-bgzip "NNN total: NNN passed"
   tallies are not re-measured one-by-one; they are annotated as
   historical/illustrative and the whole-suite figure is the
   authoritative signal. The preflight guard (Task 0) prevents this
   from silently regressing again.

2. **[DONE 2026-05-18] `seqs_per_slice` plumbed through the CRAM
   writer.** Was latent: `view -O cram,seqs_per_slice=N` produced one
   container/slice regardless of N (the `-O cram,<opt>` suffix loop
   only honored `embed_ref` and the value was parsed as a no-op).
   Fixed end-to-end via the bottom-up roll-up:
   - noodles PR #9 (merged): `records_per_slice` /
     `slices_per_container` on writer `Options` + `Builder` setters,
     default-preserving; sync+async writers and the container builder
     size from them.
   - htslib-rs PR #10 (merged, `main` `140aa4f`): `CramWriteOptions`
     + `write_cram_from_{sam,bam,path}_with_reference_and_options`
     entry points; noodles pin bumped.
   - samtools-rs PR #46 (green, both gates): `Opts.records_per_slice`
     / `slices_per_container`; `apply_output_fmt_option` parses
     `seqs_per_slice=` / `slices_per_slice=`; the `-O cram,<opt>`
     suffix loop now routes through it (still lenient on unknowns);
     CRAM-write call sites use the `_and_options` variants.
   - Verified: `ce#1000.sam` → 1 container by default, >1 with
     `seqs_per_slice=100` (via `-O cram,...` and
     `--output-fmt-option`). Bgzip-honest `test_view` parity
     unchanged (445 / 427 passed / 0 failed / 18 xfail).
   - Follow-up (not blocking): filter/region CRAM-output combinations
     with `seqs_per_slice` still use defaults (matches the prior
     `embed_ref` scoping to core full-file write paths); plumb
     `_and_options` into the filter/region CRAM writers when needed.

3. **[DONE 2026-05-18] Phase 4/5 non-fixture hardening — tracked
   follow-ups.** All items that were enumerated here landed:
   - [x] `view` paired-aware filters / reference-dependent expression
     counts — already worked; locked in count mode (#50).
   - [x] CRAM `seqs_per_slice` / `slices_per_slice` for the filter and
     region CRAM-output paths (htslib-rs #11, samtools-rs #48; full-
     file paths were Task 2).
   - [x] `cat` BAM container/BGZF block-level copy fast path (noodles
     #10 `read_raw_frame`/`EOF` → htslib-rs #12
     `append_bam_alignment_frames` → samtools-rs #53), all-or-nothing
     with the record-level path as universal fallback.
   - [x] `reheader` BAM BGZF block-level fast path (samtools-rs #54,
     same primitive/fallback).
   - [x] Threads: worker-pool byte-identity invariant extended to
     `stats` (#51); view/sort/index already covered.
   - [x] Exit codes: broadened C-vs-Rust error-class coverage (#49).

   **Remaining (open-ended, no fixture failures — scope with user
   before starting):**
   - `cat`/`reheader` **CRAM** container-level fast path (BAM done;
     CRAM still re-encodes — correct but slow).
   - Deeper thread *parallelism* parity vs C samtools (the correctness
     invariant is locked; actual internal parallelization of remaining
     readers/writers is perf-only).
   - `view` multi-file inputs; deeper CRAM streaming/perf parity.
   - `rmdup`, broader pileup/COV/GCD edge cases, `stats` CRAM
     no-region per-cycle/quality parity.

4. After every commit keep the full local gate green (with `bgzip`/
   `tabix` on `PATH`): `cargo fmt --all --check`, `cargo clippy
   --workspace --all-targets -- -D warnings`, `cargo test --workspace`,
   plus the required parity/regression subset. Update `TODO.md`,
   `docs/subcommand-coverage.md`, and `docs/test-status.md` whenever a
   group is promoted, deferred, or newly understood.

Completed samtools-rs-only batch (no htslib-rs / noodles changes required):
- ~~**SAM aux float formatting — remaining commands.**~~ **Done.** The shared `samtools_rs::sam_render` module (`format_aux_float`, `format_htslib_exponent`, `fix_sam_aux_floats`, `fix_sam_text`, `write_record`, `write_header`) now backs every noodles-`sam::io::Writer` SAM-output path: `view`, `split`, plus `reheader` SAM→SAM, `sort`, `merge`, `collate`, `addreplacerg`, `reset`, `fixmate`, `rmdup`, `markdup`, and `cat` (their `SamFile`/`SamStdout`/`Sam*Sink` sinks now wrap a plain `File`/`Stdout` and render through `sam_render`). So every SAM-text output path emits htslib `%g`-style float aux spelling. Regression covered by `sort_sam_output_uses_htslib_float_aux_spelling` plus the existing `view`/`split` fixtures.
- ~~**`fastq` index extraction × name-grouping interaction.**~~ **Done.** `emit_index_files` dedupes to **one index record per adjacent qname-group** (matching upstream `flush_rec` → `output_index`). The CASAVA barcode field is normalized exactly like htslib `fastq_format1` (`casava_barcode_field`: absent/non-sequence-first → `0`; otherwise non-alpha → `+`, lowercase → upper, e.g. `ac-gt` → `AC+GT`), `-i` index records carry the ` <rnum>:<filt>:0:<barcode>` CASAVA comment, and `GroupedSplitWriter` now **propagates the group barcode across mates** so an R2 (or other) record lacking its own `BC` gets the R1 mate's barcode in its CASAVA comment (`fill_casava_barcode`; upstream `bam_fastq.c:952`). **`bam2fq/{5,8,10,12}` now match byte-for-byte on every output (`1.fq`/`2.fq`/`s.fq` and the index/`bc` files).** New tests: `fastq_index_emits_one_record_per_qname_group_with_casava_comment`, `fastq_casava_barcode_propagates_from_r1_to_r2_mate` (plus the 45 existing fastq tests still green).
- ~~**`view --library` (`-l`)** library filter via `@RG LB:` aux lookup.~~ **Done.** `view -l STR` / `--library STR` resolves the requested library to the set of `@RG` IDs whose `LB:` equals STR (scanned from the input header for path, SAM/BAM/CRAM stdin), then a record passes iff its `RG:Z:` value is in that set (no-RG / non-matching RG excluded, matching upstream `bam_get_library`). Unlike explicit `-r` / `-R`, `-l` preserves all `@RG` header lines. Regression: `view_dash_l_filters_by_read_group_library`.
- ~~**`view -X` legacy custom-index synopsis**~~ **Done.** `view -X` / `--customized-index` accepts the legacy synopsis where the second positional is the explicit index path (`view -X in.bam in.bam.bai [region…]`); accepted as a no-op (our region queries build/find the index themselves), matching `idxstats -X`. Regression: `view_dash_cap_x_accepts_legacy_custom_index_synopsis`.
- ~~**`merge -s SEED`** seeded header reconciliation.~~ **Done.** The later merge parity batch ports upstream-style seeded `@RG`/`@PG` collision suffixing, plus `-r` filename read groups and `-c`/`-p` combine modes. The upstream `merge/{2,4,5,6,7}` fixtures are byte-exact modulo `@PG`; regression: `merge_reconciles_rg_pg_byte_exact_vs_upstream`.
- ~~**`samples` BAM index path verification**~~ **Done.** `samples -i` with a custom `-X` index path now mirrors `sam_index_load3`: an exact index file, a *directory* holding the index (`<dir>/<data-name>.bai`), or a suffix-less prefix all resolve via the shared `locate_associated_index` resolver, so index files at non-default locations register `Y`. Regression: `samples_custom_index_directory_reports_index_presence` (and the existing exact-file/pair test still passes).
- ~~**`addreplacerg --output-fmt=cram`** with a `-T` reference~~ **Done.** `addreplacerg` accepts `-O cram` / `--output-fmt cram` / `--output-fmt=cram` and `-T`/`--reference[=]FILE`; SAM/BAM input → CRAM output spools rewritten records to a temp BAM and converts via the shared `write_cram_from_bam_path_with_reference` (the `.fai` is built if missing). CRAM output without `-T` errors. Regression: `addreplacerg_writes_cram_output_with_reference`.
- ~~**`stats -d` / `--remove-dups` edge cases**~~ **Done.** The CRAM *region* path iterates real records through the same `update_record_with_targets` → `update` chokepoint as SAM/BAM, which already gates all histogram/seq/quality accumulation on `self.total` increasing (and `--remove-dups` filters `BAM_FDUP` before `total` is bumped). Verified end-to-end by `stats_remove_dups_excludes_duplicates_on_cram_region_path` (SAM→CRAM→indexed→region stats with/without `-d`). The CRAM no-region path was later unblocked by the htslib-rs all-record iterator; only optional CRAM NM recompute for exact mismatch/error-rate parity remains as a non-fixture-blocking nicety.

Items previously blocked on htslib-rs / noodles extensions — **ALL
RESOLVED** (see **Completed Library / Infra Batch** below; the
"htslib-rs Extensions Needed" list is fully checked off):
- ~~pileup-dependent commands~~ — pileup iterator done (#1); `mpileup`/
  `consensus`/`coverage`/`bedcov`/`depth` byte-exact (`consensus`
  77/77); `ampliconstats`, `targetcut`, and `phase` done (the latter
  two have no upstream fixtures).
- ~~`stats`/`checksum` CRAM no-region~~ — CRAM all-record iterator
  done (#2); wired.
- ~~`cram-size`~~ — done (#3), all 3 `cram_size.reg` byte-exact.
- ~~binary `@PG`~~ — done (#4), `view` SAM→BAM/CRAM + BAM→BAM.
- ~~`flagstat`/`idxstats` CRAM-no-ref~~ — done (#5), byte-exact.
- ~~large-reference CSI~~ — done (#12).
- ~~SAM aux float formatting~~ — resolved via `sam_render`.

## Progress Snapshot

> ⚠️ **Counts caution:** specific per-group tallies in this snapshot and
> the sections below (e.g. *"test_view 445: 427 passed"*) were measured
> **before** the bgzip trust reset and are illustrative only. The
> authoritative signal is the bgzip-honest whole-suite run: **998 total
> / 966 passed / 0 failed / 32 expected failure** on `main` `81b4d87`
> (see "Authoritative whole-suite parity"; Task 1 done 2026-05-18).
> `test.pl` emits no per-group breakdown, so per-group splits are not
> separately re-verified — cite the whole-suite figure as fact, the
> per-group numbers as historical.

**Phases 0–2 complete; Waves A/B/C/D substantially complete and
byte-exact for the upstream-fixtured subcommands. The completed
library/infra batch delivered all 12 numbered blockers
byte/fixture-verified and committed: pileup iterator, CRAM all-record
iterator, CRAM container/codec inventory for `cram-size`, binary `@PG`,
aux mutation, threads, write-index, SO-less index, region grammar, BAQ,
large-ref CSI, and embed_ref read+write.** Stable upstream groups/regression
files include: `consensus` (all 77 `consensus.reg`), `sort` (all
`test_sort`), `cram-size` (all 3), `reference` (all `test_reference`),
`dict` (all `test_dict`), `faidx`, `fqidx`, `head` (all `test_head`),
`collate`, `calmd`, `idxstats`, `quickcheck`, `index`, `cat`, `reheader`, `addreplacerg` (all
`test_addrprg`), `markdup`, `bedcov`, `split` (all `test_split`),
`coverage`, `import`, `stats`, `checksum`, `depad`, `merge`, `fixmate`, `ampliconclip`, and `ampliconstats`.
Other commands have
substantial Rust coverage and full upstream fixture parity where fixtures
exist. `phase` (`phase.c`, 843 LOC) and `targetcut` (`cut_target.c`, 257
LOC) are now ported with focused unit tests because upstream ships no fixture
groups for them. **Remaining: Phase 4/5 non-fixture hardening and polish
(per-subcommand integration tests largely already in place; thread, exit-code,
and perf triage).**

Subcommands shipped:
- ✅ required upstream parity groups/regression files:
  `quickcheck`, `dict`, `faidx`, `fqidx`, `head`, `sort`, `collate`,
  `calmd`, `idxstats`, `index`, `reference`, `cat`, `reheader`,
  `addreplacerg`, `markdup`, `bedcov`, `split`, `large_positions`,
  `coverage`, `import`, `stats`, `depad`, `merge`, `fixmate`, `reset`,
  `ampliconclip`, `ampliconstats`, `checksum`, `view`, `bam2fq`/`fastq`,
  `depth`, `mpileup`, plus `consensus` and `cram-size`.
- ✅ focused/no upstream fixture group:
  `flags`, `flagstat`, `samples`, `phase`, and `targetcut`.
- 🟡 implemented with non-fixture hardening still tracked:
  `rmdup`, plus the documented non-fixture follow-ups for otherwise
  promoted commands such as `view`, `fastq`/`fasta`/`bam2fq`, `depth`,
  `mpileup`, `stats`, `sort`, `merge`, and `tview`.

Remaining subcommands and their blockers — **resolved** (kept for
history; current truth in the banner + Progress Snapshot above):
- ~~BAM aux-tag mutation~~ — done; all parity consumers byte-exact.
- ~~pileup iterator in htslib-rs~~ — done (#1); `mpileup`/`consensus`
  (77/77) / `ampliconstats` / pileup-based `bedcov`/`coverage`/`depth`
  byte-exact; `phase` is now ported with focused Rust unit tests
  because upstream has no fixtures.
- ~~CRAM all-record iterator~~ — done (#2); `stats`/`checksum`
  no-region CRAM + `reference` MD path wired.
- ~~Other complex algorithms~~ — `cram-size` ✅ (3/3), `reference` ✅
  (full `test_reference`, embed_ref read+write), `checksum`/`markdup`/
  `sort` (incl. minimiser + template-coordinate) byte-exact. All
  subcommands listed for Phase 3 now have implementations; `depad`
  BAM/CRAM and a few non-parity niceties remain
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

See **Active Goal** above: underlying-library blockers are complete and
summarized in this file. New HTSlib-shaped API gaps should still be
implemented in `htslib-rs` first and consumed from `samtools-rs` after
the dependent htslib-rs commit is merged and pinned.

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

- [~] Propagate `threads` into BGZF/noodles worker pools instead of accepting it as an API-compatible no-op in several wrappers. Native `view_region`, `view_bed`, and `extract_unmapped_pairs` now route BAM outputs through noodles' multithreaded BGZF writer when `threads` is nonzero; `faidx`/`fqidx` BGZF input/output paths consume local/global thread counts; `index` now routes local/global thread counts into multithreaded BGZF readers for BAM and BGZF-SAM BAI/CSI construction. Remaining wrappers still need alignment reader/writer worker-pool propagation.
- [ ] Replace in-memory `sort`/`merge` implementations with streaming/external algorithms for large BAMs.
- [ ] Deepen `@PG` parity for native-generated BAM headers.
- [ ] Add broader real-world CRAM fixtures when available from VNtyper/BioScript workflows.

Current parity directions, each substantial:

1. **Close non-fixture follow-ups.** Work leaf-first from the remaining
   `**Pending:**` notes below and either implement the behavior or document
   a deliberate deferral when the full upstream harness does not cover it.
2. **Keep required CI representative.** The full harness is now required;
   keep the faster subset helpers in sync as smoke gates for stable groups
   and regression files.
3. ~~**Remove the advisory full-harness escape hatch.**~~ **Done.** The full
   `samtools/test/test.pl` run is now a required CI gate; `test_bgzip` is
   satisfied by the external htslib `bgzip` tool on `PATH`.

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

- `tview` curses/HTML and interactive viewers (`bam_tview*.c`). A narrow
  noninteractive `-d T -p REGION` text renderer is implemented for the
  large-position fixture.
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
- [~] **@PG add helper** (`samtools-rs/src/pg.rs`): shared helper now builds raw-header `@PG` lines with HTSlib-style argv stringification, generated unique IDs, and upstream `sam_hdr_add_pg` field order `ID, PN, PP, VN, CL` (PP precedes VN/CL so the upstream harness' `s/\tVN:.*//` normalization keeps `PP`), with `PP` links for terminal program chains. `cat`, `split`, `reheader`, `sort`, `merge`, `collate`, `addreplacerg`, `reset`, `fixmate`, `rmdup`, and `view`'s SAM-output paths (file-input header-only, SAM output with `-h`, plus BAM/CRAM stdin SAM/header-only) use it for default output headers and honor `--no-PG`. The upstream `reheader/{1,4}` header section now matches after harness reordering. Binary `@PG` is wired for every `view` input→binary path that noodles can decode (SAM→BAM/CRAM, BAM→BAM/CRAM including plain/filter/region/sanitizer). **Pending:** CRAM-input→binary header transform remains on the direct-copy path until the remaining noodles CRAM decode limitation is addressed; continue verifying complex merge/split/reheader `@PG` chains.
- [x] **Aux-tag list parser** (`samtools-rs/src/aux_list.rs`): port `parse_aux_list` from `sam_utils.c`. Used by `view`, `reset`, `fastq`, and future aux-aware commands.
- [~] **BED index** (`samtools-rs/src/bedidx.rs`): shared BED parser/index now stores 0-based half-open intervals by reference, skips comments/UCSC metadata, emits HTSlib-style 1-based inclusive region strings, supports overlap queries, and is used by `view -L`, `depth -b`, `bedcov`, and native `view_bed`. **Pending:** interval-tree acceleration/parity with `bedidx.c`, stricter upstream diagnostics where needed, and integration into `ampliconclip` and future `mpileup`.
- [~] **Reference helpers** (`samtools-rs/src/reference.rs`): shared FASTA helper now derives associated `.fai` paths, builds missing FASTA indexes through `htslib-rs::faidx_compat`, loads `(SN, LN)` dictionaries, and matches candidate FASTA references against BAM/CRAM `@SQ` dictionaries for `samples -f/-F`. **Pending:** mmap/FASTA sequence cache, common `--reference` option plumbing, CRAM reference resolution, and integration into `calmd`, `consensus`, `mpileup`, and `import`.
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

- [x] `view` (`sam_view.c`, 68k) — the full upstream `test_view` group passes in CI (PR #43 parity-gate, with `tabix`). The "445 total: 427 passed, 18 expected failures" tally is pre-bgzip and unverified; the trustworthy statement is "green in the PR #43 parity-gate". Implemented surface includes SAM↔SAM passthrough, SAM→BAM/CRAM, SAM/BAM/reference-backed CRAM stdin count/text/BAM/CRAM paths, stdin CRAM reference discovery from `-T`, `@SQ UR:`, or M5/`REF_PATH`, reference-backed CRAM→SAM text/count paths including flag/MAPQ filtered count mode, reference-backed CRAM→BAM full-file and region output, reference-backed BAM→CRAM and CRAM→CRAM full-file and region output, header-only / count modes, `-h` `-H` `-c` `-b` `-C` `-T` `-t` `-o` `--no-PG`, filter flags `-f`/`-F`/`-G`/`--exclude-flags`/`-q`, count-only no-reference CRAM `--save-counts` for simple summary-backed filters, reference-independent MAPQ/flag expressions, and read-group/library/aux-tag filters, region queries including `chr:pos` open-ended semantics, `-L FILE` BED restriction with two-column point semantics and positional-region intersection, `-U FILE`/`-p` for SAM/BAM/CRAM input to SAM/BAM/reference-backed CRAM output, `-e EXPR`, `--save-counts FILE`, `-m INT`, `-x`/`--remove-tag`, `--keep-tag` for SAM/BAM/CRAM input to SAM/BAM/reference-backed CRAM output, htslib-style aux float spelling for filtered SAM-input, binary→SAM, and CRAM-stdin SAM output, `-O FORMAT`, default `@PG`, `-N`, `-r`/`-R`, `-n`, `-d`/`-D`, `-B`, `-s`, `--remove-flags`, and `--fetch-pairs`. **Non-fixture follow-up:** reference-dependent expression counts, multi-file inputs, paired filters, and deeper CRAM performance/streaming parity.
- [x] `head` (`sam_view.c` shared) — SAM and BAM input; SAM/BAM/CRAM stdin header/record output; CRAM header-only modes; reference-backed CRAM record extraction for `-n N`; `-h N`, `-n N`, all-default.
- [x] `quickcheck` (`bam_quickcheck.c`) — passes byte-for-byte against `quickcheck/all.expected`.
- [x] `index` (`bam_index.c`) — full upstream `test_index` group passes (26/26), including threaded duplicates: large-reference CSI indexing/query, default BAI output, explicit `-o` and legacy `<in> <out.idx>` index destinations, `-M` multi-file indexing, `view -X` custom-index queries, `merge -X -R` with one custom index per BAM input, and `view --write-index` auto-indexing for BAM/CRAM/BGZF-SAM outputs. Local `-@ N`, attached `-@N`, `--threads N`, `--threads=N`, and top-level `--threads` now route BAM/BGZF-SAM BAI/CSI index construction through noodles multithreaded BGZF readers when nonzero. **Pending (non-fixture):** CRAI thread parity and broader throughput comparison against C samtools.
- [x] `idxstats` (`bam_stat.c`) — index-based per-reference counts for BAM, with streaming slow-path counts for SAM, unindexed BAM, and CRAM. CRAM works with or without an explicit reference via the htslib-rs synthesizing-reference summary path because idxstats only needs reference-independent reference ids and flags; tests cover both explicit-reference and no-reference CRAM, and C-vs-Rust smoke covers `dat/test_input_1_a.cram`. **Pending (non-fixture):** index-derived CRAM counting fast path instead of streaming decode.
- [x] `faidx` / `fqidx` (`faidx.c`) — index-build mode works (`samtools faidx file.fa` produces `file.fa.fai`); BGZF FASTA/FASTQ input now writes `.gzi` and can be indexed/retrieved; local region extraction works for positional regions, `-r` region files, `-o`, `.gz`/`.bgz`/`.bgzf` BGZF output, `--length`, `--write-index` for file outputs, FASTQ mode via `fqidx` and `faidx -f`, reverse-complement `-i` with mark-strand modes, `--continue`-style missing-region tolerance, and upstream-style zero/truncated region warning keywords. Local `-@` / `--threads` and top-level `--threads` now route BGZF FASTA/FASTQ read/write paths through the noodles multithreaded BGZF worker APIs when nonzero. The full upstream `test_faidx` (8/8) and `test_fqidx` (16/16) groups now pass and are promoted into the required parity subset. Non-fixture follow-up: exact warning text parity, compression-level option effects, and broader BGZI edge cases.
- [x] `dict` (`dict.c`) — sequence dictionary builder. Passes byte-for-byte against `dict.out`, `dict.alias.out`, `dict.alt.out` (run via test.pl-style stdin/file invocations).
- [x] `flagstat` / `flagstats` (`bam_stat.c`) — SAM, BAM, and CRAM input. Default + `-O json` + `-O tsv` output modes. CRAM works with or without an explicit reference via the htslib-rs synthesizing-reference summary path because flagstat only needs reference-independent flags; tests cover both explicit-reference and no-reference CRAM, and C-vs-Rust smoke covers BAM text/JSON/TSV plus `dat/test_input_1_a.cram` no-reference text. Required extending `htslib-rs::alignment_compat::AlignmentRecordSummary` with `flags_u16` / `reference_sequence_id` / `mate_reference_sequence_id` / `mapping_quality` accessors, plus BAM and CRAM summary paths.

### Wave B — File Ops

- [x] `sort` (`bam_sort.c`, 138k) — in-memory coordinate/`-n` natural-name/`-N` lexicographic-name/`-t TAG` sort for BAM, SAM, and reference-backed CRAM input. `-o`/`-o -` (stdout)/`-O sam|bam|cram`/`--output-fmt=cram`/`--write-index`/`--no-PG`. CRAM output requires a top-level `--reference` and is encoded via a temporary BAM conversion. Emits the **raw input header** (preserving @SQ/@RG field order & @CO) with @HD `SO`/`GO`/`SS` applied + text @PG; name sort uses the exact `bam_sort.c` comparator (`strnum_cmp` natural order + the `flag&0xc0/0x100/0x800` READ1/READ2/supp/sec tiebreak). Tolerates SAM `c/C/s/S/I` scalar aux integer synonyms (htslib-compatible). `--template-coordinate` ports `template_coordinate_key` / `bam1_cmp_template_coordinate` with RG→library lookup, unclipped coordinates, mate CIGAR handling, and `@HD GO:query`. `-M` minimiser sort ports `worker_minhash`, `bam1_cmp_by_minhash`, `build_minhash_index`, `-K`, `-H`, `-R`, and `-I` indexed-reference variants. **Every upstream `test_sort` fixture is byte-exact** (modulo the harness' @PG/VN normalization): coordinate, natural/lexicographic name, aux-tag, minimiser `{basic,indexed,indexed-poly}`, and template-coordinate. Tests `sort_matches_upstream_test_sort_fixtures`, `sort_minimiser_all_variants_match_upstream`, and `sort_cram_input_uses_top_level_reference` cover the CRAM input/output path. **Pending (non-fixture):** on-disk external merge for very large inputs and thread/memory caps.
- [x] `merge` (`bam_sort.c` shared) — full upstream `test_merge` group passes (28/28, including threaded duplicates). In-memory multi-input merge (BAM/SAM) supports coordinate/`-n` natural-name/`-t TAG`/`--template-coordinate` order and CRAM output with top-level `--reference` via temporary BAM conversion. **Upstream `-s SEED` @RG/@PG reconciliation implemented**: seeded `gen_unique_id` (`crate::rand48` glibc LCG) suffixes colliding IDs in header-line order per file; raw merged header (input[0] @HD verbatim for coordinate; SO/GO/SS for name/tag/template; @SQ unioned by SN plus `@SQ AN` alias resolution for SAM-to-BAM and `-L BED` regions; @RG/@PG ID/PG:/PP: remapped; @CO appended); records' `RG:Z:`/`PG:Z:` remapped, delete-then-append ordered like HTSlib, dropped with upstream `bam_translate` warning when unresolved. `-r` attaches a filename-stem @RG to every record; `-c`/`-p` combine identical @RG/@PG IDs (grouped short opts `-cp`/`-rp`). `-f`/`-o`/`-o -`/`-O sam|bam|cram`/`--output-fmt=cram`/`-b`/`-R`/`-L`/`--write-index`/`--no-PG`. Tests include `merge_reconciles_rg_pg_byte_exact_vs_upstream`, `merge_input_list_precedes_remaining_positionals_like_upstream`, `merge_tag_sort_ties_match_upstream_secondary_order`, `merge_template_coordinate_matches_upstream_fixture`, `merge_l_bed_resolves_reference_aliases`, and `merge_writes_cram_output_with_top_level_reference`. **Pending (non-fixture):** k-way streaming merge.
- [x] `collate` / `bamshuf` (`bamshuf.c`) — in-memory grouping for BAM, SAM (tolerant `c/C/s/S/I` aux reader), and reference-backed CRAM. **Non-fast order is the exact upstream `bamshuf` order**: bucket by `hash_X31_Wang(qname) % 64`, then sort each bucket by `(hash, qname, flag>>6&3)` (ported `hash_Wang` bit-mix). Fast `-f` mode mirrors the ring buffer: evict-the-oldest-after-insert so a read whose mate is further than `-r` away is deferred. SAM output emits the **raw input header** with `@HD SO:unsorted GO:query` applied (preserving input `@RG`/`@SQ` field order) + records via `sam_render`. Output format inferred from the `-o` filename extension when `--output-fmt` is absent, including `.cram`. `-o`/`-O`/`-n`/positional prefix/`--output-fmt=cram`/`--no-PG`; bare `collate <input>` now rejects with usage like upstream rather than inventing a default output prefix. CRAM output requires a top-level `--reference` and is encoded via a temporary BAM conversion. **Byte-exact vs the ENTIRE upstream `test_collate` harness (6/6)**; tests `collate_*` + `collate_matches_upstream_test_collate_fixtures` in `sort_merge.rs`, with `collate_cram_input_uses_top_level_reference` covering CRAM input/output and `collate_requires_explicit_output_destination` locking the no-implicit-output argument parity. **Pending:** on-disk hash-bucket for inputs larger than memory.
- [x] `cat` (`bam_cat.c`) — full upstream `test_cat` group passes (26/26): BAM, recompressed BAM, CRAM visible concatenation, CRAM region paths, stdout redirection, `-p 1/2` + `-p 2/2`, and `-h` replacement headers. SAM input is rejected like upstream; BAM paths remain record-level decompress + re-encode with `-o`, `-h`, `-b FILE`, default `@PG`, `--no-PG`, and BAM `-r region`. **Pending (non-fixture):** true CRAM-preserving concatenation and BGZF block-level BAM fast path; current CRAM fixture support writes the SAM-visible stream.
- [~] `split` (`bam_split.c`) — basic BAM/SAM/whole-file-CRAM-by-`@RG` splitting with per-output `@RG` header filtering and default `@PG` insertion; explicit `-d TAG` string/integer aux-tag splitting with on-demand outputs; explicit `-d RG` unknown-read-group header insertion; `-M`/`--max-split`, `-f` template (`%*`, `%!`, `%#`, `%.`) with extension-based output-format inference when no explicit format is given, `-u` unaccounted, `-h` unaccounted SAM header override, `--output-fmt sam|bam|cram`, `--no-PG`, `--write-index` BAI generation for BAM outputs, and `-p N` padding. The full upstream `test_split` group passes (18/18), CRAM input covers the `test_checksum` split/merge path, direct C-vs-Rust smoke locks `.sam` template inference on a stable RG-only BAM, and `split_sam_input_by_rg_to_cram_outputs_with_reference` covers CRAM output with top-level `--reference`. **Pending:** sorted-by-tag streaming mode and deeper upstream `@PG` byte-parity for complex chains.
- [x] `reheader` (`bam_reheader.c`) — full upstream `test_reheader` group passes (7/7): BAM header replacement, CRAM v2.1/v3.0 visible reheader, in-place harness paths, and `-c <command>` external header filtering. **Pending (non-fixture):** true CRAM-preserving in-place/binary rewrite and BAM BGZF block-level fast path; current CRAM fixture support rewrites the SAM-visible stream.
- [~] `addreplacerg` (`bam_addrprg.c`) — SAM/BAM/reference-backed CRAM add/replace `@RG` + `RG:Z`. `-r` now unescapes `\t`/`\n` so a full `@RG\tID:..\tCN:..` spec works; incremental `-r KEY:VAL`, `-R ID` (rejected if absent), default-first-`@RG`-ID, `-m overwrite_all|orphan_only`, `-w` edit, `-O sam|bam|cram` (`cram` requires `-T`/`--reference`), `-o`, `@PG`/`--no-PG`. **Byte-exact vs the whole upstream `test_addrprg` group** (`addrprg/{1,2,4,5}` + `-R` overwrite, modulo `@PG`; `addrprg/3` = expected `-R` failure); integration tests `addreplacerg_matches_upstream_group`, `addreplacerg_writes_cram_output_with_reference`, and `addreplacerg_accepts_reference_backed_cram_input`. **Pending:** mate-aware updates, full orphan-first semantics.
- [x] `fastq` / `fasta` / `bam2fq` (`bam_fastq.c`) — the full upstream `test_bam2fq` group passes (84/84), including threaded duplicates. Single-stream output works for SAM, BAM, and CRAM (records written to stdout, `-o FILE`, or `-0 FILE`), with `-f`/`--require-flags`, `--rf`/`--include-flags`, `-F`/`--exclude-flags`, `-G`, the upstream default `0x900` secondary/supplementary exclusion, read-name suffix controls (`-n`/`-N`), `-O` original-quality `OQ` tag output, `-v INT` missing-quality defaults for FASTQ, `-U`/`--UMI-tag` UMI read-name suffixes, `-i`/`--barcode-tag` CASAVA barcode fields, upstream-style name-grouped paired split outputs (`-1`/`-2`/`-s`/`-0`) that pick the best per-readpart record per qname-group and route R1+R2 to `-1`/`-2`, R1-only or R2-only singletons to `-s` (falling back to `-1`/`-2` when `-s` is absent), READ_OTHER to `-0` (falling back to `-s` when `-0` is absent), paired stdout interleaving with upstream default `/1`/`/2` suffixing when only discard files are specified, per-record interleaved output when `-1` and `-2` paths alias, SAM/BAM/CRAM selected aux comments via `-T` / compact `-Tfoo` in single and split output modes, all-tag comments via `-T ''` / `-T '*'`, `B` array aux comment formatting, FASTQ tag filtering via `-d`/`--tag TAG[:VALUE]` and `-D`/`--tag-file TAG:FILE`, accumulating `-t` (`RG,BC,QT`) and `-T TAG,...`, repeated `-d` / `-D` union semantics with mismatched-tag rejection, FASTA/FASTQ reverse-complement of reverse-strand records, headerless `.sam` index extraction, per-qname-group `--i1`/`--i2` index FASTQ extraction with `--index-format`, `--quality-tag`, `--barcode-tag`, CASAVA paired-end barcode propagation, and `--no-sc` / `--no-sc-bkp` / `--sc-aux` soft-clip trimming with backup aux fields. CRAM input uses explicit top-level `--reference` or discoverable references; generated `view -C` setup roundtrip parity depends on the vendored noodles CRAM mate/MAPQ preservation fix. **Pending (non-fixture):** broader CRAM reference-discovery edge cases and worker-thread propagation beyond accepted `-@`.
- [x] `import` (`bam_import.c`) — basic single FASTA/FASTQ and paired FASTQ (`-1`/`-2`, `--r1`/`--r2`, `-s` interleaved, plus two positional inputs) → SAM/BAM/CRAM (`-O bam` / `--bam`, `.bam` inference, `-O cram` / `--cram` / `--output-fmt=cram`, or `.cram` inference), including positional single input plus `-0` single-read alias, `-0` singleton input alongside paired `-1`/`-2`, positional interleaved FASTQ detection from `/1`/`/2` read names, no-op `--no-PG`, CASAVA parsing (`-i`) with upstream-style reverse comments, SRA name2 (`-N`), UMI extraction (`-U`/`--UMI-tag`) with reverse comments, CASAVA barcode sequence tags (`--barcode-tag`), FASTQ definition aux tags (`-T`) including upstream-style float exponent spelling, explicit index reads (`--i1`/`--i2`) for `-0`, `-s`, positional interleaved, and paired `-1`/`-2` inputs with barcode sequence/quality tags (`--barcode-tag`/`--quality-tag`) and `-b`, and read-group header/tag support (`-R`/`-r`) with repeated `-r` accumulation, `-r` precedence over `-R`, and `-r` ID validation. CRAM output encodes the unmapped import stream through a temporary BAM and empty temporary FASTA. **Byte-exact vs the full upstream `test_import` group (21/21)**, including the roundtrip cases that pipe through `fastq`; Rust regression `import_fastq_to_cram` covers CRAM output.

### Wave C — Editing / Mate-aware

- [x] `fixmate` (`bam_mate.c`) — full upstream `test_fixmate` group passes (42 total: 40 passed, 2 expected failures, including threaded duplicates). Name-grouped BAM/SAM/reference-backed CRAM input mate fixup rejects `@HD SO:coordinate`, recalculates TLEN from mate 5-prime positions, updates default MC/MQ tags, supports `-m` mate-score tags, `-c` lowercase template-CIGAR `ct:Z` tags, `-z`/`--sanitize`, `-r`, `-` stdout output, default `@PG`/`--no-PG`, raw-header SAM output for fixture byte parity, and CRAM output with a top-level `--reference` (`-O cram` / `--output-fmt=cram`) via a temporary BAM conversion. Aux updates are order-preserving (`bam_aux_del` + append-to-tail semantics, MQ before MC, `MC:Z:*` when either read is mapped). `-M` base-modification parity covers `MM`/`ML`/`MN`, draft `Mm`/`Ml` normalization, hard-clipped secondary/supplementary trimming from the primary sequence, invalid ML/MN deletion, and no-sequence cases. Integration test `fixmate_matches_upstream_group` covers the locked fixture set; Rust regression `fixmate_accepts_reference_backed_cram_input` covers CRAM input via top-level `--reference` and CRAM output.
- [x] `markdup` (`bam_markdup.c`, 89k) — **faithful upstream key/score port**. SE/PE duplicate marking for SAM/BAM/reference-backed CRAM input, with CRAM output via `-O cram` / `--output-fmt=cram` and `-T` / `--reference` or top-level `--reference` using the shared temporary-BAM conversion path. PE reads build the upstream `make_pair_key` (template default + `--mode s` sequence; unclipped coords from CIGAR & `MC` tag; `R_LE`/`R_RI` left/right discriminator so a template's two mates get distinct keys and only corresponding mates of duplicate templates collide) plus a shared `make_single_key`; the kept read of a colliding key is the one with the higher `calc_score` = Σ(base qual ≥ 15) + `ms` mate-score tag, with the QCFAIL-asymmetry override and qname `strcmp` tie-break. `-S` seeds a qname `dup_hash` from marked-duplicate reads carrying `SA`/`XA` or an unmapped mate and flags matching supplementary/secondary/unmapped records (gated on `-S`, as upstream). `-b`/`--barcode-tag`, `-c`, `-t` `do` tags, `-d` `dt:Z:SQ|LB` with the full `find_duplicate_chains` optical re-tagging (per-read `original`/`duplicate` chain links + `check_chain_against_original` + `check_duplicate_chain`), `get_coordinates_colons` optical-name parse, `--use-read-groups` (rg-keyed), `--duplicate-count` (`dc:i`), `--include-fails`, `-m`/`--mode t|s`, `-r`, `-s`, `-O`, `-o`, regex `--read-coords`/`--coords-order`/`--barcode-rgx`/`--barcode-name` (via the `regex` crate; capture-span-bounded coord parse), raw-header SAM output (preserves input `@RG`/`@SQ` order), upstream-shaped expect-fail exits for queryname sort, bad coordinate order, missing `MC`, and missing `ms`, `@PG`/`--no-PG`. **Byte-exact vs the ENTIRE upstream `test_markdup` SAM harness — `markdup/{5..18}.expected.sam` (all 14 passing fixtures) plus `1..4` expect-fail partial stdout/stderr cases**; `-s` stats counts match the promoted fixture matrix (`5..18`) including optical/barcode/read-group/duplicate-count cases. Tests `markdup_matches_upstream_test_markdup_fixtures`, `markdup_upstream_expect_fail_cases_return_one_with_expected_partial_output`, `markdup_stats_match_upstream_fixture_counts`, and `markdup_accepts_reference_backed_cram_input_and_output`.
- [~] `rmdup` (`bam_rmdup.c` + `bam_rmdupse.c`) — single-end and paired-end duplicate removal for BAM, SAM, and reference-backed CRAM inputs. SE records are keyed by `(tid, pos, reverse-flag)`; PE records pair by qname and are keyed by the canonical pair of `(tid, pos, strand)` triples, retaining the highest MAPQ/combined MAPQ record or pair. `-s`/`-S` force single-end treatment. `-O sam|bam|cram` / `--output-fmt[=]FMT` selects output format, `.cram` output paths infer CRAM, CRAM output spools kept records through a temp BAM and the shared reference-backed CRAM writer, and local `-T`/`--reference` or top-level `--reference` supplies CRAM input/output references. Default `@PG` insertion via `pg::add_samtools_pg_to_header` and `--no-PG` are supported. Regression `rmdup_accepts_reference_backed_cram_input_and_output` covers CRAM input to CRAM output; `scripts/run-byte-parity-smoke.py` now includes direct C-vs-Rust output-file and stderr diagnostic smoke cases for single-end SAM output and paired SAM output. **Pending:** broader deprecated-command parity for binary output, CRAM, and non-smoke diagnostic/stat output cases.
- [~] `calmd` / `fillmd` (`bam_md.c`) — SAM, BAM, and reference-backed CRAM input can emit SAM text with MD/NM tags recomputed against a FASTA reference via CIGAR/reference walking. The default text path now preserves unchanged `MD`/`NM` aux tags in place and skips unmapped records like upstream; direct C-vs-Rust smoke locks `calmd --no-PG mpileup.1.sam mpileup.ref.fa` stdout/stderr/exit byte parity. `-e` changes matching bases to `=` for both mapped records and unmapped records that still carry a reference/CIGAR, while still leaving unmapped records' MD/NM aux fields untouched. BAQ paths (`-r`, `-E`, `-A`) are wired through `htslib_rs::alignment_compat::recalculate_baq_*` for SAM input directly, and for BAM/CRAM input via a temporary SAM bridge. `-d` drops all aux tags except `RG`, `-q` bins base qualities, `-N` suppresses MD/NM aux updates and diagnostics, `-C cap` caps MAPQ through the upstream `sam_cap_mapq` port when `cap > 10`, and `-n max_nm` masks matching bases to `N` plus zeroes their qualities when the recomputed NM reaches the threshold; changed MD/NM tags emit upstream-shaped `[bam_fillmd1]` diagnostics unless `-Q` quiet mode is set. `-b`/`-u` emit BAM output, while `.cram`, `-O cram`, and `--output-fmt=cram` encode the recalculated SAM stream as CRAM with the supplied reference. Default `@PG` insertion via `pg::add_samtools_pg` (text-level) and `--no-PG` are supported. Tests include `calmd_dash_u_a_r_emits_bgzf_bam_like_upstream`, `calmd_writes_cram_output_with_reference`, `calmd_cap_mapping_quality_uses_sam_cap_mapq`, `calmd_max_nm_masks_matching_bases_and_qualities`, `calmd_dash_e_changes_matching_bases_for_mapped_and_unmapped_records`, `calmd_dash_q_bins_base_qualities`, `calmd_dash_cap_n_preserves_existing_md_nm_tags`, `calmd_dash_d_keeps_only_rg_aux_tag`, and `calmd_baq_accepts_bam_and_cram_input`. **Pending:** broader upstream MD/BAQ option byte parity.
- [x] `targetcut` (`cut_target.c`) — fosmid pool target cutting. Faithful port of the upstream pileup consensus + revised MAQ error-model scoring + two-state target HMM; supports `-Q`, `-i`, `-0`, `-1`, `-2`, `-f`/`--reference`, and `-o`. Upstream ships no fixture group, so coverage is focused Rust tests: unit tests for long supported interval emission, min-baseQ filtering, and attached scoring-option parsing, plus `tests/targetcut.rs` CLI integration coverage for dispatch, `-o`, quality filtering, and error exits.
- [x] `reset` (`reset.c`) — strip alignment fields (`tid`/`pos`/`cigar`/`mate_*`/`template_length`) for BAM, SAM, and CRAM inputs, set MAPQ to `0`, drop a default set of aligner aux tags (NM, MD, AS, XS, SA, MC, MQ, NH, HI, ms), clear `PROPER_PAIR`/`SECONDARY`/`SUPPLEMENTARY`/`REVERSE`/`MATE_REVERSE`, set `UNMAPPED`, set `MATE_UNMAPPED` for paired reads, reverse-restore reverse-strand sequence/quality, preserve duplicate flags with `--dupflag`, remove read-group headers/tags with `--no-RG`, remove program header chains with `--reject-PG`, add a new samtools `@PG` chain entry by default (via the shared `pg::add_samtools_pg_to_header` helper), suppress the new `@PG` with `--no-PG` while preserving existing entries (matching upstream's `noPGentry` semantics), accept SAM/BAM input from stdin/no positional input/`-`, and tolerate legacy SAM `@HD VN:1` headers. `-x`/`--remove-tag`, `--remove-tag ^...`, and `--keep-tag` now match upstream precedence, including unioning multiple keep sets and letting `--no-RG` take precedence over keeping `RG`. **Order-preserving aux drop** (input field order, HTSlib `bam_aux_del` semantics) + **fresh reset output header** (keep `@HD VN:1.6`, drop `@SQ`/`@CO`, `@RG` verbatim unless `--no-RG`, `--reject-PG` removes the named `@PG` + subsequent program headers) + format inferred from the `-o` extension (`sam_open_mode`: SAM unless `.bam`/`.cram`). `-O cram`, `--output-fmt=cram`, and `.cram` output paths encode the now-unmapped reset stream through a temporary BAM and empty temporary FASTA. `--reject-PG` uses the upstream positional rule (`reset.c:223`: keep `@PG` until the first matching `ID`, drop it and all subsequent `@PG`). CRAM input decodes with an explicit top-level `--reference`, embedded reference, or adjacent FASTA discovery for the upstream no-`-T` fixture. **Byte-exact vs the full upstream `test_reset` group (18/18): `basic.1.mp.1` (reset\|view, stdin, file), `basic.output.mp.1` (`-o` SAM from stdin), `basic.bam.input`, `basic.cram.input`, `output.nRG.*`, `output.keep.*`, `output.flg.*`, `reject.1`, `reject.2`** (harness `hskip=1` + `ignore_pg_header` where applicable); test `reset_matches_upstream_test_reset_fixtures`. Rust regression `reset_writes_cram_output` covers CRAM output. **Pending:** broader CRAM reference-discovery parity.
- [~] `ampliconclip` (`bam_ampliconclip.c`, 40k) — **faithful port**. Per-reference BED primer sites (sorted by `right`), `matching_clip_site` (binary-search + `--tolerance`/`--strand` overlap pick), `bam_trim_left`/`bam_trim_right` soft/hard clip (CIGAR/POS/SEQ/QUAL rewrite, hardclip merge, full-consume→empty), `active_query_len`-gated `--filter-len`/`--fail-len`/`--unmap-len`, `--both-ends`, `--original` (`OA` tag), `--keep-tag` (default deletes `NM`/`MD`, order-preserving), `--clipped`, `--no-excluded`, `--rejects-file`, `--primer-counts` TSV, `-f` stats, `-o`/`-O sam|bam|cram`, `-b`, raw-header `@HD SO:coordinate→unknown`, default `@PG`/`--no-PG`. SAM/BAM/reference-backed CRAM input works, and CRAM output uses `-T` / `--reference` or top-level `--reference` via the shared temporary-BAM conversion path. **Byte-exact vs the entire upstream `test_ampliconclip` harness** (10 SAM fixtures + 3 primer-counts TSVs), plus the dormant current-upstream `3_multi_ref_both_clip` `--both-ends` multi-reference edge; tests `ampliconclip_matches_upstream_test_ampliconclip_fixtures` and `ampliconclip_accepts_reference_backed_cram_input_and_output`. **Pending:** BGZF block fast path.

### Wave D — Stats / Pileup

- [x] `depth` (`bam2depth.c`) — per-position depth via sparse CIGAR walks for SAM, BAM, and indexed CRAM. `-a`/`-aa`/`-d`/`-q`/`-o`, `-H` header output, `-f` input file lists, flag filters (`-g`, `-G`/`--excl-flags`, `--incl-flags`, `--require-flags`), `-l` minimum read length filtering, `-r` region restriction, `-b` BED restriction, `-J` deletion-depth counting, `-s` paired-overlap removal, and multi-input columnar output are supported. CRAM can decode with an explicit top-level reference or, for this CIGAR-only metric path, through the htslib-rs synthetic-reference indexed query. BED order is preserved for depth output while `-a`/`-aa` match upstream's empty-region/reference-coverage behavior; deletion-only sites are emitted as zero-depth touches when `-J` is not set. **Byte-exact** vs upstream `test/mpileup:depth.reg` (55 expected passes, 1 expected failure), including the companion `mpileup -ABQ0` depth-simulation commands, plus the large-reference fixtures `depth.expected.out` and `depth_bed.expected.out`. Tests cover SAM/CRAM region depth, no-reference CRAM region depth matching the reference-backed path, multi-input columns/header/list files, flag filters, read-length filtering, and large-reference depth + BED. **Pending (non-fixture):** broader pileup-edge parity.
- [~] `coverage` (`coverage.c`) — per-reference/`-r` region `numreads`, `covbases`, `coverage`, `meandepth`, `meanbaseq`, and `meanmapq` via CIGAR walks for SAM, BAM, and indexed CRAM. CRAM can decode with an explicit top-level reference or, for the current CIGAR/quality metric path, through the htslib-rs synthetic-reference indexed query. `--min-depth` thresholds covered-base counts, `-Q`/`--min-BQ` filters low-quality bases, `-q`/`--min-MQ` filters reads, `-b`/`--bam-list` expands input filename lists, `--ff`/`--excl-flags` replaces the default filter-out flags, `--rf`/`--incl-flags` requires at least one selected flag, `-l`/`--min-read-len` filters short alignments, `-d` caps per-position depth for reported coverage/depth metrics, and multiple inputs aggregate into one row per reference/region. `-m`/`--histogram`, `-A`/`--ascii`, and `-D`/`--plot-depth` emit upstream-shaped 10-row histogram/depth-plot output with UTF-8 or ASCII glyphs, sidebars, x-axis labels, and `-w`/`--n-bins` column control; `coverage_depth_plot_uses_depth_bins_not_breadth_histogram` locks `-D` as mean-depth bins rather than the breadth histogram. **Byte-exact vs the full upstream `test_coverage` group (6/6)**, now promoted into the required parity subset: `coverage/{1..5}.expected` plus multi-input and `-Q`/`-q`; C `printf %g`/`%.3g` formatting (`c_printf_g`, `coverage.c:211`), `min_depth`-gated `meandepth`/`meanbaseq` accumulators (per-position baseq vecs), missing quality scores contributing HTSlib's `0xff` base-quality sentinel, and pileup-arrival reference row ordering. Direct C-vs-Rust smoke now includes `coverage` over `test_input_1_a.bam` without harness normalization. Tests cover no-reference CRAM region coverage matching the reference-backed path. **Pending:** broader byte-parity for non-fixture histogram edge cases such as terminal-width default bin sizing and uneven bin tails.
- [x] `bedcov` (`bedcov.c`) — total aligned-base coverage per BED region, walking each record's CIGAR for SAM, BAM, and indexed CRAM. CRAM can decode with an explicit top-level reference or, for this CIGAR-only metric path, through the htslib-rs synthetic-reference indexed query. `-Q` mapq filter, `-g`/`-G` filter-mask controls, `-j` deletion/refskip skipping, `-H` output headers, `-c` read-count columns, and `-d` depth-threshold columns are supported. **Byte-exact vs the full upstream `test_bedcov` group** (8/8): `bedcov/bedcov{,_j,_gG,_c}.expected`, attached `-g512 -G2048`, and all `-H` header cases including empty source header fields and BED12-derived placeholder columns. Direct C-vs-Rust smoke locks default, `-j`, `-g512 -G2048`, and `-c -H` output. Tests cover no-reference CRAM bedcov success.
- [x] `stats` (`stats.c` + `stats_isize.c`, 123k + 8k) — `SN` (Summary Numbers) section plus FFQ/LFQ first/last fragment quality histograms, GCF/GCL first/last fragment GC histograms, and approximate CIGAR-walk COV coverage histograms for SAM, BAM, and reference-backed CRAM region paths, including record-backed `-I`/`--id` read-group/sample filtering, `-f`/`--required-flag`, `-F`/`--filtering-flag`, `-i`/`--insert-size` insert-size cap, `-m`/`--most-inserts` insert-size bulk selection, `-l`/`--read-length`, `-q`/`--trim-quality` BWA trim counting, `-c`/`--coverage MIN,MAX,STEP` COV binning, and `-g`/`--cov-threshold` target percentage lines with target-region validation. SAM and BAM iterate records directly to populate sequence-length, quality, GC, CIGAR, NM, COV, no-reference and reference-backed GC-depth bins, and runtime coordinate-order accumulators. The emitted lines cover: raw total / filtered / sequences / runtime is sorted / 1st & last fragments / mapped / mapped+paired / unmapped / properly paired / paired / duplicated / MQ0 / QC-failed / non-primary / supplementary / total length / total first fragment length / total last fragment length / bases mapped / bases mapped (cigar) / mismatches (NM aux) / error rate / average length & per-fragment / maximum length & per-fragment / bases trimmed / average quality / singletons / insert size mean & stddev / inward, outward, and other oriented pair counts / pairs on different chromosomes / percentage of properly paired reads / target bases / target genome coverage above threshold. SAM, indexed BAM, and reference-backed CRAM positional region arguments and `-t` target files restrict the summary and COV positions, with overlapping BAM/CRAM regions de-duplicated. `--ref-stats` RFS output now matches the upstream target-file and positional-region fixture matrix, including the `-t targets ref1` case where RFS uses target-file intervals. `-d` / `--remove-dups` filters duplicate-marked primary records and their quality/GC/COV histogram contributions. Missing CRAM references fail cleanly. READ_OTHER records now match upstream for average-quality accumulation, the properly-paired percentage denominator, no-reference GCD zero-GC contribution behavior, and reference-backed MPC/GCD handling of ambiguous reference bases. **Byte-exact vs the full upstream `test_stats` group (42 total: 38 passing, 4 expected failures)**; tests include `stats_matches_upstream_stat_fixtures`, `stats_no_reference_gcd_matches_upstream_multibin_shape`, and `stats_filters_required_and_filtering_flags`. Direct C-vs-Rust smoke now locks raw default BAM output plus representative `-d`, `-f`, `-F`, `-i`, `-m`, `-l`, `-c`, `-g`, `-t`, `-q`, reference-backed, `--ref-stats`, reference-backed `--ref-stats`, and positional-region paths, including upstream-compatible produced-by and command-line headers, without harness normalization. **Pending (not fixture-covered):** broader CRAM no-region per-cycle/quality/COV parity and remaining pileup COV/GCD edge cases.
- [~] `mpileup` (`bam_plcmd.c`, 49k) — default **text** pileup implemented on the `htslib-rs::alignment_compat` pileup iterator (`pileup_from_alignment_paths_with[_reference][_and_options]` + `PileupColumn`/`PileupRead`). Supports multi-input + `-b` list (incl. `file://`), `-f` reference (plain or bgzipped, for the ref-base column, BAQ, and CRAM decode), `-r region` (incl. attached `-r17:1-2` form), `-l` BED positions, `-a`/`-aa` zero-depth text columns, glued short clusters such as `-ABQ0`, `-Q`/`--min-BQ` (default 13), `-q`/`--min-MQ`, `--ff`/`--excl-flags` (mask **replace**, default `0x704`), `--rf`/`--incl-flags`, `-A`/`--count-orphans`, `-B`/`--no-BAQ`, `-x`/`--ignore-overlaps`, `-o`. Faithful `pileup_seq` byte encoding (`.`/`,`, upper/lower mismatch, `^`+mapq head, `$` tail, `*` deletion, `<`/`>` ref-skip, `+`/`-` indels), default BAQ-adjusted qualities via htslib-rs `probaln`, HTSlib `MPLP_SMART_OVERLAPS` overlap removal and `MPLP_NO_ORPHAN` orphan filter, and the `[mpileup] N samples in M input files` stderr line. The pileup iterator tolerates SAM records with `SEQ=*` for upstream depth-simulation fixtures while still rejecting truncated non-empty sequences. **Byte-exact vs the full upstream `test_mpileup` group (7/7)** and the `test/mpileup:depth.reg` mpileup-simulation lines; direct C-vs-Rust smoke locks both stdout and `-o FILE` for a stable `-B --ff 0x14` region. Integration tests `mpileup_default_baq_matches_upstream_out1`, `mpileup_minus_b_ff_matches_upstream_out3`, `mpileup_overlap_removal_matches_upstream_out5`, plus htslib-rs BAQ/BGZF-reference tests. **Pending (non-fixture):** `@RG`-`SM` sample grouping (currently one sample per file), VCF/BCF (`-g`/`-v`) output, base-modification columns, per-position mods/qpos/qname extra columns, CRAM-without-index-via-region efficiency.
- [x] `consensus` (`bam_consensus.c`, 126k, + `consensus_pileup.c`) — **byte-exact vs ALL 77 `test/consensus/consensus.reg` cases**: both `--mode simple` (freq/score) and the default `bayesian`/`recall` Gap5 model (`calculate_consensus_gap5` + the `consensus_init` probability tables / `fast_exp`/`fast_log2`/`ph_log` math, fed by the htslib-rs pileup `nm_init` precompute `PileupRead::bayes_poly`/`bayes_nm_local`), fasta/fastq/pileup, `-a`/`-aa`, `-r`, `-T`/`--ref-qual`, `--min-MQ`/`--min-BQ`, show-del/ins, glued short options. In-process harness test `consensus_matches_upstream_consensus_reg` (77/77).
- [x] `phase` (`phase.c`) — heterozygote phasing. Faithful port of the upstream pileup-driven heterozygote discovery, revised MAQ error-model genotype likelihoods, local haplotype dynamic program, fragment phasing/ambiguity masking, optional chimera fixing, site-list controls (`-l`/`-e`), and split-BAM output prefix (`-b`, writing `*.0.bam`, `*.1.bam`, `*.chimera.bam` with `ZP:A:Y` tags for confidently phased reads). Supports `-Q`/`--min-BQ`, `-q`, `-k`, `-D`, `-F`, `-A`, `--no-PG`, and reference-backed CRAM via `-f`/`--reference`. Upstream ships no `test_phase` fixture group, so coverage is focused Rust tests: unit tests for phase-set marker/evidence emission, min-baseQ filtering, and split-BAM creation, plus `tests/phase.rs` CLI integration coverage for dispatch, split-BAM output, and error exits.
- [x] `depad` / `pad2unpad` (`padding.c`) — SAM, BAM, and reference-backed CRAM input with `-T` padded FASTA reference converts padded reference columns to unpadded coordinates/CIGAR (`I`/`P`). SAM output (`-s`), BAM output (default, `-u`, `-1`, or `.bam` output path), and CRAM output (`-O cram`, `--output-fmt=cram`, or `.cram` output path) are supported; BAM paths reuse the depadded SAM stream and CRAM paths encode through a derived unpadded reference FASTA. **Byte-exact vs the full upstream `test_depad` group (9/9)** against the `depad.001` fixture with `--no-PG`, covering SAM/BAM input and SAM/BAM output modes. Non-fixture Rust coverage adds `depad_cram_input_and_output_roundtrip`.
- [x] `ampliconstats` (`amplicon_stats.c`, 1776 LOC) — **faithful port**. Per-ref BED (file order, strand), `count_amplicon`/`bed2amplicon`, ±`pos-margin` position→amplicon lookup, `accumulate_stats` (flag filter, qname read-pair overlap removal, `depth_all`/`depth_valid`, `nreads`/`nbases`/`coverage`, `amp_dist` via TLEN±`tlen-adjust`, `tcoord` freq map), `append_lstats` (sum+sum² s.d.), full `dump_stats` (`SS`/`AMPLICON`/`F·`/`C·` incl. `depth_bin` RLE, `TCOORD` ≥ `tcoord-min-count`, `--tcoord-bin` aggregation, `FAMP`, `COMBINED` MEAN/STDDEV). `-S`/`-s`/`--use-sample-name`/`-c`/`-t`/`-d d1,d2,d3`/`-m`/`-D`/`-b`/`-a`/`-l`/`-f`/`-F`/`-o`; SAM, BAM, and reference-backed CRAM input are supported. **Byte-exact vs the entire upstream `test_ampliconstats` harness** (`stats`, `stats_mixed`, `stats_partial`, modulo the harness-stripped version/command-line lines); direct C-vs-Rust smoke now locks all three output-file shapes with those metadata lines stripped. Tests `ampliconstats_matches_upstream_test_ampliconstats_fixtures`, `ampliconstats_use_sample_name_uses_first_read_group_sample`, `ampliconstats_accepts_bam_and_reference_backed_cram_input`, and `aggregate_tcoord_merges_nearby_same_status_into_most_frequent_site`.
- [x] `cram-size` (`cram_size.c`) — **byte-exact vs the entire upstream `test/cram_size/cram_size.reg`** (`normal.out`, `verbose.out`, `encodings.out`). Faithful `cram_expand_method`/`comp_method2expanded` method decoder + verbatim tables, `Container::blocks()` walk, `cram_cid2ds` content_id→data-series map, normal (aggregate-by-cid) / `-v` (by cid+method) reports with the exact `BLOCK …%6.2f%% %-Ns` formatting + ratio/`>999%`/summary, and `-e` `cram_describe_encodings`/`cram_codec_describe` with htslib's exact DS + `tag_encoding_map` ordering. Built on the vendored-noodles CRAM inventory surface (`CompressionHeader` encodings public, `Container::blocks()`, ordered `TagEncodings`). Test `cram_size_matches_upstream_cram_size_reg` (all 3); direct C-vs-Rust smoke now also locks normal, verbose, and encodings text output.
- [x] `checksum` (`bam_checksum.c`, 47k) — the full upstream `test_checksum` group passes (14/14), including threaded duplicates. Default order-agnostic checksum output for SAM/BAM/CRAM input is implemented, including `-o`, `-f`/`-F`/`-b`, `-c`, `-N`, `-q`, `-v`, `-T`, `-O`, `-P`, `-C`, `-M`, `-B`, `-a` field-selection shorthand with upstream-style sanitizer defaults, `-z`/`--sanitize` record mutation, `-m` checksum-output merging for default/position/CIGAR/mate-column reports, selected and wildcard/exclusion aux-tag hashing for scalar/string/array tags with canonical integer encoding, read-group grouping, CRAM whole-file input via explicit-reference and no-reference htslib-rs all-record iterator paths, and split/merge checksum reports for BAM and CRAM read-group partitions.
- [x] `samples` (`bam_samples.c`) — list `@RG SM:` samples across inputs. Header-driven dedup, `-T TAG`, `-o`, `-h`, `-i` index-presence column, `-f`/`-F` FASTA dictionary matching, stdin path lists, `-X` custom index pairs, and CRAM headers are implemented.
- [x] `reference` (`reference.c`) — SAM/BAM/CRAM MD-tag reconstruction to FASTA, indexed BAM `-r`, `-o`/`-q`, embedded-reference CRAM read/write support, and `reference -e` embedded extraction are implemented. **Byte-exact vs the entire upstream `test_reference` suite**, including MD path, embed_ref extraction, and both `-r 17:1000-1500` variants; direct C-vs-Rust smoke also locks MD `-o FILE`. CRAM MD mode now discovers `@SQ UR:` references for C-generated embedded-reference CRAMs and emits region FASTA headers without a non-upstream `length:` suffix. Tests `reference_embed_ref_full_test_reference_byte_exact` and `reference_cram_md_path_with_reference_matches_upstream`.
- [x] `flags` (`bam_flags.c`) — explain a numeric BAM flag. Byte-for-byte parity with upstream.

## Phase 4: Test Harness Integration

- [x] **Parity gate setup**: confirmed the pinned upstream harness does not honor `-e samtools=<rust-binary-path>` for most commands because it constructs commands from `$$opts{bin}/samtools` after option parsing. CI stages the Rust binary at the ignored `samtools/samtools` path, runs a required filtered `test.pl` copy for the stable subset maintained in `scripts/run-passing-parity-subset.py`, runs a required `regression.sh` subset for `consensus.reg` and `cram_size.reg` via `scripts/run-passing-regression-subset.py`, and runs the full `cd samtools && perl test/test.pl` harness as a required gate. `test_bgzip` is provided by the external htslib `bgzip` tool from the CI `tabix` package.
- [x] **`@PG` handling**: upstream expected-output comparisons that add a new samtools program header are normalized in the harness rather than by editing expected fixtures. `test.pl` uses `ignore_pg_header => 1` for merge, sort, collate, fixmate, addreplacerg, split, and reset comparisons whose generated `@PG` lines carry version/path-specific `VN` / `CL` values; reheader pipelines strip `VN` before comparison and still use header reordering. The shared Rust `pg` helper emits upstream field order (`ID, PN, PP, VN, CL`) so the harness' `VN` stripping preserves `PP` links.
- [x] **Status ledger**: `docs/test-status.md` tracks the upstream `test.pl` groups as `passing` / `external`, plus legacy `partial` / `not-yet-ported` / `blocked` values for future regressions.
- [x] **Run progressively**: stable groups are enabled in `scripts/run-passing-parity-subset.py`, stable regression files run through `scripts/run-passing-regression-subset.py`, and the full upstream harness is now required.
- [~] **Rust integration tests per subcommand**: under `crates/samtools-rs/tests/<name>.rs`, write at least: happy path, error path, region/`-L`/format-flag variants, threaded variant where applicable. These run on every PR independently of the Perl gate. Current audit: every registered subcommand now has integration-test hits. `phase`, `targetcut`, and `tview` have dedicated CLI coverage in `tests/phase.rs`, `tests/targetcut.rs`, and `tests/tview.rs`; many older commands still keep their integration coverage in shared files (`misc.rs`, `sort_merge.rs`, `stats_wave_d.rs`, `view.rs`) and should be split/deepened into per-command files over time.
- [x] **Compile-side test binaries**: `samtools/test/merge/test_bam_translate`, `test_rtrans_build`, `test_trans_tbl_init`, `samtools/test/split/test_*`, `samtools/test/vcf-miniview.c` — ported to Rust test coverage under the relevant crate. Merge helper coverage lives in `commands::merge` unit tests: RG/PG translation, unknown-tag dropping/once-only warning state, forced filename-derived RG behavior, `@SQ` unioning, input-header preservation, colliding `@RG`/`@PG` translation with `@RG.PG`/`@PG.PP` field remapping, and `-r` filename-derived RG setup; `test_rtrans_build` has no meaningful assertions. Split helper coverage lives in `commands::split` unit tests and integration tests: filename-template expansion including CRAM extension, bad percent-escape rejection, read-group header filtering/synthesis, string/integer split-tag extraction, early parse/error exits, and reference-backed split CRAM outputs that read back successfully. `vcf-miniview.c` is a VCF/BCF helper, covered in `htslib-rs` `variant_io.rs` by BCF-to-VCF viewing plus the `vcf-miniview -f`-style filtered-output test.

## Phase 5: Parity Polishing

- [~] **Diff every `test_<name>` output byte-for-byte** against the C samtools outputs on a known fixture corpus (locally, dev-only). `scripts/run-byte-parity-smoke.py` now compares C samtools vs Rust samtools stdout/stderr/exit status for stable direct-smoke cases (`view` SAM text, BAM `-h` text, count-only output, MAPQ/required-flag/filtering-flag count filters, indexed BAM region text, MAPQ, `flag.proper_pair`, CIGAR-derived (`endpos`/`rlen`/`sclen`/`hclen`), and mate/reference (`mpos`/`pnext`/`tlen`/`rnext`/`mrname`/`refid`/`mrefid`/`ncigar`) expression filters, qname allow/deny filters, read-group filters including upstream's no-`RG` pass-through unless `-n`, library filters, aux-tag value/file filters, rendered BAM/CRAM binary aux-tag rewrite output, and `view --save-counts` JSON including count-only no-reference CRAM MAPQ/flag expression and read-group/library/aux-tag filter counts, plus no-reference CRAM non-region expression-filtered SAM output and rendered BAM output, `sort`, `flagstat` BAM text/JSON/TSV and CRAM no-reference text, `idxstats` BAM/CRAM, raw default `stats` BAM output plus `-d`/`-f`/`-F`/`-i`/`-m`/`-l`/`-c`/`-g`/`-t`/`-q`/reference-backed/`--ref-stats`/reference-backed `--ref-stats`/region, explicit `index -o` BAI output, `checksum` stdout and `-o`, `faidx`/`fqidx` retrieval and `.fai` index output, `faidx`/`fqidx` zero-length, truncated-region, and `--continue` missing-region warnings, `bedcov` default/`-j`/`-g512 -G2048`/`-c -H`, `coverage` default plus `-o`, fixed-width ASCII `-m` histogram, `-D` depth plot, uneven-bin-tail histogram, and `$COLUMNS`-driven default-bin histogram, `depth` stdout and `-o`, `head`, noninteractive `tview -d T -p` over indexed BGZF-SAM, successful `quickcheck`, full-fixture `quickcheck -v` output/diagnostics, `dict` stdout and `-o`, `reference` MD full/region output, MD `-o`, and `reference -e` embedded-CRAM extraction over a C-generated CRAM, `mpileup` BAQ multi-BAM list, `mpileup -B --ff` stdout and `-o`, `mpileup` overlap-removal output, `merge` seeded SAM output, `merge -r` filename-RG output, `merge --template-coordinate`, `samples` stdout and `-o`, `flags`, `addreplacerg`, `reset`, `depad -s`, `split` `.sam` template output inference and missing-RG/no-`-u` error behavior, `collate`, `consensus`, `ampliconclip -o -`, `markdup -O sam`, `fixmate -O sam`, `calmd` default SAM text, `-e` matching-base conversion, `-r -e` BAQ plus matching-base conversion, `-d` tag dropping, `-n` masking diagnostics, `-Q` quiet suppression, `-C` MAPQ capping, `-q` quality binning, and `-N` no-MD/NM update, `fastq` stdout and `-0 FILE` split-routing shape, and `import` paired read-group plus interleaved aux-output SAM) plus `rmdup` single-end/paired SAM output files and stderr diagnostics, and representative stable error cases (unknown command, missing `view`/`index`/`sort`/`flagstat`/`idxstats`/`head`/`quickcheck`/`dict`/`faidx`/`fqidx`/`checksum`/`coverage`/`depth`/`samples`/`addreplacerg`/`reset`/`consensus`/`ampliconclip`/`ampliconstats`/`cat`/`split`/`stats`/`fixmate`/`markdup`/`rmdup`/`calmd`/`fastq`/`fasta`/`bam2fq`/`reheader`/`import`/`merge`/`mpileup`/`reference`/`cram-size`/`phase`/`targetcut`/`tview`/`depad` inputs, `cat`/`reheader` SAM-input rejection, missing `bedcov` BAM input, `collate -O` missing input, mapped SAM without `@SQ`, bad-header `quickcheck`, missing-EOF quickcheck, truncated-CRAM quickcheck) without harness normalization. Validation uses an actual upstream C build, e.g. `make -C samtools HTSDIR=../htslib-rs/htslib samtools` followed by `python3 scripts/run-byte-parity-smoke.py --c-samtools ./samtools/samtools --rust-samtools ./target/debug/samtools`. **Pending:** expand from smoke cases to every promoted `test_<name>` group output and classify any diffs (real bug / acceptable cosmetic / `@PG` only). The former raw default `stats` header miss is now closed by upstream-compatible banner lines, and representative stats option output is covered by strict smoke cases.
- [~] **Threads**: `htslib-rs::bgzf_compat` exposes worker-count variants for automatic BGZF read/write, and `faidx` / `fqidx` now parse local `-@ N`, attached `-@N`, `--threads N`, `--threads=N`, plus top-level global `--threads`, routing BGZF input/output through noodles multithreaded BGZF worker APIs when nonzero. `htslib-rs::index_compat` exposes multithreaded BGZF reader variants for BAM and BGZF-SAM BAI/CSI construction, and `samtools index` routes local/global thread counts into them when nonzero. `htslib-rs::alignment_compat` also exposes multithreaded BGZF BAM-output helpers for region slicing and required-flag extraction, and native `view_region`, `view_bed`, and `extract_unmapped_pairs` use them when `threads` is nonzero. Tests cover local/global threaded FASTA and FASTQ BGZF retrieval output, threaded index BAI construction, htslib-rs worker-count index builders, and threaded native BAM output paths. **Pending:** propagate worker pools through the remaining BAM/CRAM alignment readers/writers and subcommands, then compare parallelism behavior against C samtools.
- [~] **Exit codes**: representative invalid-input exits are now locked in `tests/exit_codes.rs`: missing input paths for `view`/`index`/`sort` return 1, malformed SAM CIGAR in `view` returns 1 instead of silently passing through, mapped SAM records without `@SQ` headers return 1 like upstream, `quickcheck` preserves upstream-style bitmask exits for not-sequence-data / missing EOF / truncated CRAM, and unknown top-level commands return 1. `crates/samtools-rs-cli/tests/exit_codes.rs` locks exact upstream-shaped stderr for missing `view`/`index`/`sort`/`flagstat`/`idxstats`/`head`/`checksum`/`coverage`/`depth`/`samples`/`addreplacerg`/`reset`/`consensus`/`ampliconclip`/`ampliconstats`/`cat`/`split`/`stats`/`fixmate`/`markdup`/`rmdup`/`calmd`/`fastq`/`fasta`/`bam2fq`/`reheader`/`import`/`merge`/`mpileup`/`reference`/`cram-size`/`phase`/`targetcut`/`tview`/`depad` inputs, missing `bedcov` BAM input with exit 2, and `collate -O` missing input, including htslib's `[E::hts_open_format]` line, exact `dict` and `faidx`/`fqidx` build-missing stderr, and mapped SAM records without `@SQ` headers, including htslib's parse/read context. `scripts/run-byte-parity-smoke.py` now directly compares C-vs-Rust exit status and stderr for unknown command, missing `view`/`index`/`sort`/`flagstat`/`idxstats`/`head`/`quickcheck`/`dict`/`faidx`/`fqidx`/`checksum`/`coverage`/`depth`/`samples`/`addreplacerg`/`reset`/`consensus`/`ampliconclip`/`ampliconstats`/`cat`/`split`/`stats`/`fixmate`/`markdup`/`rmdup`/`calmd`/`fastq`/`fasta`/`bam2fq`/`reheader`/`import`/`merge`/`mpileup`/`reference`/`cram-size`/`phase`/`targetcut`/`tview`/`depad` inputs, missing `bedcov` BAM input, `collate -O` missing input, no-`@SQ` `view`, not-sequence-data `quickcheck`, missing-EOF quickcheck, and truncated-CRAM quickcheck cases. **Pending:** broaden direct C-samtools comparison across all subcommands and error classes.
- [~] **Performance triage**: the smoke bench can now compare the Rust in-process path against an external samtools binary via `SAMTOOLS_RS_BENCH_COMPARE=/path/to/samtools`, reporting per-case mean/min ratios for `view`, `sort`, `markdup`, `stats`, `mpileup`, `coverage`, `depth`, and `checksum`. Current one-iteration comparison against a locally rebuilt C samtools (`make -C samtools HTSDIR=../htslib-rs/htslib samtools`; `SAMTOOLS_RS_BENCH_ITERS=1 SAMTOOLS_RS_BENCH_COMPARE=./samtools/samtools cargo bench -p samtools-rs --bench smoke`) shows: `view` 0.11x, `sort` 0.10x, `markdup` 0.05x, `stats` 3.46x, `mpileup` 4.37x, `coverage` 0.06x, `depth` 0.09x, `checksum` 0.16x. **Pending:** broaden to every subcommand on representative datasets and optimize the over-2x `stats` / `mpileup` paths after parity.
- [x] **Bench harness**: stable custom timing harness added as `crates/samtools-rs/benches/smoke.rs` with a Cargo `[[bench]]` target. It times `view`, `sort`, `markdup`, `stats`, `mpileup`, `coverage`, `depth`, and `checksum` against existing upstream fixtures, writes outputs under a temp directory, accepts `SAMTOOLS_RS_BENCH_ITERS=N` (default 3), and optionally compares with `SAMTOOLS_RS_BENCH_COMPARE=/path/to/samtools`. The comparison path handles C-samtools CLI differences (`markdup` positional output and `stats` stdout output) while keeping Rust outputs file-backed. Smoke validation: `SAMTOOLS_RS_BENCH_ITERS=1 cargo bench -p samtools-rs --bench smoke`; comparison-path validation: `SAMTOOLS_RS_BENCH_ITERS=1 SAMTOOLS_RS_BENCH_COMPARE=./samtools/samtools cargo bench -p samtools-rs --bench smoke`.

## Completed Library / Infra Batch — ✅ ALL RESOLVED

All 12 library/infra blockers are complete, byte/fixture-verified, and
included in the current htslib-rs pin. The CRAM internals and
embed_ref gaps were closed with minimal patches to the **owned vendored
noodles fork** (`madhavajay/noodles`, an htslib-rs submodule), then
consumed from htslib-rs and samtools-rs. This section is the rolled-up
record of the deleted planning file.

- [x] **#1 pileup iterator** — `htslib_rs::alignment_compat` exposes `PileupColumn` / `PileupRead` plus `PileupOptions` with flag/mapq filters, smart overlap removal, and orphan filtering. It unlocked byte-exact `consensus` (77/77), `coverage`, `bedcov`, `depth`, `mpileup` core text cases, `ampliconstats`, and the fixtureless `targetcut`/`phase` ports.
- [x] **#2 CRAM all-record iterator + embed_ref read/write** — `query_cram_records_all_from_path[_with_reference]` supports no-region CRAM iteration. `stats` / `checksum` no-region CRAM, `reference` MD path, embedded-reference reads, and `view -O cram,embed_ref=1` writes are wired; the full upstream `test_reference` suite is byte-exact.
- [x] **#3 CRAM container/block/codec inventory** — vendored noodles exposes the needed compression-header/block surfaces and ordered tag encodings. `cram-size`, `cram-size -v`, and `cram-size -e` are byte-exact vs all three `cram_size.reg` fixtures.
- [x] **#4 binary `@PG`** — `view` injects binary `@PG` for every decodable input→binary path: SAM→BAM/CRAM and BAM→BAM/CRAM, including plain/filter/region/sanitizer paths. CRAM-input direct-copy paths remain intentionally direct until the remaining noodles CRAM limitation is removed.
- [x] **#5 CRAM `idxstats` / `flagstat` without explicit reference** — htslib-rs synthesizes an all-`N` reference repository from CRAM header `@SQ` lines so reference-independent flag/tid summaries decode without user `-T`; outputs match BAM equivalents.
- [x] **#6 index BAMs without `@HD SO:coordinate`** — htslib-rs BAI building no longer rejects missing sort-order headers; `samtools index` works on the affected fixtures.
- [x] **#7 aux mutation** — mutable `RecordBuf` aux access/removal/insertion covers all parity consumers: `addreplacerg`, `ampliconclip`, `fixmate`, and `calmd` paths are byte-exact where fixture-backed. Optional true in-place binary-resize primitives are performance-only.
- [x] **#8 threads** — `-@` / `--threads` are accepted through the relevant command paths and produce byte-identical output independent of count. Actual worker-pool performance wiring remains Phase 5/perf work.
- [x] **#9 write-index** — `sort --write-index`, `view --write-index`, and `merge --write-index` produce BAI files byte-identical to a post-pass `samtools index` build. CSI/CRAI auto-write and true inline index emission are perf/streaming follow-ups.
- [x] **#10 region grammar** — `.` and `*` region semantics are implemented for the upstream-tested SAM/count paths (`view_dot_region_means_whole_file`, `view_star_region_selects_unplaced_reads`). Binary `*` output shares the remaining binary-filter limitations.
- [x] **#11 BAQ / `probaln_glocal`** — htslib-rs BAQ/probaln surfaces are verified against realn fixtures; `calmd -uAr` BAM output and `mpileup` default BAQ-adjusted qualities are wired and accepted by their upstream harness paths.
- [x] **#12 large-reference CSI robustness** — htslib-rs auto-sizes CSI depth from the largest reference, fixing `large_chrom.bam ref2` / `ref2:1-541556283` without a noodles patch; output is byte-exact vs `dat/large_chrom.out`.

Resolved noodles-adjacent items:

- [x] **CSI query robustness for very large references/regions** — fixed entirely in htslib-rs; noodles unpatched.
- [x] **SAM aux float formatting (`f:` scalars and `B:f` arrays)** — resolved in samtools-rs via `sam_render` helpers now used by every SAM-text output path; no noodles-side patch required.

## Submodule Pinning

- [x] Pin `samtools/` to a specific upstream release tag once Phase 0 lands (record tag + commit in `README.md` and `version.rs`). Current pin: upstream tag `1.23.1`, commit `6efb9b6da35224cf804921dedecf9fb8f411365d`.
- [x] Pin `htslib-rs/` to a known-green commit when Phase 0 lands. Current pin: `61e6e72f14f251e0849e8fe87a420eff374892af` (merged all 12 completed library/infra blockers, with vendored noodles PR #6 merged; prior: `f61801c`, `3fadffc`, `8372873`, `530b27c`, `ca812dd`, `9cf30b3`, `e25f392`, `5b25622`, `da4d331`, `6bd6fb0`, `88bd29f`).

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
