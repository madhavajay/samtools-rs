# samtools-rs

Pure Rust port of [samtools](https://github.com/samtools/samtools) — the SAM/BAM/CRAM toolkit.

This project mirrors the upstream `samtools` command-line interface and ports
its test suite. Format-level I/O routes through the sibling
[`htslib-rs`](https://github.com/madhavajay/htslib-rs) workspace (a pure-Rust
HTSlib compatibility layer built on [noodles](https://github.com/zaeleus/noodles)).
There is no `cc`, no `bindgen`, and no link to C HTSlib.

See [`TODO.md`](TODO.md) for the full porting plan, scope decisions, and phase
breakdown.

## Scope

- **In scope:** all upstream `samtools` subcommands except those listed below.
- **Out of scope:** `tview` (interactive curses viewer), remote I/O backends
  (`https://`, `s3://`, `ftp://`, `gs://`), the `misc/` utilities and Perl scripts,
  the `lz4/` vendored library, and C ABI exposure.

## Layout

```
samtools-rs/
├── crates/
│   ├── samtools-rs/        # library: shared infra + one module per subcommand
│   └── samtools-rs-cli/    # binary: `samtools`
├── samtools/               # upstream C source + test suite (submodule, reference only)
├── htslib-rs/              # pure-Rust HTSlib compatibility layer (submodule)
└── TODO.md
```

## Building

Clone with submodules:

```sh
git clone --recurse-submodules git@github.com:madhavajay/samtools-rs.git
```

Build and run the Rust gate:

```sh
cargo fmt --all --check
cargo clippy --workspace --all-targets -- -D warnings
cargo test --workspace
```

Run the binary:

```sh
cargo run -p samtools-rs-cli -- --version
cargo run -p samtools-rs-cli -- <subcommand> [args]
```

## Parity testing

The upstream `samtools/test/test.pl` Perl harness is used as the parity gate.
Once a subcommand is implemented, its tests are enabled. Until then, the
harness fails on `not yet implemented` exits — see [`TODO.md`](TODO.md)
Phase 3 and [`docs/test-status.md`](docs/test-status.md) for the rolling
status.

```sh
cargo build --release -p samtools-rs-cli
cp target/release/samtools samtools/samtools
cd samtools
perl test/test.pl
```

The pinned upstream harness constructs most commands from `./samtools`, so the
Rust binary is staged at the ignored `samtools/samtools` path before running it.

`bgzip` and `tabix` **must** be on `PATH` before running the parity or
regression subsets. Upstream `test.pl` *silently skips* (does not fail)
any group whose input generation needs them, so a missing `bgzip` turns a
run into a false green. The subset runners
(`scripts/run-passing-{parity,regression}-subset.py`) now hard-error
unless both tools are found. They are prebuilt in the vendored htslib:

```sh
make -C htslib-rs/htslib tabix bgzip
export PATH="$PWD/htslib-rs/htslib:$PATH"
```

(CI installs the `tabix` apt package, which provides both.)

Local developers may regenerate expected outputs from the C samtools
(requires `autoconf`, a C toolchain, and the `htslib` C library):

```sh
cd samtools && autoreconf -i && ./configure && make
cd test && perl test.pl --redo-outputs
```

CI does not build C samtools.

## Versioning

`samtools --version` reports the upstream samtools version tracked by the
`samtools/` submodule (currently `1.23.1`, see
[`crates/samtools-rs/src/version.rs`](crates/samtools-rs/src/version.rs)). The
`@PG VN:` field emitted into output headers uses the same string so that
upstream's byte-comparison tests pass.

Current submodule pins:

- `samtools/`: upstream tag `1.23.1`, commit `6efb9b6da35224cf804921dedecf9fb8f411365d`.
- `htslib-rs/`: commit `88bd29f5f0d5e87d3f5d28da1f106a4b518e3926`.

## License

MIT, matching upstream samtools and htslib.
