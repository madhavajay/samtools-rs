//! `samtools sort` — sort alignment records.
//!
//! Mirrors `main_sort` in `bam_sort.c`. The upstream implementation is the
//! largest single file in samtools (138k LOC) and supports external k-way
//! merge with temp files, name/coordinate/tag/template-coordinate sort,
//! and many auxiliary flags.
//!
//! This initial Rust port supports **in-memory coordinate, name, or tag sort
//! for BAM/SAM/reference-backed CRAM**, which is sufficient for small/medium inputs. Records are
//! sorted by `(reference_sequence_id, alignment_start)` for coordinate mode,
//! by `qname` for name mode, or by `TAG` with coordinate/name secondary keys
//! for tag mode, then written to the output.
//!
//! Supported flags:
//!  - `-n` — name sort (default is coordinate sort).
//!  - `-t TAG` — sort by auxiliary tag, using coordinate/name as secondary key.
//!  - `-o FILE` — output file (default stdout).
//!  - `-O sam|bam`, `--output-fmt sam|bam` — output format (default: bam).
//!  - `-@`/`--threads`, `-m`/`--max-mem`, `-T`/`--temp` — accepted but ignored.
//!  - `--no-PG` — accepted, silently ignored.
//!  - `--write-index` — write a BAI next to coordinate-sorted BAM output.
//!  - `-M` — minimiser sort: faithful `bam_sort.c` `worker_minhash` + `bam1_cmp_by_minhash` + `build_minhash_index` + `minhash_with_idx[_squash]` port, with `-K` kmer (default 20, clamped 1..=31), `-H` homopolymer squash, `-R` (no reverse-strand minimiser), and `-I FILE` indexed reference. Byte-identical to all three upstream `sort/minimiser-{basic,indexed,indexed-poly}.sam` fixtures.
//!
//! Not yet supported: external merge (large inputs spill to disk),
//! template-coordinate sort, CRAM output.

use std::cmp::Ordering;
use std::ffi::OsString;
use std::fs::File;
use std::io::{self, BufRead, BufReader, Cursor, Write};
use std::path::{Path, PathBuf};
use std::process::ExitCode;

use htslib_rs::bam;
use htslib_rs::bgzf;
use htslib_rs::format::Exact;
use htslib_rs::sam::{self, alignment::RecordBuf};

use crate::diagnostics::{print_error, print_error_errno};
use crate::io as sam_io;
use crate::sam_global::current_global_args;

/// Entry point for `samtools sort`.
pub fn main(args: &[OsString]) -> ExitCode {
    let mut name_sort = false;
    let mut natural_sort = true;
    let mut output: Option<PathBuf> = None;
    let mut output_fmt = OutFmt::Bam;
    let mut input: Option<PathBuf> = None;
    let mut local_write_index = false;
    let mut no_pg = false;
    let mut tag_sort: Option<[u8; 2]> = None;
    // Minimiser (`-M`) sort state. `-K` kmer (default 20, clamped
    // 1..=31), `-H` enables homopolymer squash (default off / no_squash),
    // `-R` disables the reverse-strand minimiser, `-I` indexed reference.
    let mut minhash_mode = false;
    let mut minhash_kmer: i32 = 20;
    let mut minhash_try_rev = true;
    let mut minhash_no_squash = true;
    let mut minhash_indexed: Option<PathBuf> = None;

    let mut iter = args.iter().skip(1).peekable();
    while let Some(arg) = iter.next() {
        let s = arg.to_str().unwrap_or("");
        match s {
            "-n" | "--name" => {
                name_sort = true;
            }
            // `bam_sort.c` `-N`: name sort with lexicographical
            // (byte `strcmp`) collation instead of natural order.
            "-N" => {
                name_sort = true;
                natural_sort = false;
            }
            "-o" | "--output" => {
                output = match iter.next().and_then(|a| a.to_str()) {
                    // `-o -` means stdout (output stays None).
                    Some("-") | None => None,
                    Some(p) => Some(PathBuf::from(p)),
                };
            }
            "-t" => {
                let Some(v) = iter.next().and_then(|a| a.to_str()) else {
                    print_error("sort", "missing value for -t");
                    return ExitCode::from(1);
                };
                tag_sort = match parse_tag(v) {
                    Ok(tag) => Some(tag),
                    Err(e) => {
                        print_error("sort", e);
                        return ExitCode::from(1);
                    }
                };
            }
            "-O" | "--output-fmt" => {
                let Some(v) = iter.next().and_then(|a| a.to_str()) else {
                    print_error("sort", format!("missing value for {}", s));
                    return ExitCode::from(1);
                };
                output_fmt = match parse_output_format(v) {
                    Ok(fmt) => fmt,
                    Err(e) => {
                        print_error("sort", e);
                        return ExitCode::from(1);
                    }
                };
            }
            "-@"
            | "--threads"
            | "-m"
            | "--max-mem"
            | "-T"
            | "--temp"
            | "-l"
            | "--compression-level" => {
                let _ = iter.next();
            }
            "-M" => minhash_mode = true,
            "-H" => minhash_no_squash = false,
            "-R" => minhash_try_rev = false,
            "-K" => {
                minhash_kmer = iter
                    .next()
                    .and_then(|a| a.to_str())
                    .and_then(|v| v.parse().ok())
                    .unwrap_or(minhash_kmer);
            }
            "-I" => {
                minhash_indexed = iter.next().map(PathBuf::from);
            }
            // Glued minimiser clusters: `-MH`, `-MR`, `-MHR` (getopt
            // string `…MI:K:uRw:H`; only H/R follow M without a value).
            _ if s.starts_with("-M")
                && s.len() > 2
                && !s.starts_with("--")
                && s[2..].chars().all(|c| c == 'H' || c == 'R') =>
            {
                minhash_mode = true;
                for c in s[2..].chars() {
                    match c {
                        'H' => minhash_no_squash = false,
                        'R' => minhash_try_rev = false,
                        _ => {}
                    }
                }
            }
            // Attached `-K10` / `-Iref.fa`.
            _ if s.starts_with("-K") && s.len() > 2 && !s.starts_with("--") => {
                minhash_kmer = s[2..].parse().unwrap_or(minhash_kmer);
            }
            _ if s.starts_with("-I") && s.len() > 2 && !s.starts_with("--") => {
                minhash_indexed = Some(PathBuf::from(&s[2..]));
            }
            "--write-index" => {
                local_write_index = true;
            }
            "--no-PG" => {
                no_pg = true;
            }
            "-u" => {
                // Accepted but currently ignored (controls uncompressed output).
            }
            "--help" => {
                let _ = print_usage();
                return ExitCode::SUCCESS;
            }
            // Attached-value forms of the accepted-but-ignored options
            // (`-@4`, `-m768M`, `-l6`, `-Kprefix`, `-Tprefix`, `--threads=4`).
            _ if (s.starts_with("-@")
                || s.starts_with("-m")
                || s.starts_with("-l")
                || s.starts_with("-K")
                || s.starts_with("-T"))
                && s.len() > 2
                && !s.starts_with("--") =>
            {
                // value is in the same token; nothing to consume.
            }
            _ if s.starts_with("--threads=")
                || s.starts_with("--max-mem=")
                || s.starts_with("--compression-level=")
                || s.starts_with("--temp=") =>
            {
                // value embedded; ignored.
            }
            _ if s.starts_with('-') && s != "-" => {
                print_error(
                    "sort",
                    format!("option `{}` is not yet supported in samtools-rs sort", s),
                );
                return ExitCode::from(1);
            }
            _ => {
                if input.is_none() {
                    input = Some(PathBuf::from(arg));
                }
            }
        }
    }

    let Some(input) = input else {
        let _ = print_usage();
        return ExitCode::from(1);
    };

    let format = match sam_io::sam_open_format(&input) {
        Ok(f) => f,
        Err(e) => {
            print_error("sort", e.to_string());
            return ExitCode::from(1);
        }
    };
    if !matches!(format.exact, Exact::Sam | Exact::Bam | Exact::Cram) {
        print_error(
            "sort",
            "only SAM, BAM, and reference-backed CRAM input are currently supported",
        );
        return ExitCode::from(1);
    }

    let write_index = local_write_index || current_global_args().write_index;
    if write_index {
        if output.is_none() {
            print_error("sort", "--write-index requires -o FILE");
            return ExitCode::from(1);
        }
        if name_sort {
            print_error("sort", "--write-index requires coordinate sort output");
            return ExitCode::from(1);
        }
        if !matches!(output_fmt, OutFmt::Bam) {
            print_error("sort", "--write-index is only supported for BAM output");
            return ExitCode::from(1);
        }
    }

    let minhash = if minhash_mode {
        Some(MinhashOpts {
            kmer: minhash_kmer.clamp(1, 31),
            try_rev: minhash_try_rev,
            no_squash: minhash_no_squash,
            indexed: minhash_indexed,
        })
    } else {
        None
    };

    match run_sort(
        &input,
        output.as_deref(),
        name_sort,
        tag_sort,
        output_fmt,
        write_index,
        if no_pg { None } else { Some(args) },
        minhash,
        natural_sort,
    ) {
        Ok(()) => ExitCode::SUCCESS,
        Err(e) => {
            print_error_errno("sort", "sort failed", &e);
            ExitCode::from(1)
        }
    }
}

fn parse_output_format(raw: &str) -> Result<OutFmt, String> {
    match raw.to_ascii_lowercase().as_str() {
        "sam" => Ok(OutFmt::Sam),
        "bam" => Ok(OutFmt::Bam),
        _ => Err(format!("unsupported output format \"{}\"", raw)),
    }
}

fn parse_tag(raw: &str) -> Result<[u8; 2], String> {
    let bytes = raw.as_bytes();
    if bytes.len() == 2 {
        Ok([bytes[0], bytes[1]])
    } else {
        Err(format!("sort tag must be exactly two bytes, got {:?}", raw))
    }
}

#[derive(Clone, Copy)]
pub(crate) enum OutFmt {
    Sam,
    Bam,
}

/// `samtools sort -M` minimiser parameters.
#[derive(Clone)]
pub(crate) struct MinhashOpts {
    pub kmer: i32,
    pub try_rev: bool,
    pub no_squash: bool,
    /// `-I FILE` indexed reference (built into a kmer→refpos map).
    pub indexed: Option<PathBuf>,
}

#[allow(clippy::too_many_arguments)]
pub(crate) fn run_sort(
    input: &Path,
    output: Option<&Path>,
    name_sort: bool,
    tag_sort: Option<[u8; 2]>,
    fmt: OutFmt,
    write_index: bool,
    pg_argv: Option<&[OsString]>,
    minhash: Option<MinhashOpts>,
    natural_sort: bool,
) -> io::Result<()> {
    let format = sam_io::sam_open_format(input)?;
    let (mut header, mut records) = match format.exact {
        Exact::Sam => read_sam_records(input)?,
        Exact::Bam => read_bam_records(input)?,
        Exact::Cram => read_cram_records(input)?,
        _ => {
            return Err(io::Error::new(
                io::ErrorKind::InvalidInput,
                "only SAM, BAM, and reference-backed CRAM input are currently supported",
            ));
        }
    };

    let mut minhash_mapped = false;
    if let Some(opts) = minhash.as_ref() {
        // Faithful port of `bam_sort.c` `worker_minhash` +
        // `bam1_cmp_by_minhash`. Per unmapped record compute the
        // minimiser hash over its sequence (against the `-I` index when
        // given), reverse-complement it when the reverse-strand
        // minimiser wins, then sort unmapped records by that 62-bit key
        // (mapped records keep coordinate order).
        minhash_mapped = records.iter().any(|r| r.reference_sequence_id().is_some());
        let kmer_h = match opts.indexed.as_deref() {
            // `-w` defaults to 100 in `bam_sort.c`.
            Some(p) => Some(build_minhash_index(p, opts.kmer, 100, opts.no_squash)?),
            None => None,
        };
        let mut keyed: Vec<(MinhashKey, RecordBuf)> = records
            .drain(..)
            .map(|mut r| {
                let key = minhash_prepare(&mut r, opts, kmer_h.as_ref());
                (key, r)
            })
            .collect();
        keyed.sort_by(|a, b| minhash_cmp(&a.0, &b.0));
        records = keyed.into_iter().map(|(_, r)| r).collect();
    } else if let Some(tag) = tag_sort {
        records.sort_by(|a, b| compare_by_tag(a, b, tag, name_sort, natural_sort));
    } else if name_sort {
        records.sort_by(|a, b| name_cmp(a, b, natural_sort));
    } else {
        // Coordinate sort: by (reference_sequence_id, alignment_start).
        records.sort_by(|a, b| {
            // Records with no reference (unmapped) sort to the end.
            coordinate_key(a).cmp(&coordinate_key(b))
        });
    }

    // Sort-order tags for @HD.
    let (so, ss): (String, Option<String>) = if minhash.is_some() {
        if minhash_mapped {
            (
                "coordinate".to_string(),
                Some("coordinate:minhash".to_string()),
            )
        } else {
            ("unsorted".to_string(), Some("unsorted:minhash".to_string()))
        }
    } else if let Some(tag) = tag_sort {
        (
            "unsorted".to_string(),
            Some(format!(
                "unsorted:{}{}:{}",
                tag[0] as char,
                tag[1] as char,
                if name_sort {
                    if natural_sort {
                        "queryname:natural"
                    } else {
                        "queryname:lexicographical"
                    }
                } else {
                    "coordinate"
                }
            )),
        )
    } else if name_sort {
        (
            "queryname".to_string(),
            Some(
                if natural_sort {
                    "queryname:natural"
                } else {
                    "queryname:lexicographical"
                }
                .to_string(),
            ),
        )
    } else {
        ("coordinate".to_string(), None)
    };
    set_sort_order(&mut header, &so, ss.as_deref());

    // Emit the *raw* input header (preserving @SQ/@RG field order, @CO,
    // etc. — noodles' canonical writer reorders @RG fields) with the @HD
    // SO/SS applied and the samtools @PG appended.
    let mut header_text = apply_hd_sort_order(
        &crate::header_text::read_raw_header_text_with_format(input, format.exact)?,
        &so,
        ss.as_deref(),
    );
    if let Some(argv) = pg_argv {
        header_text = crate::pg::add_samtools_pg(&header_text, argv).map_err(io::Error::other)?;
    }

    {
        let mut writer = open_output(output, fmt, &header, &header_text)?;
        for rec in &records {
            writer.write_record(&header, rec)?;
        }
    }

    if write_index {
        let Some(path) = output else {
            return Err(io::Error::new(
                io::ErrorKind::InvalidInput,
                "--write-index requires -o FILE",
            ));
        };
        write_bam_index(path)?;
    }
    Ok(())
}

fn read_bam_records(input: &Path) -> io::Result<(sam::Header, Vec<RecordBuf>)> {
    let mut reader = bam::io::Reader::new(File::open(input)?);
    let header = reader.read_header()?;
    let mut records = Vec::new();
    loop {
        let mut record = RecordBuf::default();
        if reader.read_record_buf(&header, &mut record)? == 0 {
            break;
        }
        records.push(record);
    }
    Ok((header, records))
}

fn read_sam_records(input: &Path) -> io::Result<(sam::Header, Vec<RecordBuf>)> {
    crate::sam_compat::read_sam_records_tolerant(input)
}

fn read_cram_records(input: &Path) -> io::Result<(sam::Header, Vec<RecordBuf>)> {
    let reference = current_global_args().reference.ok_or_else(|| {
        io::Error::new(
            io::ErrorKind::InvalidInput,
            "CRAM input requires top-level --reference FILE",
        )
    })?;
    let text =
        htslib_rs::alignment_compat::view_cram_as_sam_text_from_path_with_reference_and_limit(
            input, reference, None,
        )?;
    let mut reader = sam::io::Reader::new(BufReader::new(Cursor::new(text)));
    read_sam_records_from_reader(&mut reader)
}

fn read_sam_records_from_reader<R>(
    reader: &mut sam::io::Reader<R>,
) -> io::Result<(sam::Header, Vec<RecordBuf>)>
where
    R: BufRead,
{
    let header = reader.read_header()?;
    let mut records = Vec::new();
    loop {
        let mut record = RecordBuf::default();
        if reader.read_record_buf(&header, &mut record)? == 0 {
            break;
        }
        records.push(record);
    }
    Ok((header, records))
}

fn coordinate_key(r: &RecordBuf) -> (i32, i64) {
    let tid = r
        .reference_sequence_id()
        .map(|t| t as i32)
        .unwrap_or(i32::MAX);
    let pos = r.alignment_start().map(usize::from).unwrap_or(0) as i64;
    (tid, pos)
}

fn name_key(r: &RecordBuf) -> Vec<u8> {
    r.name().map(|s| s.to_vec()).unwrap_or_default()
}

/// Port of `bam_sort.c`'s `strnum_cmp` natural-order comparison: runs of
/// digits compare numerically (leading zeros skipped, then by length,
/// then by first differing digit); everything else byte-wise.
fn strnum_cmp(a: &[u8], b: &[u8]) -> Ordering {
    let (mut ia, mut ib) = (0usize, 0usize);
    let is_digit = |c: u8| c.is_ascii_digit();
    while ia < a.len() && ib < b.len() {
        let (ca, cb) = (a[ia], b[ib]);
        if !is_digit(ca) || !is_digit(cb) {
            if ca != cb {
                return ca.cmp(&cb);
            }
            ia += 1;
            ib += 1;
        } else {
            while ia < a.len() && a[ia] == b'0' {
                ia += 1;
            }
            while ib < b.len() && b[ib] == b'0' {
                ib += 1;
            }
            while ia < a.len() && ib < b.len() && is_digit(a[ia]) && a.get(ia) == b.get(ib) {
                ia += 1;
                ib += 1;
            }
            let diff =
                a.get(ia).copied().unwrap_or(0) as i32 - b.get(ib).copied().unwrap_or(0) as i32;
            while ia < a.len() && ib < b.len() && is_digit(a[ia]) && is_digit(b[ib]) {
                ia += 1;
                ib += 1;
            }
            if ia < a.len() && is_digit(a[ia]) {
                return Ordering::Greater;
            } else if ib < b.len() && is_digit(b[ib]) {
                return Ordering::Less;
            } else if diff != 0 {
                return diff.cmp(&0);
            }
        }
    }
    let ra = ia < a.len();
    let rb = ib < b.len();
    if ra {
        Ordering::Greater
    } else if rb {
        Ordering::Less
    } else {
        Ordering::Equal
    }
}

/// `bam_sort.c` QueryName secondary key:
/// `((f&0xc0)<<8)|((f&0x100)<<3)|((f&0x800)>>3)` — READ1, READ2,
/// (primary), SUPPLEMENTARY, SECONDARY.
fn qname_flag_key(r: &RecordBuf) -> u32 {
    let f = u32::from(u16::from(r.flags()));
    ((f & 0xc0) << 8) | ((f & 0x100) << 3) | ((f & 0x800) >> 3)
}

// ----- `samtools sort -M` minimiser sort (non-indexed path) -----

/// `bam_sort.c` `#define XOR 0xdead7878beef7878`.
const MINHASH_XOR: u64 = 0xdead_7878_beef_7878;

/// Per-record sort key. For mapped records the coordinate triple is
/// used (`bam1_cmp_by_minhash` falls back to `bam1_cmp_core` whenever
/// either side is mapped); for unmapped records the 62-bit minimiser
/// key `m`, the `isize` tie-break, and the post-RC strand bit are used.
struct MinhashKey {
    mapped: bool,
    tid: u64,
    pos: u64,
    rev: u8,
    m: u64,
    isize_tb: i64,
}

/// 16-entry `bam_seqi` → 0..3 forward base table (`L` in `minhash`).
const MINHASH_L: [u64; 16] = [0, 0, 1, 0, 2, 0, 0, 0, 3, 0, 0, 0, 0, 0, 0, 0];
/// 16-entry reverse-complement base table (`R`, pre-shift).
const MINHASH_R: [u64; 16] = [0, 3, 2, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0];

/// ASCII base → 4-bit `bam_seqi` code (`=ACMGRSVTWYHKDBN`).
fn base_to_4bit(c: u8) -> usize {
    match c.to_ascii_uppercase() {
        b'=' => 0,
        b'A' => 1,
        b'C' => 2,
        b'M' => 3,
        b'G' => 4,
        b'R' => 5,
        b'S' => 6,
        b'V' => 7,
        b'T' => 8,
        b'W' => 9,
        b'Y' => 10,
        b'H' => 11,
        b'K' => 12,
        b'D' => 13,
        b'B' => 14,
        _ => 15,
    }
}

/// Faithful port of `bam_sort.c` `minhash` (the windowed scan).
/// `i_start` is the C `*curr_pos` in-value; returns
/// `(minhashf, curr_pos, is_rev, end)` where
/// `curr_pos = minhashpf - (kmer-1)` and `end = (i_end == len)`.
#[allow(clippy::too_many_arguments)]
fn minhash(
    seq: &[u8],
    kmer: i32,
    window: i32,
    i_start: i32,
    try_fwd: bool,
    try_rev: bool,
    no_squash: bool,
) -> (u64, i32, bool, bool) {
    let kmer = kmer as usize;
    let len = seq.len() as i32;
    let mask: u64 = if 2 * kmer >= 64 {
        u64::MAX
    } else {
        (1u64 << (2 * kmer)) - 1
    };
    let xor = MINHASH_XOR & mask;
    let shift = 2 * (kmer as u32 - 1);

    let i_start = i_start.max(0);
    let i_end = i_start + window.min(len - i_start);

    // Forward strand.
    let mut hashf: u64 = 0;
    let mut minhashf: u64 = u64::MAX;
    let mut minhashpf: i32 = i_start;
    if try_fwd {
        let mut last_base: i32 = -1;
        let mut i = i_start;
        let mut j = 0usize;
        while j < kmer - 1 && i < i_end {
            let base = base_to_4bit(seq[i as usize]);
            if no_squash || last_base != base as i32 {
                last_base = base as i32;
                hashf = (hashf << 2) | MINHASH_L[base];
                j += 1;
            }
            i += 1;
        }
        if no_squash {
            while i < i_end {
                let base = base_to_4bit(seq[i as usize]);
                hashf = (hashf << 2) | MINHASH_L[base];
                let hashfx = (hashf ^ MINHASH_XOR) & mask;
                if minhashf > hashfx {
                    minhashf = hashfx;
                    minhashpf = i;
                }
                i += 1;
            }
        } else {
            while i < i_end {
                let base = base_to_4bit(seq[i as usize]);
                if last_base != base as i32 {
                    last_base = base as i32;
                    hashf = (hashf << 2) | MINHASH_L[base];
                    let hashfx = (hashf ^ MINHASH_XOR) & mask;
                    if minhashf > hashfx {
                        minhashf = hashfx;
                        minhashpf = i;
                    }
                }
                i += 1;
            }
        }
    }

    let mut is_rev = false;
    if try_rev {
        let mut hashr: u64 = 0;
        let mut minhashr: u64 = u64::MAX;
        let mut minhashpr: i32 = i_start;
        let mut last_base: i32 = -1;
        let mut i = i_start;
        let mut j = 0usize;
        while j < kmer - 1 && i < len {
            let base = base_to_4bit(seq[i as usize]);
            if no_squash || last_base != base as i32 {
                last_base = base as i32;
                hashr = (hashr >> 2) | (MINHASH_R[base] << shift);
                j += 1;
            }
            i += 1;
        }
        if no_squash {
            while i < i_end {
                let base = base_to_4bit(seq[i as usize]);
                hashr = (hashr >> 2) | (MINHASH_R[base] << shift);
                if minhashr > (hashr ^ xor) {
                    minhashr = hashr ^ xor;
                    minhashpr = len - i + kmer as i32 - 2;
                }
                i += 1;
            }
        } else {
            while i < i_end {
                let base = base_to_4bit(seq[i as usize]);
                if last_base != base as i32 {
                    last_base = base as i32;
                    hashr = (hashr >> 2) | (MINHASH_R[base] << shift);
                    if minhashr > (hashr ^ xor) {
                        minhashr = hashr ^ xor;
                        minhashpr = len - i + kmer as i32 - 2;
                    }
                }
                i += 1;
            }
        }
        if minhashr < minhashf {
            minhashf = minhashr;
            minhashpf = minhashpr;
            is_rev = true;
        }
    }

    (
        minhashf,
        minhashpf - (kmer as i32 - 1),
        is_rev,
        i_end == len,
    )
}

/// `bam_sort.c` `build_minhash_index` (forward strand only): read each
/// reference sequence from `ref_path` (FASTA), slide a `window`-wide
/// minhash and record `hash -> tpos+pos` with the high `UNIQ_BIT` set
/// once a hash recurs (so unique vs duplicate placements are known).
fn build_minhash_index(
    ref_path: &Path,
    kmer: i32,
    window: i32,
    no_squash: bool,
) -> io::Result<std::collections::HashMap<u64, u64>> {
    use htslib_rs::fasta;
    const UNIQ_BIT: u64 = 1 << 60;
    let reader = File::open(ref_path).map(BufReader::new)?;
    let mut reader = fasta::io::Reader::new(reader);
    let mut map: std::collections::HashMap<u64, u64> = std::collections::HashMap::new();
    let mut tpos: u64 = 0;
    for result in reader.records() {
        let record = result?;
        let mut seq = record.sequence().as_ref().to_vec();
        seq.make_ascii_uppercase();
        let len = seq.len() as i32;
        if len < window {
            continue;
        }
        let mut pos: i32 = 0;
        loop {
            let last_pos = pos;
            let (hashf, out_pos, _rev, end) =
                minhash(&seq, kmer, window, pos, true, false, no_squash);
            let key_pos = tpos.wrapping_add(out_pos as i64 as u64);
            map.entry(hashf)
                .and_modify(|v| *v = key_pos | UNIQ_BIT)
                .or_insert(key_pos);
            pos = (last_pos + kmer).max(out_pos + 1);
            if end {
                break;
            }
        }
        tpos += seq.len() as u64;
    }
    Ok(map)
}

/// `bam_sort.c` `minhash_with_idx[_squash]` merged via the `squash`
/// flag: scan the whole read, preferring hashes that are uniquely (then
/// non-uniquely) placed in `kmer_h` over absent ones; on a hit the key
/// becomes the reference position. Returns `(minhashf, pos, dir)`.
fn minhash_with_idx(
    seq: &[u8],
    kmer: i32,
    kmer_h: &std::collections::HashMap<u64, u64>,
    try_rev: bool,
    squash: bool,
) -> (u64, i32, bool) {
    const UNIQ_BIT: u64 = 1 << 60;
    const UNIQ_MASK: u64 = UNIQ_BIT - 1;
    let uniq_test = |x: u64| (x & UNIQ_BIT) == 0;
    let kmer = kmer as usize;
    let len = seq.len() as i32;
    let mask: u64 = if 2 * kmer >= 64 {
        u64::MAX
    } else {
        (1u64 << (2 * kmer)) - 1
    };
    let xor = MINHASH_XOR & mask;
    let shift = 2 * (kmer as u32 - 1);

    // Forward.
    let mut hashf: u64 = 0;
    let mut minhashf = u64::MAX;
    let mut minhashfi = u64::MAX;
    let mut minhashfd = u64::MAX;
    let (mut minhashpf, mut minhashpfi, mut minhashpfd) = (0i32, 0i32, 0i32);
    let mut last_base: i32 = -1;
    let mut i = 0i32;
    let mut j = 0usize;
    while j < kmer - 1 && i < len {
        let base = base_to_4bit(seq[i as usize]);
        if squash && base as i32 == last_base {
            i += 1;
            continue;
        }
        last_base = base as i32;
        j += 1;
        hashf = (hashf << 2) | MINHASH_L[base];
        i += 1;
    }
    let mut found_f = 0i32;
    while i < len {
        let base = base_to_4bit(seq[i as usize]);
        if squash && base as i32 == last_base {
            i += 1;
            continue;
        }
        last_base = base as i32;
        hashf = ((hashf << 2) | MINHASH_L[base]) & mask;
        let hashfx = hashf ^ xor;
        let mut index = 0;
        if (minhashfi > hashfx || (found_f < 2 && minhashfd > hashfx))
            && let Some(&v) = kmer_h.get(&hashfx)
        {
            index = if uniq_test(v) { 2 } else { 1 };
        }
        found_f |= index;
        match index {
            2 => {
                minhashfi = hashfx;
                minhashpfi = i;
            }
            1 => {
                minhashfd = hashfx;
                minhashpfd = i;
            }
            _ => {
                if minhashf > hashfx {
                    minhashf = hashfx;
                    minhashpf = i;
                }
            }
        }
        i += 1;
    }
    if minhashfi != u64::MAX {
        minhashf = minhashfi;
        minhashpf = minhashpfi;
    } else if minhashfd != u64::MAX {
        minhashf = minhashfd;
        minhashpf = minhashpfd;
    }

    let mut dir = false;
    if try_rev {
        let mut hashr: u64 = 0;
        let mut minhashr = u64::MAX;
        let mut minhashri = u64::MAX;
        let mut minhashrd = u64::MAX;
        let (mut minhashpr, mut minhashpri, mut minhashprd) = (0i32, 0i32, 0i32);
        let mut last_base: i32 = -1;
        let mut i = 0i32;
        let mut j = 0usize;
        while j < kmer - 1 && i < len {
            let base = base_to_4bit(seq[i as usize]);
            if squash && base as i32 == last_base {
                i += 1;
                continue;
            }
            last_base = base as i32;
            j += 1;
            hashr = (hashr >> 2) | (MINHASH_R[base] << shift);
            i += 1;
        }
        let mut found_r = 0i32;
        while i < len {
            let base = base_to_4bit(seq[i as usize]);
            if squash && base as i32 == last_base {
                i += 1;
                continue;
            }
            last_base = base as i32;
            hashr = (hashr >> 2) | (MINHASH_R[base] << shift);
            let hashrx = hashr ^ xor;
            let mut index = 0;
            if (minhashri > hashrx || (found_r < 2 && minhashrd > hashrx))
                && let Some(&v) = kmer_h.get(&hashrx)
            {
                index = if uniq_test(v) { 2 } else { 1 };
            }
            found_r |= index;
            match index {
                2 => {
                    minhashri = hashrx;
                    minhashpri = i;
                }
                1 => {
                    minhashrd = hashrx;
                    minhashprd = i;
                }
                _ => {
                    if minhashr > hashrx {
                        minhashr = hashrx;
                        minhashpr = i;
                    }
                }
            }
            i += 1;
        }
        if minhashri != u64::MAX {
            minhashr = minhashri;
            minhashpr = minhashpri;
        } else if minhashrd != u64::MAX {
            minhashr = minhashrd;
            minhashpr = minhashprd;
        }
        if ((minhashf > minhashr) || (found_f == 0 && found_r != 0))
            && (found_f == 0 || found_r != 0)
        {
            minhashf = minhashr;
            minhashpf = len - minhashpr + kmer as i32 - 2;
            dir = true;
        }
    }

    // Indexed kmer → its reference position (mask off the uniq bit).
    if let Some(&v) = kmer_h.get(&minhashf) {
        minhashf = v & UNIQ_MASK;
    }
    let out = if minhashf != u64::MAX { minhashf } else { 0 };
    (out, minhashpf, dir)
}

/// `bam_sort.c` `reverse_complement`: reverse-complement the sequence,
/// reverse the quality scores, and toggle the REVERSE flag bit.
fn minhash_reverse_complement(r: &mut RecordBuf) {
    let seq = r.sequence_mut().as_mut();
    seq.reverse();
    for b in seq.iter_mut() {
        *b = minhash_complement(*b);
    }
    r.quality_scores_mut().as_mut().reverse();
    let f = r.flags_mut();
    *f ^= htslib_rs::sam::alignment::record::Flags::REVERSE_COMPLEMENTED;
}

fn minhash_complement(b: u8) -> u8 {
    match b {
        b'A' => b'T',
        b'C' => b'G',
        b'G' => b'C',
        b'T' => b'A',
        b'M' => b'K',
        b'R' => b'Y',
        b'W' => b'W',
        b'S' => b'S',
        b'Y' => b'R',
        b'K' => b'M',
        b'V' => b'B',
        b'H' => b'D',
        b'D' => b'H',
        b'B' => b'V',
        b'a' => b't',
        b'c' => b'g',
        b'g' => b'c',
        b't' => b'a',
        b'm' => b'k',
        b'r' => b'y',
        b'w' => b'w',
        b's' => b's',
        b'y' => b'r',
        b'k' => b'm',
        b'v' => b'b',
        b'h' => b'd',
        b'd' => b'h',
        b'b' => b'v',
        b'N' | b'n' => b,
        _ => b'N',
    }
}

/// `worker_minhash` for one record + key extraction. With `kmer_h` the
/// indexed path (`minhash_with_idx[_squash]`, `mh -= pos`) is taken;
/// otherwise the non-indexed path (`mh += 1<<30`, `isize = 65535-pos`).
fn minhash_prepare(
    r: &mut RecordBuf,
    opts: &MinhashOpts,
    kmer_h: Option<&std::collections::HashMap<u64, u64>>,
) -> MinhashKey {
    if let Some(tid) = r.reference_sequence_id() {
        // Mapped: keep coordinate order (bam1_cmp_core path).
        let pos = r.alignment_start().map(usize::from).unwrap_or(0) as u64; // core.pos+1
        let rev = u8::from(r.flags().is_reverse_complemented());
        return MinhashKey {
            mapped: true,
            tid: tid as u64,
            pos,
            rev,
            m: 0,
            isize_tb: 0,
        };
    }

    let seq: Vec<u8> = r.sequence().as_ref().to_vec();
    let len = seq.len() as i32;
    let (mh, isize_tb, is_rev) = if let Some(kmer_h) = kmer_h {
        let (minhashf, pos, dir) =
            minhash_with_idx(&seq, opts.kmer, kmer_h, opts.try_rev, !opts.no_squash);
        if dir {
            minhash_reverse_complement(r);
        }
        // worker_minhash indexed branch: `mh -= pos; pos = 0;`.
        (minhashf.wrapping_sub(pos as i64 as u64), 0i64, dir)
    } else {
        let (minhashf, curr, is_rev, _end) =
            minhash(&seq, opts.kmer, len, 0, true, opts.try_rev, opts.no_squash);
        if is_rev {
            minhash_reverse_complement(r);
        }
        let mh = minhashf.wrapping_add(1 << 30);
        let isize_tb = if 65535 - curr >= 0 { 65535 - curr } else { 0 } as i64;
        (mh, isize_tb, is_rev)
    };
    let m = (((mh >> 31) & 0x7fff_ffff) << 31) | (mh & 0x7fff_ffff);
    MinhashKey {
        mapped: false,
        tid: u64::MAX, // core.tid == -1 cast to uint64
        pos: 0,
        rev: u8::from(is_rev),
        m,
        isize_tb,
    }
}

/// `bam1_cmp_by_minhash` → `Ordering` (stable sort preserves input
/// order on `Equal`, mirroring `ks_mergesort`).
fn minhash_cmp(a: &MinhashKey, b: &MinhashKey) -> Ordering {
    if a.mapped || b.mapped {
        // bam1_cmp_core, MinHash (non-QueryName) branch: tid, then
        // pos+1, then reverse bit.
        return a
            .tid
            .cmp(&b.tid)
            .then_with(|| a.pos.cmp(&b.pos))
            .then_with(|| a.rev.cmp(&b.rev));
    }
    a.m.cmp(&b.m)
        // bigger isize sorts first (A.isize > B.isize → A before B).
        .then_with(|| b.isize_tb.cmp(&a.isize_tb))
        // bam1_cmp_core tail: m equal ⇒ pos equal ⇒ reverse bit asc.
        .then_with(|| a.rev.cmp(&b.rev))
}

/// Full `bam_sort.c` QueryName comparator. `bam_sort.c`'s `strnum_cmp`
/// is natural-order unless `-N` (`!natural_sort`), where it is a plain
/// byte `strcmp`.
fn name_cmp(a: &RecordBuf, b: &RecordBuf, natural_sort: bool) -> Ordering {
    let (ka, kb) = (name_key(a), name_key(b));
    let primary = if natural_sort {
        strnum_cmp(&ka, &kb)
    } else {
        ka.cmp(&kb)
    };
    primary.then_with(|| qname_flag_key(a).cmp(&qname_flag_key(b)))
}

fn compare_by_tag(
    a: &RecordBuf,
    b: &RecordBuf,
    tag: [u8; 2],
    name_sort: bool,
    natural_sort: bool,
) -> Ordering {
    tag_sort_value(a, tag)
        .cmp(&tag_sort_value(b, tag))
        .then_with(|| {
            if name_sort {
                name_cmp(a, b, natural_sort)
            } else {
                coordinate_key(a).cmp(&coordinate_key(b))
            }
        })
        .then_with(|| name_key(a).cmp(&name_key(b)))
}

#[derive(Clone, Debug)]
enum TagSortValue {
    Missing,
    Character(u8),
    Array(String),
    Text(Vec<u8>),
    Int(i64),
    Float(f32),
}

impl TagSortValue {
    fn rank(&self) -> u8 {
        match self {
            Self::Missing => 0,
            Self::Character(_) => b'A',
            Self::Array(_) => b'B',
            Self::Text(_) => b'H',
            Self::Int(_) => b'c',
            Self::Float(_) => b'f',
        }
    }
}

impl PartialEq for TagSortValue {
    fn eq(&self, other: &Self) -> bool {
        self.cmp(other) == Ordering::Equal
    }
}

impl Eq for TagSortValue {}

impl PartialOrd for TagSortValue {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        Some(self.cmp(other))
    }
}

impl Ord for TagSortValue {
    fn cmp(&self, other: &Self) -> Ordering {
        use TagSortValue::*;

        match (self, other) {
            (Missing, Missing) => Ordering::Equal,
            (Missing, _) => Ordering::Less,
            (_, Missing) => Ordering::Greater,
            (Int(a), Int(b)) => a.cmp(b),
            (Float(a), Float(b)) => a.total_cmp(b),
            (Int(a), Float(b)) => (*a as f32).total_cmp(b),
            (Float(a), Int(b)) => a.total_cmp(&(*b as f32)),
            (Character(a), Character(b)) => a.cmp(b),
            (Text(a), Text(b)) => a.cmp(b),
            (Array(a), Array(b)) => a.cmp(b),
            _ => self.rank().cmp(&other.rank()),
        }
    }
}

fn tag_sort_value(record: &RecordBuf, tag: [u8; 2]) -> TagSortValue {
    use sam::alignment::record_buf::data::field::Value;

    match record.data().get(&tag) {
        None => TagSortValue::Missing,
        Some(Value::Character(c)) => TagSortValue::Character(*c),
        Some(Value::Int8(n)) => TagSortValue::Int(i64::from(*n)),
        Some(Value::UInt8(n)) => TagSortValue::Int(i64::from(*n)),
        Some(Value::Int16(n)) => TagSortValue::Int(i64::from(*n)),
        Some(Value::UInt16(n)) => TagSortValue::Int(i64::from(*n)),
        Some(Value::Int32(n)) => TagSortValue::Int(i64::from(*n)),
        Some(Value::UInt32(n)) => TagSortValue::Int(i64::from(*n)),
        Some(Value::Float(n)) => TagSortValue::Float(*n),
        Some(Value::String(s)) | Some(Value::Hex(s)) => TagSortValue::Text(s.to_vec()),
        Some(Value::Array(array)) => TagSortValue::Array(format!("{:?}", array)),
    }
}

fn set_sort_order(header: &mut sam::Header, so: &str, ss: Option<&str>) {
    use bstr::BString;
    use sam::header::record::value::map::{self, Map};
    if let Some(hd) = header.header_mut() {
        hd.other_fields_mut()
            .insert(map::header::tag::SORT_ORDER, BString::from(so));
        match ss {
            Some(ss) => {
                hd.other_fields_mut()
                    .insert(map::header::tag::SUBSORT_ORDER, BString::from(ss));
            }
            None => {
                hd.other_fields_mut()
                    .shift_remove(&map::header::tag::SUBSORT_ORDER);
            }
        }
    } else {
        let mut hd: Map<map::Header> = Map::default();
        hd.other_fields_mut()
            .insert(map::header::tag::SORT_ORDER, BString::from(so));
        if let Some(ss) = ss {
            hd.other_fields_mut()
                .insert(map::header::tag::SUBSORT_ORDER, BString::from(ss));
        }
        *header.header_mut() = Some(hd);
    }
}

trait SortSink {
    fn write_record(&mut self, header: &sam::Header, record: &RecordBuf) -> io::Result<()>;
}

struct BamFile(bam::io::Writer<bgzf::io::Writer<File>>);
struct BamStdout(bam::io::Writer<bgzf::io::Writer<io::Stdout>>);
struct SamFile(File);
struct SamStdout(io::Stdout);

impl SortSink for BamFile {
    fn write_record(&mut self, header: &sam::Header, record: &RecordBuf) -> io::Result<()> {
        use sam::alignment::io::Write as _;
        self.0.write_alignment_record(header, record)
    }
}
impl SortSink for BamStdout {
    fn write_record(&mut self, header: &sam::Header, record: &RecordBuf) -> io::Result<()> {
        use sam::alignment::io::Write as _;
        self.0.write_alignment_record(header, record)
    }
}
impl SortSink for SamFile {
    fn write_record(&mut self, header: &sam::Header, record: &RecordBuf) -> io::Result<()> {
        // Shared renderer: htslib `%g` float aux spelling.
        crate::sam_render::write_record(&mut self.0, header, record)
    }
}
impl SortSink for SamStdout {
    fn write_record(&mut self, header: &sam::Header, record: &RecordBuf) -> io::Result<()> {
        crate::sam_render::write_record(&mut self.0, header, record)
    }
}

/// Replaces/sets the `@HD` line's `SO:`/`SS:` fields in raw header text,
/// preserving every other line and field verbatim (so `@RG`/`@SQ`/`@CO`
/// keep their original byte form). Inserts an `@HD` if absent.
fn apply_hd_sort_order(raw: &str, so: &str, ss: Option<&str>) -> String {
    let mut lines: Vec<String> = Vec::new();
    let mut had_hd = false;
    for line in raw.lines() {
        if line.starts_with("@HD") {
            had_hd = true;
            let mut fields: Vec<&str> = line
                .split('\t')
                .filter(|f| !f.starts_with("SO:") && !f.starts_with("SS:"))
                .collect();
            let mut nl = fields.join("\t");
            if fields.is_empty() {
                nl.push_str("@HD");
            }
            nl.push_str(&format!("\tSO:{so}"));
            if let Some(ss) = ss {
                nl.push_str(&format!("\tSS:{ss}"));
            }
            lines.push(nl);
            let _ = &mut fields;
        } else {
            lines.push(line.to_string());
        }
    }
    if !had_hd {
        let hd = match ss {
            Some(ss) => format!("@HD\tVN:1.6\tSO:{so}\tSS:{ss}"),
            None => format!("@HD\tVN:1.6\tSO:{so}"),
        };
        lines.insert(0, hd);
    }
    let mut s = lines.join("\n");
    s.push('\n');
    s
}

fn open_output(
    out: Option<&Path>,
    fmt: OutFmt,
    header: &sam::Header,
    header_text: &str,
) -> io::Result<Box<dyn SortSink>> {
    match (out, fmt) {
        (Some(p), OutFmt::Sam) => {
            let mut file = File::create(p)?;
            file.write_all(header_text.as_bytes())?;
            Ok(Box::new(SamFile(file)))
        }
        (Some(p), OutFmt::Bam) => {
            let file = File::create(p)?;
            let mut writer = bam::io::Writer::new(file);
            writer.write_header(header)?;
            Ok(Box::new(BamFile(writer)))
        }
        (None, OutFmt::Sam) => {
            let mut stdout = io::stdout();
            stdout.write_all(header_text.as_bytes())?;
            Ok(Box::new(SamStdout(stdout)))
        }
        (None, OutFmt::Bam) => {
            let mut writer = bam::io::Writer::new(io::stdout());
            writer.write_header(header)?;
            Ok(Box::new(BamStdout(writer)))
        }
    }
}

fn write_bam_index(path: &Path) -> io::Result<()> {
    let index = htslib_rs::index_compat::build_bai(path)?;
    htslib_rs::index_compat::write_bai(append_extension(path, "bai"), &index)
}

fn append_extension(path: &Path, ext: &str) -> PathBuf {
    let mut s = path.as_os_str().to_os_string();
    s.push(".");
    s.push(ext);
    PathBuf::from(s)
}

fn print_usage() -> io::Result<()> {
    let mut w = io::stderr().lock();
    writeln!(w, "Usage: samtools sort [options] <in.bam|in.sam|in.cram>")?;
    writeln!(
        w,
        "  -n              sort by read name (default: coordinate)"
    )?;
    writeln!(
        w,
        "  -t TAG          sort by auxiliary tag, then coordinate/name"
    )?;
    writeln!(w, "  -o FILE         write output to FILE (default stdout)")?;
    writeln!(w, "  --output-fmt sam|bam")?;
    writeln!(
        w,
        "  -@/-m/-T/-K     accepted but currently ignored (in-memory sort only)"
    )?;
    writeln!(w, "  --write-index   write BAI index for BAM file output")?;
    Ok(())
}
