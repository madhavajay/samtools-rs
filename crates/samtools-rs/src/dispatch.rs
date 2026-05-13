//! Top-level `samtools <subcommand>` dispatcher.
//!
//! Mirrors the upstream `main()` and `usage()` in
//! `samtools/bamtk.c`, including subcommand aliases.

use std::ffi::OsString;
use std::io::{self, Write};
use std::process::ExitCode;

use crate::commands;
use crate::version::{HTSLIB_VERSION, SAMTOOLS_VERSION};

/// Subcommand entry-point signature. Each subcommand receives its own
/// argv slice (the first element is the subcommand name itself, matching
/// the convention upstream uses when calling `main_xxx(argc-1, argv+1)`
/// from `bamtk.c`).
pub type SubcommandMain = fn(&[OsString]) -> ExitCode;

/// One subcommand registration row. `aliases` are alternate names
/// (e.g. `idxstats` and `idxstat`) that all dispatch to `entry`.
pub struct Subcommand {
    /// Canonical name, used in help output.
    pub name: &'static str,
    /// Alternative names the user may type.
    pub aliases: &'static [&'static str],
    /// Entry point.
    pub entry: SubcommandMain,
}

/// Top-level entry point used by the binary crate. Parses
/// `argv[1]` and dispatches to the matching subcommand, or prints help.
pub fn run(args: Vec<OsString>) -> ExitCode {
    if args.len() < 2 {
        let _ = print_usage(&mut io::stderr());
        return ExitCode::from(1);
    }

    let sub = args[1].clone();
    let sub_str = match sub.to_str() {
        Some(s) => s,
        None => {
            let _ = writeln!(
                io::stderr(),
                "[main] unrecognized command '{}'",
                sub.to_string_lossy()
            );
            return ExitCode::from(1);
        }
    };

    match sub_str {
        "help" | "--help" => {
            if args.len() == 2 {
                let _ = print_usage(&mut io::stdout());
                return ExitCode::SUCCESS;
            }
            // "samtools help CMD [...]" -> "samtools CMD" with no args
            let cmd = args[2].clone();
            let new_args = vec![OsString::from("samtools"), cmd];
            return run(new_args);
        }
        "version" | "--version" => {
            let _ = print_long_version(&mut io::stdout());
            return ExitCode::SUCCESS;
        }
        "--version-only" => {
            let _ = writeln!(
                io::stdout(),
                "{}+htslib-{}",
                SAMTOOLS_VERSION,
                HTSLIB_VERSION
            );
            return ExitCode::SUCCESS;
        }
        "pileup" => {
            let _ = writeln!(
                io::stderr(),
                "[main] The `pileup' command has been removed. Please use `mpileup' instead."
            );
            return ExitCode::from(1);
        }
        _ => {}
    }

    for entry in SUBCOMMANDS {
        if entry.name == sub_str || entry.aliases.contains(&sub_str) {
            // Pass argv from the subcommand name onwards (mirrors C's argv+1).
            let sub_args: Vec<OsString> = args[1..].to_vec();
            let code = (entry.entry)(&sub_args);
            // Flush stdout; mirror bamtk.c's final fclose(stdout) check.
            let mut out = io::stdout();
            let _ = out.flush();
            return code;
        }
    }

    let _ = writeln!(io::stderr(), "[main] unrecognized command '{}'", sub_str);
    ExitCode::from(1)
}

/// Registry of every supported subcommand. Order follows upstream
/// `bamtk.c`'s usage groupings.
pub const SUBCOMMANDS: &[Subcommand] = &[
    // -- Indexing
    Subcommand {
        name: "dict",
        aliases: &[],
        entry: commands::dict::main,
    },
    Subcommand {
        name: "faidx",
        aliases: &[],
        entry: commands::faidx::main,
    },
    Subcommand {
        name: "fqidx",
        aliases: &[],
        entry: commands::fqidx::main,
    },
    Subcommand {
        name: "index",
        aliases: &[],
        entry: commands::index::main,
    },
    // -- Editing
    Subcommand {
        name: "calmd",
        aliases: &["fillmd"],
        entry: commands::calmd::main,
    },
    Subcommand {
        name: "fixmate",
        aliases: &[],
        entry: commands::fixmate::main,
    },
    Subcommand {
        name: "reheader",
        aliases: &[],
        entry: commands::reheader::main,
    },
    Subcommand {
        name: "targetcut",
        aliases: &[],
        entry: commands::targetcut::main,
    },
    Subcommand {
        name: "addreplacerg",
        aliases: &[],
        entry: commands::addreplacerg::main,
    },
    Subcommand {
        name: "markdup",
        aliases: &[],
        entry: commands::markdup::main,
    },
    Subcommand {
        name: "ampliconclip",
        aliases: &[],
        entry: commands::ampliconclip::main,
    },
    // -- File operations
    Subcommand {
        name: "collate",
        aliases: &["bamshuf"],
        entry: commands::collate::main,
    },
    Subcommand {
        name: "cat",
        aliases: &[],
        entry: commands::cat::main,
    },
    Subcommand {
        name: "consensus",
        aliases: &[],
        entry: commands::consensus::main,
    },
    Subcommand {
        name: "merge",
        aliases: &[],
        entry: commands::merge::main,
    },
    Subcommand {
        name: "mpileup",
        aliases: &[],
        entry: commands::mpileup::main,
    },
    Subcommand {
        name: "sort",
        aliases: &[],
        entry: commands::sort::main,
    },
    Subcommand {
        name: "split",
        aliases: &[],
        entry: commands::split::main,
    },
    Subcommand {
        name: "quickcheck",
        aliases: &[],
        entry: commands::quickcheck::main,
    },
    Subcommand {
        name: "fastq",
        aliases: &["fasta", "bam2fq"],
        entry: commands::fastq::main,
    },
    Subcommand {
        name: "import",
        aliases: &[],
        entry: commands::import::main,
    },
    Subcommand {
        name: "reference",
        aliases: &[],
        entry: commands::reference::main,
    },
    Subcommand {
        name: "reset",
        aliases: &[],
        entry: commands::reset::main,
    },
    Subcommand {
        name: "rmdup",
        aliases: &[],
        entry: commands::rmdup::main,
    },
    // -- Statistics
    Subcommand {
        name: "bedcov",
        aliases: &[],
        entry: commands::bedcov::main,
    },
    Subcommand {
        name: "coverage",
        aliases: &[],
        entry: commands::coverage::main,
    },
    Subcommand {
        name: "depth",
        aliases: &[],
        entry: commands::depth::main,
    },
    Subcommand {
        name: "flagstat",
        aliases: &["flagstats"],
        entry: commands::flagstat::main,
    },
    Subcommand {
        name: "idxstats",
        aliases: &["idxstat"],
        entry: commands::idxstats::main,
    },
    Subcommand {
        name: "cram-size",
        aliases: &[],
        entry: commands::cram_size::main,
    },
    Subcommand {
        name: "phase",
        aliases: &[],
        entry: commands::phase::main,
    },
    Subcommand {
        name: "stats",
        aliases: &["stat"],
        entry: commands::stats::main,
    },
    Subcommand {
        name: "ampliconstats",
        aliases: &[],
        entry: commands::ampliconstats::main,
    },
    Subcommand {
        name: "checksum",
        aliases: &[],
        entry: commands::checksum::main,
    },
    // -- Viewing
    Subcommand {
        name: "flags",
        aliases: &["flag"],
        entry: commands::flags::main,
    },
    Subcommand {
        name: "head",
        aliases: &[],
        entry: commands::head::main,
    },
    Subcommand {
        name: "view",
        aliases: &[],
        entry: commands::view::main,
    },
    Subcommand {
        name: "depad",
        aliases: &["pad2unpad"],
        entry: commands::depad::main,
    },
    Subcommand {
        name: "samples",
        aliases: &[],
        entry: commands::samples::main,
    },
];

/// Write the top-level usage banner. Matches upstream's `usage()` in
/// `bamtk.c` so the help text is byte-identical modulo the version line.
pub fn print_usage<W: Write>(w: &mut W) -> io::Result<()> {
    writeln!(w)?;
    writeln!(
        w,
        "Program: samtools (Tools for alignments in the SAM format)"
    )?;
    writeln!(
        w,
        "Version: {} (using htslib {})",
        SAMTOOLS_VERSION, HTSLIB_VERSION
    )?;
    writeln!(w)?;
    writeln!(w, "Usage:   samtools <command> [options]")?;
    writeln!(w)?;
    writeln!(w, "Commands:")?;
    writeln!(w, "  -- Indexing")?;
    writeln!(w, "     dict           create a sequence dictionary file")?;
    writeln!(w, "     faidx          index/extract FASTA")?;
    writeln!(w, "     fqidx          index/extract FASTQ")?;
    writeln!(w, "     index          index alignment")?;
    writeln!(w)?;
    writeln!(w, "  -- Editing")?;
    writeln!(
        w,
        "     calmd          recalculate MD/NM tags and '=' bases"
    )?;
    writeln!(w, "     fixmate        fix mate information")?;
    writeln!(w, "     reheader       replace BAM header")?;
    writeln!(
        w,
        "     targetcut      cut fosmid regions (for fosmid pool only)"
    )?;
    writeln!(w, "     addreplacerg   adds or replaces RG tags")?;
    writeln!(w, "     markdup        mark duplicates")?;
    writeln!(w, "     ampliconclip   clip oligos from the end of reads")?;
    writeln!(w)?;
    writeln!(w, "  -- File operations")?;
    writeln!(
        w,
        "     collate        shuffle and group alignments by name"
    )?;
    writeln!(w, "     cat            concatenate BAMs")?;
    writeln!(
        w,
        "     consensus      produce a consensus Pileup/FASTA/FASTQ"
    )?;
    writeln!(w, "     merge          merge sorted alignments")?;
    writeln!(w, "     mpileup        multi-way pileup")?;
    writeln!(w, "     sort           sort alignment file")?;
    writeln!(w, "     split          splits a file by read group")?;
    writeln!(
        w,
        "     quickcheck     quickly check if SAM/BAM/CRAM file appears intact"
    )?;
    writeln!(w, "     fastq          converts a BAM to a FASTQ")?;
    writeln!(w, "     fasta          converts a BAM to a FASTA")?;
    writeln!(
        w,
        "     import         Converts FASTA or FASTQ files to SAM/BAM/CRAM"
    )?;
    writeln!(
        w,
        "     reference      Generates a reference from aligned data"
    )?;
    writeln!(w, "     reset          Reverts aligner changes in reads")?;
    writeln!(w)?;
    writeln!(w, "  -- Statistics")?;
    writeln!(w, "     bedcov         read depth per BED region")?;
    writeln!(
        w,
        "     coverage       alignment depth and percent coverage"
    )?;
    writeln!(w, "     depth          compute the depth")?;
    writeln!(w, "     flagstat       simple stats")?;
    writeln!(w, "     idxstats       BAM index stats")?;
    writeln!(
        w,
        "     cram-size      list CRAM Content-ID and Data-Series sizes"
    )?;
    writeln!(w, "     phase          phase heterozygotes")?;
    writeln!(w, "     stats          generate stats (former bamcheck)")?;
    writeln!(w, "     ampliconstats  generate amplicon specific stats")?;
    writeln!(
        w,
        "     checksum       produce order-agnostic checksums of sequence content"
    )?;
    writeln!(w)?;
    writeln!(w, "  -- Viewing")?;
    writeln!(w, "     flags          explain BAM flags")?;
    writeln!(w, "     head           header viewer")?;
    writeln!(w, "     tview          text alignment viewer")?;
    writeln!(w, "     view           SAM<->BAM<->CRAM conversion")?;
    writeln!(w, "     depad          convert padded BAM to unpadded BAM")?;
    writeln!(
        w,
        "     samples        list the samples in a set of SAM/BAM/CRAM files"
    )?;
    writeln!(w)?;
    writeln!(w, "  -- Misc")?;
    writeln!(
        w,
        "     help [cmd]     display this help message or help for [cmd]"
    )?;
    writeln!(w, "     version        detailed version information")?;
    writeln!(w)?;
    Ok(())
}

/// Write the verbose `--version` / `version` banner.
pub fn print_long_version<W: Write>(w: &mut W) -> io::Result<()> {
    writeln!(w, "samtools {}", SAMTOOLS_VERSION)?;
    writeln!(w, "Using htslib {}", HTSLIB_VERSION)?;
    writeln!(w, "Copyright (C) 2025 Genome Research Ltd.")?;
    writeln!(w)?;
    writeln!(w, "Samtools compilation details:")?;
    writeln!(w, "    Features:       build=cargo curses=no")?;
    writeln!(w, "    HTSDIR:         htslib-rs (pure Rust)")?;
    Ok(())
}
