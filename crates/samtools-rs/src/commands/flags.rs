//! `samtools flags FLAGS...` — convert between textual and numeric flag
//! representations.
//!
//! Mirrors `bam_flags.c` in upstream samtools. Output must match byte-for-byte:
//! `0x<hex>\t<dec>\t<NAMES>\n` per FLAGS argument, with the usage banner using
//! HTSlib's exact wording.

use std::ffi::OsString;
use std::io::{self, Write};
use std::process::ExitCode;

use crate::bam_flag::{FLAG_NAMES, flag_to_str, str_to_flag};
use crate::diagnostics::print_error;

/// Entry point for `samtools flags`.
pub fn main(args: &[OsString]) -> ExitCode {
    if args.len() < 2 {
        let _ = write_usage(&mut io::stdout());
        return ExitCode::SUCCESS;
    }

    let mut stdout = io::stdout().lock();
    for arg in &args[1..] {
        let Some(s) = arg.to_str() else {
            print_error(
                "flags",
                format!("Could not parse \"{}\"", arg.to_string_lossy()),
            );
            let _ = write_usage(&mut io::stderr());
            return ExitCode::from(1);
        };
        let Some(mask) = str_to_flag(s) else {
            print_error("flags", format!("Could not parse \"{}\"", s));
            let _ = write_usage(&mut io::stderr());
            return ExitCode::from(1);
        };
        let _ = writeln!(
            stdout,
            "0x{:x}\t{}\t{}",
            mask,
            mask,
            flag_to_str(mask as u32)
        );
    }
    ExitCode::SUCCESS
}

fn write_usage<W: Write>(w: &mut W) -> io::Result<()> {
    writeln!(
        w,
        "About: Convert between textual and numeric flag representation"
    )?;
    writeln!(w, "Usage: samtools flags FLAGS...")?;
    writeln!(w)?;
    writeln!(
        w,
        "Each FLAGS argument is either an INT (in decimal/hexadecimal/octal) representing"
    )?;
    writeln!(
        w,
        "a combination of the following numeric flag values, or a comma-separated string"
    )?;
    writeln!(
        w,
        "NAME,...,NAME representing a combination of the following flag names:"
    )?;
    writeln!(w)?;
    for (bit, name) in FLAG_NAMES {
        writeln!(
            w,
            "{:#6x} {:>5}  {:<15}{}",
            bit,
            bit,
            name,
            flag_description(*bit)
        )?;
    }
    Ok(())
}

fn flag_description(bit: u32) -> &'static str {
    match bit {
        0x1 => "paired-end / multiple-segment sequencing technology",
        0x2 => "each segment properly aligned according to aligner",
        0x4 => "segment unmapped",
        0x8 => "next segment in the template unmapped",
        0x10 => "SEQ is reverse complemented",
        0x20 => "SEQ of next segment in template is rev.complemented",
        0x40 => "the first segment in the template",
        0x80 => "the last segment in the template",
        0x100 => "secondary alignment",
        0x200 => "not passing quality controls or other filters",
        0x400 => "PCR or optical duplicate",
        0x800 => "supplementary alignment",
        _ => "",
    }
}
