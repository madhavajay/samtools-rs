//! Helpers for constructing `@PG` header lines.

use std::collections::HashSet;
use std::ffi::OsString;
use std::io;

use htslib_rs::sam;

use crate::version::SAMTOOLS_VERSION;

/// Joins argv for a `@PG CL:` tag using HTSlib's `stringify_argv` behavior.
pub fn stringify_argv(args: &[OsString]) -> String {
    args.iter()
        .map(|arg| arg.to_string_lossy().replace('\t', " "))
        .collect::<Vec<_>>()
        .join(" ")
}

/// Adds samtools' standard `@PG` line(s) to raw SAM header text.
pub fn add_samtools_pg(header_text: &str, argv: &[OsString]) -> Result<String, String> {
    let command_line = stringify_argv(argv);
    add_pg(
        header_text,
        PgOptions {
            name: "samtools",
            version: Some(SAMTOOLS_VERSION),
            command_line: Some(&command_line),
        },
    )
}

/// Adds samtools' standard `@PG` line(s) to a typed SAM header.
pub fn add_samtools_pg_to_header(
    header: &sam::Header,
    argv: &[OsString],
) -> io::Result<sam::Header> {
    let mut bytes = Vec::new();
    {
        let mut writer = sam::io::Writer::new(&mut bytes);
        writer.write_header(header)?;
    }
    let header_text =
        String::from_utf8(bytes).map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))?;
    let header_text = add_samtools_pg(&header_text, argv)
        .map_err(|e| io::Error::new(io::ErrorKind::InvalidInput, e))?;

    header_text
        .parse()
        .map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))
}

/// Program-line fields used when adding a `@PG` record.
#[derive(Clone, Debug, Eq, PartialEq)]
pub struct PgOptions<'a> {
    pub name: &'a str,
    pub version: Option<&'a str>,
    pub command_line: Option<&'a str>,
}

/// Adds `@PG` line(s) to raw SAM header text, preserving existing line order.
pub fn add_pg(header_text: &str, options: PgOptions<'_>) -> Result<String, String> {
    validate_field("PN", options.name)?;
    if let Some(version) = options.version {
        validate_field("VN", version)?;
    }
    if let Some(command_line) = options.command_line {
        validate_field("CL", command_line)?;
    }

    let chains = PgChains::parse(header_text);
    let terminals = chains.terminals();
    let mut used_ids = chains.ids;
    let mut output = String::from(header_text);

    if !output.is_empty() && !output.ends_with('\n') {
        output.push('\n');
    }

    if terminals.is_empty() {
        let id = unique_pg_id(options.name, &mut used_ids);
        push_pg_line(&mut output, &id, None, &options);
    } else {
        for terminal in terminals {
            let id = unique_pg_id(options.name, &mut used_ids);
            push_pg_line(&mut output, &id, Some(&terminal), &options);
        }
    }

    Ok(output)
}

fn push_pg_line(output: &mut String, id: &str, pp: Option<&str>, options: &PgOptions<'_>) {
    // Field order matches upstream `sam_hdr_add_pg`: ID, PN, PP, VN, CL.
    // (Upstream places PP before VN/CL; the test harness strips `\tVN:.*`
    // for version-independent comparison, so PP must precede VN.)
    output.push_str("@PG\tID:");
    output.push_str(id);
    output.push_str("\tPN:");
    output.push_str(options.name);

    if let Some(pp) = pp {
        output.push_str("\tPP:");
        output.push_str(pp);
    }

    if let Some(version) = options.version {
        output.push_str("\tVN:");
        output.push_str(version);
    }

    if let Some(command_line) = options.command_line {
        output.push_str("\tCL:");
        output.push_str(command_line);
    }

    output.push('\n');
}

fn unique_pg_id(base: &str, used_ids: &mut HashSet<String>) -> String {
    if used_ids.insert(base.to_owned()) {
        return base.to_owned();
    }

    for i in 1.. {
        let candidate = format!("{base}.{i}");
        if used_ids.insert(candidate.clone()) {
            return candidate;
        }
    }

    unreachable!("unbounded integer suffix search should always find an ID")
}

fn validate_field(tag: &str, value: &str) -> Result<(), String> {
    if value.is_empty() {
        return Err(format!("@PG {tag} value must not be empty"));
    }

    if value.contains('\n') || value.contains('\r') || value.contains('\t') {
        return Err(format!(
            "@PG {tag} value contains an invalid control character"
        ));
    }

    Ok(())
}

#[derive(Debug, Default)]
struct PgChains {
    ids: HashSet<String>,
    id_order: Vec<String>,
    parents: HashSet<String>,
}

impl PgChains {
    fn parse(header_text: &str) -> Self {
        let mut chains = Self::default();

        for line in header_text.lines() {
            if !line.starts_with("@PG\t") {
                continue;
            }

            let mut id = None;
            let mut parent = None;
            for field in line.split('\t').skip(1) {
                if let Some(value) = field.strip_prefix("ID:") {
                    id = Some(value.to_owned());
                } else if let Some(value) = field.strip_prefix("PP:") {
                    parent = Some(value.to_owned());
                }
            }

            if let Some(id) = id
                && chains.ids.insert(id.clone())
            {
                chains.id_order.push(id);
            }
            if let Some(parent) = parent {
                chains.parents.insert(parent);
            }
        }

        chains
    }

    fn terminals(&self) -> Vec<String> {
        self.id_order
            .iter()
            .filter(|id| !self.parents.contains(*id))
            .cloned()
            .collect()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn stringify_argv_joins_with_spaces_and_replaces_tabs() {
        let args = vec![
            OsString::from("samtools"),
            OsString::from("view"),
            OsString::from("read\tname.bam"),
        ];

        assert_eq!(stringify_argv(&args), "samtools view read name.bam");
    }

    #[test]
    fn adds_single_pg_line_without_existing_programs() {
        let header = "@HD\tVN:1.6\n@SQ\tSN:sq0\tLN:8\n";
        let out = add_pg(
            header,
            PgOptions {
                name: "samtools",
                version: Some("1.23.1"),
                command_line: Some("samtools view in.bam"),
            },
        )
        .unwrap();

        assert_eq!(
            out,
            concat!(
                "@HD\tVN:1.6\n",
                "@SQ\tSN:sq0\tLN:8\n",
                "@PG\tID:samtools\tPN:samtools\tVN:1.23.1\tCL:samtools view in.bam\n"
            )
        );
    }

    #[test]
    fn links_to_each_terminal_pg_chain_with_unique_ids() {
        let header = concat!(
            "@HD\tVN:1.6\n",
            "@PG\tID:prog1\tPN:prog1\n",
            "@PG\tID:prog2\tPN:prog2\tPP:prog1\n",
            "@PG\tID:samtools\tPN:samtools\n",
        );

        let out = add_pg(
            header,
            PgOptions {
                name: "samtools",
                version: None,
                command_line: None,
            },
        )
        .unwrap();

        assert_eq!(
            out,
            concat!(
                "@HD\tVN:1.6\n",
                "@PG\tID:prog1\tPN:prog1\n",
                "@PG\tID:prog2\tPN:prog2\tPP:prog1\n",
                "@PG\tID:samtools\tPN:samtools\n",
                "@PG\tID:samtools.1\tPN:samtools\tPP:prog2\n",
                "@PG\tID:samtools.2\tPN:samtools\tPP:samtools\n",
            )
        );
    }

    #[test]
    fn rejects_invalid_pg_field_text() {
        assert_eq!(
            add_pg(
                "",
                PgOptions {
                    name: "sam\ttools",
                    version: None,
                    command_line: None,
                },
            )
            .unwrap_err(),
            "@PG PN value contains an invalid control character"
        );
    }
}
