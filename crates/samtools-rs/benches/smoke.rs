use std::ffi::OsString;
use std::path::{Path, PathBuf};
use std::process::{Command, ExitCode, Stdio};
use std::time::{Duration, Instant};

use samtools_rs::run as samtools_run;

struct Case {
    name: &'static str,
    build_args: fn(&Path, usize) -> Vec<OsString>,
    build_external_args: Option<fn(&Path, usize) -> Vec<OsString>>,
    external_stdout_file: bool,
}

fn main() {
    let iterations = std::env::var("SAMTOOLS_RS_BENCH_ITERS")
        .ok()
        .and_then(|value| value.parse::<usize>().ok())
        .unwrap_or(3)
        .max(1);
    let tmp = std::env::temp_dir().join(format!("samtools-rs-bench-smoke-{}", std::process::id()));
    let _ = std::fs::remove_dir_all(&tmp);
    let rust_tmp = tmp.join("rust");
    let compare_tmp = tmp.join("compare");
    std::fs::create_dir_all(&rust_tmp).unwrap();
    let compare = std::env::var_os("SAMTOOLS_RS_BENCH_COMPARE")
        .filter(|value| !value.is_empty())
        .map(PathBuf::from)
        .map(resolve_compare_path);
    if compare.is_some() {
        std::fs::create_dir_all(&compare_tmp).unwrap();
    }

    let cases = [
        Case {
            name: "view",
            build_args: view_args,
            build_external_args: None,
            external_stdout_file: false,
        },
        Case {
            name: "sort",
            build_args: sort_args,
            build_external_args: None,
            external_stdout_file: false,
        },
        Case {
            name: "markdup",
            build_args: markdup_args,
            build_external_args: None,
            external_stdout_file: false,
        },
        Case {
            name: "stats",
            build_args: stats_args,
            build_external_args: Some(stats_external_args),
            external_stdout_file: true,
        },
        Case {
            name: "mpileup",
            build_args: mpileup_args,
            build_external_args: None,
            external_stdout_file: false,
        },
        Case {
            name: "coverage",
            build_args: coverage_args,
            build_external_args: None,
            external_stdout_file: false,
        },
        Case {
            name: "depth",
            build_args: depth_args,
            build_external_args: None,
            external_stdout_file: false,
        },
        Case {
            name: "checksum",
            build_args: checksum_args,
            build_external_args: None,
            external_stdout_file: false,
        },
    ];

    println!("samtools-rs smoke bench: {iterations} iteration(s)");
    if let Some(path) = compare.as_deref() {
        println!("external samtools comparison: {}", path.display());
    }
    for case in cases {
        let durations = run_case(&case, &rust_tmp, iterations);
        let rust_stats = summarize(&durations);
        if let Some(path) = compare.as_deref() {
            let compare_durations = run_external_case(&case, &compare_tmp, iterations, path);
            let compare_stats = summarize(&compare_durations);
            let ratio = rust_stats.mean.as_secs_f64() / compare_stats.mean.as_secs_f64();
            println!(
                "{:<8} rust min={:>8?} mean={:>8?}  external min={:>8?} mean={:>8?}  ratio={:.2}x",
                case.name,
                rust_stats.min,
                rust_stats.mean,
                compare_stats.min,
                compare_stats.mean,
                ratio
            );
        } else {
            println!(
                "{:<8} min={:>8?} mean={:>8?}",
                case.name, rust_stats.min, rust_stats.mean
            );
        }
    }
}

struct Stats {
    min: Duration,
    mean: Duration,
}

fn summarize(durations: &[Duration]) -> Stats {
    let total: Duration = durations.iter().copied().sum();
    Stats {
        mean: total / durations.len() as u32,
        min: durations.iter().copied().min().unwrap(),
    }
}

fn run_case(case: &Case, tmp: &Path, iterations: usize) -> Vec<Duration> {
    let mut durations = Vec::with_capacity(iterations);
    for i in 0..iterations {
        let args = (case.build_args)(tmp, i);
        let start = Instant::now();
        let code = samtools_run(args);
        let elapsed = start.elapsed();
        let status = exit_to_u8(code);
        assert_eq!(status, 0, "{} failed with exit code {}", case.name, status);
        durations.push(elapsed);
    }
    durations
}

fn run_external_case(case: &Case, tmp: &Path, iterations: usize, samtools: &Path) -> Vec<Duration> {
    let mut durations = Vec::with_capacity(iterations);
    for i in 0..iterations {
        let build_args = case.build_external_args.unwrap_or(case.build_args);
        let args = build_args(tmp, i);
        let stdout = if case.external_stdout_file {
            let path = tmp.join(format!("{}-{i}.external.out", case.name));
            Stdio::from(std::fs::File::create(path).unwrap())
        } else {
            Stdio::inherit()
        };
        let start = Instant::now();
        let status = Command::new(samtools)
            .args(args.iter().skip(1))
            .stdout(stdout)
            .status()
            .unwrap_or_else(|e| panic!("failed to run {}: {e}", samtools.display()));
        let elapsed = start.elapsed();
        assert!(
            status.success(),
            "{} external samtools failed with status {}",
            case.name,
            status
        );
        durations.push(elapsed);
    }
    durations
}

fn argv(rest: &[String]) -> Vec<OsString> {
    std::iter::once(OsString::from("samtools"))
        .chain(rest.iter().map(OsString::from))
        .collect()
}

fn fixtures_dir() -> PathBuf {
    let manifest = env!("CARGO_MANIFEST_DIR");
    PathBuf::from(manifest)
        .parent()
        .unwrap()
        .parent()
        .unwrap()
        .join("samtools")
        .join("test")
}

fn resolve_compare_path(path: PathBuf) -> PathBuf {
    if path.is_absolute() || path.exists() {
        return path;
    }

    let workspace = fixtures_dir()
        .parent()
        .and_then(Path::parent)
        .map(Path::to_path_buf);
    if let Some(workspace) = workspace {
        let candidate = workspace.join(&path);
        if candidate.exists() {
            return candidate;
        }
    }

    path
}

fn exit_to_u8(code: ExitCode) -> u8 {
    format!("{:?}", code)
        .chars()
        .filter(|c| c.is_ascii_digit())
        .collect::<String>()
        .parse()
        .unwrap_or(255)
}

fn view_args(tmp: &Path, i: usize) -> Vec<OsString> {
    let d = fixtures_dir();
    argv(&[
        "view".into(),
        "-h".into(),
        "-o".into(),
        tmp.join(format!("view-{i}.sam")).display().to_string(),
        d.join("dat/test_input_1_a.bam").display().to_string(),
    ])
}

fn sort_args(tmp: &Path, i: usize) -> Vec<OsString> {
    let d = fixtures_dir();
    argv(&[
        "sort".into(),
        "-O".into(),
        "bam".into(),
        "-o".into(),
        tmp.join(format!("sort-{i}.bam")).display().to_string(),
        d.join("dat/test_input_1_a.bam").display().to_string(),
    ])
}

fn markdup_args(tmp: &Path, i: usize) -> Vec<OsString> {
    let d = fixtures_dir();
    argv(&[
        "markdup".into(),
        "--no-PG".into(),
        "-O".into(),
        "sam".into(),
        d.join("markdup/5_markdup.sam").display().to_string(),
        tmp.join(format!("markdup-{i}.sam")).display().to_string(),
    ])
}

fn stats_args(tmp: &Path, i: usize) -> Vec<OsString> {
    let d = fixtures_dir();
    argv(&[
        "stats".into(),
        "-o".into(),
        tmp.join(format!("stats-{i}.txt")).display().to_string(),
        d.join("dat/mpileup.1.sam").display().to_string(),
    ])
}

fn stats_external_args(_tmp: &Path, _i: usize) -> Vec<OsString> {
    let d = fixtures_dir();
    argv(&[
        "stats".into(),
        d.join("dat/mpileup.1.sam").display().to_string(),
    ])
}

fn mpileup_args(tmp: &Path, i: usize) -> Vec<OsString> {
    let d = fixtures_dir();
    argv(&[
        "mpileup".into(),
        "-f".into(),
        d.join("dat/mpileup.ref.fa").display().to_string(),
        "-o".into(),
        tmp.join(format!("mpileup-{i}.txt")).display().to_string(),
        d.join("dat/mpileup.1.sam").display().to_string(),
    ])
}

fn coverage_args(tmp: &Path, i: usize) -> Vec<OsString> {
    let d = fixtures_dir();
    argv(&[
        "coverage".into(),
        "-o".into(),
        tmp.join(format!("coverage-{i}.txt")).display().to_string(),
        d.join("dat/test_input_1_a.bam").display().to_string(),
    ])
}

fn depth_args(tmp: &Path, i: usize) -> Vec<OsString> {
    let d = fixtures_dir();
    argv(&[
        "depth".into(),
        "-o".into(),
        tmp.join(format!("depth-{i}.txt")).display().to_string(),
        "-r".into(),
        "ref1:1-200".into(),
        d.join("dat/test_input_1_a.bam").display().to_string(),
    ])
}

fn checksum_args(tmp: &Path, i: usize) -> Vec<OsString> {
    let d = fixtures_dir();
    argv(&[
        "checksum".into(),
        "-o".into(),
        tmp.join(format!("checksum-{i}.txt")).display().to_string(),
        d.join("checksum/chk1.bam").display().to_string(),
    ])
}
