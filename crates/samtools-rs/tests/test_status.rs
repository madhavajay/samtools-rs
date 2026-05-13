use std::collections::BTreeSet;
use std::path::PathBuf;

fn repo_root() -> PathBuf {
    let manifest = env!("CARGO_MANIFEST_DIR");
    PathBuf::from(manifest)
        .parent()
        .unwrap()
        .parent()
        .unwrap()
        .to_path_buf()
}

#[test]
fn test_status_lists_all_top_level_test_pl_groups() {
    let root = repo_root();
    let test_pl =
        std::fs::read_to_string(root.join("samtools").join("test").join("test.pl")).unwrap();
    let status = std::fs::read_to_string(root.join("docs").join("test-status.md")).unwrap();

    let mut groups = BTreeSet::new();

    for line in test_pl.lines() {
        if line.starts_with("print \"\\nNumber of tests:") {
            break;
        }

        let Some(rest) = line.strip_prefix("test_") else {
            continue;
        };
        let Some((name, _)) = rest.split_once('(') else {
            continue;
        };
        groups.insert(format!("test_{name}"));
    }

    let missing: Vec<_> = groups
        .iter()
        .filter(|group| !status.contains(&format!("| `{group}` |")))
        .cloned()
        .collect();

    assert!(
        missing.is_empty(),
        "docs/test-status.md is missing upstream test.pl groups: {missing:?}"
    );
}
