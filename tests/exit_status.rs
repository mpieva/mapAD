use std::process::Command;

use tempfile::tempdir;

#[test]
fn index_fails_on_missing_reference() {
    let temp_dir = tempdir().unwrap();
    let reference_path = temp_dir.path().join("missing.fa");

    let status = Command::new(env!("CARGO_BIN_EXE_mapad"))
        .arg("index")
        .arg("--reference")
        .arg(&reference_path)
        .status()
        .unwrap();

    assert!(!status.success());
}

#[test]
fn map_fails_on_missing_reads() {
    let temp_dir = tempdir().unwrap();
    let reads_path = temp_dir.path().join("missing.bam");
    let reference_path = temp_dir.path().join("missing.fa");
    let output_path = temp_dir.path().join("out.bam");

    let status = Command::new(env!("CARGO_BIN_EXE_mapad"))
        .arg("map")
        .arg("--reads")
        .arg(&reads_path)
        .arg("--reference")
        .arg(&reference_path)
        .arg("--output")
        .arg(&output_path)
        .args(["--library", "double_stranded"])
        .args([
            "-f", "0.5", "-d", "0.01", "-s", "0.5", "-i", "0.001", "-p", "0.03",
        ])
        .status()
        .unwrap();

    assert!(!status.success());
}
