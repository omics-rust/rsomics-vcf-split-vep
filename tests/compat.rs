//! Byte differential vs `bcftools +split-vep` across the implemented modes:
//! `-l` (list), `-f` (extract, with `-d` and `-s worst`), and `-c` (annotate VCF).

use std::path::PathBuf;
use std::process::{Command, Output};

fn golden() -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("tests/golden/annotated.vcf.gz")
}

fn ours_bin() -> &'static str {
    env!("CARGO_BIN_EXE_rsomics-vcf-split-vep")
}

fn stdout_of(mut cmd: Command, args: &[&str], file: &PathBuf) -> Vec<u8> {
    cmd.args(args).arg(file).output().unwrap().stdout
}

fn run_ours(args: &[&str], file: &PathBuf) -> Output {
    Command::new(ours_bin())
        .args(args)
        .arg(file)
        .output()
        .unwrap()
}

/// Drop bcftools' provenance/version header lines that we deliberately don't emit.
fn strip_provenance(bytes: &[u8]) -> String {
    String::from_utf8_lossy(bytes)
        .lines()
        .filter(|l| !l.starts_with("##bcftools") && !l.starts_with("##source"))
        .collect::<Vec<_>>()
        .join("\n")
}

#[test]
fn matches_bcftools_split_vep() {
    if Command::new("bcftools").arg("--version").output().is_err() {
        eprintln!("SKIP: bcftools not found");
        return;
    }
    let g = golden();

    // Extraction / list modes: byte-identical stdout.
    let exact: &[&[&str]] = &[
        &["-a", "BCSQ", "-l"],
        &["-a", "BCSQ", "-f", r"%POS\t%Consequence\n"],
        // The @-reference records (dna_change empty) are dropped by both.
        &["-a", "BCSQ", "-f", r"%POS\t%dna_change\n"],
        &[
            "-a",
            "BCSQ",
            "-d",
            "-f",
            r"%CHROM\t%POS\t%Consequence\t%gene\t%amino_acid_change\n",
        ],
        &[
            "-a",
            "BCSQ",
            "-s",
            "worst",
            "-f",
            r"%CHROM\t%POS\t%Consequence\n",
        ],
    ];
    for args in exact {
        let o = stdout_of(Command::new(ours_bin()), args, &g);
        let mut bcf = Command::new("bcftools");
        bcf.arg("+split-vep");
        let b = stdout_of(bcf, args, &g);
        assert_eq!(
            o,
            b,
            "diverged from bcftools +split-vep for {args:?}:\nours={}\nbcftools={}",
            String::from_utf8_lossy(&o),
            String::from_utf8_lossy(&b),
        );
    }

    // -c (VCF output): identical after dropping bcftools' provenance header lines.
    let annotate: &[&[&str]] = &[
        &["-a", "BCSQ", "-c", "Consequence,gene"],
        &["-a", "BCSQ", "-c", "-"],
        &["-a", "BCSQ", "-c", "0,1"],
    ];
    for args in annotate {
        let o = stdout_of(Command::new(ours_bin()), args, &g);
        let mut bcf = Command::new("bcftools");
        bcf.arg("+split-vep");
        let b = stdout_of(bcf, args, &g);
        assert_eq!(
            strip_provenance(&o),
            strip_provenance(&b),
            "diverged from bcftools +split-vep for {args:?}"
        );
    }
}

/// bcftools drops a record whose every requested CSQ subfield is empty/missing.
/// The golden has 12 `BCSQ=@POS` reference records (a single column), so
/// requesting `%dna_change` (schema index 6) drops exactly those 12, leaving 188
/// of the 200 annotated records. Expected count is bcftools 1.23.1-verified.
#[test]
fn empty_subfield_records_dropped() {
    let g = golden();
    let out = run_ours(&["-a", "BCSQ", "-f", r"%POS\t%dna_change\n"], &g);
    assert!(out.status.success());
    let text = String::from_utf8(out.stdout).unwrap();
    assert_eq!(
        text.lines().count(),
        188,
        "expected the 12 @-ref rows dropped"
    );
    for at_ref_pos in ["1578", "2341", "11304"] {
        assert!(
            !text
                .lines()
                .any(|l| l.starts_with(&format!("{at_ref_pos}\t"))),
            "@-ref record {at_ref_pos} should have been dropped"
        );
    }

    // -u does not rescue a defined-but-missing CSQ subfield: still dropped.
    let out_u = run_ours(&["-a", "BCSQ", "-u", "-f", r"%POS\t%dna_change\n"], &g);
    assert!(out_u.status.success());
    assert_eq!(
        String::from_utf8(out_u.stdout).unwrap().lines().count(),
        188
    );
}

/// A `-f` tag that is neither a CSQ subfield nor a header-defined INFO tag is a
/// fatal error without `-u`; with `-u` it renders as "." and the @-ref records
/// survive (no CSQ subfield is requested, so nothing is dropped).
#[test]
fn undefined_tag_fails_loud() {
    let g = golden();
    let bad = run_ours(&["-a", "BCSQ", "-f", r"%POS\t%nonexistent_tag\n"], &g);
    assert!(!bad.status.success(), "undefined tag must fail loud");
    let err = String::from_utf8(bad.stderr).unwrap();
    assert!(
        err.contains("no such tag defined in the VCF header: INFO/nonexistent_tag"),
        "unexpected stderr: {err}"
    );

    let ok = run_ours(&["-a", "BCSQ", "-u", "-f", r"%POS\t%nonexistent_tag\n"], &g);
    assert!(ok.status.success(), "-u should accept the undefined tag");
    let text = String::from_utf8(ok.stdout).unwrap();
    assert_eq!(
        text.lines().count(),
        200,
        "no CSQ subfield requested: no drop"
    );
    assert!(text.lines().all(|l| l.ends_with("\t.")));
    assert!(text.lines().any(|l| l.starts_with("1578\t")));
}

/// A CSQ entry with fewer columns than the Consequence field index is fatal in
/// bcftools ("Too few columns"), not silently padded. `malformed_csq.vcf` puts
/// Consequence at index 1 and a one-column entry at chr1:200.
#[test]
fn malformed_csq_too_few_columns() {
    let f = PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("tests/golden/malformed_csq.vcf");
    let out = run_ours(&["-a", "CSQ", "-f", r"%POS\t%Gene\n"], &f);
    assert!(!out.status.success(), "too-few-columns CSQ must fail loud");
    let err = String::from_utf8(out.stderr).unwrap();
    assert!(
        err.contains("Too few columns at chr1:200") && err.contains("(Consequence)"),
        "unexpected stderr: {err}"
    );
}
