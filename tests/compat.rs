//! Byte differential vs `bcftools +split-vep` across the implemented modes:
//! `-l` (list), `-f` (extract, with `-d` and `-s worst`), and `-c` (annotate VCF).

use std::path::PathBuf;
use std::process::Command;

fn golden() -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("tests/golden/annotated.vcf.gz")
}

fn ours_bin() -> &'static str {
    env!("CARGO_BIN_EXE_rsomics-vcf-split-vep")
}

fn stdout_of(mut cmd: Command, args: &[&str], file: &PathBuf) -> Vec<u8> {
    cmd.args(args).arg(file).output().unwrap().stdout
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
