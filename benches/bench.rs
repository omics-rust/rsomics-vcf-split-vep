use criterion::{Criterion, criterion_group, criterion_main};
use std::hint::black_box;
use std::path::PathBuf;
use std::process::Command;

fn bench_split_vep(c: &mut Criterion) {
    let bin = env!("CARGO_BIN_EXE_rsomics-vcf-split-vep");
    let manifest = PathBuf::from(env!("CARGO_MANIFEST_DIR"));
    let vcf = manifest.join("tests/golden/annotated.vcf.gz");
    c.bench_function("rsomics-vcf-split-vep -d -f golden", |b| {
        b.iter(|| {
            let out = Command::new(black_box(bin))
                .args([
                    "-a",
                    "BCSQ",
                    "-d",
                    "-f",
                    "%CHROM\t%POS\t%Consequence\t%gene\n",
                ])
                .arg(&vcf)
                .output()
                .unwrap();
            assert!(out.status.success());
        });
    });
}

criterion_group!(benches, bench_split_vep);
criterion_main!(benches);
