# rsomics-vcf-split-vep

Query and extract structured consequence annotations (VEP `CSQ` or bcftools/csq
`BCSQ`) from a VCF — a Rust port of `bcftools +split-vep`. Parses the `Format:
A|B|C` field list from the annotation's INFO header, then either lists the
fields, extracts them to a TSV with a query-style format string, or writes them
back into the VCF as new INFO tags.

## Install

```sh
cargo install rsomics-vcf-split-vep
```

## Usage

```sh
rsomics-vcf-split-vep -l in.vcf.gz                                  # list CSQ/BCSQ subfields
rsomics-vcf-split-vep -f '%CHROM\t%POS\t%Consequence\t%gene\n' in.vcf.gz   # extract to TSV
rsomics-vcf-split-vep -d -f '%POS\t%Consequence\n' in.vcf.gz        # one line per transcript
rsomics-vcf-split-vep -s worst -f '%POS\t%Consequence\n' in.vcf.gz  # worst transcript only
rsomics-vcf-split-vep -c Consequence,gene in.vcf.gz                 # add INFO tags (VCF out)
rsomics-vcf-split-vep -a BCSQ -c - in.vcf.gz                        # all fields → INFO
```

| flag | meaning |
|---|---|
| `-a, --annotation` | INFO tag to parse (default: CSQ, falling back to BCSQ) |
| `-l, --list` | print `<index>\t<name>` per subfield and exit |
| `-f, --format` | query-style output (`%CHROM`, `%POS`, `%<csq-field>`, `%INFO/TAG`) |
| `-d, --duplicate` | one output line per transcript/allele consequence |
| `-s, --select` | transcript selection: `worst` (most-severe) or `all` (default) |
| `-c, --columns` | VCF output: extract fields (names or 0-based indexes, `:Type`, `-` for all) into INFO |
| `-u, --allow-undef-tags` | print `.` for undefined tags |

Plain and bgzipped VCF input are auto-detected. The default consequence severity
scale (used by `-s worst`) matches bcftools'.

Not yet implemented (rare/advanced): `-g` gene lists, `-S FILE` custom severity
scales, regex/expression transcript selection, `-A` all-fields-in-format, and
`-H` header for `-f`.

## Origin

Independent Rust reimplementation of `bcftools +split-vep`. bcftools is
MIT-licensed; the `Format:`-header parsing, the `0xidx<TAB>name` list output, the
query-format semantics (missing-field padding with `.`, comma-joining across
transcripts, drop-without-`-f`), the `-c` INFO-injection shape, and the default
severity scale were taken from bcftools' MIT-licensed source. Output is verified
byte-for-byte against `bcftools +split-vep` across `-l`/`-f`/`-d`/`-s`/`-c` in
`tests/compat.rs`.

License: MIT OR Apache-2.0.
Upstream credit: [bcftools](https://github.com/samtools/bcftools) (MIT).
