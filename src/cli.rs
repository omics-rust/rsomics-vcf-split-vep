use std::fs::File;
use std::io::{self, BufRead, BufReader, BufWriter, Read, Write};
use std::path::PathBuf;

use clap::Parser;
use flate2::read::MultiGzDecoder;
use rsomics_common::{CommonFlags, Result, RsomicsError, Tool, ToolMeta};
use rsomics_vcf_split_vep::{Annotator, Extractor, Schema};

pub const META: ToolMeta = ToolMeta {
    name: env!("CARGO_PKG_NAME"),
    version: env!("CARGO_PKG_VERSION"),
};

#[derive(Parser, Debug)]
#[command(
    name = "rsomics-vcf-split-vep",
    version,
    about = "Query and extract CSQ/BCSQ annotations — port of bcftools +split-vep"
)]
pub struct Cli {
    /// Input VCF (plain or bgzipped); "-" or omitted reads stdin.
    #[arg(default_value = "-")]
    pub input: PathBuf,

    /// INFO annotation tag to parse (default: CSQ, falling back to BCSQ).
    #[arg(short = 'a', long = "annotation")]
    pub annotation: Option<String>,

    /// List the annotation subfields as "<index>\t<name>" and exit.
    #[arg(short = 'l', long = "list")]
    pub list: bool,

    /// Output format, bcftools query-style: %CHROM, %POS, %<csq-field>, %INFO/TAG.
    #[arg(short = 'f', long = "format")]
    pub format: Option<String>,

    /// Emit one output line per transcript/allele consequence.
    #[arg(short = 'd', long = "duplicate")]
    pub duplicate: bool,

    /// Transcript selection: "worst" (most-severe transcript only) or "all" (default).
    /// Accepts a trailing :CSQ:PRN as in bcftools but only the TR part is honored.
    #[arg(short = 's', long = "select")]
    pub select: Option<String>,

    /// VCF-output mode: extract these CSQ fields into INFO tags. Comma list of
    /// names or 0-based indexes, each optionally `:Type`; `-` extracts all fields.
    #[arg(short = 'c', long = "columns")]
    pub columns: Option<String>,

    /// Print "." for undefined tags instead of an empty string.
    #[arg(short = 'u', long = "allow-undef-tags")]
    pub allow_undef: bool,

    /// Output file (default: stdout).
    #[arg(short = 'o', long = "output")]
    pub output: Option<PathBuf>,

    #[command(flatten)]
    pub common: CommonFlags,
}

impl Tool for Cli {
    fn meta() -> ToolMeta {
        META
    }

    fn common(&self) -> &CommonFlags {
        &self.common
    }

    fn execute(self) -> Result<()> {
        let mut reader = open_reader(&self.input)?;
        let sink: Box<dyn Write> = match &self.output {
            Some(p) => Box::new(File::create(p).map_err(RsomicsError::Io)?),
            None if self.common.json => Box::new(io::sink()),
            None => Box::new(io::stdout()),
        };
        let mut out = BufWriter::new(sink);

        let mut header: Vec<String> = Vec::new();
        let mut extractor: Option<Extractor> = None;
        let mut annotator: Option<Annotator> = None;
        let mut header_done = false;
        let mut line = String::new();
        loop {
            line.clear();
            if reader.read_line(&mut line).map_err(RsomicsError::Io)? == 0 {
                break;
            }
            let trimmed = line.trim_end_matches(['\n', '\r']);

            if !header_done && trimmed.starts_with('#') {
                if !trimmed.starts_with("#CHROM") {
                    header.push(trimmed.to_string());
                    continue;
                }
                // #CHROM line: header is complete — resolve the output mode.
                let schema = Schema::from_header(&header, self.annotation.as_deref())?;
                if self.list {
                    Extractor::list_fields(&schema, &mut out)?;
                    return out.flush().map_err(RsomicsError::Io);
                }
                if let Some(spec) = &self.columns {
                    let ann = Annotator::new(schema, spec, self.duplicate)?;
                    for h in &header {
                        writeln!(out, "{h}").map_err(RsomicsError::Io)?;
                    }
                    for hl in ann.header_lines() {
                        writeln!(out, "{hl}").map_err(RsomicsError::Io)?;
                    }
                    writeln!(out, "{trimmed}").map_err(RsomicsError::Io)?;
                    annotator = Some(ann);
                } else {
                    let format = self.format.as_deref().ok_or_else(|| {
                        RsomicsError::InvalidInput(
                            "one of -l/--list, -f/--format, or -c/--columns is required".into(),
                        )
                    })?;
                    let select_worst = match self
                        .select
                        .as_deref()
                        .map(|s| s.split(':').next().unwrap_or(""))
                    {
                        None | Some("" | "all") => false,
                        Some("worst") => true,
                        Some(other) => {
                            return Err(RsomicsError::InvalidInput(format!(
                                "-s transcript selection '{other}' not supported (only 'worst' or 'all')"
                            )));
                        }
                    };
                    let info_tags = rsomics_vcf_split_vep::info_tag_ids(&header);
                    extractor = Some(Extractor::new(
                        format,
                        schema,
                        &info_tags,
                        self.duplicate,
                        self.allow_undef,
                        select_worst,
                    )?);
                }
                header_done = true;
                continue;
            }

            if let Some(ann) = &annotator {
                for rec in ann.annotate(trimmed) {
                    writeln!(out, "{rec}").map_err(RsomicsError::Io)?;
                }
            } else if let Some(ex) = &extractor {
                ex.emit(trimmed, &mut out)?;
            }
        }

        if !header_done && self.list {
            let schema = Schema::from_header(&header, self.annotation.as_deref())?;
            Extractor::list_fields(&schema, &mut out)?;
        }
        out.flush().map_err(RsomicsError::Io)
    }
}

fn open_reader(path: &PathBuf) -> Result<Box<dyn BufRead>> {
    let raw: Box<dyn Read> = if path.as_os_str() == "-" {
        Box::new(io::stdin())
    } else {
        Box::new(File::open(path).map_err(RsomicsError::Io)?)
    };
    let mut br = BufReader::new(raw);
    let gz = {
        let head = br.fill_buf().map_err(RsomicsError::Io)?;
        head.len() >= 2 && head[0] == 0x1f && head[1] == 0x8b
    };
    if gz {
        Ok(Box::new(BufReader::new(MultiGzDecoder::new(br))))
    } else {
        Ok(Box::new(br))
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use clap::CommandFactory;

    #[test]
    fn cli_debug_assert() {
        Cli::command().debug_assert();
    }
}
