//! Core of `bcftools +split-vep`: parse the CSQ/BCSQ schema from the VCF header
//! and extract per-transcript subfields via a query-style format string.

use std::io::Write;

use rsomics_common::{Result, RsomicsError};

/// The annotation schema: the tag (CSQ/BCSQ) and its `Format: A|B|C` field names.
pub struct Schema {
    pub tag: String,
    pub fields: Vec<String>,
}

impl Schema {
    /// Locate the schema in the VCF header. With `prefer`, use that tag; otherwise
    /// try CSQ, then BCSQ (matching bcftools' default).
    pub fn from_header(header: &[String], prefer: Option<&str>) -> Result<Schema> {
        let candidates: Vec<&str> = match prefer {
            Some(t) => vec![t],
            None => vec!["CSQ", "BCSQ"],
        };
        for tag in candidates {
            if let Some(s) = parse_one(header, tag) {
                return Ok(s);
            }
        }
        Err(RsomicsError::InvalidInput(
            "no CSQ/BCSQ INFO header with a 'Format:' description found".into(),
        ))
    }

    fn index_of(&self, name: &str) -> Option<usize> {
        self.fields.iter().position(|f| f == name)
    }
}

fn parse_one(header: &[String], tag: &str) -> Option<Schema> {
    let needle = format!("##INFO=<ID={tag},");
    let line = header.iter().find(|l| l.starts_with(&needle))?;
    let start = line.find("Format: ")? + "Format: ".len();
    let rest = &line[start..];
    let end = rest.find('"').unwrap_or(rest.len());
    let fields = rest[..end]
        .split('|')
        .map(|s| s.trim().to_string())
        .collect();
    Some(Schema {
        tag: tag.to_string(),
        fields,
    })
}

enum Token {
    Lit(String),
    Field(String),
}

fn parse_format(fmt: &str) -> Vec<Token> {
    let mut tokens = Vec::new();
    let mut lit = String::new();
    let bytes = fmt.as_bytes();
    let mut i = 0;
    while i < bytes.len() {
        match bytes[i] {
            b'\\' if i + 1 < bytes.len() => {
                lit.push(match bytes[i + 1] {
                    b't' => '\t',
                    b'n' => '\n',
                    other => other as char,
                });
                i += 2;
            }
            b'%' => {
                if !lit.is_empty() {
                    tokens.push(Token::Lit(std::mem::take(&mut lit)));
                }
                i += 1;
                let s = i;
                while i < bytes.len()
                    && (bytes[i].is_ascii_alphanumeric() || matches!(bytes[i], b'_' | b'/'))
                {
                    i += 1;
                }
                tokens.push(Token::Field(fmt[s..i].to_string()));
            }
            other => {
                lit.push(other as char);
                i += 1;
            }
        }
    }
    if !lit.is_empty() {
        tokens.push(Token::Lit(lit));
    }
    tokens
}

const VCF_COLS: [&str; 7] = ["CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER"];

/// bcftools' default consequence severity scale, ascending. A consequence's
/// rank is the highest index any of these substrings matches within it.
const SEVERITY: [&[&str]; 20] = [
    &["intergenic"],
    &["feature_truncation", "feature_elongation"],
    &["regulatory"],
    &["TF_binding_site", "TFBS"],
    &["downstream", "upstream"],
    &["non_coding_transcript", "non_coding"],
    &["intron", "NMD_transcript"],
    &["non_coding_transcript_exon"],
    &["5_prime_utr", "3_prime_utr"],
    &["coding_sequence", "mature_miRNA"],
    &["stop_retained", "start_retained", "synonymous"],
    &["incomplete_terminal_codon"],
    &["splice_region"],
    &["missense", "inframe", "protein_altering"],
    &["transcript_amplification"],
    &["exon_loss"],
    &["disruptive"],
    &["start_lost", "stop_lost", "stop_gained", "frameshift"],
    &["splice_acceptor", "splice_donor"],
    &["transcript_ablation"],
];

fn severity_rank(consequence: &str) -> i32 {
    let mut best = -1i32;
    for (rank, terms) in SEVERITY.iter().enumerate() {
        if terms.iter().any(|t| consequence.contains(t)) {
            best = rank as i32;
        }
    }
    best
}

/// Renders records to a query-style stream, extracting CSQ subfields.
pub struct Extractor {
    tokens: Vec<Token>,
    schema: Schema,
    duplicate: bool,
    allow_undef: bool,
    /// `-s worst`: keep only the single most-severe transcript per record.
    select_worst: bool,
    consequence_idx: Option<usize>,
}

impl Extractor {
    pub fn new(
        format: &str,
        schema: Schema,
        duplicate: bool,
        allow_undef: bool,
        select_worst: bool,
    ) -> Self {
        let consequence_idx = schema.index_of("Consequence");
        Self {
            tokens: parse_format(format),
            schema,
            duplicate,
            allow_undef,
            select_worst,
            consequence_idx,
        }
    }

    /// List the schema fields as `<index>\t<name>` (the `-l` mode).
    pub fn list_fields<W: Write>(schema: &Schema, out: &mut W) -> Result<()> {
        for (i, f) in schema.fields.iter().enumerate() {
            writeln!(out, "{i}\t{f}").map_err(RsomicsError::Io)?;
        }
        Ok(())
    }

    /// Emit output lines for one VCF data record. Records lacking the annotation
    /// produce no output (bcftools' default with `-f`).
    pub fn emit<W: Write>(&self, line: &str, out: &mut W) -> Result<()> {
        let cols: Vec<&str> = line.split('\t').collect();
        if cols.len() < 8 {
            return Ok(());
        }
        let entries = annotation_entries(cols[7], &self.schema.tag);
        if entries.is_empty() {
            return Ok(());
        }
        // Each entry's pipe-split subfields, padded to schema width with ".".
        let mut split: Vec<Vec<&str>> = entries
            .iter()
            .map(|e| e.split('|').collect::<Vec<_>>())
            .collect();

        // -s worst: keep only the single most-severe transcript (first on ties).
        if self.select_worst
            && split.len() > 1
            && let Some(ci) = self.consequence_idx
        {
            let mut best_i = 0;
            let mut best_rank = i32::MIN;
            for (i, e) in split.iter().enumerate() {
                let r = severity_rank(e.get(ci).copied().unwrap_or(""));
                if r > best_rank {
                    best_rank = r;
                    best_i = i;
                }
            }
            split = vec![split.swap_remove(best_i)];
        }

        if self.duplicate {
            for entry in &split {
                self.write_line(&cols, Some(entry), &split, out)?;
            }
        } else {
            self.write_line(&cols, None, &split, out)?;
        }
        Ok(())
    }

    fn write_line<W: Write>(
        &self,
        cols: &[&str],
        current: Option<&Vec<&str>>,
        all: &[Vec<&str>],
        out: &mut W,
    ) -> Result<()> {
        let mut buf = String::new();
        for tok in &self.tokens {
            match tok {
                Token::Lit(s) => buf.push_str(s),
                Token::Field(name) => buf.push_str(&self.resolve(cols, current, all, name)),
            }
        }
        out.write_all(buf.as_bytes()).map_err(RsomicsError::Io)
    }

    fn resolve(
        &self,
        cols: &[&str],
        current: Option<&Vec<&str>>,
        all: &[Vec<&str>],
        name: &str,
    ) -> String {
        if let Some(idx) = VCF_COLS.iter().position(|c| *c == name) {
            return cols.get(idx).copied().unwrap_or(".").to_string();
        }
        if let Some(key) = name.strip_prefix("INFO/") {
            return info_value(cols[7], key).unwrap_or_else(|| ".".to_string());
        }
        if let Some(fi) = self.schema.index_of(name) {
            return self.csq_value(current, all, fi);
        }
        if let Some(v) = info_value(cols[7], name) {
            return v;
        }
        if self.allow_undef {
            ".".to_string()
        } else {
            String::new()
        }
    }

    /// A CSQ subfield: in `-d` mode the current entry's value, otherwise the
    /// comma-joined values across all entries (matching bcftools).
    fn csq_value(&self, current: Option<&Vec<&str>>, all: &[Vec<&str>], fi: usize) -> String {
        fn field_at<'a>(e: &[&'a str], fi: usize) -> &'a str {
            match e.get(fi) {
                Some(s) if !s.is_empty() => s,
                _ => ".",
            }
        }
        match current {
            Some(e) => field_at(e, fi).to_string(),
            None => all
                .iter()
                .map(|e| field_at(e, fi))
                .collect::<Vec<_>>()
                .join(","),
        }
    }
}

/// Split an INFO column's `TAG=a|b,c|d` value into per-transcript entries.
fn annotation_entries(info: &str, tag: &str) -> Vec<String> {
    let prefix = format!("{tag}=");
    for field in info.split(';') {
        if let Some(v) = field.strip_prefix(&prefix) {
            return v.split(',').map(|s| s.to_string()).collect();
        }
    }
    Vec::new()
}

fn info_value(info: &str, key: &str) -> Option<String> {
    let prefix = format!("{key}=");
    for field in info.split(';') {
        if field == key {
            return Some(".".to_string());
        }
        if let Some(v) = field.strip_prefix(&prefix) {
            return Some(v.to_string());
        }
    }
    None
}

/// VCF-output mode (`-c`): copy the VCF through, injecting an INFO header and
/// appending a `name=value` INFO tag per selected CSQ subfield.
pub struct Annotator {
    schema: Schema,
    /// (csq subfield index, output INFO tag name, VCF Type).
    fields: Vec<(usize, String, String)>,
    duplicate: bool,
}

impl Annotator {
    /// `spec` is a comma list of `NAME[:TYPE]` / `INDEX[:TYPE]`, or `-` for all fields.
    pub fn new(schema: Schema, spec: &str, duplicate: bool) -> Result<Self> {
        let mut fields = Vec::new();
        if spec == "-" {
            for (i, name) in schema.fields.iter().enumerate() {
                fields.push((i, name.clone(), "String".to_string()));
            }
        } else {
            for item in spec.split(',') {
                let mut parts = item.splitn(2, ':');
                let key = parts.next().unwrap_or("");
                let ty = match parts.next() {
                    Some("Integer" | "Int") => "Integer",
                    Some("Float" | "Real") => "Float",
                    _ => "String",
                }
                .to_string();
                let idx = if key.bytes().all(|b| b.is_ascii_digit()) && !key.is_empty() {
                    key.parse::<usize>()
                        .ok()
                        .filter(|i| *i < schema.fields.len())
                } else {
                    schema.index_of(key)
                };
                let idx = idx.ok_or_else(|| {
                    RsomicsError::InvalidInput(format!("unknown CSQ field for -c: {key}"))
                })?;
                fields.push((idx, schema.fields[idx].clone(), ty));
            }
        }
        Ok(Self {
            schema,
            fields,
            duplicate,
        })
    }

    /// The `##INFO` header lines to inject (one per selected field).
    pub fn header_lines(&self) -> Vec<String> {
        self.fields
            .iter()
            .map(|(_, name, ty)| {
                format!(
                    "##INFO=<ID={name},Number=.,Type={ty},Description=\"The {name} field from INFO/{}\">",
                    self.schema.tag
                )
            })
            .collect()
    }

    /// Annotate one VCF data record, returning the output line(s).
    pub fn annotate(&self, line: &str) -> Vec<String> {
        let cols: Vec<&str> = line.split('\t').collect();
        if cols.len() < 8 {
            return vec![line.to_string()];
        }
        let entries = annotation_entries(cols[7], &self.schema.tag);
        if entries.is_empty() {
            return vec![line.to_string()];
        }
        let split: Vec<Vec<&str>> = entries.iter().map(|e| e.split('|').collect()).collect();
        let field_at = |e: &[&str], fi: usize| -> String {
            e.get(fi)
                .copied()
                .filter(|s| !s.is_empty())
                .unwrap_or(".")
                .to_string()
        };

        if self.duplicate {
            split
                .iter()
                .map(|e| {
                    let adds: Vec<String> = self
                        .fields
                        .iter()
                        .map(|(fi, name, _)| format!("{name}={}", field_at(e, *fi)))
                        .collect();
                    rebuild(&cols, &adds)
                })
                .collect()
        } else {
            let adds: Vec<String> = self
                .fields
                .iter()
                .map(|(fi, name, _)| {
                    let joined = split
                        .iter()
                        .map(|e| field_at(e, *fi))
                        .collect::<Vec<_>>()
                        .join(",");
                    format!("{name}={joined}")
                })
                .collect();
            vec![rebuild(&cols, &adds)]
        }
    }
}

/// Re-emit a VCF record with extra INFO tags appended.
fn rebuild(cols: &[&str], adds: &[String]) -> String {
    let mut out: Vec<String> = cols.iter().map(|s| s.to_string()).collect();
    let info = &mut out[7];
    let extra = adds.join(";");
    if info == "." {
        *info = extra;
    } else {
        info.push(';');
        info.push_str(&extra);
    }
    out.join("\t")
}
