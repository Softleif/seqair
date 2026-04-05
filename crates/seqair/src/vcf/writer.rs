//! VCF text writer. Serializes [`VcfRecord`]s as tab-delimited VCF lines,
//! optionally BGZF-compressed with TBI index co-production.

use super::error::VcfError;
use super::header::VcfHeader;
use super::index_builder::IndexBuilder;
use super::record::{Filters, InfoValue, SampleValue, VcfRecord};
use crate::bam::bgzf_writer::BgzfWriter;
use std::io::Write;
use std::sync::Arc;

/// Output mode for the VCF writer.
enum Output<W: Write> {
    Plain(W),
    Bgzf { bgzf: BgzfWriter<W>, index: IndexBuilder },
}

// r[impl vcf_writer.header_first]
// r[impl vcf_writer.tab_delimited]
// r[impl vcf_writer.buffer_reuse]
/// VCF text writer with optional BGZF compression and TBI index co-production.
pub struct VcfWriter<W: Write> {
    output: Output<W>,
    header: Arc<VcfHeader>,
    buf: Vec<u8>,
    header_written: bool,
}

impl<W: Write> VcfWriter<W> {
    /// Create a plain (uncompressed) VCF writer.
    pub fn new(inner: W, header: Arc<VcfHeader>) -> Self {
        Self {
            output: Output::Plain(inner),
            header,
            buf: Vec::with_capacity(4096),
            header_written: false,
        }
    }

    /// Create a BGZF-compressed VCF writer with TBI index co-production.
    pub fn bgzf(inner: W, header: Arc<VcfHeader>) -> Self {
        let n_refs = header.contigs().len();
        let bgzf = BgzfWriter::new(inner);
        Self {
            output: Output::Bgzf {
                bgzf,
                index: IndexBuilder::tbi(n_refs, crate::bam::bgzf::VirtualOffset(0)),
            },
            header,
            buf: Vec::with_capacity(4096),
            header_written: false,
        }
    }

    // r[impl vcf_writer.header_first]
    /// Write the VCF header. Must be called exactly once before any `write_record`.
    pub fn write_header(&mut self) -> Result<(), VcfError> {
        let text = self.header.to_vcf_text();
        match &mut self.output {
            Output::Plain(w) => w.write_all(text.as_bytes()).map_err(VcfError::Io)?,
            Output::Bgzf { bgzf, .. } => bgzf.write_all(text.as_bytes())?,
        }
        self.header_written = true;

        // Update index with header end offset for BGZF mode
        if let Output::Bgzf { bgzf, index } = &mut self.output {
            *index = IndexBuilder::tbi(self.header.contigs().len(), bgzf.virtual_offset());
        }
        Ok(())
    }

    // r[impl vcf_writer.tab_delimited]
    // r[impl vcf_writer.buffer_reuse]
    // r[impl vcf_writer.validation]
    /// Write a single VCF record.
    pub fn write_record(&mut self, record: &VcfRecord) -> Result<(), VcfError> {
        if !self.header_written {
            return Err(VcfError::HeaderNotWritten);
        }

        // r[impl vcf_writer.buffer_reuse]
        self.buf.clear();

        // CHROM
        self.buf.extend_from_slice(record.contig.as_bytes());
        self.buf.push(b'\t');

        // POS (1-based)
        let mut itoa_buf = itoa::Buffer::new();
        self.buf.extend_from_slice(itoa_buf.format(record.pos.get()).as_bytes());
        self.buf.push(b'\t');

        // ID
        match &record.id {
            Some(id) => self.buf.extend_from_slice(id.as_bytes()),
            None => self.buf.push(b'.'),
        }
        self.buf.push(b'\t');

        // REF — zero-alloc, writes directly into buffer
        record.alleles.write_ref_into(&mut self.buf);
        self.buf.push(b'\t');

        // ALT — zero-alloc, writes directly into buffer
        let n_alts = record.alleles.write_alts_into(&mut self.buf);
        if n_alts == 0 {
            self.buf.push(b'.');
        }
        self.buf.push(b'\t');

        // QUAL
        // r[impl vcf_writer.float_precision]
        match record.qual {
            Some(q) => {
                let mut ryu_buf = ryu::Buffer::new();
                self.buf.extend_from_slice(ryu_buf.format(q).as_bytes());
            }
            None => self.buf.push(b'.'),
        }
        self.buf.push(b'\t');

        // FILTER
        // r[impl vcf_writer.missing_dot]
        match &record.filters {
            Filters::Pass => self.buf.extend_from_slice(b"PASS"),
            Filters::Failed(ids) => {
                for (i, id) in ids.iter().enumerate() {
                    if i > 0 {
                        self.buf.push(b';');
                    }
                    self.buf.extend_from_slice(id.as_bytes());
                }
            }
            Filters::NotApplied => self.buf.push(b'.'),
        }
        self.buf.push(b'\t');

        // INFO
        // r[impl vcf_writer.info_serialization]
        if record.info.is_empty() {
            self.buf.push(b'.');
        } else {
            for (i, (key, value)) in record.info.iter().enumerate() {
                if i > 0 {
                    self.buf.push(b';');
                }
                self.buf.extend_from_slice(key.as_bytes());
                match value {
                    InfoValue::Flag => {} // key only, no =
                    _ => {
                        self.buf.push(b'=');
                        write_info_value(&mut self.buf, value);
                    }
                }
            }
        }

        // FORMAT + samples
        // r[impl vcf_writer.format_serialization]
        if !record.samples.format_keys.is_empty() && !record.samples.values.is_empty() {
            self.buf.push(b'\t');
            for (i, key) in record.samples.format_keys.iter().enumerate() {
                if i > 0 {
                    self.buf.push(b':');
                }
                self.buf.extend_from_slice(key.as_bytes());
            }

            for sample_values in &record.samples.values {
                self.buf.push(b'\t');
                for (i, val) in sample_values.iter().enumerate() {
                    if i > 0 {
                        self.buf.push(b':');
                    }
                    write_sample_value(&mut self.buf, val);
                }
            }
        }

        self.buf.push(b'\n');

        // Write the serialized record to the output
        match &mut self.output {
            Output::Plain(w) => {
                w.write_all(&self.buf).map_err(VcfError::Io)?;
            }
            Output::Bgzf { bgzf, index } => {
                bgzf.flush_if_needed(self.buf.len())?;
                bgzf.write_all(&self.buf)?;

                let tid = self.header.contig_id(&record.contig)? as i32;
                let beg = record.pos.to_zero_based().get() as u64;
                let end = beg.saturating_add(record.alleles.rlen() as u64);
                index
                    .push(tid, beg, end, bgzf.virtual_offset())
                    .map_err(|e| VcfError::Io(std::io::Error::other(e.to_string())))?;
            }
        }

        Ok(())
    }

    // r[impl vcf_writer.finish]
    /// Flush and finalize the writer.
    /// For BGZF: writes EOF marker. Returns the index builder if present.
    pub fn finish(self) -> Result<Option<IndexBuilder>, VcfError> {
        match self.output {
            Output::Plain(mut w) => {
                w.flush().map_err(VcfError::Io)?;
                Ok(None)
            }
            Output::Bgzf { bgzf, mut index } => {
                let voff = bgzf.virtual_offset();
                index
                    .finish(voff)
                    .map_err(|e| VcfError::Io(std::io::Error::other(e.to_string())))?;
                bgzf.finish()?;
                Ok(Some(index))
            }
        }
    }
}

// r[impl vcf_writer.integer_format]
fn write_info_value(buf: &mut Vec<u8>, value: &InfoValue) {
    let mut itoa_buf = itoa::Buffer::new();
    let mut ryu_buf = ryu::Buffer::new();
    match value {
        InfoValue::Integer(v) => buf.extend_from_slice(itoa_buf.format(*v).as_bytes()),
        InfoValue::Float(v) => buf.extend_from_slice(ryu_buf.format(*v).as_bytes()),
        InfoValue::Flag => {} // handled at call site
        InfoValue::String(s) => percent_encode_into(buf, s.as_bytes()),
        InfoValue::IntegerArray(arr) => {
            for (i, v) in arr.iter().enumerate() {
                if i > 0 {
                    buf.push(b',');
                }
                match v {
                    Some(n) => buf.extend_from_slice(itoa_buf.format(*n).as_bytes()),
                    None => buf.push(b'.'),
                }
            }
        }
        InfoValue::FloatArray(arr) => {
            for (i, v) in arr.iter().enumerate() {
                if i > 0 {
                    buf.push(b',');
                }
                match v {
                    Some(f) => buf.extend_from_slice(ryu_buf.format(*f).as_bytes()),
                    None => buf.push(b'.'),
                }
            }
        }
        InfoValue::StringArray(arr) => {
            for (i, v) in arr.iter().enumerate() {
                if i > 0 {
                    buf.push(b',');
                }
                match v {
                    Some(s) => percent_encode_into(buf, s.as_bytes()),
                    None => buf.push(b'.'),
                }
            }
        }
    }
}

// r[impl vcf_writer.percent_encoding]
/// Percent-encode special characters in VCF field values.
fn percent_encode_into(buf: &mut Vec<u8>, data: &[u8]) {
    for &b in data {
        match b {
            b':' => buf.extend_from_slice(b"%3A"),
            b';' => buf.extend_from_slice(b"%3B"),
            b'=' => buf.extend_from_slice(b"%3D"),
            b'%' => buf.extend_from_slice(b"%25"),
            b',' => buf.extend_from_slice(b"%2C"),
            b'\t' => buf.extend_from_slice(b"%09"),
            b'\n' => buf.extend_from_slice(b"%0A"),
            b'\r' => buf.extend_from_slice(b"%0D"),
            _ => buf.push(b),
        }
    }
}

// r[impl vcf_writer.genotype_serialization]
fn write_sample_value(buf: &mut Vec<u8>, value: &SampleValue) {
    let mut itoa_buf = itoa::Buffer::new();
    let mut ryu_buf = ryu::Buffer::new();
    match value {
        SampleValue::Missing => buf.push(b'.'),
        SampleValue::Integer(v) => buf.extend_from_slice(itoa_buf.format(*v).as_bytes()),
        SampleValue::Float(v) => buf.extend_from_slice(ryu_buf.format(*v).as_bytes()),
        SampleValue::String(s) => percent_encode_into(buf, s.as_bytes()),
        SampleValue::Genotype(gt) => {
            for (i, allele) in gt.alleles.iter().enumerate() {
                if i > 0 {
                    // Separator: phased[i-1] determines | or /
                    let phased = gt.phased.get(i.saturating_sub(1)).copied().unwrap_or(false);
                    buf.push(if phased { b'|' } else { b'/' });
                }
                match allele {
                    Some(idx) => buf.extend_from_slice(itoa_buf.format(*idx).as_bytes()),
                    None => buf.push(b'.'),
                }
            }
        }
        SampleValue::IntegerArray(arr) => {
            for (i, v) in arr.iter().enumerate() {
                if i > 0 {
                    buf.push(b',');
                }
                match v {
                    Some(n) => buf.extend_from_slice(itoa_buf.format(*n).as_bytes()),
                    None => buf.push(b'.'),
                }
            }
        }
        SampleValue::FloatArray(arr) => {
            for (i, v) in arr.iter().enumerate() {
                if i > 0 {
                    buf.push(b',');
                }
                match v {
                    Some(f) => buf.extend_from_slice(ryu_buf.format(*f).as_bytes()),
                    None => buf.push(b'.'),
                }
            }
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::vcf::alleles::Alleles;
    use crate::vcf::header::{ContigDef, FormatDef, InfoDef, Number, ValueType};
    use crate::vcf::record::{Genotype, SampleValue, VcfRecordBuilder};
    use seqair_types::{Base, One, Pos, SmolStr};

    fn test_header() -> Arc<VcfHeader> {
        Arc::new(
            VcfHeader::builder()
                .add_contig("chr1", ContigDef { length: Some(1000) })
                .unwrap()
                .add_info(
                    "DP",
                    InfoDef {
                        number: Number::Count(1),
                        typ: ValueType::Integer,
                        description: SmolStr::from("Depth"),
                    },
                )
                .unwrap()
                .add_info(
                    "DB",
                    InfoDef {
                        number: Number::Count(0),
                        typ: ValueType::Flag,
                        description: SmolStr::from("dbSNP"),
                    },
                )
                .unwrap()
                .add_format(
                    "GT",
                    FormatDef {
                        number: Number::Count(1),
                        typ: ValueType::String,
                        description: SmolStr::from("Genotype"),
                    },
                )
                .unwrap()
                .add_format(
                    "DP",
                    FormatDef {
                        number: Number::Count(1),
                        typ: ValueType::Integer,
                        description: SmolStr::from("Read Depth"),
                    },
                )
                .unwrap()
                .add_sample("S1")
                .unwrap()
                .build()
                .unwrap(),
        )
    }

    /// Write records and return the output as a String.
    fn write_vcf(header: &Arc<VcfHeader>, records: &[VcfRecord]) -> String {
        let mut output = Vec::new();
        {
            let mut writer = VcfWriter::new(&mut output, header.clone());
            writer.write_header().unwrap();
            for record in records {
                writer.write_record(record).unwrap();
            }
            writer.finish().unwrap();
        }
        String::from_utf8(output).unwrap()
    }

    // r[verify vcf_writer.header_first]
    #[test]
    fn write_header_before_records() {
        let header = test_header();
        let text = write_vcf(&header, &[]);
        assert!(text.starts_with("##fileformat=VCFv4.3\n"));
        assert!(text.contains("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO"));
    }

    /// Get the last data line fields from VCF text output.
    fn last_line_fields(text: &str) -> Vec<String> {
        text.lines().last().unwrap().split('\t').map(String::from).collect()
    }

    // r[verify vcf_writer.tab_delimited]
    // r[verify vcf_writer.missing_dot]
    #[test]
    fn write_simple_snv() {
        let header = test_header();
        let record = VcfRecordBuilder::new(
            "chr1",
            Pos::<One>::new(100).unwrap(),
            Alleles::snv(Base::A, Base::T).unwrap(),
        )
        .qual(30.0)
        .filter_pass()
        .info_integer("DP", 50)
        .build(&header)
        .unwrap();

        let text = write_vcf(&header, &[record]);
        let f = last_line_fields(&text);
        assert_eq!(f[0], "chr1");
        assert_eq!(f[1], "100");
        assert_eq!(f[2], ".");
        assert_eq!(f[3], "A");
        assert_eq!(f[4], "T");
        assert_eq!(f[5], "30.0");
        assert_eq!(f[6], "PASS");
        assert_eq!(f[7], "DP=50");
    }

    // r[verify vcf_writer.info_serialization]
    #[test]
    fn write_info_flag() {
        let header = test_header();
        let record =
            VcfRecordBuilder::new("chr1", Pos::<One>::new(1).unwrap(), Alleles::reference(Base::A))
                .info_flag("DB")
                .build(&header)
                .unwrap();

        let text = write_vcf(&header, &[record]);
        let f = last_line_fields(&text);
        assert_eq!(f[7], "DB");
    }

    // r[verify vcf_writer.genotype_serialization]
    #[test]
    fn write_genotype() {
        let header = test_header();
        let record = VcfRecordBuilder::new(
            "chr1",
            Pos::<One>::new(100).unwrap(),
            Alleles::snv(Base::A, Base::T).unwrap(),
        )
        .filter_pass()
        .format_keys(&["GT", "DP"])
        .add_sample(vec![SampleValue::Genotype(Genotype::unphased(0, 1)), SampleValue::Integer(30)])
        .build(&header)
        .unwrap();

        let text = write_vcf(&header, &[record]);
        let f = last_line_fields(&text);
        assert_eq!(f[8], "GT:DP");
        assert_eq!(f[9], "0/1:30");
    }

    // r[verify vcf_writer.genotype_serialization]
    #[test]
    fn write_phased_genotype() {
        let header = test_header();
        let record = VcfRecordBuilder::new(
            "chr1",
            Pos::<One>::new(100).unwrap(),
            Alleles::snv(Base::A, Base::T).unwrap(),
        )
        .filter_pass()
        .format_keys(&["GT"])
        .add_sample(vec![SampleValue::Genotype(Genotype::phased_diploid(0, 1))])
        .build(&header)
        .unwrap();

        let text = write_vcf(&header, &[record]);
        let f = last_line_fields(&text);
        assert_eq!(f[9], "0|1");
    }

    // r[verify vcf_writer.missing_dot]
    #[test]
    fn write_reference_site() {
        let header = test_header();
        let record =
            VcfRecordBuilder::new("chr1", Pos::<One>::new(1).unwrap(), Alleles::reference(Base::G))
                .build(&header)
                .unwrap();

        let text = write_vcf(&header, &[record]);
        let f = last_line_fields(&text);
        assert_eq!(f[2], ".");
        assert_eq!(f[4], ".");
        assert_eq!(f[5], ".");
        assert_eq!(f[6], ".");
        assert_eq!(f[7], ".");
    }

    // r[verify vcf_writer.output_formats]
    #[test]
    fn bgzf_writer_produces_valid_output() {
        let header = test_header();
        let record = VcfRecordBuilder::new(
            "chr1",
            Pos::<One>::new(100).unwrap(),
            Alleles::snv(Base::A, Base::T).unwrap(),
        )
        .filter_pass()
        .build(&header)
        .unwrap();

        let mut output = Vec::new();
        let mut writer = VcfWriter::bgzf(&mut output, header);
        writer.write_header().unwrap();
        writer.write_record(&record).unwrap();
        let index = writer.finish().unwrap();

        // Should have produced BGZF output
        assert!(output.len() > 28);
        assert_eq!(output[0], 0x1f);
        assert_eq!(output[1], 0x8b);

        // Should have an index
        assert!(index.is_some());
    }

    // r[verify vcf_writer.percent_encoding]
    #[test]
    fn percent_encoding_special_chars() {
        let mut buf = Vec::new();
        percent_encode_into(&mut buf, b"hello:world;key=val%done,x\ty\nz\rend");
        let result = String::from_utf8(buf).unwrap();
        assert_eq!(result, "hello%3Aworld%3Bkey%3Dval%25done%2Cx%09y%0Az%0Dend");
    }

    // r[verify vcf_writer.percent_encoding]
    #[test]
    fn info_string_value_is_percent_encoded() {
        let header = Arc::new(
            VcfHeader::builder()
                .add_contig("chr1", ContigDef { length: Some(1000) })
                .unwrap()
                .add_info(
                    "ANN",
                    InfoDef {
                        number: Number::Count(1),
                        typ: ValueType::String,
                        description: SmolStr::from("Annotation"),
                    },
                )
                .unwrap()
                .build()
                .unwrap(),
        );
        let record =
            VcfRecordBuilder::new("chr1", Pos::<One>::new(1).unwrap(), Alleles::reference(Base::A))
                .info_string("ANN", "val:with;special=chars")
                .build(&header)
                .unwrap();

        let text = write_vcf(&header, &[record]);
        let f = last_line_fields(&text);
        assert_eq!(f[7], "ANN=val%3Awith%3Bspecial%3Dchars");
    }
}
