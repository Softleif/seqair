use crate::{
    bam::{BamHeader, IndexedBamReader, PileupEngine, RecordStore},
    cram::reader::IndexedCramReader,
    fasta::IndexedFastaReader,
    reader::{ReaderError, indexed::CursorReader},
    sam::reader::IndexedSamReader,
};
use seqair_types::{Base, Pos, Zero};
use std::rc::Rc;

/// Alignment reader backend — indexed (BGZF+index) or plain (linear scan).
enum AlignmentBackend {
    Indexed(CursorReader),
    PlainSam(IndexedSamReader<std::io::Cursor<Vec<u8>>>),
}

/// In-memory reader bundle for fuzzing. No file I/O — all data from byte slices.
/// Uses the same `IndexedReader` enum as production code, just with `Cursor` I/O.
pub struct FuzzReaders {
    alignment: AlignmentBackend,
    fasta: IndexedFastaReader<std::io::Cursor<Vec<u8>>>,
    store: RecordStore,
    fasta_buf: Vec<u8>,
}

impl FuzzReaders {
    /// Build BAM-based readers from raw bytes: BAM + BAI + FASTA.gz + FAI + GZI.
    pub fn from_bam_bytes(
        bam_data: Vec<u8>,
        bai_data: &[u8],
        fasta_gz_data: Vec<u8>,
        fai_contents: &str,
        gzi_data: &[u8],
    ) -> Result<Self, ReaderError> {
        let bam = IndexedBamReader::from_bytes(bam_data, bai_data)?;
        let fasta = IndexedFastaReader::from_bytes(fasta_gz_data, fai_contents, gzi_data)
            .map_err(|source| ReaderError::FastaOpen { source })?;
        Ok(FuzzReaders {
            alignment: AlignmentBackend::Indexed(CursorReader::Bam(bam)),
            fasta,
            store: RecordStore::default(),
            fasta_buf: Vec::new(),
        })
    }

    /// Build CRAM-based readers from raw bytes.
    pub fn from_cram_bytes(
        cram_data: Vec<u8>,
        crai_text: &[u8],
        fasta_gz_data: Vec<u8>,
        fai_contents: &str,
        gzi_data: &[u8],
    ) -> Result<Self, ReaderError> {
        let cram = IndexedCramReader::from_bytes(
            cram_data,
            crai_text,
            fasta_gz_data.clone(),
            fai_contents,
            gzi_data,
        )?;
        let fasta = IndexedFastaReader::from_bytes(fasta_gz_data, fai_contents, gzi_data)
            .map_err(|source| ReaderError::FastaOpen { source })?;
        Ok(FuzzReaders {
            alignment: AlignmentBackend::Indexed(CursorReader::Cram(Box::new(cram))),
            fasta,
            store: RecordStore::default(),
            fasta_buf: Vec::new(),
        })
    }

    /// Build BGZF-compressed SAM reader from raw bytes: BGZF SAM + TBI index + FASTA.gz + FAI + GZI.
    pub fn from_sam_bytes(
        sam_data: Vec<u8>,
        tbi_data: &[u8],
        fasta_gz_data: Vec<u8>,
        fai_contents: &str,
        gzi_data: &[u8],
    ) -> Result<Self, ReaderError> {
        let sam = IndexedSamReader::from_bytes(sam_data, tbi_data)?;
        let fasta = IndexedFastaReader::from_bytes(fasta_gz_data, fai_contents, gzi_data)
            .map_err(|source| ReaderError::FastaOpen { source })?;
        Ok(FuzzReaders {
            alignment: AlignmentBackend::Indexed(CursorReader::Sam(sam)),
            fasta,
            store: RecordStore::default(),
            fasta_buf: Vec::new(),
        })
    }

    /// Build plain (uncompressed) SAM reader — no BGZF, no index, linear scan only.
    /// Maximizes SAM parsing coverage by removing the BGZF framing requirement.
    pub fn from_plain_sam_bytes(
        sam_data: Vec<u8>,
        fasta_gz_data: Vec<u8>,
        fai_contents: &str,
        gzi_data: &[u8],
    ) -> Result<Self, ReaderError> {
        let sam = IndexedSamReader::from_plain_bytes(sam_data)?;
        let fasta = IndexedFastaReader::from_bytes(fasta_gz_data, fai_contents, gzi_data)
            .map_err(|source| ReaderError::FastaOpen { source })?;
        Ok(FuzzReaders {
            alignment: AlignmentBackend::PlainSam(sam),
            fasta,
            store: RecordStore::default(),
            fasta_buf: Vec::new(),
        })
    }

    pub fn header(&self) -> &BamHeader {
        match &self.alignment {
            AlignmentBackend::Indexed(r) => r.header(),
            AlignmentBackend::PlainSam(r) => r.header(),
        }
    }

    /// Full pileup pipeline: fetch_into → PileupEngine with reference sequence.
    pub fn pileup(
        &mut self,
        tid: u32,
        start: Pos<Zero>,
        end: Pos<Zero>,
    ) -> Result<PileupEngine, ReaderError> {
        match &mut self.alignment {
            AlignmentBackend::Indexed(r) => r.fetch_into(tid, start, end, &mut self.store)?,
            AlignmentBackend::PlainSam(r) => {
                r.fetch_plain_into(tid, start, end, &mut self.store).map_err(ReaderError::from)?
            }
        };
        let store = std::mem::take(&mut self.store);
        let mut engine = PileupEngine::new(store, start, end);

        // Try to set reference sequence
        if let Some(name) = self.header().target_name(tid) {
            let name = name.to_owned();
            if self.fasta.fetch_seq_into(&name, start, end, &mut self.fasta_buf).is_ok() {
                let buf = std::mem::take(&mut self.fasta_buf);
                let bases = Base::from_ascii_vec(buf);
                let ref_seq = crate::bam::pileup::RefSeq::new(Rc::from(bases), start);
                engine.set_reference_seq(ref_seq);
            }
        }

        Ok(engine)
    }

    fn fetch_records(
        &mut self,
        tid: u32,
        start: Pos<Zero>,
        end: Pos<Zero>,
        store: &mut RecordStore,
    ) -> Result<usize, ReaderError> {
        match &mut self.alignment {
            AlignmentBackend::Indexed(r) => r.fetch_into(tid, start, end, store),
            AlignmentBackend::PlainSam(r) => {
                r.fetch_plain_into(tid, start, end, store).map_err(ReaderError::from)
            }
        }
    }

    pub fn fetch_into(
        &mut self,
        tid: u32,
        start: Pos<Zero>,
        end: Pos<Zero>,
        store: &mut RecordStore,
    ) -> Result<usize, ReaderError> {
        self.fetch_records(tid, start, end, store)
    }

    pub fn recover_store(&mut self, engine: &mut PileupEngine) {
        if let Some(store) = engine.take_store() {
            self.store = store;
        }
    }
}
