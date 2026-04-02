use crate::{
    ReaderError,
    bam::{BamHeader, IndexedBamReader, PileupEngine, RecordStore},
    cram::reader::IndexedCramReader,
    fasta::IndexedFastaReader,
    reader::indexed::CursorReader,
};
use seqair_types::{Base, Pos, Zero};
use std::rc::Rc;

/// In-memory reader bundle for fuzzing. No file I/O — all data from byte slices.
/// Uses the same `IndexedReader` enum as production code, just with `Cursor` I/O.
pub struct FuzzReaders {
    alignment: CursorReader,
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
            alignment: CursorReader::Bam(bam),
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
            alignment: CursorReader::Cram(Box::new(cram)),
            fasta,
            store: RecordStore::default(),
            fasta_buf: Vec::new(),
        })
    }

    pub fn header(&self) -> &BamHeader {
        self.alignment.header()
    }

    /// Full pileup pipeline: fetch_into → PileupEngine with reference sequence.
    pub fn pileup(
        &mut self,
        tid: u32,
        start: Pos<Zero>,
        end: Pos<Zero>,
    ) -> Result<PileupEngine, ReaderError> {
        self.alignment.fetch_into(tid, start, end, &mut self.store)?;
        let store = std::mem::take(&mut self.store);
        let mut engine = PileupEngine::new(store, start, end);

        // Try to set reference sequence
        if let Some(name) = self.alignment.header().target_name(tid) {
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

    pub fn fetch_into(
        &mut self,
        tid: u32,
        start: Pos<Zero>,
        end: Pos<Zero>,
        store: &mut RecordStore,
    ) -> Result<usize, ReaderError> {
        self.alignment.fetch_into(tid, start, end, store)
    }

    pub fn recover_store(&mut self, engine: &mut PileupEngine) {
        if let Some(store) = engine.take_store() {
            self.store = store;
        }
    }
}
