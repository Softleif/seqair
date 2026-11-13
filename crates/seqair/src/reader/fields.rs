//! Which per-record fields a reader decodes ([`RecordFields`]), and the
//! record type a reader built with them yields ([`PartialRecord`]), which
//! refuses every field outside them. See `r[unified.record_fields]`.

use crate::bam::{
    CigarOp, MateFields, RawBamRecord,
    aux::{Aux, AuxLayoutError, AuxTag, AuxValue, FromAuxValue, GetAuxError},
};
use crate::reader::Tid;
use seqair_types::{BamFlags, Base, BaseQuality, IupacBase, MappingQuality, Pos0, QPos, SmallVec};

// r[impl unified.record_fields]
/// The per-record fields a consumer reads. A reader built with it yields
/// [`PartialRecord`]s, which refuse every field outside it, the same way on
/// every format; the reader may skip decoding them (CRAM, SAM tags).
///
/// Start from [`Self::core`] and add what you read:
///
/// ```
/// use seqair::reader::{AuxTags, RecordFields};
///
/// let fields = RecordFields::core().with_quality().with_aux(AuxTags::only([*b"MC"]));
/// assert!(!fields.sequence());
/// ```
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct RecordFields {
    names: bool,
    seq: bool,
    qual: bool,
    aux: AuxTags,
}

impl Default for RecordFields {
    fn default() -> Self {
        Self::all()
    }
}

impl RecordFields {
    /// Every field: what a reader decodes unless told otherwise.
    pub fn all() -> Self {
        Self { names: true, seq: true, qual: true, aux: AuxTags::All }
    }

    /// Only the fields always delivered: flags, reference, position, MAPQ,
    /// CIGAR, mate reference and position, TLEN, and the read length. No
    /// name, sequence, qualities or tags.
    pub fn core() -> Self {
        Self { names: false, seq: false, qual: false, aux: AuxTags::None }
    }

    #[must_use]
    pub fn with_names(mut self) -> Self {
        self.names = true;
        self
    }

    #[must_use]
    pub fn with_sequence(mut self) -> Self {
        self.seq = true;
        self
    }

    #[must_use]
    pub fn with_quality(mut self) -> Self {
        self.qual = true;
        self
    }

    /// Replace the tags kept with `aux`.
    #[must_use]
    pub fn with_aux(mut self, aux: AuxTags) -> Self {
        self.aux = aux;
        self
    }

    /// The fields of both, for a consumer serving several tools.
    #[must_use]
    pub fn union(self, other: Self) -> Self {
        Self {
            names: self.names || other.names,
            seq: self.seq || other.seq,
            qual: self.qual || other.qual,
            aux: self.aux.union(other.aux),
        }
    }

    pub fn names(&self) -> bool {
        self.names
    }

    pub fn sequence(&self) -> bool {
        self.seq
    }

    pub fn quality(&self) -> bool {
        self.qual
    }

    pub fn aux(&self) -> &AuxTags {
        &self.aux
    }

    /// Whether this is [`Self::all`].
    pub fn is_all(&self) -> bool {
        *self == Self::all()
    }
}

/// The aux tags a [`RecordFields`] keeps.
#[derive(Debug, Clone, Default, PartialEq, Eq)]
pub enum AuxTags {
    #[default]
    All,
    None,
    /// These tags; build with [`Self::only`].
    Only(SmallVec<AuxTag, 4>),
}

impl AuxTags {
    /// The tags listed, without duplicates, in the order given.
    pub fn only(tags: impl IntoIterator<Item = impl Into<AuxTag>>) -> Self {
        let mut kept = SmallVec::new();
        for tag in tags {
            let tag = tag.into();
            if !kept.contains(&tag) {
                kept.push(tag);
            }
        }
        Self::Only(kept)
    }

    /// Whether `tag` is kept.
    pub fn contains(&self, tag: [u8; 2]) -> bool {
        match self {
            Self::All => true,
            Self::None => false,
            Self::Only(tags) => tags.contains(&AuxTag(tag)),
        }
    }

    #[must_use]
    pub fn union(self, other: Self) -> Self {
        match (self, other) {
            (Self::All, _) | (_, Self::All) => Self::All,
            (Self::None, x) | (x, Self::None) => x,
            (Self::Only(a), Self::Only(b)) => Self::only(a.into_iter().chain(b)),
        }
    }
}

// r[impl unified.partial_record.skipped]
/// A field a [`PartialRecord`] was not decoded with.
#[derive(Debug, Clone, Copy, PartialEq, Eq, thiserror::Error)]
#[error("{field:?} was not requested")]
pub struct Skipped {
    pub field: Field,
}

/// A field [`RecordFields`] can leave out.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
#[non_exhaustive]
pub enum Field {
    Name,
    Seq,
    Qual,
    Aux(AuxTag),
}

// r[impl unified.partial_record]
// r[impl unified.partial_record.no_payload]
/// A record read with [`RecordFields`]. The core fields are always there;
/// every other field is a `Result` that is `Err(Skipped)` exactly when it
/// was not requested, whatever the format and whatever it decoded.
///
/// It has no payload accessor and no writer takes it: a skipped field holds
/// filler. [`Self::into_complete`] gives the [`RawBamRecord`] back for a
/// record read with every field.
///
/// ```compile_fail
/// # use seqair::{bam::BamWriter, reader::PartialRecord};
/// # fn f<W: std::io::Write>(writer: &mut BamWriter<W>, record: &PartialRecord) {
/// writer.write_raw(record);
/// # }
/// ```
///
/// ```compile_fail
/// # fn f(record: &seqair::reader::PartialRecord) {
/// let _ = record.as_payload();
/// # }
/// ```
#[derive(Debug, Clone, Default, PartialEq, Eq)]
pub struct PartialRecord {
    raw: RawBamRecord,
    fields: RecordFields,
}

impl PartialRecord {
    /// Wrap a record decoded with `fields`, for tests and other crates'
    /// readers; `raw` must hold every field `fields` names.
    pub fn new(raw: RawBamRecord, fields: RecordFields) -> Self {
        Self { raw, fields }
    }

    /// The fields this record was read with.
    pub fn fields(&self) -> &RecordFields {
        &self.fields
    }

    /// The full record, if it was read with every field.
    ///
    /// # Errors
    /// The record itself, unchanged, when a field was skipped.
    pub fn into_complete(self) -> Result<RawBamRecord, Self> {
        if self.fields.is_all() { Ok(self.raw) } else { Err(self) }
    }

    fn gate(&self, requested: bool, field: Field) -> Result<&RawBamRecord, Skipped> {
        if requested { Ok(&self.raw) } else { Err(Skipped { field }) }
    }

    pub fn tid(&self) -> Option<Tid> {
        self.raw.tid()
    }

    pub fn pos(&self) -> Option<Pos0> {
        self.raw.pos()
    }

    pub fn mate_fields(&self) -> MateFields {
        self.raw.mate_fields()
    }

    pub fn flags(&self) -> BamFlags {
        self.raw.flags()
    }

    pub fn mapq(&self) -> MappingQuality {
        self.raw.mapq()
    }

    pub fn cigar_len(&self) -> usize {
        self.raw.cigar_len()
    }

    /// As [`RawBamRecord::cigar`]: a long CIGAR is resolved whatever tags
    /// were requested.
    pub fn cigar(&self) -> impl ExactSizeIterator<Item = CigarOp> + '_ {
        self.raw.cigar()
    }

    pub fn end_pos(&self) -> Option<Pos0> {
        self.raw.end_pos()
    }

    // r[impl unified.partial_record.l_seq]
    /// The read length, whether or not the sequence was requested.
    pub fn seq_len(&self) -> usize {
        self.raw.seq_len()
    }

    /// The read name, as [`RawBamRecord::qname`].
    ///
    /// # Errors
    /// Names were not requested.
    pub fn qname(&self) -> Result<&[u8], Skipped> {
        Ok(self.gate(self.fields.names, Field::Name)?.qname())
    }

    /// The sequence.
    ///
    /// # Errors
    /// The sequence was not requested.
    pub fn seq(&self) -> Result<SeqRef<'_>, Skipped> {
        Ok(SeqRef { raw: self.gate(self.fields.seq, Field::Seq)? })
    }

    /// The base qualities, `Ok(None)` when the record has none, as
    /// [`RawBamRecord::qual`].
    ///
    /// # Errors
    /// Qualities were not requested.
    pub fn qual(&self) -> Result<Option<&[BaseQuality]>, Skipped> {
        Ok(self.gate(self.fields.qual, Field::Qual)?.qual())
    }

    /// The requested tags.
    pub fn aux(&self) -> RequestedAux<'_> {
        RequestedAux { aux: self.raw.aux(), tags: &self.fields.aux }
    }

    /// The record for a reader to fill, now marked as read with `fields`;
    /// reassigned only when they differ, so refilling stays cheap.
    pub(crate) fn raw_mut(&mut self, fields: &RecordFields) -> &mut RawBamRecord {
        if self.fields != *fields {
            self.fields.clone_from(fields);
        }
        &mut self.raw
    }

    #[cfg(test)]
    pub(crate) fn from_payload(
        payload: Vec<u8>,
        header: &crate::bam::BamHeader,
        fields: RecordFields,
    ) -> Result<Self, crate::bam::DecodeError> {
        Ok(Self { raw: RawBamRecord::from_payload(payload, header)?, fields })
    }
}

/// A [`PartialRecord`]'s sequence, which was requested.
#[derive(Debug, Clone, Copy)]
pub struct SeqRef<'a> {
    raw: &'a RawBamRecord,
}

impl<'a> SeqRef<'a> {
    pub fn len(&self) -> usize {
        self.raw.seq_len()
    }

    pub fn is_empty(&self) -> bool {
        self.raw.seq_len() == 0
    }

    /// As [`RawBamRecord::packed_seq`].
    pub fn packed(&self) -> &'a [u8] {
        self.raw.packed_seq()
    }

    /// As [`RawBamRecord::iupac_bases`].
    pub fn iupac_bases(&self) -> impl ExactSizeIterator<Item = IupacBase> + 'a {
        self.raw.iupac_bases()
    }

    /// As [`RawBamRecord::iupac_bases_into`].
    pub fn iupac_bases_into(&self, out: &mut Vec<IupacBase>) {
        self.raw.iupac_bases_into(out);
    }

    /// As [`RawBamRecord::base_at`].
    pub fn base_at(&self, qpos: QPos) -> Option<IupacBase> {
        self.raw.base_at(qpos)
    }

    /// As [`RawBamRecord::decode_seq_into`].
    pub fn decode_into(&self, out: &mut Vec<u8>) {
        self.raw.decode_seq_into(out);
    }

    /// As [`RawBamRecord::bases`].
    pub fn bases(&self) -> Vec<Base> {
        self.raw.bases()
    }
}

/// A [`PartialRecord`]'s tags, of which only the requested ones are seen.
#[derive(Debug, Clone, Copy)]
pub struct RequestedAux<'a> {
    aux: Aux<'a>,
    tags: &'a AuxTags,
}

impl<'a> RequestedAux<'a> {
    /// As [`Aux::get`].
    ///
    /// # Errors
    /// [`GetAuxError::NotRequested`] for a tag that was not requested,
    /// whether or not the record has it; otherwise as [`Aux::get`].
    #[must_use = "get() returns a fallible Result — ignoring it silently drops errors"]
    pub fn get<T: FromAuxValue<'a>>(&self, tag: impl AsRef<[u8]>) -> Result<T, GetAuxError> {
        let bytes = tag.as_ref();
        if let Ok(tag) = <[u8; 2]>::try_from(bytes)
            && !self.tags.contains(tag)
        {
            return Err(GetAuxError::NotRequested { tag: AuxTag(tag) });
        }
        self.aux.get(bytes)
    }

    /// The requested tags the record has, as [`Aux::iter_tags`].
    pub fn iter_tags(&self) -> impl Iterator<Item = ([u8; 2], AuxValue<'a>)> + 'a {
        let tags = self.tags;
        self.aux.iter_tags().filter(move |(tag, _)| tags.contains(*tag))
    }

    /// The requested tags the record has, as [`Aux::iter_checked`].
    pub fn iter_checked(
        &self,
    ) -> impl Iterator<Item = Result<([u8; 2], AuxValue<'a>), AuxLayoutError>> + 'a {
        let tags = self.tags;
        self.aux
            .iter_checked()
            .filter(move |item| item.as_ref().map_or(true, |(t, _)| tags.contains(*t)))
    }
}

mod sealed {
    pub trait Sealed {}
    impl Sealed for super::Complete {}
    impl Sealed for super::Partial {}
}

/// What a reader yields: [`Complete`] records ([`RawBamRecord`]), or
/// [`Partial`] ones ([`PartialRecord`]) once given [`RecordFields`].
pub trait RecordKind: sealed::Sealed + Send + Sync + 'static {
    type Record: Default + std::fmt::Debug + Send;

    #[doc(hidden)]
    fn wrap(raw: RawBamRecord, fields: &RecordFields) -> Self::Record;

    #[doc(hidden)]
    fn raw_mut<'r>(record: &'r mut Self::Record, fields: &RecordFields) -> &'r mut RawBamRecord;
}

/// Readers yield [`RawBamRecord`]s with every field.
#[derive(Debug)]
pub enum Complete {}

/// Readers yield [`PartialRecord`]s.
#[derive(Debug)]
pub enum Partial {}

impl RecordKind for Complete {
    type Record = RawBamRecord;

    fn wrap(raw: RawBamRecord, _: &RecordFields) -> RawBamRecord {
        raw
    }

    fn raw_mut<'r>(record: &'r mut RawBamRecord, _: &RecordFields) -> &'r mut RawBamRecord {
        record
    }
}

impl RecordKind for Partial {
    type Record = PartialRecord;

    fn wrap(raw: RawBamRecord, fields: &RecordFields) -> PartialRecord {
        PartialRecord { raw, fields: fields.clone() }
    }

    fn raw_mut<'r>(record: &'r mut PartialRecord, fields: &RecordFields) -> &'r mut RawBamRecord {
        record.raw_mut(fields)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn union_keeps_every_field_of_either() {
        let a = RecordFields::core().with_quality().with_aux(AuxTags::only([*b"NM"]));
        let b = RecordFields::core().with_names().with_aux(AuxTags::only([*b"MC", *b"NM"]));
        let u = a.union(b);
        assert!(u.names() && u.quality() && !u.sequence());
        assert_eq!(u.aux(), &AuxTags::only([*b"NM", *b"MC"]));
        assert!(RecordFields::core().union(RecordFields::all()).is_all());
    }

    fn complete() -> RawBamRecord {
        use crate::bam::{AuxData, OwnedBamRecord, test_util};
        let mut aux = AuxData::new();
        aux.set_int(*b"NM", 2).unwrap();
        aux.set_string(*b"MC", b"4M").unwrap();
        let owned = OwnedBamRecord::builder(Some(test_util::tid(0)), Pos0::new(5), b"r1".to_vec())
            .cigar(vec![CigarOp::new(crate::bam::cigar::CigarOpType::Match, 4)])
            .seq(vec![IupacBase::A; 4])
            .qual(vec![BaseQuality::from_byte(30); 4])
            .aux(aux)
            .build()
            .unwrap();
        let mut payload = Vec::new();
        owned.to_bam_bytes(&mut payload).unwrap();
        RawBamRecord::from_payload(payload, test_util::header()).unwrap()
    }

    #[test]
    fn a_partial_record_refuses_what_was_not_requested() {
        let header = crate::bam::test_util::header();
        let complete = complete();
        let payload = complete.as_payload().to_vec();
        let fields = RecordFields::core().with_quality().with_aux(AuxTags::only([*b"MC"]));
        let partial = PartialRecord::from_payload(payload, header, fields).unwrap();

        assert_eq!(partial.qname(), Err(Skipped { field: Field::Name }));
        assert!(partial.seq().is_err());
        assert!(partial.qual().is_ok());
        assert_eq!(partial.seq_len(), complete.seq_len());
        assert!(matches!(
            partial.aux().get::<i64>("NM"),
            Err(GetAuxError::NotRequested { tag: AuxTag(t) }) if t == *b"NM"
        ));
        assert_eq!(partial.aux().get::<&str>("MC").unwrap(), "4M");
        assert!(matches!(partial.aux().get::<i64>("XX"), Err(GetAuxError::NotRequested { .. })));
        let seen: Vec<_> = partial.aux().iter_tags().map(|(t, _)| t).collect();
        assert_eq!(seen, [*b"MC"]);
        assert!(partial.into_complete().is_err());
        let full = PartialRecord::new(complete.clone(), RecordFields::all());
        assert_eq!(full.into_complete(), Ok(complete));
    }
}
