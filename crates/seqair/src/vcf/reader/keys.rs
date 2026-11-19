//! A file's field table and resolving keys through it
//! (`r[vcf_reader.field_table]`, `r[vcf_values.file_keys+3]`,
//! `r[vcf_values.optional_keys+2]`, `r[vcf_values.lenient_keys]`).

use super::fields::optional;
use super::file::{FileShared, VcfFile};
use super::index::DictId;
use super::reserved;
use super::text_keys::TextKeys;
use crate::vcf::error::{KeyError, KeyKind};
use crate::vcf::header::{FormatDef, HeaderId, InfoDef, Number, ValueType, VcfHeader};
use crate::vcf::record_encoder::{
    FieldId, FieldKind, FilterId, FormatKey, InfoKey, Section, check_genotype, check_kind,
};
use seqair_types::SmolStr;
use std::collections::HashMap;
use std::marker::PhantomData;
use std::sync::atomic::{AtomicBool, AtomicUsize, Ordering};
use std::sync::{OnceLock, PoisonError, RwLock};
use tracing::{debug, warn};

/// Where an id's definition in one section comes from.
#[derive(Debug)]
enum Decl<D: 'static> {
    /// The header declares it: its place in the header's INFO or FORMAT
    /// map, so a lookup by id reads it without hashing the name.
    Header(usize),
    /// VCF 4.5 reserves it and the header leaves it out (VCF text only).
    Reserved(&'static str, &'static D),
}

// Not derived: those would ask `D: Copy`.
impl<D> Clone for Decl<D> {
    fn clone(&self) -> Self {
        *self
    }
}
impl<D> Copy for Decl<D> {}

/// A definition in one section, and the `Type` it gives.
#[derive(Debug)]
struct Declared<D: 'static> {
    typ: ValueType,
    decl: Decl<D>,
}

impl<D> Clone for Declared<D> {
    fn clone(&self) -> Self {
        *self
    }
}
impl<D> Copy for Declared<D> {}

// r[impl vcf_reader.field_table]
// r[impl vcf_values.text_resolution+2]
/// The INFO/FORMAT keys a file's records can name: the header's dictionary
/// and, for VCF text, the reserved keys it leaves out (ids past the
/// header's), then keys records use or callers resolve that neither
/// declares, read as `Type=String, Number=.` (ids after those).
#[derive(Debug)]
pub(crate) struct FieldTable {
    /// Header and reserved names to ids: what every record looks up first.
    keys: TextKeys,
    /// By id: the INFO / FORMAT definition, `None` where there is none.
    info: Vec<Option<Declared<InfoDef>>>,
    format: Vec<Option<Declared<FormatDef>>>,
    /// The header dictionary's size: ids from here are reserved, then extra.
    header_len: usize,
    /// Names of the reserved ids, from `header_len`.
    reserved_names: Vec<&'static str>,
    /// The id the table's keys carry when the header does not declare them:
    /// the file's readers accept it, every writer refuses it.
    pub(crate) id: HeaderId,
    /// Keys neither declares, by name; their names by id after the reserved.
    extra: RwLock<HashMap<Box<[u8]>, DictId>>,
    extra_names: Arena,
    warned: AtomicBool,
    /// The file's path, for the warning.
    pub(crate) path: Option<std::path::PathBuf>,
}

impl FieldTable {
    pub(crate) fn new(header: &VcfHeader, text: bool) -> Self {
        let dict = header.string_map();
        let header_len = dict.len();
        let mut info = vec![None; header_len];
        let mut format = vec![None; header_len];
        for (at, (name, def)) in header.infos().iter().enumerate() {
            if let Some(i) = dict.get(name) {
                set(&mut info, i, Declared { typ: def.typ, decl: Decl::Header(at) });
            }
        }
        for (at, (name, def)) in header.formats().iter().enumerate() {
            if let Some(i) = dict.get(name) {
                set(&mut format, i, Declared { typ: def.typ, decl: Decl::Header(at) });
            }
        }
        let mut reserved_names: Vec<&'static str> = Vec::new();
        if text {
            // The id of reserved `name`: the header's, or a new one past it.
            let mut id_of = |name: &'static str| -> Option<usize> {
                if let Some(i) = dict.get(name) {
                    return Some(i);
                }
                if let Some(at) = reserved_names.iter().position(|&n| n == name) {
                    return Some(header_len.saturating_add(at));
                }
                let id = header_len.saturating_add(reserved_names.len());
                // r[impl vcf_reader.limits+2]
                DictId::new(u32::try_from(id).ok()?)?;
                reserved_names.push(name);
                Some(id)
            };
            for (name, def) in reserved::INFO {
                if !header.infos().contains_key(*name)
                    && let Some(id) = id_of(name)
                {
                    set(&mut info, id, Declared { typ: def.typ, decl: Decl::Reserved(name, def) });
                }
            }
            for (name, def) in reserved::FORMAT {
                if !header.formats().contains_key(*name)
                    && let Some(id) = id_of(name)
                {
                    let declared = Declared { typ: def.typ, decl: Decl::Reserved(name, def) };
                    set(&mut format, id, declared);
                }
            }
        }
        let base = header_len.saturating_add(reserved_names.len());
        info.resize(base, None);
        format.resize(base, None);
        let reserved_ids =
            reserved_names.iter().enumerate().map(|(i, &n)| (header_len.saturating_add(i), n));
        Self {
            keys: TextKeys::new(dict.iter().chain(reserved_ids)),
            info,
            format,
            header_len,
            reserved_names,
            id: HeaderId::next(),
            extra: RwLock::default(),
            extra_names: Arena::default(),
            warned: AtomicBool::new(false),
            path: None,
        }
    }

    /// Ids below this one are the header's or reserved.
    pub(crate) fn base(&self) -> usize {
        self.header_len.saturating_add(self.reserved_names.len())
    }

    /// The `Type` id `id` has in `kind`'s section, if it has a definition
    /// there: the header's or a reserved one.
    pub(crate) fn declared(&self, kind: KeyKind, id: DictId) -> Option<ValueType> {
        let i = id.get() as usize;
        match kind {
            KeyKind::Info => self.info.get(i).copied().flatten().map(|d| d.typ),
            KeyKind::Format => self.format.get(i).copied().flatten().map(|d| d.typ),
            KeyKind::Filter => None,
        }
    }

    /// The `Type` field `name` of `kind`'s section has, if it has one.
    pub(crate) fn declared_name(&self, kind: KeyKind, name: &str) -> Option<ValueType> {
        self.declared(kind, self.keys.get(name.as_bytes())?)
    }

    /// The header's or a reserved id of `name`, whatever its section.
    pub(crate) fn get(&self, name: &[u8]) -> Option<DictId> {
        self.keys.get(name)
    }

    /// The name of id `id`.
    pub(crate) fn name<'a>(&'a self, header: &'a VcfHeader, id: DictId) -> Option<&'a str> {
        let i = id.get() as usize;
        if i < self.header_len {
            return header.string_map().name(i);
        }
        if let Some(&name) = self.reserved_names.get(i.saturating_sub(self.header_len)) {
            return Some(name);
        }
        self.extra_names.get(i.checked_sub(self.base())?).map(SmolStr::as_str)
    }

    /// Id `id`'s name and INFO definition: the header's, a reserved one, or
    /// the undeclared one.
    pub(crate) fn info_def<'a>(
        &'a self,
        header: &'a VcfHeader,
        id: DictId,
    ) -> Option<(&'a str, &'a InfoDef)> {
        Some(match self.info.get(id.get() as usize).copied().flatten().map(|d| d.decl) {
            Some(Decl::Header(at)) => header.infos().get_index(at).map(|(n, d)| (n.as_str(), d))?,
            Some(Decl::Reserved(name, def)) => (name, def),
            None => (self.name(header, id)?, &reserved::UNDECLARED_INFO),
        })
    }

    /// Id `id`'s name and FORMAT definition (see [`Self::info_def`]).
    pub(crate) fn format_def<'a>(
        &'a self,
        header: &'a VcfHeader,
        id: DictId,
    ) -> Option<(&'a str, &'a FormatDef)> {
        Some(match self.format.get(id.get() as usize).copied().flatten().map(|d| d.decl) {
            Some(Decl::Header(at)) => {
                header.formats().get_index(at).map(|(n, d)| (n.as_str(), d))?
            }
            Some(Decl::Reserved(name, def)) => (name, def),
            None => (self.name(header, id)?, &reserved::UNDECLARED_FORMAT),
        })
    }

    // r[impl vcf_reader.undeclared+5]
    /// The id of `key`, a VCF text record's `kind` key, appending it to the
    /// table if neither the header nor VCF defines it. A key the header does
    /// not declare, reserved or not, is warned of once per file. `None` only
    /// for a key no id can be given: not UTF-8, or past 2^16 ids.
    #[inline]
    pub(crate) fn text_id(&self, kind: KeyKind, key: &[u8]) -> Option<DictId> {
        match self.keys.get(key) {
            Some(id) if self.declared(kind, id).is_some() => {
                if id.get() as usize >= self.header_len {
                    self.warn_record(kind, key);
                }
                Some(id)
            }
            _ => self.undeclared(kind, key, true),
        }
    }

    /// The id of an undeclared key, appended on first use; `record`: a
    /// record uses it, which is warned of once per file.
    #[cold]
    fn undeclared(&self, kind: KeyKind, key: &[u8], record: bool) -> Option<DictId> {
        if record {
            self.warn_record(kind, key);
        }
        if let Some(&id) = self.extra.read().unwrap_or_else(PoisonError::into_inner).get(key) {
            return Some(id);
        }
        let name = std::str::from_utf8(key).ok()?;
        let mut extras = self.extra.write().unwrap_or_else(PoisonError::into_inner);
        if let Some(&id) = extras.get(key) {
            return Some(id);
        }
        let id = match self.keys.get(key) {
            Some(id) => id,
            None => {
                let index = self.extra_names.len();
                // r[impl vcf_reader.limits+2]
                let id =
                    u32::try_from(self.base().saturating_add(index)).ok().and_then(DictId::new)?;
                self.extra_names.push(index, SmolStr::from(name));
                id
            }
        };
        extras.insert(key.into(), id);
        Some(id)
    }

    /// Log, once per file, that a record uses a key the header does not declare.
    #[inline]
    pub(crate) fn warn_record(&self, kind: KeyKind, key: &[u8]) {
        if !self.warned.load(Ordering::Relaxed) {
            self.warn_record_cold(kind, key);
        }
    }

    #[cold]
    fn warn_record_cold(&self, kind: KeyKind, key: &[u8]) {
        if !self.warned.swap(true, Ordering::Relaxed) {
            warn!(
                %kind,
                key = std::str::from_utf8(key).ok(),
                path = ?self.path,
                "a record uses a key the header does not declare; plain keys do not read it, \
                 lenient keys do (further undeclared keys in this file are not reported)"
            );
        }
    }

    /// The field `name` of `kind`'s section, as marker `V` reads it: the
    /// header's, or for a `lenient` key also one it does not declare.
    fn key<V: FieldKind>(
        &self,
        header: &VcfHeader,
        kind: KeyKind,
        name: &str,
        lenient: bool,
    ) -> Result<FieldId, KeyError> {
        let section = if kind == KeyKind::Format { Section::Format } else { Section::Info };
        let declared = match section {
            Section::Info => header.info_key::<V>(name).map(|k| k.0),
            Section::Format => header.format_key::<V>(name).map(|k| k.0),
        };
        match declared {
            Err(KeyError::UnknownField { .. }) if lenient => {}
            declared => return declared,
        }
        let reserved = match section {
            Section::Info => {
                reserved::INFO.iter().find(|(n, _)| *n == name).map(|(_, d)| (d.typ, d.number))
            }
            Section::Format => {
                reserved::FORMAT.iter().find(|(n, _)| *n == name).map(|(_, d)| (d.typ, d.number))
            }
        };
        let unknown = || KeyError::UnknownField { name: name.into(), kind };
        let (typ, number) = match reserved {
            Some((typ, number)) => {
                check_kind::<V>(section, name, typ, number)?;
                (typ, number)
            }
            None => {
                check_genotype::<V>(section, name)?;
                if !V::reads_type(ValueType::String) || !V::holds_number(Number::Unknown) {
                    return Err(unknown());
                }
                (ValueType::String, Number::Unknown)
            }
        };
        // A VCF text table has the reserved ids; a BCF one gives the name an
        // id past the dictionary, which no record can name.
        let id = self.undeclared(kind, name.as_bytes(), false).ok_or_else(unknown)?;
        // r[impl vcf_reader.undeclared+5]
        // Warned of where a record uses it: most lenient keys find nothing.
        debug!(%kind, key = name, %typ, %number, "a lenient key for a field the header does not declare");
        Ok(FieldId::new(id, name.into(), typ, number, section, self.id))
    }
}

fn set<T: Copy>(table: &mut Vec<Option<T>>, id: usize, value: T) {
    if table.len() <= id {
        table.resize(id.saturating_add(1), None);
    }
    if let Some(slot) = table.get_mut(id) {
        *slot = Some(value);
    }
}

/// Names that live as long as the table, appended under its write lock and
/// read without one: chunks of doubling size, each allocated once.
#[derive(Debug)]
struct Arena {
    chunks: [OnceLock<Box<[OnceLock<SmolStr>]>>; ARENA_CHUNKS],
    len: AtomicUsize,
}

/// The first chunk's size; chunk `k` holds `ARENA_FIRST << k`.
const ARENA_FIRST: usize = 16;
/// Enough chunks for every id a dictionary can hold (2^16).
const ARENA_CHUNKS: usize = 13;

impl Default for Arena {
    fn default() -> Self {
        Self { chunks: std::array::from_fn(|_| OnceLock::new()), len: AtomicUsize::new(0) }
    }
}

impl Arena {
    /// The chunk and offset of entry `i`.
    fn locate(i: usize) -> (usize, usize) {
        let k = (i / ARENA_FIRST).saturating_add(1).ilog2() as usize;
        let start = ARENA_FIRST.saturating_mul((1usize << k).saturating_sub(1));
        (k, i.saturating_sub(start))
    }

    fn len(&self) -> usize {
        self.len.load(Ordering::Acquire)
    }

    fn get(&self, i: usize) -> Option<&SmolStr> {
        let (k, at) = Self::locate(i);
        self.chunks.get(k)?.get()?.get(at)?.get()
    }

    /// Store entry `i`, the next one, under the table's write lock.
    fn push(&self, i: usize, name: SmolStr) {
        let (k, at) = Self::locate(i);
        let Some(chunk) = self.chunks.get(k) else { return };
        let chunk = chunk.get_or_init(|| (0..ARENA_FIRST << k).map(|_| OnceLock::new()).collect());
        if let Some(slot) = chunk.get(at)
            && slot.set(name).is_ok()
        {
            self.len.store(i.saturating_add(1), Ordering::Release);
        }
    }
}

impl FileShared {
    pub(crate) fn info_key<V: FieldKind>(
        &self,
        name: &str,
        lenient: bool,
    ) -> Result<InfoKey<V>, KeyError> {
        Ok(InfoKey(self.fields.key::<V>(&self.header, KeyKind::Info, name, lenient)?, PhantomData))
    }

    pub(crate) fn format_key<V: FieldKind>(
        &self,
        name: &str,
        lenient: bool,
    ) -> Result<FormatKey<V>, KeyError> {
        let id = self.fields.key::<V>(&self.header, KeyKind::Format, name, lenient)?;
        Ok(FormatKey(id, PhantomData))
    }
}

impl VcfFile {
    // r[impl vcf_values.file_keys+3]
    /// Resolve an INFO key for reading this file, checking marker `V`
    /// against the header's definition: `let dp: InfoInt = file.info_key("DP")?;`.
    ///
    /// Resolves what [`VcfHeader::info_key`](crate::vcf::VcfHeader::info_key)
    /// resolves, the same for VCF text and BCF. A key the header does not
    /// declare needs [`Self::lenient`].
    ///
    /// # Errors
    /// [`KeyError::UnknownField`] for a field the header does not declare,
    /// `TypeMismatch`/`NumberMismatch`/`NotGenotype` for a marker that
    /// cannot read it.
    pub fn info_key<V: FieldKind>(&self, name: &str) -> Result<InfoKey<V>, KeyError> {
        self.0.info_key(name, false)
    }

    // r[impl vcf_values.file_keys+3]
    /// Resolve a FORMAT key for reading this file (see [`Self::info_key`]).
    ///
    /// # Errors
    /// As [`Self::info_key`].
    pub fn format_key<V: FieldKind>(&self, name: &str) -> Result<FormatKey<V>, KeyError> {
        self.0.format_key(name, false)
    }

    // r[impl vcf_values.optional_keys+2]
    /// An INFO key, or `None` if the header does not declare the field: for
    /// an optional annotation, one `?` and one `Option`.
    ///
    /// # Errors
    /// As [`Self::info_key`], but for `UnknownField`.
    pub fn info_key_opt<V: FieldKind>(&self, name: &str) -> Result<Option<InfoKey<V>>, KeyError> {
        optional(self.info_key(name))
    }

    // r[impl vcf_values.optional_keys+2]
    /// A FORMAT key, or `None` if the header does not declare the field.
    ///
    /// # Errors
    /// As [`Self::format_key`], but for `UnknownField`.
    pub fn format_key_opt<V: FieldKind>(
        &self,
        name: &str,
    ) -> Result<Option<FormatKey<V>>, KeyError> {
        optional(self.format_key(name))
    }

    // r[impl vcf_values.lenient_keys]
    /// Keys that also resolve fields the header does not declare, for
    /// reading files whose headers forgot them.
    pub fn lenient(&self) -> LenientKeys<'_> {
        LenientKeys(&self.0)
    }

    // r[impl vcf_values.file_keys+3]
    /// A FILTER the header declares, as [`VcfHeader::filter`](crate::vcf::VcfHeader::filter)
    /// and [`PlanContext::filter`](super::PlanContext::filter) resolve it.
    ///
    /// # Errors
    /// [`KeyError::UnknownFilter`].
    pub fn filter(&self, name: &str) -> Result<FilterId, KeyError> {
        self.0.header.filter(name)
    }
}

// r[impl vcf_values.reader_keys+2]
/// A reader's key resolution, from [`Reader::keys`](super::Reader::keys):
/// what [`VcfFile::info_key`] and its siblings resolve, for a reader
/// whose file is not at hand (a [`VcfStream`](super::VcfStream)'s).
#[derive(Debug, Clone, Copy)]
pub struct Keys<'a>(pub(crate) &'a FileShared);

impl<'a> Keys<'a> {
    /// As [`VcfFile::info_key`].
    ///
    /// # Errors
    /// As [`VcfFile::info_key`].
    pub fn info_key<V: FieldKind>(&self, name: &str) -> Result<InfoKey<V>, KeyError> {
        self.0.info_key(name, false)
    }

    /// As [`VcfFile::format_key`].
    ///
    /// # Errors
    /// As [`VcfFile::format_key`].
    pub fn format_key<V: FieldKind>(&self, name: &str) -> Result<FormatKey<V>, KeyError> {
        self.0.format_key(name, false)
    }

    /// As [`VcfFile::info_key_opt`].
    ///
    /// # Errors
    /// As [`Self::info_key`], but for `UnknownField`.
    pub fn info_key_opt<V: FieldKind>(&self, name: &str) -> Result<Option<InfoKey<V>>, KeyError> {
        optional(self.info_key(name))
    }

    /// As [`VcfFile::format_key_opt`].
    ///
    /// # Errors
    /// As [`Self::format_key`], but for `UnknownField`.
    pub fn format_key_opt<V: FieldKind>(
        &self,
        name: &str,
    ) -> Result<Option<FormatKey<V>>, KeyError> {
        optional(self.format_key(name))
    }

    /// As [`VcfFile::lenient`].
    pub fn lenient(&self) -> LenientKeys<'a> {
        LenientKeys(self.0)
    }

    /// As [`VcfFile::filter`].
    ///
    /// # Errors
    /// As [`VcfFile::filter`].
    pub fn filter(&self, name: &str) -> Result<FilterId, KeyError> {
        self.0.header.filter(name)
    }
}

// r[impl vcf_values.lenient_keys]
/// Key resolution that accepts fields the header does not declare, from
/// [`VcfFile::lenient`] or [`PlanContext::lenient`](super::PlanContext::lenient).
///
/// A declared field resolves as [`VcfFile::info_key`] resolves it. An
/// undeclared one VCF 4.5 reserves resolves with its reserved `Type` and
/// `Number`; any other undeclared one resolves as `Type=String, Number=.`
/// for a [`Str`](crate::vcf::Str) or [`Any`](crate::vcf::Any) marker only.
/// A record using an undeclared key is warned of once per file. VCF text and BCF resolve
/// the same names to the same definitions; a BCF record can only name
/// declared fields, so in BCF such a key finds nothing. A writer refuses
/// the key, as its header has no id for the field.
#[derive(Debug, Clone, Copy)]
pub struct LenientKeys<'a>(pub(crate) &'a FileShared);

impl LenientKeys<'_> {
    /// An INFO key, declared or not.
    ///
    /// # Errors
    /// [`KeyError::UnknownField`] for an undeclared field that is not
    /// reserved, unless `V` reads a String; as [`VcfFile::info_key`] else.
    pub fn info_key<V: FieldKind>(&self, name: &str) -> Result<InfoKey<V>, KeyError> {
        self.0.info_key(name, true)
    }

    /// A FORMAT key, declared or not (see [`Self::info_key`]).
    ///
    /// # Errors
    /// As [`Self::info_key`].
    pub fn format_key<V: FieldKind>(&self, name: &str) -> Result<FormatKey<V>, KeyError> {
        self.0.format_key(name, true)
    }

    /// An INFO key, or `None` where [`Self::info_key`] fails with `UnknownField`.
    ///
    /// # Errors
    /// As [`Self::info_key`], but for `UnknownField`.
    pub fn info_key_opt<V: FieldKind>(&self, name: &str) -> Result<Option<InfoKey<V>>, KeyError> {
        optional(self.info_key(name))
    }

    /// A FORMAT key, or `None` where [`Self::format_key`] fails with `UnknownField`.
    ///
    /// # Errors
    /// As [`Self::format_key`], but for `UnknownField`.
    pub fn format_key_opt<V: FieldKind>(
        &self,
        name: &str,
    ) -> Result<Option<FormatKey<V>>, KeyError> {
        optional(self.format_key(name))
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::vcf::record_encoder::Str;
    use std::sync::Arc;

    /// Counts `warn!` events.
    #[derive(Default)]
    struct Count(AtomicUsize);

    impl tracing::Subscriber for Count {
        fn enabled(&self, _: &tracing::Metadata<'_>) -> bool {
            true
        }
        fn new_span(&self, _: &tracing::span::Attributes<'_>) -> tracing::span::Id {
            tracing::span::Id::from_u64(1)
        }
        fn record(&self, _: &tracing::span::Id, _: &tracing::span::Record<'_>) {}
        fn record_follows_from(&self, _: &tracing::span::Id, _: &tracing::span::Id) {}
        fn event(&self, event: &tracing::Event<'_>) {
            if *event.metadata().level() == tracing::Level::WARN {
                self.0.fetch_add(1, Ordering::Relaxed);
            }
        }
        fn enter(&self, _: &tracing::span::Id) {}
        fn exit(&self, _: &tracing::span::Id) {}
    }

    // r[verify vcf_reader.undeclared+5]
    // r[verify vcf_values.lenient_keys]
    /// Resolving keys, strict or lenient, warns of nothing; the first record
    /// using an undeclared key, reserved or not, warns once per file.
    #[test]
    fn undeclared_keys_are_logged_once_where_records_use_them() {
        let header = VcfHeader::parse(
            "##fileformat=VCFv4.3\n##INFO=<ID=DP,Number=1,Type=Integer,Description=\"d\">\n\
             #CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n",
        )
        .unwrap();
        let count = Arc::new(Count::default());
        let seen = || count.0.load(Ordering::Relaxed);
        tracing::subscriber::with_default(Arc::clone(&count), || {
            let fields = FieldTable::new(&header, true);
            assert!(matches!(
                fields.key::<Str>(&header, KeyKind::Info, "SVTYPE", false),
                Err(KeyError::UnknownField { .. })
            ));
            fields.key::<Str>(&header, KeyKind::Info, "SVTYPE", true).unwrap();
            let trid = fields.key::<Str>(&header, KeyKind::Info, "TRID", true).unwrap();
            assert_eq!(seen(), 0, "resolving warns of nothing");
            fields.text_id(KeyKind::Info, b"DP").unwrap();
            assert_eq!(seen(), 0, "a declared key");
            assert_eq!(fields.text_id(KeyKind::Info, b"TRID"), Some(trid.id));
            assert_eq!(seen(), 1, "a record's undeclared key, resolved or not");
            let xx = fields.text_id(KeyKind::Info, b"XX").unwrap();
            fields.text_id(KeyKind::Info, b"SVLEN").unwrap();
            assert_eq!(seen(), 1, "once per file");
            assert_eq!(fields.text_id(KeyKind::Format, b"XX"), Some(xx), "one id per name");
            assert_eq!(fields.name(&header, xx), Some("XX"));

            let fields = FieldTable::new(&header, true);
            fields.text_id(KeyKind::Info, b"SVLEN").unwrap();
            assert_eq!(seen(), 2, "a reserved key the header leaves out");
        });
    }

    #[test]
    fn the_arena_finds_what_it_was_given() {
        let arena = Arena::default();
        for i in 0..1000 {
            arena.push(i, SmolStr::from(i.to_string()));
        }
        assert_eq!(arena.len(), 1000);
        for i in 0..1000 {
            assert_eq!(arena.get(i).map(SmolStr::as_str), Some(i.to_string().as_str()));
        }
        assert_eq!(arena.get(1000), None);
        let last = ARENA_FIRST * ((1 << ARENA_CHUNKS) - 1) - 1;
        assert_eq!(
            Arena::locate(last),
            (ARENA_CHUNKS - 1, (ARENA_FIRST << (ARENA_CHUNKS - 1)) - 1)
        );
        assert!(last >= 1 << 16, "every id fits");
    }
}
