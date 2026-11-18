//! The value views a reader hands out, the same for BCF and VCF text
//! (`r[vcf_values.access]`): [`OptValues`] yields each element as `Option<T>`,
//! [`Values`] promises none is missing. Both are lazy, allocation-free,
//! `ExactSizeIterator`s with `get(i)`.

use super::bcf_value::{BcfNumber, BcfVec};

/// A number a VCF or BCF field holds: `i32` or `f32`.
pub(crate) trait Number: BcfNumber {}
impl Number for i32 {}
impl Number for f32 {}

/// A vector with an element missing where [`Values`] promised none.
#[derive(Debug, Clone, Copy, PartialEq, Eq, thiserror::Error)]
#[error("element {index} of a vector that must not have missing elements is missing")]
pub struct MissingElement {
    pub index: usize,
}

/// A number vector, missing elements as `None`.
#[derive(Debug, Clone)]
pub(crate) enum OptValues<'a, T> {
    Bcf(BcfVec<'a, T>),
}

impl<'a, T: Number> OptValues<'a, T> {
    /// Element `i` of those not yet yielded; `None` past the end.
    #[expect(
        clippy::option_option,
        reason = "in range, then the element as the iterator yields it"
    )]
    pub(crate) fn get(&self, i: usize) -> Option<Option<T>> {
        match self {
            Self::Bcf(v) => v.get(i),
        }
    }

    /// No elements, or a lone missing one: how both formats spell a field
    /// that is `.` as a whole.
    pub(crate) fn is_missing(&self) -> bool {
        match self.len() {
            0 => true,
            1 => matches!(self.get(0), Some(None)),
            _ => false,
        }
    }

    // r[impl vcf_values.access]
    /// The same elements, checked to have none missing.
    pub(crate) fn all_present(self) -> Result<Values<'a, T>, MissingElement> {
        match self.clone().position(|v| v.is_none()) {
            Some(index) => Err(MissingElement { index }),
            None => Ok(Values(self)),
        }
    }
}

impl<T: Number> Iterator for OptValues<'_, T> {
    type Item = Option<T>;

    #[inline]
    fn next(&mut self) -> Option<Option<T>> {
        match self {
            Self::Bcf(v) => v.next(),
        }
    }

    fn size_hint(&self) -> (usize, Option<usize>) {
        match self {
            Self::Bcf(v) => v.size_hint(),
        }
    }
}

impl<T: Number> ExactSizeIterator for OptValues<'_, T> {}

/// A number vector with no missing element (checked by [`OptValues::all_present`]).
#[derive(Debug, Clone)]
pub(crate) struct Values<'a, T>(OptValues<'a, T>);

impl<T: Number> Values<'_, T> {
    /// Element `i` of those not yet yielded; `None` past the end.
    pub(crate) fn get(&self, i: usize) -> Option<T> {
        self.0.get(i).flatten()
    }
}

impl<T: Number> Iterator for Values<'_, T> {
    type Item = T;

    #[inline]
    fn next(&mut self) -> Option<T> {
        let v = self.0.next()?;
        debug_assert!(v.is_some(), "all_present checked every element");
        v
    }

    fn size_hint(&self) -> (usize, Option<usize>) {
        self.0.size_hint()
    }
}

impl<T: Number> ExactSizeIterator for Values<'_, T> {}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::vcf::reader::bcf_value::TypedValue;

    fn bcf(bytes: &[u8]) -> OptValues<'_, i32> {
        OptValues::Bcf(TypedValue::parse(bytes).unwrap().0.ints().unwrap())
    }

    // r[verify vcf_values.access]
    #[test]
    fn all_present_rejects_the_first_missing_element() {
        assert_eq!(
            bcf(&[0x31, 1, 0x80, 0x80]).all_present().unwrap_err(),
            MissingElement { index: 1 }
        );
        let v = bcf(&[0x31, 1, 2, 0x81]).all_present().unwrap();
        assert_eq!((v.len(), v.get(1), v.get(2)), (2, Some(2), None));
        assert_eq!(v.collect::<Vec<_>>(), [1, 2]);
    }

    #[test]
    fn a_lone_missing_element_or_nothing_is_missing() {
        assert!(bcf(&[0x01]).is_missing());
        assert!(bcf(&[0x11, 0x80]).is_missing());
        assert!(bcf(&[0x21, 0x80, 0x81]).is_missing());
        assert!(!bcf(&[0x21, 0x80, 0x80]).is_missing());
        assert!(!bcf(&[0x11, 0]).is_missing());
    }
}
