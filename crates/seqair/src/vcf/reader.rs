//! VCF/BCF reading (`docs/spec/2-vcf-1-reader.md`).
//!
//! So far only the value decoders a reader is built on: BCF typed values
//! ([`bcf_value`]), VCF text values ([`text_value`]), and the views both
//! backends share ([`values`]). None of them knows a `VcfHeader`.

#![cfg_attr(not(test), expect(dead_code, reason = "used by the VCF reader in a later commit"))]

pub(crate) mod bcf_value;
pub(crate) mod text_value;
pub(crate) mod values;
