//! The banded pair-HMM on a GPU, through wgpu compute shaders (Metal, Vulkan,
//! DX12). Behind the `gpu` feature.
//!
//! One invocation scores one (read, haplotype) pair, sweeping the read's rows
//! over the band with the batch kernel's arithmetic: the same recurrence in
//! the same order, the same flush of subnormal cells, the same power-of-two
//! renormalisation every eight rows, and the same `f64` finish on the CPU. So
//! a score is the CPU kernels' score **up to what the GPU compiler does to
//! `f32` arithmetic**, which is not nothing: Metal compiles with fast-math and
//! SPIR-V permits contracting `a * b + c` into a fused multiply-add, so bit
//! parity with [`align_strips`](crate::align_strips) is not promised, and the
//! oracle for this kernel is a tolerance, not equality. Measured, it is
//! bit-identical to the strip kernel on an Apple M4 Pro and an AMD RX 5700 XT
//! except for pairs scoring below log10 −40, which a GPU that flushes
//! subnormal intermediates can score lower.
//!
//! The one thing the CPU entry points do after their kernels -- rescoring a
//! pair that finished below the floor `f32` can vouch for, in `f64` -- needs
//! the pair's inputs, which a launch does not keep. So a score comes back as
//! a [`GpuScore`], carrying the floor [`GpuPairs`] worked out when the pair was
//! pushed: [`GpuScore::trusted`] is the score wherever a CPU entry point would
//! return the kernel's score unchanged (all but about 1% of real pairs), and
//! [`GpuScore::or_rescore`] takes the pair to score the rest.
//!
//! ```no_run
//! use compair::gpu::{GpuAligner, GpuContext, GpuPairs};
//! # fn main() -> Result<(), Box<dyn std::error::Error>> {
//! # let (haplotypes, reads): (Vec<compair::Haplotype>, Vec<(compair::Read, compair::Band)>) = (vec![], vec![]);
//! let emission = compair::StandardEmission::default();
//! let mut aligner = GpuAligner::new(GpuContext::new()?);
//! let mut pairs = GpuPairs::new();
//! // Every read against every haplotype, read-major as `Workspace::align_reads`
//! // lays them out; `range` is where this locus's scores land in the launch.
//! let range = pairs.push_reads(&haplotypes, &reads, &emission)?;
//! // In push order; `submit` returns once the launch is queued.
//! let scores = aligner.submit(&pairs)?.collect()?;
//! let per_read = scores.get(range).unwrap_or_default().chunks_exact(haplotypes.len());
//! for ((read, band), scores) in reads.iter().zip(per_read) {
//!     for (haplotype, score) in haplotypes.iter().zip(scores) {
//!         // What `Workspace::align_reads` returns for this pair.
//!         let score = score.or_rescore(haplotype, read, &emission, *band);
//!     }
//! }
//! # Ok(()) }
//! ```
//!
//! `examples/gpu.rs` runs this against the CPU, falling back to it when there
//! is no adapter.

mod device;
mod emulate;
mod plan;
mod shader;

use std::time::Duration;

use seqair_types::Strand;

use crate::{Band, Emission, Haplotype, Log10Likelihood, Read};

pub use device::{GpuAligner, GpuContext, Handle, KernelOptions, Wait};
#[doc(hidden)]
pub use emulate::Subnormals;
pub use plan::{GpuPairs, HaplotypeSlot, ReadSlot};
#[doc(hidden)]
pub use shader::{Contraction, Style, Variant};

/// One pair's score from the GPU kernel, and the floor below which its `f32`
/// arithmetic cannot vouch for it.
///
/// The floor is the one every CPU entry point applies after its kernel (see
/// the crate docs' **Precision**), worked out from the pair when it was
/// pushed, so a score passes or fails it exactly where the CPU's would.
#[derive(Debug, Clone, Copy, PartialEq)]
#[must_use = "a GPU score is not the pair's score until it is trusted or rescored"]
pub struct GpuScore {
    score: Log10Likelihood,
    floor: Option<f64>,
}

impl GpuScore {
    pub(crate) fn new(score: Log10Likelihood, floor: Option<f64>) -> Self {
        Self { score, floor }
    }

    /// The score, where a CPU entry point would return the kernel's score as
    /// it is; `None` where it would score the pair again in `f64`.
    #[must_use]
    pub fn trusted(self) -> Option<Log10Likelihood> {
        self.floor.is_some_and(|floor| self.score.get() >= floor).then_some(self.score)
    }

    /// The pair's score as the CPU entry points return it: [`Self::trusted`],
    /// or the pair scored again by [`align_banded_f64`](crate::align_banded_f64).
    /// The inputs must be the ones the pair was pushed with.
    pub fn or_rescore<E: Emission + ?Sized>(
        self,
        haplotype: &Haplotype,
        read: &Read,
        emission: &E,
        band: Band,
    ) -> Log10Likelihood {
        self.trusted().unwrap_or_else(|| crate::align_banded_f64(haplotype, read, emission, band))
    }

    /// The kernel's own score, before the floor: for comparing it with the
    /// CPU kernels, not for calling with.
    #[must_use]
    pub const fn untrusted(self) -> Log10Likelihood {
        self.score
    }
}

/// Everything the GPU path can fail with.
#[derive(Debug, thiserror::Error)]
pub enum GpuError {
    #[error("no GPU adapter found")]
    NoAdapter(#[source] wgpu::RequestAdapterError),

    #[error("creating the GPU device failed")]
    Device(#[from] wgpu::RequestDeviceError),

    #[error("the kernel for {variant:?} did not compile")]
    Shader {
        variant: Variant,
        #[source]
        source: Box<wgpu::Error>,
    },

    #[error("writing the kernel's source failed")]
    ShaderSource(#[from] core::fmt::Error),

    #[error("the GPU did not finish {step} within {timeout:?}")]
    Timeout { step: &'static str, timeout: Duration },

    #[error("waiting for the GPU failed")]
    Poll(#[source] wgpu::PollError),

    #[error("mapping the scores for reading failed")]
    Map(#[from] wgpu::BufferAsyncError),

    #[error("the scores were not mapped after the GPU finished")]
    MapPending,

    #[error("reading the mapped scores failed")]
    MapRange(#[from] wgpu::MapRangeError),

    #[error("the {what} do not fit the kernel's 32-bit indices")]
    TooLarge { what: &'static str },

    #[error("the {what} buffer needs {bytes} bytes, over the device's binding limit of {limit}")]
    BufferTooLarge { what: &'static str, bytes: u64, limit: u64 },

    #[error("a slot that this `GpuPairs` did not hand out")]
    UnknownSlot,

    #[error("a read on {read:?} paired with a haplotype pushed for {haplotype:?}")]
    StrandMismatch { read: Strand, haplotype: Strand },

    #[error("the shader cache's lock was poisoned by a panic")]
    Poisoned,
}
