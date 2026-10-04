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
//! Scores come back without their pairs' inputs, so the one thing the CPU
//! entry points do after their kernels -- rescoring a pair that finished
//! below `f32`'s floor in `f64` -- is the caller's: pass each score through
//! [`trusted`](crate::trusted).
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
//! // Then each through `compair::trusted` with its pair, as the CPU does.
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

pub use device::{GpuAligner, GpuContext, Handle, KernelOptions, Wait};
#[doc(hidden)]
pub use emulate::Subnormals;
pub use plan::{GpuPairs, HaplotypeSlot, ReadSlot};
#[doc(hidden)]
pub use shader::{Contraction, Style, Variant};

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
