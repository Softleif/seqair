//! The device, the compiled shaders, and one launch at a time.

use std::{
    collections::HashMap,
    sync::{Arc, Mutex, mpsc},
    time::{Duration, Instant},
};

use bytemuck::{Pod, Zeroable};

use super::{
    GpuError,
    plan::{GpuPairs, PairRecord, TRANSITIONS},
    shader::{self, Contraction, Style, Variant},
};
use crate::types::Log10Likelihood;

/// A GPU device and the shaders compiled for it, shared by every
/// [`GpuAligner`] built on it.
///
/// Shaders are compiled on first use, one per band class and kernel variant,
/// and kept; see [`GpuAligner`] for what a class is.
pub struct GpuContext {
    device: wgpu::Device,
    queue: wgpu::Queue,
    info: wgpu::AdapterInfo,
    /// The adapter maps storage buffers into host memory (Apple Silicon,
    /// integrated GPUs), so results are read in place rather than copied to a
    /// staging buffer first.
    uma: bool,
    /// The device can time a pass; see [`KernelOptions::timestamps`].
    timestamps: bool,
    limits: wgpu::Limits,
    layout: wgpu::BindGroupLayout,
    pipeline_layout: wgpu::PipelineLayout,
    pipelines: Mutex<HashMap<Variant, wgpu::ComputePipeline>>,
    /// `plan::TRANSITIONS`, uploaded once.
    transitions: wgpu::Buffer,
}

impl core::fmt::Debug for GpuContext {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        f.debug_struct("GpuContext")
            .field("adapter", &self.info.name)
            .field("uma", &self.uma)
            .finish()
    }
}

impl GpuContext {
    /// The highest-performance adapter wgpu finds, with the limits raised to
    /// what it supports for storage buffers.
    pub fn new() -> Result<Arc<Self>, GpuError> {
        pollster::block_on(Self::new_async())
    }

    async fn new_async() -> Result<Arc<Self>, GpuError> {
        let instance = wgpu::Instance::new(&wgpu::InstanceDescriptor::default());
        let adapter = instance
            .request_adapter(&wgpu::RequestAdapterOptions {
                power_preference: wgpu::PowerPreference::HighPerformance,
                compatible_surface: None,
                force_fallback_adapter: false,
            })
            .await
            .map_err(GpuError::NoAdapter)?;
        let info = adapter.get_info();
        // Vulkan reports mappable primary buffers on a discrete card too, as
        // host-visible device memory; reading results out of VRAM over the bus
        // is not what this path is for, so it is taken on shared memory only.
        let uma = adapter.features().contains(wgpu::Features::MAPPABLE_PRIMARY_BUFFERS)
            && matches!(info.device_type, wgpu::DeviceType::IntegratedGpu | wgpu::DeviceType::Cpu);
        let timestamps = adapter.features().contains(wgpu::Features::TIMESTAMP_QUERY);
        let supported = adapter.limits();
        let required_limits = wgpu::Limits {
            max_storage_buffer_binding_size: supported.max_storage_buffer_binding_size,
            max_buffer_size: supported.max_buffer_size,
            ..wgpu::Limits::default()
        };
        let (device, queue) = adapter
            .request_device(&wgpu::DeviceDescriptor {
                label: Some("compair::gpu"),
                required_features: if uma {
                    wgpu::Features::MAPPABLE_PRIMARY_BUFFERS
                } else {
                    wgpu::Features::empty()
                } | if timestamps {
                    wgpu::Features::TIMESTAMP_QUERY
                } else {
                    wgpu::Features::empty()
                },
                required_limits,
                ..Default::default()
            })
            .await?;
        // wgpu's default for an error outside a scope is to panic; a caller
        // scoring reads would rather have it logged, and every call this crate
        // makes that can fail on the device is inside a scope anyway.
        device.on_uncaptured_error(Arc::new(|error| {
            tracing::error!(%error, "uncaptured wgpu error");
        }));
        tracing::debug!(adapter = %info.name, backend = ?info.backend, uma, "compair GPU context");
        let limits = device.limits();

        let storage = |binding, read_only| wgpu::BindGroupLayoutEntry {
            binding,
            visibility: wgpu::ShaderStages::COMPUTE,
            ty: wgpu::BindingType::Buffer {
                ty: wgpu::BufferBindingType::Storage { read_only },
                has_dynamic_offset: false,
                min_binding_size: None,
            },
            count: None,
        };
        let layout = device.create_bind_group_layout(&wgpu::BindGroupLayoutDescriptor {
            label: Some("compair::gpu::layout"),
            entries: &[
                storage(0, true),
                storage(1, true),
                storage(2, true),
                storage(3, false),
                storage(5, true),
                wgpu::BindGroupLayoutEntry {
                    binding: 4,
                    visibility: wgpu::ShaderStages::COMPUTE,
                    ty: wgpu::BindingType::Buffer {
                        ty: wgpu::BufferBindingType::Uniform,
                        has_dynamic_offset: true,
                        min_binding_size: wgpu::BufferSize::new(PARAMS_BYTES),
                    },
                    count: None,
                },
            ],
        });
        let pipeline_layout = device.create_pipeline_layout(&wgpu::PipelineLayoutDescriptor {
            label: Some("compair::gpu::pipeline_layout"),
            bind_group_layouts: &[&layout],
            immediate_size: 0,
        });
        let transitions = wgpu::util::DeviceExt::create_buffer_init(
            &device,
            &wgpu::util::BufferInitDescriptor {
                label: Some("compair::gpu::transitions"),
                contents: bytemuck::cast_slice(&TRANSITIONS),
                usage: wgpu::BufferUsages::STORAGE,
            },
        );
        Ok(Arc::new(Self {
            transitions,
            device,
            queue,
            info,
            uma,
            timestamps,
            limits,
            layout,
            pipeline_layout,
            pipelines: Mutex::new(HashMap::new()),
        }))
    }

    /// What wgpu reports about the adapter.
    #[must_use]
    pub fn adapter_info(&self) -> &wgpu::AdapterInfo {
        &self.info
    }

    /// Whether results are read in place rather than through a staging copy.
    #[must_use]
    pub fn is_uma(&self) -> bool {
        self.uma
    }

    /// A four-byte copy and its readback, waited for as `wait` says: the
    /// floor under any launch's latency.
    #[doc(hidden)]
    pub fn round_trip(&self, wait: Wait) -> Result<(), GpuError> {
        let buffer = |usage| {
            self.device.create_buffer(&wgpu::BufferDescriptor {
                label: None,
                size: 4,
                usage,
                mapped_at_creation: false,
            })
        };
        let source = buffer(wgpu::BufferUsages::COPY_SRC);
        let target = buffer(wgpu::BufferUsages::COPY_DST | wgpu::BufferUsages::MAP_READ);
        let mut encoder =
            self.device.create_command_encoder(&wgpu::CommandEncoderDescriptor { label: None });
        encoder.copy_buffer_to_buffer(&source, 0, &target, 0, 4);
        let submission = self.queue.submit(Some(encoder.finish()));
        map_and_wait(self, wait, target.slice(..), submission, Duration::from_secs(10))
    }

    /// The compiled pipeline for `variant`, compiling it on first use.
    fn pipeline(&self, variant: Variant) -> Result<wgpu::ComputePipeline, GpuError> {
        let mut cache = self.pipelines.lock().map_err(|_| GpuError::Poisoned)?;
        if let Some(pipeline) = cache.get(&variant) {
            return Ok(pipeline.clone());
        }
        let source = shader::source(variant)?;
        let scope = self.device.push_error_scope(wgpu::ErrorFilter::Validation);
        let module = self.device.create_shader_module(wgpu::ShaderModuleDescriptor {
            label: Some("compair::gpu::kernel"),
            source: wgpu::ShaderSource::Wgsl(source.into()),
        });
        let pipeline = self.device.create_compute_pipeline(&wgpu::ComputePipelineDescriptor {
            label: Some("compair::gpu::kernel"),
            layout: Some(&self.pipeline_layout),
            module: &module,
            entry_point: Some("main"),
            compilation_options: wgpu::PipelineCompilationOptions::default(),
            cache: None,
        });
        if let Some(error) = pollster::block_on(scope.pop()) {
            return Err(GpuError::Shader { variant, source: Box::new(error) });
        }
        tracing::debug!(?variant, "compiled compair GPU kernel");
        cache.insert(variant, pipeline.clone());
        Ok(pipeline)
    }
}

/// How the kernel is compiled. The default is the fastest variant that is
/// bit-identical to the CPU kernels on both GPUs measured; the others exist
/// for the notes' experiments.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
#[doc(hidden)]
pub struct KernelOptions {
    /// `None` picks per backend and band class; see [`auto_style`].
    pub style: Option<Style>,
    pub contraction: Contraction,
    pub workgroup_size: u32,
    pub wait: Wait,
    /// Time each launch's compute pass on the GPU, where the device can;
    /// [`GpuAligner::last_kernel_time`] reports it.
    pub timestamps: bool,
}

impl Default for KernelOptions {
    fn default() -> Self {
        Self {
            style: None,
            contraction: Contraction::Blocked,
            workgroup_size: 64,
            wait: Wait::Spin,
            timestamps: false,
        }
    }
}

/// How [`Handle::collect`] waits for the GPU.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
#[doc(hidden)]
pub enum Wait {
    /// wgpu's blocking wait. On Metal (wgpu-hal 28) that checks the command
    /// buffer's status and sleeps a millisecond between checks, which puts a
    /// 1 ms floor under every launch.
    Block,
    /// Non-blocking polls with a `yield_now` between them, until the scores
    /// are mapped: the launch's own latency, for a core's worth of polling.
    Spin,
}

/// The widest band class the unrolled style is generated for. Past it the
/// shader's size and the driver's compile time grow with the width while the
/// row buffer spills out of registers anyway (on RDNA from a half-width of
/// ~40), so the loop style is no worse.
const UNROLL_LIMIT: u32 = 64;

/// The style a class runs in when the options leave it open: the loop on
/// Metal, whose compiler keeps its dynamically indexed arrays cheap, and the
/// unrolled band elsewhere, because RADV lowers a dynamic index into the
/// register file to a chain of branches and runs the loop 3-5x slower. See
/// `docs/notes/gpu.md` for the measurements.
fn auto_style(backend: wgpu::Backend, half_width: u32) -> Style {
    match backend {
        wgpu::Backend::Metal => Style::Loop,
        _ if half_width <= UNROLL_LIMIT => Style::Unrolled,
        _ => Style::Loop,
    }
}

/// The uniform each dispatch reads: which slice of the pairs it runs.
#[repr(C)]
#[derive(Debug, Clone, Copy, Pod, Zeroable)]
struct Params {
    first: u32,
    count: u32,
    /// Always zero, but a uniform, so the compiler cannot fold it away; see
    /// [`Contraction::Blocked`].
    zero: u32,
    pad: u32,
}

const PARAMS_BYTES: u64 = size_of::<Params>() as u64;

/// A buffer that grows to the largest launch it has held.
#[derive(Debug)]
struct Grown {
    buffer: wgpu::Buffer,
    capacity: u64,
}

/// The per-aligner buffers, in the shader's binding order.
#[derive(Debug, Default)]
struct Buffers {
    pairs: Option<Grown>,
    rows: Option<Grown>,
    weights: Option<Grown>,
    scores: Option<Grown>,
    params: Option<Grown>,
    staging: Option<Grown>,
}

/// Scores [`GpuPairs`] on a [`GpuContext`], one launch at a time.
///
/// A launch groups its pairs by **band class** -- the half-width rounded up to
/// a multiple of four -- and dispatches each class through a shader whose
/// row buffer is sized for it, because WGSL has no private arrays of run-time
/// size. A pair narrower than its class computes the extra positions masked.
///
/// The aligner owns its buffers, so a launch is [`submit`](Self::submit) then
/// [`Handle::collect`], and the handle borrows the aligner: a second launch
/// cannot overwrite the buffers of one still in flight. Build a second aligner
/// on the same context to overlap filling one launch with running another.
#[derive(Debug)]
pub struct GpuAligner {
    context: Arc<GpuContext>,
    options: KernelOptions,
    timeout: Duration,
    buffers: Buffers,
    /// Launch slot to pushed-pair index, and the records in launch order.
    order: Vec<u32>,
    records: Vec<PairRecord>,
    params: Vec<u8>,
    dispatches: Vec<(Variant, u32, u32)>,
    timing: Option<Timing>,
    last_kernel_time: Option<Duration>,
}

/// The query set and buffers a timed launch writes its two timestamps into.
#[derive(Debug)]
struct Timing {
    queries: wgpu::QuerySet,
    resolve: wgpu::Buffer,
    readback: wgpu::Buffer,
}

impl Timing {
    fn new(device: &wgpu::Device) -> Self {
        let size = 2 * size_of::<u64>() as u64;
        Self {
            queries: device.create_query_set(&wgpu::QuerySetDescriptor {
                label: Some("compair::gpu::timing"),
                ty: wgpu::QueryType::Timestamp,
                count: 2,
            }),
            resolve: device.create_buffer(&wgpu::BufferDescriptor {
                label: Some("compair::gpu::timing_resolve"),
                size,
                usage: wgpu::BufferUsages::QUERY_RESOLVE | wgpu::BufferUsages::COPY_SRC,
                mapped_at_creation: false,
            }),
            readback: device.create_buffer(&wgpu::BufferDescriptor {
                label: Some("compair::gpu::timing_readback"),
                size,
                usage: wgpu::BufferUsages::MAP_READ | wgpu::BufferUsages::COPY_DST,
                mapped_at_creation: false,
            }),
        }
    }
}

impl GpuAligner {
    #[must_use]
    pub fn new(context: Arc<GpuContext>) -> Self {
        Self::with_options(context, KernelOptions::default())
    }

    #[doc(hidden)]
    #[must_use]
    pub fn with_options(context: Arc<GpuContext>, options: KernelOptions) -> Self {
        let timing =
            (options.timestamps && context.timestamps).then(|| Timing::new(&context.device));
        Self {
            context,
            options,
            timeout: Duration::from_secs(10),
            buffers: Buffers::default(),
            order: Vec::new(),
            records: Vec::new(),
            params: Vec::new(),
            dispatches: Vec::new(),
            timing,
            last_kernel_time: None,
        }
    }

    /// The GPU's own time for the last collected launch's compute pass, when
    /// [`KernelOptions::timestamps`] was set and the device can time one.
    #[doc(hidden)]
    #[must_use]
    pub fn last_kernel_time(&self) -> Option<Duration> {
        self.last_kernel_time
    }

    /// How long [`Handle::collect`] waits for the GPU before giving up.
    #[must_use]
    pub fn with_timeout(mut self, timeout: Duration) -> Self {
        self.timeout = timeout;
        self
    }

    #[must_use]
    pub fn context(&self) -> &Arc<GpuContext> {
        &self.context
    }

    /// Scores every pair and waits for the result.
    pub fn align(&mut self, pairs: &GpuPairs) -> Result<Vec<Log10Likelihood>, GpuError> {
        self.submit(pairs)?.collect()
    }

    /// Uploads the pairs and starts the launch; the scores come from the
    /// returned handle's [`collect`](Handle::collect).
    pub fn submit(&mut self, pairs: &GpuPairs) -> Result<Handle<'_>, GpuError> {
        let context = Arc::clone(&self.context);
        let workgroup = self.options.workgroup_size;

        // Launch order: grouped by class, stable within one, so each class is
        // one contiguous slice of the pair buffer.
        self.order.clear();
        self.order.extend(
            pairs
                .pairs
                .iter()
                .enumerate()
                .filter(|(_, pending)| pending.class.is_some())
                .filter_map(|(index, _)| u32::try_from(index).ok()),
        );
        self.order.sort_by_key(|&index| pairs.pairs.get(index as usize).and_then(|p| p.class));
        self.records.clear();
        self.records.extend(
            self.order
                .iter()
                .filter_map(|&index| pairs.pairs.get(index as usize))
                .map(|p| p.record),
        );

        // One dispatch per class, split where a class outgrows one dispatch's
        // workgroup count.
        self.dispatches.clear();
        let per_dispatch =
            context.limits.max_compute_workgroups_per_dimension.saturating_mul(workgroup);
        let mut start = 0usize;
        while let Some(&first) = self.order.get(start) {
            let class = pairs.pairs.get(first as usize).and_then(|p| p.class).unwrap_or(0);
            let run = self
                .order
                .get(start..)
                .unwrap_or_default()
                .iter()
                .take_while(|&&index| {
                    pairs.pairs.get(index as usize).and_then(|p| p.class) == Some(class)
                })
                .count();
            let variant = Variant {
                half_width: class,
                style: self
                    .options
                    .style
                    .unwrap_or_else(|| auto_style(context.info.backend, class)),
                contraction: self.options.contraction,
                workgroup_size: workgroup,
            };
            let mut done = 0usize;
            while done < run {
                let count = (run - done).min(per_dispatch as usize);
                let first_slot = u32::try_from(start + done)
                    .map_err(|_| GpuError::TooLarge { what: "pairs" })?;
                let count_u32 =
                    u32::try_from(count).map_err(|_| GpuError::TooLarge { what: "pairs" })?;
                self.dispatches.push((variant, first_slot, count_u32));
                done += count;
            }
            start += run;
        }

        let submission = if self.dispatches.is_empty() { None } else { Some(self.encode(pairs)?) };
        Ok(Handle { aligner: self, submission, len: pairs.pairs.len() })
    }

    /// Uploads, binds and dispatches `self.dispatches`, and queues the
    /// readback copy on a discrete GPU.
    fn encode(&mut self, pairs: &GpuPairs) -> Result<wgpu::SubmissionIndex, GpuError> {
        let context = Arc::clone(&self.context);
        let device = &context.device;
        let queue = &context.queue;
        let alignment =
            u64::from(context.limits.min_uniform_buffer_offset_alignment).max(PARAMS_BYTES);
        let stride =
            usize::try_from(alignment).map_err(|_| GpuError::TooLarge { what: "params" })?;

        self.params.clear();
        self.params.resize(self.dispatches.len() * stride, 0);
        for (slot, &(_, first, count)) in self.dispatches.iter().enumerate() {
            let params = Params { first, count, zero: 0, pad: 0 };
            if let Some(target) =
                self.params.get_mut(slot * stride..slot * stride + size_of::<Params>())
            {
                target.copy_from_slice(bytemuck::bytes_of(&params));
            }
        }

        let scores_bytes = bytes(self.records.len(), size_of::<[u32; 2]>());
        let storage_limit = u64::from(context.limits.max_storage_buffer_binding_size);
        let copy = wgpu::BufferUsages::COPY_DST;
        let storage = wgpu::BufferUsages::STORAGE;
        let uploads: [(&str, &[u8]); 3] = [
            ("pairs", bytemuck::cast_slice(&self.records)),
            ("rows", bytemuck::cast_slice(&pairs.rows)),
            ("weights", bytemuck::cast_slice(&pairs.weights)),
        ];
        for (what, data) in uploads {
            let bytes = data.len() as u64;
            if bytes > storage_limit {
                return Err(GpuError::BufferTooLarge { what, bytes, limit: storage_limit });
            }
        }
        if scores_bytes > storage_limit {
            return Err(GpuError::BufferTooLarge {
                what: "scores",
                bytes: scores_bytes,
                limit: storage_limit,
            });
        }

        let [(_, pair_bytes), (_, row_bytes), (_, weight_bytes)] = uploads;
        let pairs_buffer =
            grow(device, &mut self.buffers.pairs, "pairs", pair_bytes.len() as u64, storage | copy);
        queue.write_buffer(pairs_buffer, 0, pair_bytes);
        let rows_buffer =
            grow(device, &mut self.buffers.rows, "rows", row_bytes.len() as u64, storage | copy);
        queue.write_buffer(rows_buffer, 0, row_bytes);
        let weights_buffer = grow(
            device,
            &mut self.buffers.weights,
            "weights",
            weight_bytes.len() as u64,
            storage | copy,
        );
        queue.write_buffer(weights_buffer, 0, weight_bytes);
        let params_buffer = grow(
            device,
            &mut self.buffers.params,
            "params",
            self.params.len() as u64,
            wgpu::BufferUsages::UNIFORM | copy,
        );
        queue.write_buffer(params_buffer, 0, &self.params);
        let scores_usage = if context.uma {
            storage | wgpu::BufferUsages::MAP_READ
        } else {
            storage | wgpu::BufferUsages::COPY_SRC
        };
        let scores_buffer =
            grow(device, &mut self.buffers.scores, "scores", scores_bytes, scores_usage);

        let bind_group = device.create_bind_group(&wgpu::BindGroupDescriptor {
            label: Some("compair::gpu::bind_group"),
            layout: &context.layout,
            entries: &[
                wgpu::BindGroupEntry {
                    binding: 0,
                    resource: whole(pairs_buffer, pair_bytes.len() as u64),
                },
                wgpu::BindGroupEntry {
                    binding: 1,
                    resource: whole(rows_buffer, row_bytes.len() as u64),
                },
                wgpu::BindGroupEntry {
                    binding: 2,
                    resource: whole(weights_buffer, weight_bytes.len() as u64),
                },
                wgpu::BindGroupEntry {
                    binding: 3,
                    resource: whole(scores_buffer, scores_bytes.max(4)),
                },
                wgpu::BindGroupEntry {
                    binding: 5,
                    resource: context.transitions.as_entire_binding(),
                },
                wgpu::BindGroupEntry {
                    binding: 4,
                    resource: wgpu::BindingResource::Buffer(wgpu::BufferBinding {
                        buffer: params_buffer,
                        offset: 0,
                        size: wgpu::BufferSize::new(PARAMS_BYTES),
                    }),
                },
            ],
        });

        let mut encoder = device.create_command_encoder(&wgpu::CommandEncoderDescriptor {
            label: Some("compair::gpu"),
        });
        {
            let mut pass = encoder.begin_compute_pass(&wgpu::ComputePassDescriptor {
                label: Some("compair::gpu::pass"),
                timestamp_writes: self.timing.as_ref().map(|timing| {
                    wgpu::ComputePassTimestampWrites {
                        query_set: &timing.queries,
                        beginning_of_pass_write_index: Some(0),
                        end_of_pass_write_index: Some(1),
                    }
                }),
            });
            for (slot, &(variant, _, count)) in self.dispatches.iter().enumerate() {
                let pipeline = context.pipeline(variant)?;
                pass.set_pipeline(&pipeline);
                let offset = u32::try_from(slot * stride)
                    .map_err(|_| GpuError::TooLarge { what: "params" })?;
                pass.set_bind_group(0, &bind_group, &[offset]);
                pass.dispatch_workgroups(count.div_ceil(variant.workgroup_size), 1, 1);
            }
        }
        if !context.uma {
            let staging = grow(
                device,
                &mut self.buffers.staging,
                "staging",
                scores_bytes,
                wgpu::BufferUsages::MAP_READ | copy,
            );
            encoder.copy_buffer_to_buffer(scores_buffer, 0, staging, 0, scores_bytes);
        }
        if let Some(timing) = &self.timing {
            encoder.resolve_query_set(&timing.queries, 0..2, &timing.resolve, 0);
            encoder.copy_buffer_to_buffer(
                &timing.resolve,
                0,
                &timing.readback,
                0,
                2 * size_of::<u64>() as u64,
            );
        }
        Ok(queue.submit(Some(encoder.finish())))
    }
}

/// The first `size` bytes of `buffer` as a binding.
fn whole(buffer: &wgpu::Buffer, size: u64) -> wgpu::BindingResource<'_> {
    wgpu::BindingResource::Buffer(wgpu::BufferBinding {
        buffer,
        offset: 0,
        size: wgpu::BufferSize::new(size),
    })
}

/// `count` items of `size` bytes, as a buffer size. A slice's length in bytes
/// always fits a `u64`.
fn bytes(count: usize, size: usize) -> u64 {
    (count.saturating_mul(size)) as u64
}

/// The buffer in `slot`, replaced by one of at least `size` bytes (and never
/// less than four, since wgpu rejects empty bindings) if it is smaller.
fn grow<'a>(
    device: &wgpu::Device,
    slot: &'a mut Option<Grown>,
    label: &'static str,
    size: u64,
    usage: wgpu::BufferUsages,
) -> &'a wgpu::Buffer {
    let size = size.max(4).next_multiple_of(4);
    let grown = match slot.take() {
        Some(grown) if grown.capacity >= size => grown,
        _ => {
            // Half again as much, so a caller whose launches creep upwards
            // does not reallocate every time.
            let capacity = size.saturating_add(size / 2).next_multiple_of(4);
            let buffer = device.create_buffer(&wgpu::BufferDescriptor {
                label: Some(label),
                size: capacity,
                usage,
                mapped_at_creation: false,
            });
            Grown { buffer, capacity }
        }
    };
    &slot.insert(grown).buffer
}

/// A launch in flight. [`collect`](Self::collect) waits for it.
#[must_use = "a launch's scores are only read by `collect`"]
pub struct Handle<'a> {
    aligner: &'a mut GpuAligner,
    submission: Option<wgpu::SubmissionIndex>,
    len: usize,
}

impl Handle<'_> {
    /// Waits for the GPU and returns one score per pushed pair, in push order.
    pub fn collect(self) -> Result<Vec<Log10Likelihood>, GpuError> {
        let mut out = vec![Log10Likelihood::IMPOSSIBLE; self.len];
        let Some(submission) = self.submission else {
            return Ok(out);
        };
        let aligner = self.aligner;
        let context = &aligner.context;
        let timeout = aligner.timeout;
        let bytes = bytes(aligner.records.len(), size_of::<[u32; 2]>());
        let buffer = if context.uma { &aligner.buffers.scores } else { &aligner.buffers.staging };
        let buffer = buffer.as_ref().map(|grown| &grown.buffer).ok_or(GpuError::MapPending)?;
        let slice = buffer.slice(..bytes);
        map_and_wait(context, aligner.options.wait, slice, submission.clone(), timeout)?;
        {
            let view = slice.get_mapped_range();
            let scores: &[[u32; 2]] = bytemuck::cast_slice(&view);
            for (&index, &[sum, exponent]) in aligner.order.iter().zip(scores) {
                if let Some(slot) = out.get_mut(index as usize) {
                    *slot = finish(f32::from_bits(sum), i32::from_ne_bytes(exponent.to_ne_bytes()));
                }
            }
        }
        buffer.unmap();
        aligner.last_kernel_time = None;
        if let Some(timing) = &aligner.timing {
            let slice = timing.readback.slice(..);
            map_and_wait(context, aligner.options.wait, slice, submission, timeout)?;
            {
                let view = slice.get_mapped_range();
                if let [begin, end] = bytemuck::cast_slice::<u8, u64>(&view) {
                    let ticks = end.saturating_sub(*begin);
                    #[allow(
                        clippy::cast_precision_loss,
                        reason = "a kernel's tick count, for a diagnostic"
                    )]
                    let nanos = ticks as f64 * f64::from(context.queue.get_timestamp_period());
                    aligner.last_kernel_time = Some(Duration::from_secs_f64(nanos / 1e9));
                }
            }
            timing.readback.unmap();
        }
        Ok(out)
    }
}

/// Maps `slice` for reading and waits until it is, which is after
/// `submission` has finished.
fn map_and_wait(
    context: &GpuContext,
    wait: Wait,
    slice: wgpu::BufferSlice<'_>,
    submission: wgpu::SubmissionIndex,
    timeout: Duration,
) -> Result<(), GpuError> {
    let (sender, receiver) = mpsc::channel();
    slice.map_async(wgpu::MapMode::Read, move |result| {
        // The receiver outlives the wait below; if it is gone, nobody is
        // waiting for this result.
        if sender.send(result).is_err() {
            tracing::debug!("a compair GPU readback finished after its handle was dropped");
        }
    });
    match wait {
        Wait::Block => {
            context
                .device
                .poll(wgpu::PollType::Wait {
                    submission_index: Some(submission),
                    timeout: Some(timeout),
                })
                .map_err(|error| match error {
                    wgpu::PollError::Timeout => GpuError::Timeout { step: "scoring", timeout },
                    other => GpuError::Poll(other),
                })?;
            receiver.try_recv().map_err(|_| GpuError::MapPending)??;
        }
        Wait::Spin => {
            let started = Instant::now();
            loop {
                context.device.poll(wgpu::PollType::Poll).map_err(GpuError::Poll)?;
                if let Ok(result) = receiver.try_recv() {
                    result?;
                    break;
                }
                if started.elapsed() > timeout {
                    return Err(GpuError::Timeout { step: "scoring", timeout });
                }
                std::thread::yield_now();
            }
        }
    }
    Ok(())
}

/// The batch kernel's last step, on the CPU in `f64` exactly as there: the
/// shader returns the last row's `f32` sum and the power of two it is scaled
/// by, and a `log10` in `f64` is not something WGSL has.
pub(super) fn finish(sum: f32, exponent: i32) -> Log10Likelihood {
    if sum > 0.0 {
        Log10Likelihood::new(
            f64::from(sum).log10() - f64::from(exponent) * core::f64::consts::LOG10_2,
        )
    } else {
        Log10Likelihood::IMPOSSIBLE
    }
}
