//! Scoring on the GPU, with the CPU as the fallback and the check.
//!
//! A launch has a fixed cost of a millisecond or so, so the GPU pays for
//! itself at thousands of pairs per launch, not tens: collect many loci into
//! one `GpuPairs` rather than launching per locus. This scores 4,096 simulated
//! reads against a reference and an alternate haplotype in one launch.
//!
//! `cargo run -p compair --release --features gpu --example gpu`

use std::time::Instant;

use compair::{
    Band, Base, BaseQuality, Haplotype, Read, StandardEmission, Strand, Workspace,
    gpu::{GpuAligner, GpuContext, GpuError, GpuPairs},
};

const READS: usize = 4096;
const READ_LEN: usize = 100;

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let (haplotypes, reads) = simulate()?;
    let haplotypes: Vec<&Haplotype> = haplotypes.iter().collect();
    let reads: Vec<(&Read, Band)> = reads.iter().map(|(read, band)| (read, *band)).collect();
    let emission = StandardEmission::default();

    let mut workspace = Workspace::new();
    let mut on_cpu = Vec::new();
    let started = Instant::now();
    workspace.align_reads(&haplotypes, &reads, &emission, &mut on_cpu);
    println!("CPU, one thread: {} pairs in {:?}", on_cpu.len(), started.elapsed());

    let context = match GpuContext::new() {
        Ok(context) => context,
        // No GPU is not an error for a caller that has the CPU path.
        Err(GpuError::NoAdapter(error)) => {
            println!("no GPU adapter ({error}); the CPU scores above are the answer");
            return Ok(());
        }
        Err(error) => return Err(error.into()),
    };
    println!("GPU: {}", context.adapter_info().name);

    // The emission is folded in as reads and haplotypes are pushed, and a
    // haplotype's terms depend on the strand of the reads scored against it,
    // so push each haplotype once per strand and reuse the slots.
    let mut pairs = GpuPairs::new();
    let mut slots = Vec::with_capacity(2 * haplotypes.len());
    for strand in [Strand::OT, Strand::OB] {
        for haplotype in &haplotypes {
            slots.push((strand, pairs.push_haplotype(haplotype, strand, &emission)?));
        }
    }
    for &(read, band) in &reads {
        let read_slot = pairs.push_read(read, &emission)?;
        for &(strand, haplotype_slot) in &slots {
            if strand == read.strand() {
                pairs.push_pair(read_slot, haplotype_slot, band)?;
            }
        }
    }

    // The first launch compiles the shaders; time the second.
    let mut aligner = GpuAligner::new(context);
    aligner.align(&pairs)?;
    let started = Instant::now();
    // `submit` returns as soon as the work is queued, so a caller can fill
    // the next `GpuPairs` before `collect` waits for these scores.
    let on_gpu = aligner.submit(&pairs)?.collect()?;
    println!("GPU, one launch: {} pairs in {:?}", on_gpu.len(), started.elapsed());

    // Scores come back in push order, which is `align_reads`' order here.
    // They are not promised bit-identical to the CPU's -- GPU compilers may
    // reorder `f32` arithmetic -- so compare with a tolerance.
    let worst = on_cpu
        .iter()
        .zip(&on_gpu)
        .map(|(cpu, gpu)| (cpu.get() - gpu.get()).abs())
        .filter(|difference| difference.is_finite())
        .fold(0.0, f64::max);
    let identical = on_cpu.iter().zip(&on_gpu).filter(|(cpu, gpu)| cpu == gpu).count();
    println!(
        "{identical} of {} scores bit-identical, largest difference {worst:.1e} log10",
        on_gpu.len()
    );
    Ok(())
}

/// A locus's haplotypes, and its reads each with its band.
type Locus = (Vec<Haplotype>, Vec<(Read, Band)>);

/// Two 300 bp haplotypes that differ by a 3 bp deletion in the middle, and
/// reads cut from either, on both strands, with a sequencing error each.
fn simulate() -> Result<Locus, compair::Error> {
    // xorshift64, so every run scores the same reads.
    let mut state = 0x9e37_79b9_7f4a_7c15u64;
    let mut next = move |below: u64| {
        state ^= state << 13;
        state ^= state >> 7;
        state ^= state << 17;
        usize::try_from(state % below).unwrap_or(0)
    };
    let reference: Vec<Base> =
        (0..300).map(|_| Base::KNOWN.get(next(4)).copied().unwrap_or(Base::A)).collect();
    let mut deletion = reference.clone();
    deletion.drain(150..153);

    let mut reads = Vec::with_capacity(READS);
    for index in 0..READS {
        let source = if index % 2 == 0 { &reference } else { &deletion };
        let start = next(190);
        let mut bases =
            source.get(start..start + READ_LEN).map(<[Base]>::to_vec).unwrap_or_default();
        if let Some(base) = bases.get_mut(next(READ_LEN as u64)) {
            *base = base.inverse();
        }
        let strand = if next(2) == 0 { Strand::OT } else { Strand::OB };
        let read = Read::uniform(
            bases,
            &[BaseQuality::from_byte(30); READ_LEN],
            BaseQuality::from_byte(45),
            BaseQuality::from_byte(45),
            BaseQuality::from_byte(10),
            strand,
        )?;
        reads.push((read, Band::anchored(i32::try_from(start).unwrap_or(0))));
    }
    Ok((vec![Haplotype::new(reference), Haplotype::new(deletion)], reads))
}
