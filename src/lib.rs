//! Lossy archives for Oxford Nanopore and PacBio long reads that pipelines read like the original.
//!
//! Records are edited as encoded BAM and written with pasta as CRAM (3.0 by default, 3.1 with
//! [`policy::Mode::Max`]) or BAM:
//!
//! - quality scores are binned ([`quality::Bins`]), optionally keeping native scores in
//!   reference homopolymers ([`homopolymer::Runs`]), or replaced by one score per read;
//! - base-modification probabilities (`ML`) are quantized to a few levels;
//! - `MD`/`NM` that CRAM cannot regenerate are left out, and with `Max`, secondary alignments.
//!
//! [`policy::Policy::preset`] gives the validated policy for each [`policy::Platform`] and
//! [`policy::Mode`]; a TOML [`policy::Policy`] can replace it.

mod alignments;
pub mod fastq;
pub mod homopolymer;
pub mod methylation;
pub mod pack;
mod pipeline;
pub mod policy;
pub mod quality;
pub mod raw;
pub mod transform;
pub mod verify;
