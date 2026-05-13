//! MIRAGE: Mutation-encoded Inference of RNA Activity via Generative Effects
//!
//! Rust port of the R MIRAGE package. Estimates site-level RNA mutation or
//! conversion signal from matched treatment and control read-count tables.

pub mod motif;
pub mod stats;
pub mod optim;
pub mod mle;
pub mod prior;
pub mod inference;
pub mod io;

pub use inference::{
    estimate_inference_with_empirical, estimate_inference_with_prior,
    EmpiricalParams, PriorParams, EmpiricalResult, PriorResult, BgMethod, BgTarget,
    Lambda1Mode, NmerMode, SiteRecord, HomoSite, HeterSite,
};
pub use prior::{compute_prior, MotifFreqExp, MotifFreqBg, MotifPriorRow};
