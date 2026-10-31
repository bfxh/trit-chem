//! Deterministic, dependency-free thermochemistry and reaction engine.
//!
//! Scope of M0:
//!   * 118 elements with the properties the AE-ENGINE F0 design asks for;
//!   * standard-state thermochemistry at 1 bar / 1 M, with NIST Shomate
//!     temperature dependence where available and constant Cp elsewhere;
//!   * reaction construction with element and charge balance checks;
//!   * dG ranking, mass-action kinetics (Arrhenius + RK4) and a
//!     thermodynamically downhill reaction cascade.
//!
//! Deliberate omissions (see docs/DESIGN.md): activity coefficients, phase
//! equilibrium beyond pure condensed phases, and stiff ODE integration.
#![forbid(unsafe_code)]
#![warn(missing_docs)]

pub mod element;
pub mod equilibrium;
pub mod formula;
pub mod kinetics;
pub mod reaction;
pub mod redox;
pub mod rules;
pub mod species;
pub mod state;
pub mod thermo;

mod elements_data;
mod species_data;

pub use element::{Element, Phase};
pub use formula::{Formula, FormulaError};
pub use reaction::{Reaction, ReactionError, Term};
pub use species::{Db, SPECIES, Species};
pub use state::State;
pub use thermo::{Cp, R, ShomateBlock, T0};

/// Crate version, mirrored into CLI diagnostics.
pub const VERSION: &str = env!("CARGO_PKG_VERSION");
