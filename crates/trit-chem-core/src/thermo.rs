//! Standard-state thermochemistry.
//!
//! Standard state: ideal gas at 1 bar, solutes at 1 mol/L, pure condensed
//! phases at unit activity. The Shomate equations follow NIST WebBook:
//!
//! ```text
//! Cp = A + B*t + C*t^2 + D*t^3 + E/t^2
//! H  - H(298.15) = A*t + B*t^2/2 + C*t^3/3 + D*t^4/4 - E/t + F - H
//! S  = A*ln(t) + B*t + C*t^2/2 + D*t^3/3 - E/(2*t^2) + G
//! t = T / 1000
//! ```
//!
//! Note the two same-named coefficients: `F` and `H` are entries 5 and 7 of
//! the eight-number block.

/// Molar gas constant, J/(mol K).
pub const R: f64 = 8.314_462_618_153_24;
/// Thermochemical reference temperature, K.
pub const T0: f64 = 298.15;
/// Standard pressure, Pa.
pub const P0_PA: f64 = 1.0e5;

/// One fitted Shomate range `[lo, hi]` in kelvin.
#[derive(Clone, Copy, Debug)]
pub struct ShomateBlock {
    /// Lower temperature bound, K.
    pub lo: f64,
    /// Upper temperature bound, K.
    pub hi: f64,
    /// Coefficients `[A, B, C, D, E, F, G, H]`.
    pub a: [f64; 8],
}

impl ShomateBlock {
    fn t(&self, t: f64) -> f64 {
        t / 1000.0
    }

    /// Heat capacity at `t`, J/(mol K).
    pub fn cp(&self, t: f64) -> f64 {
        let x = self.t(t);
        self.a[0] + self.a[1] * x + self.a[2] * x * x + self.a[3] * x * x * x + self.a[4] / (x * x)
    }

    /// `H(T) - H(298.15)` in kJ/mol.
    fn h_rel(&self, t: f64) -> f64 {
        let x = self.t(t);
        self.a[0] * x
            + self.a[1] * x * x / 2.0
            + self.a[2] * x * x * x / 3.0
            + self.a[3] * x * x * x * x / 4.0
            - self.a[4] / x
            + self.a[5]
            - self.a[7]
    }

    /// Standard entropy in J/(mol K).
    fn s_abs(&self, t: f64) -> f64 {
        let x = self.t(t);
        self.a[0] * x.ln() + self.a[1] * x + self.a[2] * x * x / 2.0 + self.a[3] * x * x * x / 3.0
            - self.a[4] / (2.0 * x * x)
            + self.a[6]
    }
}

/// Temperature model for one species.
#[derive(Clone, Copy, Debug)]
pub enum Cp {
    /// Temperature-independent heat capacity, J/(mol K).
    Const(f64),
    /// Piecewise Shomate fits, ordered by `lo`.
    Shomate(&'static [ShomateBlock]),
}

impl Cp {
    /// Pick the fit covering `t`, or the nearest one when `t` is outside every
    /// fitted range (M0 extrapolates instead of refusing; DESIGN.md says so).
    pub fn block_for(&self, t: f64) -> Option<&'static ShomateBlock> {
        match self {
            Cp::Const(_) => None,
            Cp::Shomate(blocks) => {
                let dist = |b: &ShomateBlock| (b.lo - t).max(t - b.hi).max(0.0);
                blocks
                    .iter()
                    .find(|b| t >= b.lo && t <= b.hi)
                    .or_else(|| blocks.iter().min_by(|x, y| dist(x).total_cmp(&dist(y))))
            }
        }
    }

    /// Heat capacity at `t`, J/(mol K).
    pub fn cp(&self, t: f64) -> f64 {
        match self {
            Cp::Const(c) => *c,
            Cp::Shomate(_) => self.block_for(t).map_or(0.0, |b| b.cp(t)),
        }
    }

    /// `H(T) - H(298.15)` in kJ/mol.
    pub fn h_rel_kj(&self, t: f64) -> f64 {
        match self {
            Cp::Const(c) => c * (t - T0) / 1000.0,
            Cp::Shomate(_) => match self.block_for(t) {
                Some(b) => b.h_rel(t) - b.h_rel(T0),
                None => 0.0,
            },
        }
    }

    /// `S(T) - S(298.15)` in J/(mol K).
    pub fn s_rel(&self, t: f64) -> f64 {
        match self {
            Cp::Const(c) => c * (t / T0).ln(),
            Cp::Shomate(_) => match self.block_for(t) {
                Some(b) => b.s_abs(t) - b.s_abs(T0),
                None => 0.0,
            },
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    // NIST WebBook, H2 gas, first fit.
    const H2: [f64; 8] = [
        33.066_178,
        -11.363_417,
        11.432_816,
        -2.772_874,
        -0.158_558,
        -9.980_797,
        172.707_974,
        0.0,
    ];

    #[test]
    fn shomate_reproduces_reference_entropy() {
        let b = ShomateBlock {
            lo: 298.0,
            hi: 1000.0,
            a: H2,
        };
        assert!(
            (b.s_abs(298.15) - 130.68).abs() < 0.05,
            "{}",
            b.s_abs(298.15)
        );
        assert!(b.h_rel(298.15).abs() < 0.01, "{}", b.h_rel(298.15));
        assert!((b.cp(298.15) - 28.836).abs() < 0.05, "{}", b.cp(298.15));
    }

    #[test]
    fn constant_cp_integrates_correctly() {
        let cp = Cp::Const(75.291);
        assert!((cp.cp(400.0) - 75.291).abs() < 1e-12);
        let dh = cp.h_rel_kj(373.15);
        assert!((dh - 75.291 * 75.0 / 1000.0).abs() < 1e-9, "{dh}");
        let ds = cp.s_rel(373.15);
        assert!(
            (ds - 75.291 * (373.15_f64 / 298.15).ln()).abs() < 1e-9,
            "{ds}"
        );
    }
}
