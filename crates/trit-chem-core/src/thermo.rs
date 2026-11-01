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

    /// Absolute standard entropy in J/(mol K); used only by the unit test that
    /// checks the fits against tabulated values.
    #[cfg(test)]
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
    ///
    /// Rows whose fits do not reach 298.15 K carry a generated bridge segment,
    /// so the extrapolation that matters in practice is the one above `hi`.
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
            // Each NIST block is already referenced to H(298.15), so use it
            // directly. Subtracting a value from the *same* block at 298.15
            // would extrapolate that block outside its range and put a step in
            // H(T) at every range boundary (measured: 3.1 kJ/mol for CO2 at
            // 1200 K before this fix).
            Cp::Shomate(_) => match self.block_for(t) {
                Some(b) => b.h_rel(t),
                None => 0.0,
            },
        }
    }

    /// `S(T) - S(298.15)` in J/(mol K).
    ///
    /// Shomate blocks give absolute entropy, but the reference value for the
    /// *first* block only equals the tabulated `S(298.15)`; for a higher block
    /// `s_abs(298.15)` is an extrapolation and would put a step in S(T) at the
    /// range boundary. Integrate `Cp/T` from 298.15 instead: continuous across
    /// ranges and anchored at the tabulated entropy.
    ///
    /// A version that split the integral at every range edge was measured and
    /// rejected: against a 4000-slice reference the single pass here is out by
    /// +0.05 J/(mol K) for molten NaCl at 1738 K and the split version by
    /// -0.11, because each segment's endpoint still samples the block that
    /// claims the boundary. The kink is not what costs accuracy.
    pub fn s_rel(&self, t: f64) -> f64 {
        match self {
            Cp::Const(c) => c * (t / T0).ln(),
            Cp::Shomate(_) => integrate(|x| self.cp(x) / x, T0, t, 48),
        }
    }
}

/// Composite Simpson integration with a fixed even number of intervals.
///
/// Used for `S(T) - S(298.15)`; 48 intervals keep the error far below the 0.01
/// J/(mol K) that matters, at ~50 Cp evaluations per call. If `s()` ever shows
/// up in a hot loop, cache it per (species, temperature).
fn integrate<F: Fn(f64) -> f64>(f: F, a: f64, b: f64, n: usize) -> f64 {
    let n = n.max(2) & !1;
    let h = (b - a) / n as f64;
    let mut sum = f(a) + f(b);
    for k in 1..n {
        sum += f(a + h * k as f64) * if k % 2 == 1 { 4.0 } else { 2.0 };
    }
    sum * h / 3.0
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

    /// NIST fits some phases only from their melting or boiling point upward
    /// (Fe gas: 3133.345-6000 K, NaCl liquid: 1074-2500 K). Run the polynomial
    /// down to 298.15 K and its `E/t^2` term takes over - measured 79 J/(mol K)
    /// for iron vapour against the real 25.15, and 428 J/(mol K) for molten salt
    /// against 70.4 at the fit's own lower edge. The generator bridges those rows.
    #[test]
    fn no_row_extrapolates_a_runaway_heat_capacity_down_to_298_k() {
        let db = crate::species::Db::standard();
        let mut bridged = 0;
        for sp in db.species() {
            let Cp::Shomate(blocks) = sp.cp else { continue };
            let first = &blocks[0];
            if first.lo > T0 {
                let edge = first.cp(first.lo);
                let at298 = sp.cp.cp(T0);
                assert!(
                    edge > 0.0 && (at298 / edge - 1.0).abs() <= 0.3,
                    "{}: Cp(298.15) = {at298} against {edge} at the lower edge {} K",
                    sp.id,
                    first.lo
                );
            } else if first.a[1] == 0.0 && first.a[4] == 0.0 && blocks.len() > 1 {
                bridged += 1;
                let seam = first.hi;
                // The bridge's Cp is an effective constant, so it does not have
                // to meet the fit's own value at the seam - but the gap must
                // stay a rounding of the same physics, not a new phase.
                let below = sp.cp.cp(seam - 1e-6);
                let above = sp.cp.cp(seam + 1e-6);
                assert!(
                    above > 0.0 && (below / above - 1.0).abs() < 0.2,
                    "{}: Cp jumps {below} -> {above} at the bridge seam {seam} K",
                    sp.id
                );
                assert!(
                    (sp.cp.h_rel_kj(seam - 1e-6) - sp.cp.h_rel_kj(seam + 1e-6)).abs() < 0.1,
                    "{}: H step at the bridge seam {seam} K",
                    sp.id
                );
            }
        }
        assert!(
            bridged >= 2,
            "expected at least two bridged rows, got {bridged}"
        );
    }
}
