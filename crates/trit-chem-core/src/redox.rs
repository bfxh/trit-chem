//! Standard reduction potentials.
//!
//! A potential belongs to a *half-reaction*, not to an element. The design
//! document's element table has a "standard electrode potential" column; that
//! column cannot be filled per element without picking an arbitrary couple, so
//! the data lives here instead.
//!
//! Values: standard reduction potentials at 25 C (aqueous, 1 M / 1 bar),
//! transcribed from the compilation in data/external/thermo_table.txt (Table 3).

use crate::thermo::R;

/// Faraday constant, C/mol.
pub const F: f64 = 96_485.332_123_3;

/// One half-reaction written as a reduction.
#[derive(Clone, Copy, Debug)]
pub struct HalfReaction {
    /// Reaction as written, e.g. `"2H+ + 2e- -> H2"`.
    pub label: &'static str,
    /// Electrons transferred.
    pub n: u8,
    /// Standard reduction potential, V.
    pub e0: f64,
}

impl HalfReaction {
    /// Standard Gibbs energy of the reduction, kJ/mol.
    pub fn dg_kj(&self) -> f64 {
        -f64::from(self.n) * F * self.e0 / 1000.0
    }

    /// Nernst potential at temperature `t` and reaction quotient `q`.
    pub fn nernst(&self, t: f64, q: f64) -> f64 {
        self.e0 - R * t / (f64::from(self.n) * F) * q.max(1e-300).ln()
    }
}

/// Cell potential for `cathode - anode`.
pub fn cell_potential(cathode: &HalfReaction, anode: &HalfReaction) -> f64 {
    cathode.e0 - anode.e0
}

/// True when the cell as written (cathode reduced, anode oxidised) is
/// spontaneous under standard conditions.
pub fn is_spontaneous(cathode: &HalfReaction, anode: &HalfReaction) -> bool {
    cell_potential(cathode, anode) > 0.0
}

/// A verified subset of the electrochemical series, strongest oxidiser first.
pub const HALF_REACTIONS: &[HalfReaction] = &[
    HalfReaction {
        label: "F2 + 2e- -> 2F-",
        n: 2,
        e0: 2.87,
    },
    HalfReaction {
        label: "O3 + 2H+ + 2e- -> O2 + H2O",
        n: 2,
        e0: 2.07,
    },
    HalfReaction {
        label: "H2O2 + 2H+ + 2e- -> 2H2O",
        n: 2,
        e0: 1.78,
    },
    HalfReaction {
        label: "MnO4- + 8H+ + 5e- -> Mn2+ + 4H2O",
        n: 5,
        e0: 1.51,
    },
    HalfReaction {
        label: "Au3+ + 3e- -> Au",
        n: 3,
        e0: 1.40,
    },
    HalfReaction {
        label: "Cl2 + 2e- -> 2Cl-",
        n: 2,
        e0: 1.36,
    },
    HalfReaction {
        label: "Cr2O7 2- + 14H+ + 6e- -> 2Cr3+ + 7H2O",
        n: 6,
        e0: 1.33,
    },
    HalfReaction {
        label: "O2 + 4H+ + 4e- -> 2H2O",
        n: 4,
        e0: 1.23,
    },
    HalfReaction {
        label: "MnO2 + 4H+ + 2e- -> Mn2+ + 2H2O",
        n: 2,
        e0: 1.23,
    },
    HalfReaction {
        label: "Br2 + 2e- -> 2Br-",
        n: 2,
        e0: 1.09,
    },
    HalfReaction {
        label: "NO3- + 4H+ + 3e- -> NO + 2H2O",
        n: 3,
        e0: 0.96,
    },
    HalfReaction {
        label: "Ag+ + e- -> Ag",
        n: 1,
        e0: 0.80,
    },
    HalfReaction {
        label: "Fe3+ + e- -> Fe2+",
        n: 1,
        e0: 0.77,
    },
    HalfReaction {
        label: "O2 + 2H+ + 2e- -> H2O2",
        n: 2,
        e0: 0.68,
    },
    HalfReaction {
        label: "I2 + 2e- -> 2I-",
        n: 2,
        e0: 0.54,
    },
    HalfReaction {
        label: "Cu+ + e- -> Cu",
        n: 1,
        e0: 0.52,
    },
    HalfReaction {
        label: "O2 + 2H2O + 4e- -> 4OH-",
        n: 4,
        e0: 0.40,
    },
    HalfReaction {
        label: "Cu2+ + 2e- -> Cu",
        n: 2,
        e0: 0.34,
    },
    HalfReaction {
        label: "AgCl + e- -> Ag + Cl-",
        n: 1,
        e0: 0.22,
    },
    HalfReaction {
        label: "SO4 2- + 4H+ + 2e- -> H2SO3 + H2O",
        n: 2,
        e0: 0.17,
    },
    HalfReaction {
        label: "Cu2+ + e- -> Cu+",
        n: 1,
        e0: 0.15,
    },
    HalfReaction {
        label: "Sn4+ + 2e- -> Sn2+",
        n: 2,
        e0: 0.15,
    },
    HalfReaction {
        label: "2H+ + 2e- -> H2",
        n: 2,
        e0: 0.00,
    },
    HalfReaction {
        label: "Pb2+ + 2e- -> Pb",
        n: 2,
        e0: -0.13,
    },
    HalfReaction {
        label: "Ni2+ + 2e- -> Ni",
        n: 2,
        e0: -0.23,
    },
    HalfReaction {
        label: "Co2+ + 2e- -> Co",
        n: 2,
        e0: -0.28,
    },
    HalfReaction {
        label: "Cd2+ + 2e- -> Cd",
        n: 2,
        e0: -0.40,
    },
    HalfReaction {
        label: "Fe2+ + 2e- -> Fe",
        n: 2,
        e0: -0.44,
    },
    HalfReaction {
        label: "Cr3+ + 3e- -> Cr",
        n: 3,
        e0: -0.74,
    },
    HalfReaction {
        label: "Zn2+ + 2e- -> Zn",
        n: 2,
        e0: -0.76,
    },
    HalfReaction {
        label: "2H2O + 2e- -> H2 + 2OH-",
        n: 2,
        e0: -0.83,
    },
    HalfReaction {
        label: "Mn2+ + 2e- -> Mn",
        n: 2,
        e0: -1.18,
    },
    HalfReaction {
        label: "Al3+ + 3e- -> Al",
        n: 3,
        e0: -1.66,
    },
    HalfReaction {
        label: "Mg2+ + 2e- -> Mg",
        n: 2,
        e0: -2.36,
    },
    HalfReaction {
        label: "Na+ + e- -> Na",
        n: 1,
        e0: -2.71,
    },
    HalfReaction {
        label: "Ca2+ + 2e- -> Ca",
        n: 2,
        e0: -2.87,
    },
    HalfReaction {
        label: "K+ + e- -> K",
        n: 1,
        e0: -2.93,
    },
    HalfReaction {
        label: "Li+ + e- -> Li",
        n: 1,
        e0: -3.05,
    },
];

/// Look up a half-reaction by exact label.
pub fn by_label(label: &str) -> Option<&'static HalfReaction> {
    HALF_REACTIONS.iter().find(|h| h.label == label)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn hydrogen_electrode_is_the_zero() {
        let h = by_label("2H+ + 2e- -> H2").unwrap();
        assert_eq!(h.e0, 0.0);
        assert!(h.dg_kj().abs() < 1e-9);
    }

    #[test]
    fn sodium_metal_reduces_water() {
        let water = by_label("2H2O + 2e- -> H2 + 2OH-").unwrap();
        let sodium = by_label("Na+ + e- -> Na").unwrap();
        // Sodium is the anode (oxidised), water is reduced at the cathode.
        assert!(is_spontaneous(water, sodium));
        let e = cell_potential(water, sodium);
        assert!((e - 1.88).abs() < 1e-9, "{e}");
        assert!((water.dg_kj() - 160.2).abs() < 0.5, "{}", water.dg_kj());
    }

    #[test]
    fn nernst_reduces_to_standard_potential_at_q_one() {
        let cu = by_label("Cu2+ + 2e- -> Cu").unwrap();
        assert!((cu.nernst(298.15, 1.0) - cu.e0).abs() < 1e-12);
        // Ten-fold lower [Cu2+] (Q = 10) lowers the potential by ~29.6 mV.
        let e = cu.nernst(298.15, 10.0);
        assert!((cu.e0 - e - 0.0296).abs() < 2e-3, "{}", cu.e0 - e);
    }
}
