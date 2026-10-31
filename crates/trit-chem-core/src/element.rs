//! The 118-element table and lookups.

use crate::elements_data::ELEMENTS;

/// Aggregation state used by both elements and species.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash)]
pub enum Phase {
    /// Ideal gas.
    Gas,
    /// Pure liquid.
    Liquid,
    /// Pure solid.
    Solid,
    /// Solute at unit molarity standard state.
    Aqueous,
    /// Unclassified (only used by element tables for exotic entries).
    Unknown,
}

/// One row of the periodic table.
///
/// `first_ie_ev` is the first ionization energy in eV (NIST ASD values via the
/// source JSON, converted from kJ/mol). `outer_shell_electrons` is the
/// population of the outermost shell, not the full valence count of d-block
/// elements -- those two differ and the design doc conflated them.
#[derive(Clone, Copy, Debug)]
pub struct Element {
    /// Atomic number.
    pub z: u8,
    /// IUPAC symbol.
    pub symbol: &'static str,
    /// English name.
    pub name: &'static str,
    /// Standard atomic weight, g/mol.
    pub mass: f64,
    /// Pauling electronegativity (dimensionless), where defined.
    pub electronegativity: Option<f64>,
    /// First ionization energy, eV.
    pub first_ie_ev: Option<f64>,
    /// Empirical atomic radius, pm.
    pub atomic_radius_pm: Option<f64>,
    /// Electrons in the outermost shell.
    pub outer_shell_electrons: u8,
    /// Common oxidation states (CRC Handbook; predicted for Z >= 104).
    pub oxidation_states: &'static [i8],
    /// Phase at 298.15 K and 1 bar.
    pub phase_298: Phase,
    /// Melting point, K.
    pub melt_k: Option<f64>,
    /// Boiling point, K.
    pub boil_k: Option<f64>,
    /// Density at 298 K, g/cm3.
    pub density_g_cm3: Option<f64>,
    /// Periodic-table category string from the source table.
    pub category: &'static str,
}

impl Element {
    /// True for alkali, alkaline-earth, transition and post-transition metals.
    pub fn is_metal(&self) -> bool {
        self.category.contains("metal")
            && !self.category.contains("nonmetal")
            && !self.category.contains("metalloid")
    }

    /// True when `ox` appears in the curated common-state list.
    pub fn has_oxidation_state(&self, ox: i8) -> bool {
        self.oxidation_states.contains(&ox)
    }

    /// Smallest common positive oxidation state, if any.
    pub fn lowest_cation(&self) -> Option<i8> {
        self.oxidation_states
            .iter()
            .copied()
            .filter(|o| *o > 0)
            .min()
    }

    /// Largest common positive oxidation state, if any.
    pub fn highest_cation(&self) -> Option<i8> {
        self.oxidation_states
            .iter()
            .copied()
            .filter(|o| *o > 0)
            .max()
    }
}

/// All 118 elements, ordered by atomic number.
pub fn all() -> &'static [Element; 118] {
    &ELEMENTS
}

/// Look up by atomic number.
pub fn by_z(z: u8) -> Option<&'static Element> {
    ELEMENTS.iter().find(|e| e.z == z)
}

/// Look up by symbol; accepts exact (`Na`) and normalized (`na`, `NA`) forms.
pub fn by_symbol(symbol: &str) -> Option<&'static Element> {
    if let Some(hit) = ELEMENTS.iter().find(|e| e.symbol == symbol) {
        return Some(hit);
    }
    let mut chars = symbol.chars();
    let first = chars.next()?.to_ascii_uppercase();
    let rest: String = chars.map(|c| c.to_ascii_lowercase()).collect();
    let normalized = format!("{first}{rest}");
    ELEMENTS.iter().find(|e| e.symbol == normalized)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn table_is_complete_and_sorted() {
        assert_eq!(all().len(), 118);
        for (i, e) in all().iter().enumerate() {
            assert_eq!(usize::from(e.z), i + 1);
        }
    }

    #[test]
    fn spot_checks_against_handbook_values() {
        let na = by_symbol("Na").unwrap();
        assert!((na.mass - 22.989_769_282).abs() < 1e-9);
        assert_eq!(na.electronegativity, Some(0.93));
        assert_eq!(na.outer_shell_electrons, 1);
        assert!(na.is_metal());

        let cl = by_symbol("Cl").unwrap();
        assert!((cl.first_ie_ev.unwrap() - 12.9676).abs() < 5e-3);
        assert!(!cl.is_metal());

        assert_eq!(by_symbol("na").unwrap().z, 11);
        assert_eq!(by_z(26).unwrap().symbol, "Fe");
        assert!(by_z(119).is_none());
    }
}
