//! Chemical formula parsing.
//!
//! Supported syntax: element symbols, integer subscripts, nested parentheses
//! and case-insensitive symbols (`na2so4`). Charges and hydrate dots are not
//! part of a formula here -- charge lives on [`crate::Species`].

use core::fmt;

use crate::element;

/// A parsed formula: element counts keyed by atomic number, sorted by Z.
#[derive(Clone, Debug, Default, PartialEq)]
pub struct Formula {
    /// `(Z, count)` pairs, sorted by Z, counts > 0.
    pub counts: Vec<(u8, f64)>,
}

/// Parse failure.
#[derive(Clone, Debug, PartialEq, Eq)]
pub enum FormulaError {
    /// Unknown element symbol.
    UnknownSymbol(String),
    /// Unbalanced or misplaced parenthesis.
    UnbalancedParentheses,
    /// Unexpected character.
    Unexpected(char),
    /// Empty input.
    Empty,
}

impl fmt::Display for FormulaError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            FormulaError::UnknownSymbol(s) => write!(f, "unknown element symbol: {s}"),
            FormulaError::UnbalancedParentheses => write!(f, "unbalanced parentheses"),
            FormulaError::Unexpected(c) => write!(f, "unexpected character: {c:?}"),
            FormulaError::Empty => write!(f, "empty formula"),
        }
    }
}

impl std::error::Error for FormulaError {}

impl Formula {
    /// Parse `"Ca(OH)2"` style formulas.
    pub fn parse(input: &str) -> Result<Self, FormulaError> {
        let raw: Vec<char> = input.chars().filter(|c| !c.is_whitespace()).collect();
        let case_sensitive = raw.iter().any(|c| c.is_ascii_uppercase());
        let s: Vec<char> = normalize_case(raw);
        if s.is_empty() {
            return Err(FormulaError::Empty);
        }
        let mut stack: Vec<Vec<(u8, f64)>> = vec![Vec::new()];
        let mut i = 0;
        while i < s.len() {
            let c = s[i];
            if c == '(' {
                stack.push(Vec::new());
                i += 1;
            } else if c == ')' {
                if stack.len() == 1 {
                    return Err(FormulaError::UnbalancedParentheses);
                }
                i += 1;
                let (mult, next) = read_multiplier(&s, i);
                i = next;
                let inner = stack.pop().expect("checked len");
                let top = stack.last_mut().expect("non-empty");
                for (z, n) in inner {
                    top.push((z, n * mult));
                }
            } else if c.is_ascii_alphabetic() {
                let bad: String = s[i..]
                    .iter()
                    .take_while(|c| c.is_ascii_alphabetic())
                    .take(3)
                    .collect();
                let (element, len) =
                    match_symbol(&s, i, case_sensitive).ok_or(FormulaError::UnknownSymbol(bad))?;
                i += len;
                let (mult, next) = read_multiplier(&s, i);
                i = next;
                stack.last_mut().expect("non-empty").push((element.z, mult));
            } else {
                return Err(FormulaError::Unexpected(c));
            }
        }
        if stack.len() != 1 {
            return Err(FormulaError::UnbalancedParentheses);
        }
        let mut counts: Vec<(u8, f64)> = Vec::new();
        for (z, n) in stack.pop().expect("one level") {
            match counts.iter_mut().find(|(zz, _)| *zz == z) {
                Some(slot) => slot.1 += n,
                None => counts.push((z, n)),
            }
        }
        counts.retain(|(_, n)| *n > 0.0);
        counts.sort_by_key(|(z, _)| *z);
        Ok(Formula { counts })
    }

    /// Count of element `z`.
    pub fn count(&self, z: u8) -> f64 {
        self.counts
            .iter()
            .find(|(zz, _)| *zz == z)
            .map_or(0.0, |(_, n)| *n)
    }

    /// Molar mass in g/mol from the IUPAC standard atomic weights.
    pub fn molar_mass(&self) -> f64 {
        self.counts
            .iter()
            .map(|(z, n)| element::by_z(*z).map_or(0.0, |e| e.mass) * n)
            .sum()
    }
}

/// Longest-prefix element match.
///
/// Case-sensitive input uses exact symbol spelling (`CO` = C + O, `Co` =
/// cobalt). All-lowercase input falls back to the capitalized form, which is
/// inherently ambiguous for `co` (cobalt wins); callers that care should spell
/// formulas with proper case.
fn match_symbol(
    s: &[char],
    start: usize,
    case_sensitive: bool,
) -> Option<(&'static crate::Element, usize)> {
    let max = (s.len() - start).min(3);
    for len in (1..=max).rev() {
        if !s[start..start + len].iter().all(char::is_ascii_alphabetic) {
            continue;
        }
        let sym: String = s[start..start + len].iter().collect();
        if case_sensitive {
            if let Some(e) = element::all().iter().find(|e| e.symbol == sym) {
                return Some((e, len));
            }
        } else if let Some(e) = element::by_symbol(&sym)
            && sym.eq_ignore_ascii_case(e.symbol)
        {
            return Some((e, len));
        }
    }
    None
}

/// All-lowercase input is capitalized per element token; mixed case is kept.
fn normalize_case(chars: Vec<char>) -> Vec<char> {
    if chars.iter().any(|c| c.is_ascii_uppercase()) {
        return chars;
    }
    let mut out = Vec::with_capacity(chars.len());
    let mut at_start = true;
    for c in chars {
        if c.is_ascii_alphabetic() {
            out.push(if at_start { c.to_ascii_uppercase() } else { c });
            at_start = false;
        } else {
            out.push(c);
            at_start = true;
        }
    }
    out
}

fn read_multiplier(s: &[char], start: usize) -> (f64, usize) {
    let mut i = start;
    let mut seen = false;
    let mut value = 0.0_f64;
    while i < s.len() && s[i].is_ascii_digit() {
        value = value * 10.0 + f64::from(s[i].to_digit(10).expect("digit"));
        seen = true;
        i += 1;
    }
    if seen { (value, i) } else { (1.0, i) }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn roundtrip(f: &str) -> Vec<(u8, f64)> {
        Formula::parse(f).unwrap().counts
    }

    #[test]
    fn simple_and_nested() {
        assert_eq!(roundtrip("H2O"), vec![(1, 2.0), (8, 1.0)]);
        assert_eq!(roundtrip("Ca(OH)2"), vec![(1, 2.0), (8, 2.0), (20, 1.0)]);
        assert_eq!(
            roundtrip("Al2(SO4)3"),
            vec![(8, 12.0), (13, 2.0), (16, 3.0)]
        );
        assert_eq!(roundtrip("NaCl"), vec![(11, 1.0), (17, 1.0)]);
    }

    #[test]
    fn mixed_case_and_whitespace() {
        assert_eq!(roundtrip("h2 o"), roundtrip("H2O"));
        assert_eq!(roundtrip("na2so4"), roundtrip("Na2SO4"));
    }

    #[test]
    fn molar_mass_matches_known_values() {
        let water = Formula::parse("H2O").unwrap().molar_mass();
        assert!((water - 18.015).abs() < 0.02, "{water}");
        let salt = Formula::parse("NaCl").unwrap().molar_mass();
        assert!((salt - 58.44).abs() < 0.05, "{salt}");
    }

    #[test]
    fn rejects_bad_input() {
        assert_eq!(Formula::parse(""), Err(FormulaError::Empty));
        assert_eq!(
            Formula::parse("(H2O"),
            Err(FormulaError::UnbalancedParentheses)
        );
        assert_eq!(
            Formula::parse("Xx2"),
            Err(FormulaError::UnknownSymbol("Xx".into()))
        );
    }
}
