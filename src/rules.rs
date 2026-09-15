use crate::elements::get_element;
use crate::types::FormulaCandidate;

/// Validate a candidate formula against the Seven Golden Rules.
/// Returns true if the candidate passes all enabled checks.
pub fn validate(candidate: &FormulaCandidate, strict: bool) -> bool {
    rdbe_check(candidate) && element_ratio_check(candidate, strict) && nitrogen_rule_check(candidate)
}

/// Rule 1: RDBE must be >= -0.5 and a multiple of 0.5
fn rdbe_check(candidate: &FormulaCandidate) -> bool {
    let dbe = candidate.dbe;
    if dbe < -0.5 {
        return false;
    }
    // DBE should be a multiple of 0.5 (integers or half-integers)
    let twice_dbe = dbe * 2.0;
    (twice_dbe - twice_dbe.round()).abs() < 1e-6
}

/// Rule 2: Element ratio constraints (Kind & Fiehn 2007, Table 2)
fn element_ratio_check(candidate: &FormulaCandidate, strict: bool) -> bool {
    let get = |sym: &str| -> u32 {
        candidate
            .composition
            .iter()
            .find(|(s, _)| *s == sym)
            .map(|&(_, n)| n)
            .unwrap_or(0)
    };

    let c = get("C");
    if c == 0 {
        return true; // Can't check ratios without carbon
    }
    let cf = c as f64;

    let checks: &[(&str, f64, f64, f64, f64)] = &[
        //  element, loose_min, loose_max, strict_min, strict_max (extended / common range)
        ("H", 0.1, 6.0, 0.2, 3.1),
        ("F", 0.0, 6.0, 0.0, 1.5),
        ("Cl", 0.0, 2.0, 0.0, 0.8),
        ("Br", 0.0, 2.0, 0.0, 0.8),
        ("N", 0.0, 4.0, 0.0, 1.3),
        ("O", 0.0, 3.0, 0.0, 1.2),
        ("P", 0.0, 2.0, 0.0, 0.3),
        ("S", 0.0, 3.0, 0.0, 0.8),
        ("Si", 0.0, 1.0, 0.0, 0.5),
    ];

    for &(sym, loose_min, loose_max, strict_min, strict_max) in checks {
        let n = get(sym);
        if n == 0 {
            continue;
        }
        let ratio = n as f64 / cf;
        let (lo, hi) = if strict {
            (strict_min, strict_max)
        } else {
            (loose_min, loose_max)
        };
        if ratio < lo || ratio > hi {
            return false;
        }
    }
    true
}

/// Rule 4: Nitrogen rule (odd/even parity check)
/// For molecules with even nominal mass and even N count: H should be even
/// For molecules with odd N count: molecular mass parity flips
fn nitrogen_rule_check(candidate: &FormulaCandidate) -> bool {
    let get = |sym: &str| -> u32 {
        candidate
            .composition
            .iter()
            .find(|(s, _)| *s == sym)
            .map(|&(_, n)| n)
            .unwrap_or(0)
    };

    let n_count = get("N");
    // Nominal mass from mass numbers: rounding the exact mass flips parity once the
    // accumulated mass defect exceeds 0.5 Da (e.g. C42H82NO8P, 759.578 Da).
    let nominal_mass: u32 = candidate
        .composition
        .iter()
        .map(|&(s, n)| n * get_element(s).map_or(0, |e| e.isotopes[0].mass_number))
        .sum();

    let mass_is_odd = nominal_mass % 2 == 1;
    let n_is_odd = n_count % 2 == 1;

    // For standard organic molecules (C, H, N, O, P, S):
    // Nitrogen rule: M is odd iff N count is odd
    mass_is_odd == n_is_odd
}
