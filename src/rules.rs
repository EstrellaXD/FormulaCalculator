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

/// Rule 2: Element ratio constraints (Kind & Fiehn 2007, Table 1)
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
        //  element, loose_min, loose_max, strict_min, strict_max
        ("H", 0.1, 6.0, 0.2, 3.1),
        ("N", 0.0, 4.0, 0.0, 1.3),
        ("O", 0.0, 3.0, 0.0, 1.2),
        ("P", 0.0, 6.0, 0.0, 0.3),
        ("S", 0.0, 2.0, 0.0, 0.8),
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
    let nominal_mass = candidate.monoisotopic_mass.round() as u64;

    // If N is even, nominal mass should be even for typical organic molecules.
    // If N is odd, nominal mass should be odd.
    // This is a soft check — only flag clear violations.
    let mass_is_odd = nominal_mass % 2 == 1;
    let n_is_odd = n_count % 2 == 1;

    // For standard organic molecules (C, H, N, O, P, S):
    // Nitrogen rule: M is odd iff N count is odd
    mass_is_odd == n_is_odd
}

/// Partial RDBE check for pruning during search.
/// Given partial composition (some elements assigned, some not yet),
/// check if DBE can still land in a valid range.
///
/// `assigned` = [(symbol, valence, count)] for already-assigned elements
/// `remaining_h_max` = maximum H atoms still possible
///
/// Returns false if the branch is impossible.
pub fn partial_dbe_feasible(
    assigned: &[(&str, i32, u32)],
    remaining_h_max: u32,
) -> bool {
    // Compute partial DBE from assigned elements
    let mut dbe = 1.0_f64;
    for &(_, valence, count) in assigned {
        dbe += count as f64 * (valence as f64 - 2.0) / 2.0;
    }

    // With remaining H (valence=1), each H decreases DBE by 0.5.
    // Minimum possible final DBE = dbe - remaining_h_max * 0.5
    // Maximum possible final DBE = dbe (if no more H added, ignoring other elements)
    // For validity we need final DBE >= -0.5
    // Best case (lowest DBE): all remaining H used → dbe - remaining_h_max * 0.5 >= -0.5?
    // That's always potentially satisfiable if current dbe isn't impossibly negative.
    // Worst case: minimum DBE with max H must be achievable: dbe - max_h * 0.5 can be >= -0.5

    // Current DBE must not be so low that even adding 0 more H can't fix it.
    // Actually, the useful check is: current partial DBE shouldn't be so HIGH that
    // even max H can't bring it to a reasonable range, and shouldn't be so LOW
    // that it's already below -0.5 without any more negative contributors.

    // Simple check: current partial DBE (before H) shouldn't be impossibly negative
    // because H only makes it more negative.
    // If dbe is already < -0.5 and we haven't added H yet, adding H only lowers it more
    // — but we might have H already in `assigned`. This check is mainly useful when
    // we've assigned non-H elements and want to know if any H count will work.
    let min_final_dbe = dbe - remaining_h_max as f64 * 0.5;
    let max_final_dbe = dbe;

    // DBE must be >= -0.5 → max_final_dbe >= -0.5
    // DBE shouldn't be absurdly large (practical limit ~40 for 1000 Da organics)
    max_final_dbe >= -0.5 && min_final_dbe <= 50.0
}
