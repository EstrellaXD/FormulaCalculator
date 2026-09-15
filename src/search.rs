use crate::elements::HYDROGEN;
use crate::rules;
use crate::types::{
    build_formula_string, calculate_dbe, ElementConstraint, FormulaCandidate, SearchParams,
};

/// Maximum number of element slots (fixed-size array for zero-alloc recursion).
const MAX_ELEMENTS: usize = 12;

/// Search for all molecular formulas matching the target mass within tolerance.
pub fn search(params: &SearchParams) -> Vec<FormulaCandidate> {
    let abs_tol = params.target_mass * params.tolerance_ppm * 1e-6;
    let mass_lo = params.target_mass - abs_tol;
    let mass_hi = params.target_mass + abs_tol;

    // Separate H from other elements; H goes last for direct-calculation optimization
    let mut non_h: Vec<&ElementConstraint> = Vec::new();
    let mut h_constraint: Option<&ElementConstraint> = None;

    for c in &params.element_constraints {
        if c.element.symbol == "H" {
            h_constraint = Some(c);
        } else {
            non_h.push(c);
        }
    }

    // Sort non-H elements by monoisotopic mass descending (heavier first = better pruning)
    non_h.sort_by(|a, b| {
        b.element
            .monoisotopic_mass
            .partial_cmp(&a.element.monoisotopic_mass)
            .unwrap()
    });

    let n_non_h = non_h.len();

    // Precompute suffix max mass: suffix_max[i] = sum of max_count * mono_mass for elements i..n
    // (not including H; H is handled separately)
    let mut suffix_max = vec![0.0_f64; n_non_h + 1];
    for i in (0..n_non_h).rev() {
        suffix_max[i] = suffix_max[i + 1]
            + non_h[i].max_count as f64 * non_h[i].element.monoisotopic_mass;
    }

    // suffix_h_cap[i] = most H that elements i..n can add to the RDBE cap
    let mut suffix_h_cap = vec![0_i64; n_non_h + 1];
    for i in (0..n_non_h).rev() {
        suffix_h_cap[i] = suffix_h_cap[i + 1]
            + non_h[i].max_count as i64 * (non_h[i].element.valence as i64 - 2).max(0);
    }

    let h_mass = HYDROGEN.monoisotopic_mass;
    let h_min = h_constraint.map_or(0, |c| c.min_count);
    let h_max = h_constraint.map_or(0, |c| c.max_count);

    let mut results = Vec::new();
    let mut counts = [0u32; MAX_ELEMENTS]; // counts[i] = count of non_h[i], counts[n_non_h] = H count

    recurse(
        0,
        0.0,
        3,
        &non_h,
        &suffix_max,
        &suffix_h_cap,
        mass_lo,
        mass_hi,
        params.target_mass,
        h_mass,
        h_min,
        h_max,
        &mut counts,
        &mut results,
        params.apply_rules,
        params.strict_rules,
    );

    // Keep the max_results smallest absolute mass errors
    let by_error = |a: &FormulaCandidate, b: &FormulaCandidate| {
        a.mass_error_da.abs().total_cmp(&b.mass_error_da.abs())
    };
    if results.len() > params.max_results {
        results.select_nth_unstable_by(params.max_results, by_error);
        results.truncate(params.max_results);
    }
    results.sort_by(by_error);
    results
}

#[allow(clippy::too_many_arguments)]
fn recurse(
    depth: usize,
    acc_mass: f64,
    h_cap: i64,
    non_h: &[&ElementConstraint],
    suffix_max: &[f64],
    suffix_h_cap: &[i64],
    mass_lo: f64,
    mass_hi: f64,
    target_mass: f64,
    h_mass: f64,
    h_min: u32,
    h_max: u32,
    counts: &mut [u32; MAX_ELEMENTS],
    results: &mut Vec<FormulaCandidate>,
    apply_rules: bool,
    strict_rules: bool,
) {
    if depth == non_h.len() {
        // H layer: direct calculation instead of looping
        solve_hydrogen(
            acc_mass,
            h_cap,
            non_h,
            mass_lo,
            mass_hi,
            target_mass,
            h_mass,
            h_min,
            h_max,
            counts,
            results,
            apply_rules,
            strict_rules,
        );
        return;
    }

    let elem = non_h[depth];
    let remaining_for_elem = mass_hi - acc_mass;
    if remaining_for_elem < 0.0 {
        return;
    }
    let mass_based_max = (remaining_for_elem / elem.element.monoisotopic_mass).floor() as u32;
    let max_n = elem.max_count.min(mass_based_max);

    // Skip counts too small to reach mass_lo. Reach grows with n (every element outweighs
    // H), so solve the reach bound for n; with rules, max(cap, 0) >= cap and >= 0 makes
    // the smaller of the two linear solutions a safe start. The check in the loop stays exact.
    let m = elem.element.monoisotopic_mass;
    let short = mass_lo - acc_mass - suffix_max[depth + 1];
    let mut start = (short - h_max as f64 * h_mass) / m;
    if apply_rules {
        let cap = (h_cap + suffix_h_cap[depth + 1]) as f64;
        let dv = elem.element.valence as f64 - 2.0;
        start = start.max((short / m).min((short - cap * h_mass) / (m + dv * h_mass)));
    }
    let min_n = elem.min_count.max(start.floor().max(0.0) as u32);

    for n in min_n..=max_n {
        let new_mass = acc_mass + n as f64 * elem.element.monoisotopic_mass;

        // Upper bound pruning: already exceeded
        if new_mass > mass_hi {
            break;
        }

        // Lower bound pruning: even with max remaining elements, can't reach target
        let new_h_cap = h_cap + n as i64 * (elem.element.valence as i64 - 2);
        let max_remaining_h = h_limit(new_h_cap + suffix_h_cap[depth + 1], h_max, apply_rules);
        if new_mass + suffix_max[depth + 1] + max_remaining_h as f64 * h_mass < mass_lo {
            continue;
        }

        counts[depth] = n;
        recurse(
            depth + 1,
            new_mass,
            new_h_cap,
            non_h,
            suffix_max,
            suffix_h_cap,
            mass_lo,
            mass_hi,
            target_mass,
            h_mass,
            h_min,
            h_max,
            counts,
            results,
            apply_rules,
            strict_rules,
        );
    }
}

#[allow(clippy::too_many_arguments)]
fn solve_hydrogen(
    acc_mass: f64,
    h_cap: i64,
    non_h: &[&ElementConstraint],
    mass_lo: f64,
    mass_hi: f64,
    target_mass: f64,
    h_mass: f64,
    h_min: u32,
    h_max: u32,
    counts: &[u32; MAX_ELEMENTS],
    results: &mut Vec<FormulaCandidate>,
    apply_rules: bool,
    strict_rules: bool,
) {
    // Every H count whose total mass falls inside [mass_lo, mass_hi]
    let lo = ((mass_lo - acc_mass) / h_mass).ceil().max(h_min as f64);
    let hi = ((mass_hi - acc_mass) / h_mass)
        .floor()
        .min(h_limit(h_cap, h_max, apply_rules) as f64);
    if lo > hi {
        return;
    }

    for h_count in lo as u32..=hi as u32 {
        let total_mass = acc_mass + h_count as f64 * h_mass;

        let mut composition: Vec<(&'static str, u32)> = Vec::with_capacity(non_h.len() + 1);
        for (i, c) in non_h.iter().enumerate() {
            if counts[i] > 0 {
                composition.push((c.element.symbol, counts[i]));
            }
        }
        if h_count > 0 {
            composition.push(("H", h_count));
        }

        // Must have at least one element
        if composition.is_empty() {
            continue;
        }

        let dbe = calculate_dbe(&composition);
        let error_da = total_mass - target_mass;
        let error_ppm = error_da / target_mass * 1e6;

        let candidate = FormulaCandidate {
            formula_string: build_formula_string(&composition),
            composition,
            monoisotopic_mass: total_mass,
            mass_error_da: error_da,
            mass_error_ppm: error_ppm,
            dbe,
            isotope_pattern: None,
            adduct_name: None,
            observed_mz: None,
        };

        if apply_rules && !rules::validate(&candidate, strict_rules) {
            continue;
        }

        results.push(candidate);
    }
}

/// H upper bound. With rules on, RDBE >= -0.5 implies H <= 3 + sum(n * (valence - 2)),
/// where `h_cap` is that right-hand side.
fn h_limit(h_cap: i64, h_max: u32, apply_rules: bool) -> u32 {
    if apply_rules {
        h_cap.clamp(0, h_max as i64) as u32
    } else {
        h_max
    }
}
