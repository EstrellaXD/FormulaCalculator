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

    let h_mass = HYDROGEN.monoisotopic_mass;
    let h_min = h_constraint.map_or(0, |c| c.min_count);
    let h_max = h_constraint.map_or(
        (mass_hi / h_mass).floor() as u32,
        |c| c.max_count,
    );

    let mut results = Vec::new();
    let mut counts = [0u32; MAX_ELEMENTS]; // counts[i] = count of non_h[i], counts[n_non_h] = H count

    // Build element info for partial DBE checking
    let non_h_valences: Vec<(&str, i32)> = non_h
        .iter()
        .map(|c| (c.element.symbol, c.element.valence))
        .collect();

    recurse(
        0,
        0.0,
        &non_h,
        &non_h_valences,
        &suffix_max,
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
        params.max_results,
    );

    // Sort by absolute mass error
    results.sort_by(|a, b| {
        a.mass_error_da
            .abs()
            .partial_cmp(&b.mass_error_da.abs())
            .unwrap()
    });

    if results.len() > params.max_results {
        results.truncate(params.max_results);
    }

    results
}

#[allow(clippy::too_many_arguments)]
fn recurse(
    depth: usize,
    acc_mass: f64,
    non_h: &[&ElementConstraint],
    non_h_valences: &[(&str, i32)],
    suffix_max: &[f64],
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
    max_results: usize,
) {
    if results.len() >= max_results * 2 {
        return; // safety cap to prevent runaway
    }

    if depth == non_h.len() {
        // H layer: direct calculation instead of looping
        solve_hydrogen(
            acc_mass,
            non_h,
            non_h_valences,
            mass_lo,
            mass_hi,
            target_mass,
            h_mass,
            h_min,
            h_max,
            counts,
            depth,
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

    for n in elem.min_count..=max_n {
        let new_mass = acc_mass + n as f64 * elem.element.monoisotopic_mass;

        // Upper bound pruning: already exceeded
        if new_mass > mass_hi {
            break;
        }

        // Lower bound pruning: even with max remaining elements, can't reach target
        // Include H in the remaining capacity
        let max_remaining_non_h = suffix_max[depth + 1];
        let max_remaining_h = h_max as f64 * h_mass;
        if new_mass + max_remaining_non_h + max_remaining_h < mass_lo {
            continue;
        }

        // Partial DBE pruning
        if apply_rules && depth >= 1 {
            // Build partial assignment for DBE check
            let assigned: Vec<(&str, i32, u32)> = (0..=depth)
                .map(|i| {
                    let c = if i < depth { counts[i] } else { n };
                    (non_h_valences[i].0, non_h_valences[i].1, c)
                })
                .collect();
            if !rules::partial_dbe_feasible(&assigned, h_max) {
                continue;
            }
        }

        counts[depth] = n;
        recurse(
            depth + 1,
            new_mass,
            non_h,
            non_h_valences,
            suffix_max,
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
            max_results,
        );
    }
}

#[allow(clippy::too_many_arguments)]
fn solve_hydrogen(
    acc_mass: f64,
    non_h: &[&ElementConstraint],
    _non_h_valences: &[(&str, i32)],
    mass_lo: f64,
    mass_hi: f64,
    target_mass: f64,
    h_mass: f64,
    h_min: u32,
    h_max: u32,
    counts: &mut [u32; MAX_ELEMENTS],
    h_depth: usize,
    results: &mut Vec<FormulaCandidate>,
    apply_rules: bool,
    strict_rules: bool,
) {
    let remaining = target_mass - acc_mass;
    let h_count_f = remaining / h_mass;
    let h_count_rounded = h_count_f.round();

    if h_count_rounded < h_min as f64 || h_count_rounded > h_max as f64 || h_count_rounded < 0.0 {
        return;
    }

    let h_count = h_count_rounded as u32;
    let total_mass = acc_mass + h_count as f64 * h_mass;

    if total_mass < mass_lo || total_mass > mass_hi {
        return;
    }

    counts[h_depth] = h_count;

    // Build composition
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
        return;
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
        return;
    }

    results.push(candidate);
}

/// Parallel search: split by C count using rayon.
pub fn search_parallel(params: &SearchParams) -> Vec<FormulaCandidate> {
    use rayon::prelude::*;

    let abs_tol = params.target_mass * params.tolerance_ppm * 1e-6;
    let mass_hi = params.target_mass + abs_tol;

    // Find C constraint
    let c_idx = params
        .element_constraints
        .iter()
        .position(|c| c.element.symbol == "C");

    let c_constraint = match c_idx {
        Some(i) => &params.element_constraints[i],
        None => return search(params), // no C → fallback to single-threaded
    };

    let c_min = c_constraint.min_count;
    let c_max = c_constraint.max_count.min(
        (mass_hi / c_constraint.element.monoisotopic_mass).floor() as u32,
    );

    let mut all_results: Vec<FormulaCandidate> = (c_min..=c_max)
        .into_par_iter()
        .flat_map(|c_count| {
            // Build params with fixed C count
            let mut sub_params = params.clone();
            for c in &mut sub_params.element_constraints {
                if c.element.symbol == "C" {
                    c.min_count = c_count;
                    c.max_count = c_count;
                }
            }
            search(&sub_params)
        })
        .collect();

    all_results.sort_by(|a, b| {
        a.mass_error_da
            .abs()
            .partial_cmp(&b.mass_error_da.abs())
            .unwrap()
    });

    all_results.truncate(params.max_results);
    all_results
}
