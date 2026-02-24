use crate::types::ElementConstraint;

/// Compute the maximum count of each element based on target mass + tolerance.
/// max_count = floor((target_mass + abs_tolerance) / monoisotopic_mass)
/// Then cap by any user-specified maximum.
pub fn compute_bounds(
    target_mass: f64,
    tolerance_ppm: f64,
    constraints: &mut [ElementConstraint],
) {
    let abs_tol = target_mass * tolerance_ppm * 1e-6;
    let mass_hi = target_mass + abs_tol;

    for c in constraints.iter_mut() {
        let mass_based_max = (mass_hi / c.element.monoisotopic_mass).floor() as u32;
        c.max_count = c.max_count.min(mass_based_max);
        // Ensure min <= max
        if c.min_count > c.max_count {
            c.min_count = c.max_count;
        }
    }
}

/// Build default constraints for a list of elements with sensible max bounds.
/// The caller should then call compute_bounds() to refine based on target mass.
pub fn default_constraints(
    elements: &[&'static crate::types::Element],
) -> Vec<ElementConstraint> {
    elements
        .iter()
        .map(|&el| ElementConstraint {
            element: el,
            min_count: 0,
            max_count: u32::MAX, // will be refined by compute_bounds
        })
        .collect()
}
