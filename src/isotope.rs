use crate::types::{Element, IsotopePeak};

/// Calculate the theoretical isotope distribution for a molecular formula.
/// Returns peaks normalized so the highest peak = 1.0.
pub fn isotope_distribution(
    composition: &[(&'static str, u32)],
    elements: &[&'static Element],
) -> Vec<IsotopePeak> {
    // Start with a delta at mass=0, probability=1.0
    let mut pattern: Vec<(f64, f64)> = vec![(0.0, 1.0)];

    for &(symbol, count) in composition {
        if count == 0 {
            continue;
        }
        let elem = match elements.iter().find(|e| e.symbol == symbol) {
            Some(e) => e,
            None => continue,
        };

        // Build single-atom isotope pattern for this element
        let atom_pattern: Vec<(f64, f64)> = elem
            .isotopes
            .iter()
            .map(|iso| (iso.exact_mass, iso.abundance))
            .collect();

        // Raise to power `count` using repeated squaring convolution
        let elem_pattern = power_convolve(&atom_pattern, count);

        // Convolve with running total
        pattern = convolve(&pattern, &elem_pattern);
    }

    // Convert to IsotopePeak, normalize
    if pattern.is_empty() {
        return vec![];
    }

    let max_intensity = pattern.iter().map(|&(_, p)| p).fold(0.0_f64, f64::max);
    if max_intensity <= 0.0 {
        return vec![];
    }

    pattern
        .into_iter()
        .map(|(mass, prob)| IsotopePeak {
            mass,
            relative_intensity: prob / max_intensity,
        })
        .filter(|p| p.relative_intensity > 1e-4) // drop negligible peaks
        .collect()
}

/// Convolve two isotope distributions.
/// Each distribution is a list of (mass, probability) pairs.
/// Peaks with the same nominal mass are merged.
fn convolve(a: &[(f64, f64)], b: &[(f64, f64)]) -> Vec<(f64, f64)> {
    let mut result: Vec<(f64, f64)> = Vec::with_capacity(a.len() * b.len());

    for &(ma, pa) in a {
        for &(mb, pb) in b {
            result.push((ma + mb, pa * pb));
        }
    }

    merge_by_nominal_mass(&mut result);
    // Prune tiny peaks to keep convolution tractable
    result.retain(|&(_, p)| p > 1e-12);
    result
}

/// Raise a base pattern to power n using repeated squaring.
fn power_convolve(base: &[(f64, f64)], mut n: u32) -> Vec<(f64, f64)> {
    if n == 0 {
        return vec![(0.0, 1.0)];
    }

    let mut result: Vec<(f64, f64)> = vec![(0.0, 1.0)];
    let mut current = base.to_vec();

    while n > 0 {
        if n & 1 == 1 {
            result = convolve(&result, &current);
        }
        current = convolve(&current, &current);
        n >>= 1;
    }

    result
}

/// Merge peaks with the same nominal mass (round to nearest integer).
/// Masses are averaged weighted by probability; probabilities are summed.
fn merge_by_nominal_mass(peaks: &mut Vec<(f64, f64)>) {
    if peaks.is_empty() {
        return;
    }

    // Sort by nominal mass
    peaks.sort_by(|a, b| {
        (a.0.round() as i64)
            .cmp(&(b.0.round() as i64))
            .then(a.0.partial_cmp(&b.0).unwrap())
    });

    let mut merged: Vec<(f64, f64)> = Vec::new();
    let mut current_nominal = peaks[0].0.round() as i64;
    let mut mass_sum = 0.0_f64;
    let mut prob_sum = 0.0_f64;

    for &(mass, prob) in peaks.iter() {
        let nominal = mass.round() as i64;
        if nominal == current_nominal {
            mass_sum += mass * prob;
            prob_sum += prob;
        } else {
            if prob_sum > 0.0 {
                merged.push((mass_sum / prob_sum, prob_sum));
            }
            current_nominal = nominal;
            mass_sum = mass * prob;
            prob_sum = prob;
        }
    }
    if prob_sum > 0.0 {
        merged.push((mass_sum / prob_sum, prob_sum));
    }

    *peaks = merged;
}
