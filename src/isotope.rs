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

        // Single-atom pattern indexed by nucleon shift from the lightest isotope
        let base = elem.isotopes[0].mass_number;
        let width = elem.isotopes[elem.isotopes.len() - 1].mass_number - base + 1;
        let mut atom_pattern = vec![(0.0, 0.0); width as usize];
        for iso in elem.isotopes {
            atom_pattern[(iso.mass_number - base) as usize] = (iso.exact_mass, iso.abundance);
        }

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

/// Convolve two isotope distributions indexed by nucleon shift.
/// Each entry is (probability-weighted mean mass, probability). Grouping by the
/// integer shift, not by rounding exact masses, keeps M+1 and M+2 apart when
/// the mass defect puts a peak near a .5 Da boundary.
fn convolve(a: &[(f64, f64)], b: &[(f64, f64)]) -> Vec<(f64, f64)> {
    let mut result = vec![(0.0, 0.0); a.len() + b.len() - 1];

    for (i, &(ma, pa)) in a.iter().enumerate() {
        for (j, &(mb, pb)) in b.iter().enumerate() {
            let p = pa * pb;
            result[i + j].0 += (ma + mb) * p;
            result[i + j].1 += p;
        }
    }

    for r in &mut result {
        if r.1 > 0.0 {
            r.0 /= r.1;
        }
    }
    // Prune the negligible high-mass tail to keep convolution tractable
    while result.len() > 1 && result[result.len() - 1].1 <= 1e-12 {
        result.pop();
    }
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
