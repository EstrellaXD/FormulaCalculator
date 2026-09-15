use std::fmt;

#[derive(Debug, Clone, Copy)]
pub struct Isotope {
    pub mass_number: u32,
    pub exact_mass: f64,
    pub abundance: f64,
}

#[derive(Debug, Clone, Copy)]
pub struct Element {
    pub symbol: &'static str,
    pub atomic_number: u32,
    pub valence: i32,
    pub isotopes: &'static [Isotope],
    pub monoisotopic_mass: f64,
}

#[derive(Debug, Clone)]
pub struct ElementConstraint {
    pub element: &'static Element,
    pub min_count: u32,
    pub max_count: u32,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum IonMode {
    Positive,
    Negative,
    Neutral,
}

#[derive(Debug, Clone)]
pub struct Adduct {
    pub name: &'static str,
    pub mass_adjustment: f64,
    pub charge: i32,
    pub mode: IonMode,
}

#[derive(Debug, Clone)]
pub struct SearchParams {
    pub target_mass: f64,
    pub tolerance_ppm: f64,
    pub element_constraints: Vec<ElementConstraint>,
    pub adduct: Option<Adduct>,
    pub apply_rules: bool,
    pub strict_rules: bool,
    pub max_results: usize,
    pub calc_isotopes: bool,
}

#[derive(Debug, Clone)]
pub struct IsotopePeak {
    pub mass: f64,
    pub relative_intensity: f64,
}

#[derive(Debug, Clone)]
pub struct FormulaCandidate {
    pub composition: Vec<(&'static str, u32)>,
    pub monoisotopic_mass: f64,
    pub mass_error_da: f64,
    pub mass_error_ppm: f64,
    pub dbe: f64,
    pub formula_string: String,
    pub isotope_pattern: Option<Vec<IsotopePeak>>,
    pub adduct_name: Option<String>,
    pub observed_mz: Option<f64>,
}

impl fmt::Display for FormulaCandidate {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "{}", self.formula_string)
    }
}

pub fn build_formula_string(composition: &[(&str, u32)]) -> String {
    // Hill order: C, H, then alphabetical; alphabetical throughout when there is no C
    let has_c = composition.iter().any(|&(el, n)| el == "C" && n > 0);
    let order: &[&str] = if has_c {
        &["C", "H", "Br", "Cl", "F", "I", "N", "O", "P", "S", "Si"]
    } else {
        &["Br", "C", "Cl", "F", "H", "I", "N", "O", "P", "S", "Si"]
    };
    let mut s = String::new();
    for &sym in order {
        if let Some(&(_, count)) = composition.iter().find(|(el, _)| *el == sym) {
            if count > 0 {
                s.push_str(sym);
                if count > 1 {
                    s.push_str(&count.to_string());
                }
            }
        }
    }
    s
}

pub fn calculate_dbe(composition: &[(&str, u32)]) -> f64 {
    let mut dbe = 1.0_f64;
    for &(sym, count) in composition {
        let valence = match sym {
            "C" | "Si" => 4,
            "H" | "F" | "Cl" | "Br" | "I" | "Na" | "K" => 1,
            "N" | "P" => 3,
            "O" | "S" => 2,
            _ => 0,
        };
        dbe += count as f64 * (valence as f64 - 2.0) / 2.0;
    }
    dbe
}
