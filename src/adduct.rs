use crate::types::{Adduct, IonMode};

pub static ADDUCTS_POSITIVE: &[Adduct] = &[
    Adduct { name: "[M+H]+",   mass_adjustment: 1.007_276_47,  charge: 1, mode: IonMode::Positive },
    Adduct { name: "[M+Na]+",  mass_adjustment: 22.989_218_02,  charge: 1, mode: IonMode::Positive },
    Adduct { name: "[M+K]+",   mass_adjustment: 38.963_706_68,  charge: 1, mode: IonMode::Positive },
    Adduct { name: "[M+NH4]+", mass_adjustment: 18.034_164_07,  charge: 1, mode: IonMode::Positive },
    Adduct { name: "[M+2H]2+", mass_adjustment: 2.014_552_94,   charge: 2, mode: IonMode::Positive },
];

pub static ADDUCTS_NEGATIVE: &[Adduct] = &[
    Adduct { name: "[M-H]-",      mass_adjustment: -1.007_276_47,  charge: -1, mode: IonMode::Negative },
    Adduct { name: "[M+Cl]-",     mass_adjustment: 34.969_402_00,  charge: -1, mode: IonMode::Negative },
    Adduct { name: "[M+HCOO]-",   mass_adjustment: 44.998_201_60,  charge: -1, mode: IonMode::Negative },
    Adduct { name: "[M+CH3COO]-", mass_adjustment: 59.013_851_38,  charge: -1, mode: IonMode::Negative },
    Adduct { name: "[M-2H]2-",    mass_adjustment: -2.014_552_94,  charge: -2, mode: IonMode::Negative },
];

/// Convert observed m/z to neutral monoisotopic mass given an adduct.
pub fn mz_to_neutral(mz: f64, adduct: &Adduct) -> f64 {
    mz * adduct.charge.unsigned_abs() as f64 - adduct.mass_adjustment
}

/// Convert neutral monoisotopic mass to expected m/z given an adduct.
pub fn neutral_to_mz(neutral_mass: f64, adduct: &Adduct) -> f64 {
    (neutral_mass + adduct.mass_adjustment) / adduct.charge.unsigned_abs() as f64
}

/// Find an adduct by name (case-insensitive, flexible matching).
pub fn find_adduct(name: &str) -> Option<&'static Adduct> {
    let normalized = name
        .replace(' ', "")
        .replace('⁺', "+")
        .replace('⁻', "-")
        .to_lowercase();

    ADDUCTS_POSITIVE
        .iter()
        .chain(ADDUCTS_NEGATIVE.iter())
        .find(|a| {
            a.name
                .replace(' ', "")
                .to_lowercase()
                == normalized
        })
}

/// Get all adducts for a given ion mode.
pub fn adducts_for_mode(mode: IonMode) -> &'static [Adduct] {
    match mode {
        IonMode::Positive => ADDUCTS_POSITIVE,
        IonMode::Negative => ADDUCTS_NEGATIVE,
        IonMode::Neutral => &[],
    }
}
