use crate::types::{Element, Isotope};

// IUPAC 2021 / NIST atomic weights and isotopic compositions
// Monoisotopic mass = exact mass of the most abundant isotope

const C_ISOTOPES: &[Isotope] = &[
    Isotope { mass_number: 12, exact_mass: 12.000_000_000_0, abundance: 0.9893 },
    Isotope { mass_number: 13, exact_mass: 13.003_354_835_3, abundance: 0.0107 },
];

const H_ISOTOPES: &[Isotope] = &[
    Isotope { mass_number: 1, exact_mass: 1.007_825_032_07, abundance: 0.999_885 },
    Isotope { mass_number: 2, exact_mass: 2.014_101_778_12, abundance: 0.000_115 },
];

const N_ISOTOPES: &[Isotope] = &[
    Isotope { mass_number: 14, exact_mass: 14.003_074_004_8, abundance: 0.996_36 },
    Isotope { mass_number: 15, exact_mass: 15.000_108_898_3, abundance: 0.003_64 },
];

const O_ISOTOPES: &[Isotope] = &[
    Isotope { mass_number: 16, exact_mass: 15.994_914_619_56, abundance: 0.997_57 },
    Isotope { mass_number: 17, exact_mass: 16.999_131_757_0, abundance: 0.000_38 },
    Isotope { mass_number: 18, exact_mass: 17.999_161_001_0, abundance: 0.002_05 },
];

const P_ISOTOPES: &[Isotope] = &[
    Isotope { mass_number: 31, exact_mass: 30.973_761_63, abundance: 1.0 },
];

const S_ISOTOPES: &[Isotope] = &[
    Isotope { mass_number: 32, exact_mass: 31.972_071_00, abundance: 0.9499 },
    Isotope { mass_number: 33, exact_mass: 32.971_458_76, abundance: 0.0075 },
    Isotope { mass_number: 34, exact_mass: 33.967_867_01, abundance: 0.0425 },
    Isotope { mass_number: 36, exact_mass: 35.967_081_00, abundance: 0.0001 },
];

// Extended elements

const F_ISOTOPES: &[Isotope] = &[
    Isotope { mass_number: 19, exact_mass: 18.998_403_22, abundance: 1.0 },
];

const CL_ISOTOPES: &[Isotope] = &[
    Isotope { mass_number: 35, exact_mass: 34.968_852_68, abundance: 0.7576 },
    Isotope { mass_number: 37, exact_mass: 36.965_902_60, abundance: 0.2424 },
];

const BR_ISOTOPES: &[Isotope] = &[
    Isotope { mass_number: 79, exact_mass: 78.918_337_1, abundance: 0.5069 },
    Isotope { mass_number: 81, exact_mass: 80.916_291_0, abundance: 0.4931 },
];

const I_ISOTOPES: &[Isotope] = &[
    Isotope { mass_number: 127, exact_mass: 126.904_473, abundance: 1.0 },
];

const SI_ISOTOPES: &[Isotope] = &[
    Isotope { mass_number: 28, exact_mass: 27.976_926_532_5, abundance: 0.922_23 },
    Isotope { mass_number: 29, exact_mass: 28.976_494_700, abundance: 0.046_85 },
    Isotope { mass_number: 30, exact_mass: 29.973_770_171, abundance: 0.030_92 },
];

pub static CARBON: Element = Element {
    symbol: "C", atomic_number: 6, valence: 4,
    isotopes: C_ISOTOPES, monoisotopic_mass: 12.000_000_000_0,
};
pub static HYDROGEN: Element = Element {
    symbol: "H", atomic_number: 1, valence: 1,
    isotopes: H_ISOTOPES, monoisotopic_mass: 1.007_825_032_07,
};
pub static NITROGEN: Element = Element {
    symbol: "N", atomic_number: 7, valence: 3,
    isotopes: N_ISOTOPES, monoisotopic_mass: 14.003_074_004_8,
};
pub static OXYGEN: Element = Element {
    symbol: "O", atomic_number: 8, valence: 2,
    isotopes: O_ISOTOPES, monoisotopic_mass: 15.994_914_619_56,
};
pub static PHOSPHORUS: Element = Element {
    symbol: "P", atomic_number: 15, valence: 3, // RDBE convention (Kind & Fiehn eq. 1)
    isotopes: P_ISOTOPES, monoisotopic_mass: 30.973_761_63,
};
pub static SULFUR: Element = Element {
    symbol: "S", atomic_number: 16, valence: 2,
    isotopes: S_ISOTOPES, monoisotopic_mass: 31.972_071_00,
};
pub static FLUORINE: Element = Element {
    symbol: "F", atomic_number: 9, valence: 1,
    isotopes: F_ISOTOPES, monoisotopic_mass: 18.998_403_22,
};
pub static CHLORINE: Element = Element {
    symbol: "Cl", atomic_number: 17, valence: 1,
    isotopes: CL_ISOTOPES, monoisotopic_mass: 34.968_852_68,
};
pub static BROMINE: Element = Element {
    symbol: "Br", atomic_number: 35, valence: 1,
    isotopes: BR_ISOTOPES, monoisotopic_mass: 78.918_337_1,
};
pub static IODINE: Element = Element {
    symbol: "I", atomic_number: 53, valence: 1,
    isotopes: I_ISOTOPES, monoisotopic_mass: 126.904_473,
};
pub static SILICON: Element = Element {
    symbol: "Si", atomic_number: 14, valence: 4,
    isotopes: SI_ISOTOPES, monoisotopic_mass: 27.976_926_532_5,
};

pub fn get_element(symbol: &str) -> Option<&'static Element> {
    match symbol {
        "C" => Some(&CARBON),
        "H" => Some(&HYDROGEN),
        "N" => Some(&NITROGEN),
        "O" => Some(&OXYGEN),
        "P" => Some(&PHOSPHORUS),
        "S" => Some(&SULFUR),
        "F" => Some(&FLUORINE),
        "Cl" => Some(&CHLORINE),
        "Br" => Some(&BROMINE),
        "I" => Some(&IODINE),
        "Si" => Some(&SILICON),
        _ => None,
    }
}

/// Parse an element string like "CHNOPS" or "CHNOPSClBr" into element references.
/// Multi-character symbols (Cl, Br, Si) are recognized.
pub fn parse_element_string(s: &str) -> Vec<&'static Element> {
    let mut elements = Vec::new();
    let chars: Vec<char> = s.chars().collect();
    let mut i = 0;
    while i < chars.len() {
        // Try two-character symbol first
        if i + 1 < chars.len() && chars[i + 1].is_lowercase() {
            let sym: String = chars[i..=i + 1].iter().collect();
            if let Some(el) = get_element(&sym) {
                elements.push(el);
                i += 2;
                continue;
            }
        }
        // Single-character symbol
        let sym = chars[i].to_string();
        if let Some(el) = get_element(&sym) {
            elements.push(el);
        }
        i += 1;
    }
    elements
}

pub const ELECTRON_MASS: f64 = 0.000_548_579_909_07;
pub const PROTON_MASS: f64 = 1.007_276_47;
