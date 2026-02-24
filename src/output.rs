use crate::types::FormulaCandidate;
use comfy_table::{presets::UTF8_FULL_CONDENSED, Cell, ContentArrangement, Table};

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum OutputFormat {
    Table,
    Json,
    Csv,
}

pub fn print_results(
    candidates: &[FormulaCandidate],
    format: OutputFormat,
    show_isotopes: bool,
    show_adduct: bool,
) {
    match format {
        OutputFormat::Table => print_table(candidates, show_isotopes, show_adduct),
        OutputFormat::Json => print_json(candidates),
        OutputFormat::Csv => print_csv(candidates, show_isotopes, show_adduct),
    }
}

fn print_table(candidates: &[FormulaCandidate], show_isotopes: bool, show_adduct: bool) {
    if candidates.is_empty() {
        println!("No candidates found.");
        return;
    }

    let mut table = Table::new();
    table
        .load_preset(UTF8_FULL_CONDENSED)
        .set_content_arrangement(ContentArrangement::Dynamic);

    let mut header = vec!["#", "Formula", "Exact Mass", "\u{0394}(Da)", "\u{0394}(ppm)", "DBE"];
    if show_adduct {
        header.push("Adduct");
        header.push("m/z");
    }
    if show_isotopes {
        header.push("M+1%");
        header.push("M+2%");
    }
    table.set_header(header);

    for (i, c) in candidates.iter().enumerate() {
        let mut row: Vec<Cell> = vec![
            Cell::new(i + 1),
            Cell::new(&c.formula_string),
            Cell::new(format!("{:.6}", c.monoisotopic_mass)),
            Cell::new(format!("{:+.6}", c.mass_error_da)),
            Cell::new(format!("{:+.2}", c.mass_error_ppm)),
            Cell::new(format!("{:.1}", c.dbe)),
        ];

        if show_adduct {
            row.push(Cell::new(
                c.adduct_name.as_deref().unwrap_or("-"),
            ));
            row.push(Cell::new(
                c.observed_mz
                    .map(|mz| format!("{:.4}", mz))
                    .unwrap_or_else(|| "-".to_string()),
            ));
        }

        if show_isotopes {
            if let Some(ref pattern) = c.isotope_pattern {
                let m1 = pattern
                    .iter()
                    .find(|p| (p.mass - c.monoisotopic_mass).round() as i32 == 1)
                    .map(|p| format!("{:.2}", p.relative_intensity * 100.0))
                    .unwrap_or_else(|| "-".to_string());
                let m2 = pattern
                    .iter()
                    .find(|p| (p.mass - c.monoisotopic_mass).round() as i32 == 2)
                    .map(|p| format!("{:.2}", p.relative_intensity * 100.0))
                    .unwrap_or_else(|| "-".to_string());
                row.push(Cell::new(m1));
                row.push(Cell::new(m2));
            } else {
                row.push(Cell::new("-"));
                row.push(Cell::new("-"));
            }
        }

        table.add_row(row);
    }

    println!("{table}");
    println!("\nFound {} candidate(s).", candidates.len());
}

fn print_json(candidates: &[FormulaCandidate]) {
    #[derive(serde::Serialize)]
    struct JsonCandidate {
        formula: String,
        exact_mass: f64,
        mass_error_da: f64,
        mass_error_ppm: f64,
        dbe: f64,
        #[serde(skip_serializing_if = "Option::is_none")]
        adduct: Option<String>,
        #[serde(skip_serializing_if = "Option::is_none")]
        observed_mz: Option<f64>,
        composition: Vec<(String, u32)>,
        #[serde(skip_serializing_if = "Option::is_none")]
        isotope_pattern: Option<Vec<JsonPeak>>,
    }

    #[derive(serde::Serialize)]
    struct JsonPeak {
        mass: f64,
        relative_intensity: f64,
    }

    let output: Vec<JsonCandidate> = candidates
        .iter()
        .map(|c| JsonCandidate {
            formula: c.formula_string.clone(),
            exact_mass: c.monoisotopic_mass,
            mass_error_da: c.mass_error_da,
            mass_error_ppm: c.mass_error_ppm,
            dbe: c.dbe,
            adduct: c.adduct_name.clone(),
            observed_mz: c.observed_mz,
            composition: c.composition.iter().map(|&(s, n)| (s.to_string(), n)).collect(),
            isotope_pattern: c.isotope_pattern.as_ref().map(|pat| {
                pat.iter()
                    .map(|p| JsonPeak {
                        mass: p.mass,
                        relative_intensity: p.relative_intensity,
                    })
                    .collect()
            }),
        })
        .collect();

    println!("{}", serde_json::to_string_pretty(&output).unwrap());
}

fn print_csv(candidates: &[FormulaCandidate], show_isotopes: bool, show_adduct: bool) {
    let mut header = "formula,exact_mass,mass_error_da,mass_error_ppm,dbe".to_string();
    if show_adduct {
        header.push_str(",adduct,observed_mz");
    }
    if show_isotopes {
        header.push_str(",m_plus_1_pct,m_plus_2_pct");
    }
    println!("{header}");

    for c in candidates {
        let mut line = format!(
            "{},{:.6},{:+.6},{:+.2},{:.1}",
            c.formula_string, c.monoisotopic_mass, c.mass_error_da, c.mass_error_ppm, c.dbe
        );
        if show_adduct {
            line.push_str(&format!(
                ",{},{}",
                c.adduct_name.as_deref().unwrap_or(""),
                c.observed_mz
                    .map(|mz| format!("{:.4}", mz))
                    .unwrap_or_default()
            ));
        }
        if show_isotopes {
            if let Some(ref pattern) = c.isotope_pattern {
                let m1 = pattern
                    .iter()
                    .find(|p| (p.mass - c.monoisotopic_mass).round() as i32 == 1)
                    .map(|p| format!("{:.4}", p.relative_intensity * 100.0))
                    .unwrap_or_else(|| "0".to_string());
                let m2 = pattern
                    .iter()
                    .find(|p| (p.mass - c.monoisotopic_mass).round() as i32 == 2)
                    .map(|p| format!("{:.4}", p.relative_intensity * 100.0))
                    .unwrap_or_else(|| "0".to_string());
                line.push_str(&format!(",{m1},{m2}"));
            } else {
                line.push_str(",0,0");
            }
        }
        println!("{line}");
    }
}
