use clap::Parser;
use std::io::{self, BufRead};

use mfcalc::adduct;
use mfcalc::bounds;
use mfcalc::elements;
use mfcalc::isotope;
use mfcalc::output::{self, OutputFormat};
use mfcalc::search;
use mfcalc::types::{IonMode, SearchParams};

#[derive(Parser, Debug)]
#[command(
    name = "mfcalc",
    about = "Molecular formula calculator from exact mass",
    version
)]
struct Cli {
    /// Target exact mass or m/z value (Da)
    mass: Option<f64>,

    /// Mass tolerance in ppm
    #[arg(short, long, default_value_t = 2.0)]
    ppm: f64,

    /// Elements to consider (e.g. "CHNOPS", "CHNOPSCl")
    #[arg(short, long, default_value = "CHNOPS")]
    elements: String,

    /// Adduct type, e.g. "[M+H]+". If omitted with -m pos/neg, all adducts are tried.
    #[arg(short, long)]
    adduct: Option<String>,

    /// Ion mode: pos / neg / neutral. pos/neg automatically tries all adducts.
    #[arg(short, long, default_value = "neutral")]
    mode: String,

    /// Use strict element ratio filters
    #[arg(long)]
    strict: bool,

    /// Disable Seven Golden Rules filtering
    #[arg(long)]
    no_rules: bool,

    /// Calculate isotope patterns for results
    #[arg(long)]
    isotope: bool,

    /// Maximum results to return
    #[arg(long, default_value_t = 50)]
    max_results: usize,

    /// Output format: table / json / csv
    #[arg(short, long, default_value = "table")]
    output: String,

    /// Batch mode: read masses from file (one per line)
    #[arg(long)]
    batch: Option<String>,

    /// Use parallel search (split by C count)
    #[arg(long)]
    parallel: bool,
}

fn main() {
    let cli = Cli::parse();

    let output_format = match cli.output.as_str() {
        "json" => OutputFormat::Json,
        "csv" => OutputFormat::Csv,
        _ => OutputFormat::Table,
    };

    let ion_mode = match cli.mode.as_str() {
        "pos" | "positive" => IonMode::Positive,
        "neg" | "negative" => IonMode::Negative,
        _ => IonMode::Neutral,
    };

    let masses: Vec<f64> = if let Some(ref batch_file) = cli.batch {
        read_batch_file(batch_file)
    } else if let Some(mass) = cli.mass {
        vec![mass]
    } else {
        eprintln!("Error: provide a mass value or --batch file");
        std::process::exit(1);
    };

    // Resolve which adducts to try:
    // - Specific -a "[M+H]+": only that one
    // - -m pos/neg without -a: all adducts for that mode
    // - -m neutral (default): no adduct, treat input as neutral mass
    let adducts_to_try: Vec<Option<&mfcalc::types::Adduct>> =
        if let Some(ref adduct_name) = cli.adduct {
            match adduct::find_adduct(adduct_name) {
                Some(a) => vec![Some(a)],
                None => {
                    eprintln!("Error: unknown adduct '{adduct_name}'");
                    eprintln!("Available adducts:");
                    for a in adduct::ADDUCTS_POSITIVE.iter().chain(adduct::ADDUCTS_NEGATIVE) {
                        eprintln!("  {}", a.name);
                    }
                    std::process::exit(1);
                }
            }
        } else if ion_mode != IonMode::Neutral {
            // Auto-try all adducts for the given ion mode
            adduct::adducts_for_mode(ion_mode)
                .iter()
                .map(Some)
                .collect()
        } else {
            vec![None]
        };

    let show_adduct_column = adducts_to_try.len() > 1
        || adducts_to_try.iter().any(|a| a.is_some());

    let elem_list = elements::parse_element_string(&cli.elements);
    if elem_list.is_empty() {
        eprintln!("Error: no valid elements in '{}'", cli.elements);
        std::process::exit(1);
    }

    for (mass_idx, &input_mass) in masses.iter().enumerate() {
        // Collect candidates across all adducts for this mass, then merge + sort
        let mut all_candidates = Vec::new();

        for adduct_opt in &adducts_to_try {
            let neutral_mass = match adduct_opt {
                Some(a) => adduct::mz_to_neutral(input_mass, a),
                None => input_mass,
            };

            if neutral_mass <= 0.0 {
                continue;
            }

            let mut constraints = bounds::default_constraints(&elem_list);
            bounds::compute_bounds(neutral_mass, cli.ppm, &mut constraints);

            let params = SearchParams {
                target_mass: neutral_mass,
                tolerance_ppm: cli.ppm,
                element_constraints: constraints,
                adduct: adduct_opt.cloned(),
                apply_rules: !cli.no_rules,
                strict_rules: cli.strict,
                max_results: cli.max_results,
                calc_isotopes: cli.isotope,
            };

            let mut candidates = if cli.parallel {
                search::search_parallel(&params)
            } else {
                search::search(&params)
            };

            // Tag each candidate with its adduct and observed m/z
            for c in &mut candidates {
                if let Some(a) = adduct_opt {
                    c.adduct_name = Some(a.name.to_string());
                    c.observed_mz = Some(adduct::neutral_to_mz(c.monoisotopic_mass, a));
                }
            }

            // Calculate isotope patterns if requested
            if cli.isotope {
                let elem_refs: Vec<&'static mfcalc::types::Element> =
                    elem_list.iter().copied().collect();
                for c in &mut candidates {
                    c.isotope_pattern =
                        Some(isotope::isotope_distribution(&c.composition, &elem_refs));
                }
            }

            all_candidates.extend(candidates);
        }

        // Sort merged results by absolute mass error
        all_candidates.sort_by(|a, b| {
            a.mass_error_da
                .abs()
                .partial_cmp(&b.mass_error_da.abs())
                .unwrap()
        });
        all_candidates.truncate(cli.max_results);

        // Print header
        let abs_tol = input_mass * cli.ppm * 1e-6;
        let elem_syms: Vec<&str> = elem_list.iter().map(|e| e.symbol).collect();
        let rules_str = if cli.no_rules {
            "disabled"
        } else if cli.strict {
            "Seven Golden Rules (strict)"
        } else {
            "Seven Golden Rules (loose)"
        };

        if ion_mode == IonMode::Neutral {
            println!(
                "Target mass: {:.6} Da | Tolerance: {:.1} ppm (\u{00b1}{:.6} Da)",
                input_mass, cli.ppm, abs_tol
            );
        } else {
            let mode_str = if ion_mode == IonMode::Positive { "pos" } else { "neg" };
            println!(
                "Input m/z: {:.6} | Mode: {} | Tolerance: {:.1} ppm",
                input_mass, mode_str, cli.ppm
            );
        }
        println!("Elements: {} | Rules: {}", elem_syms.join(" "), rules_str);
        println!();

        output::print_results(&all_candidates, output_format, cli.isotope, show_adduct_column);

        if mass_idx + 1 < masses.len() {
            println!("---");
        }
    }
}

fn read_batch_file(path: &str) -> Vec<f64> {
    let file = match std::fs::File::open(path) {
        Ok(f) => f,
        Err(e) => {
            eprintln!("Error reading batch file '{path}': {e}");
            std::process::exit(1);
        }
    };

    io::BufReader::new(file)
        .lines()
        .filter_map(|line| {
            let line = line.ok()?;
            let trimmed = line.trim();
            if trimmed.is_empty() || trimmed.starts_with('#') {
                return None;
            }
            trimmed.parse::<f64>().ok()
        })
        .collect()
}
