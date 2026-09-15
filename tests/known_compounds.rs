use mfcalc::bounds;
use mfcalc::elements;
use mfcalc::isotope;
use mfcalc::search;
use mfcalc::types::SearchParams;

fn search_for(mass: f64, ppm: f64, elem_str: &str) -> Vec<mfcalc::types::FormulaCandidate> {
    let elem_list = elements::parse_element_string(elem_str);
    let mut constraints = bounds::default_constraints(&elem_list);
    bounds::compute_bounds(mass, ppm, &mut constraints);

    let params = SearchParams {
        target_mass: mass,
        tolerance_ppm: ppm,
        element_constraints: constraints,
        adduct: None,
        apply_rules: true,
        strict_rules: false,
        max_results: 100,
        calc_isotopes: false,
    };
    search::search(&params)
}

fn assert_formula_found(mass: f64, expected_formula: &str, ppm: f64) {
    let results = search_for(mass, ppm, "CHNOPS");
    let found = results.iter().any(|c| c.formula_string == expected_formula);
    assert!(
        found,
        "Expected {} at mass {:.5} within {} ppm, but not found. Got: {:?}",
        expected_formula,
        mass,
        ppm,
        results.iter().map(|c| &c.formula_string).collect::<Vec<_>>()
    );
}

#[test]
fn test_glycine() {
    assert_formula_found(75.03203, "C2H5NO2", 2.0);
}

#[test]
fn test_alanine() {
    assert_formula_found(89.04768, "C3H7NO2", 2.0);
}

#[test]
fn test_glucose() {
    assert_formula_found(180.06339, "C6H12O6", 2.0);
}

#[test]
fn test_adenine() {
    assert_formula_found(135.0545, "C5H5N5", 2.0);
}

#[test]
fn test_atp() {
    assert_formula_found(506.99576, "C10H16N5O13P3", 2.0);
}

#[test]
fn test_caffeine() {
    assert_formula_found(194.08038, "C8H10N4O2", 2.0);
}

#[test]
fn test_cholesterol() {
    assert_formula_found(386.35486, "C27H46O", 2.0);
}

#[test]
fn test_putrescine() {
    assert_formula_found(88.10005, "C4H12N2", 2.0);
}

#[test]
fn test_spermidine() {
    assert_formula_found(145.15789, "C7H19N3", 2.0);
}

#[test]
fn test_dbe_values() {
    let results = search_for(194.08038, 2.0, "CHNOPS");
    let caffeine = results.iter().find(|c| c.formula_string == "C8H10N4O2").unwrap();
    assert!((caffeine.dbe - 6.0).abs() < 0.01, "Caffeine DBE should be 6.0");
}

#[test]
fn test_mass_error_within_tolerance() {
    let results = search_for(180.06339, 2.0, "CHNOPS");
    for c in &results {
        assert!(
            c.mass_error_ppm.abs() <= 2.0,
            "Mass error {:.2} ppm exceeds tolerance for {}",
            c.mass_error_ppm,
            c.formula_string
        );
    }
}

#[test]
fn test_wider_tolerance_gives_more_results() {
    let narrow = search_for(180.06339, 0.5, "CHNOPS");
    let wide = search_for(180.06339, 5.0, "CHNOPS");
    assert!(
        wide.len() >= narrow.len(),
        "Wider tolerance should give at least as many results"
    );
}

#[test]
fn test_isotope_distribution_glucose() {
    let elem_list = elements::parse_element_string("CHNOPS");
    let composition: Vec<(&str, u32)> = vec![("C", 6), ("H", 12), ("O", 6)];
    let elem_refs: Vec<&'static mfcalc::types::Element> = elem_list.iter().copied().collect();
    let pattern = isotope::isotope_distribution(&composition, &elem_refs);

    assert!(!pattern.is_empty());
    // M+0 should be the base peak (1.0)
    assert!((pattern[0].relative_intensity - 1.0).abs() < 0.01);
    // M+1 for C6H12O6: ~6 * 1.07% (C13) + small contributions ≈ 6-7%
    let m1 = pattern.iter().find(|p| (p.mass - pattern[0].mass).round() as i32 == 1);
    assert!(m1.is_some(), "Should have M+1 peak");
    let m1_pct = m1.unwrap().relative_intensity * 100.0;
    assert!(m1_pct > 5.0 && m1_pct < 10.0, "M+1 for glucose should be ~6-7%, got {:.1}%", m1_pct);
}

#[test]
fn test_adduct_conversion() {
    use mfcalc::adduct;

    let h_adduct = adduct::find_adduct("[M+H]+").unwrap();
    // Glucose neutral mass ≈ 180.06339
    let mz = adduct::neutral_to_mz(180.06339, h_adduct);
    // [M+H]+ should be ~181.07066
    assert!((mz - 181.07066).abs() < 0.001);

    let neutral = adduct::mz_to_neutral(mz, h_adduct);
    assert!((neutral - 180.06339).abs() < 0.0001);
}

#[test]
fn test_no_rules_mode() {
    // Without rules, we should get more results
    let elem_list = elements::parse_element_string("CHNOPS");
    let mut constraints = bounds::default_constraints(&elem_list);
    bounds::compute_bounds(180.06339, 2.0, &mut constraints);

    let with_rules = search::search(&SearchParams {
        target_mass: 180.06339,
        tolerance_ppm: 2.0,
        element_constraints: constraints.clone(),
        adduct: None,
        apply_rules: true,
        strict_rules: false,
        max_results: 1000,
        calc_isotopes: false,
    });

    let without_rules = search::search(&SearchParams {
        target_mass: 180.06339,
        tolerance_ppm: 2.0,
        element_constraints: constraints,
        adduct: None,
        apply_rules: false,
        strict_rules: false,
        max_results: 1000,
        calc_isotopes: false,
    });

    assert!(
        without_rules.len() >= with_rules.len(),
        "Disabling rules should give at least as many results"
    );
}

#[test]
fn test_lipid_with_large_mass_defect_passes_rules() {
    // PC(34:1): exact-mass rounding gives nominal 760 (even) with N=1
    assert_formula_found(759.57781, "C42H82NO8P", 2.0);
}

#[test]
fn test_perfluorinated_compound_found() {
    // PFOA: 15 F drive the partial DBE negative before C is assigned
    let results = search_for(413.97370, 2.0, "CHNOPSF");
    assert!(results.iter().any(|c| c.formula_string == "C8HF15O2"));
}

#[test]
fn test_polychlorinated_compound_found() {
    // Tetrachloroethylene: Cl/C = 2 is the extended-range limit
    let results = search_for(163.87541, 2.0, "CHNOPSCl");
    assert!(results.iter().any(|c| c.formula_string == "C2Cl4"));
}

#[test]
fn test_no_hydrogen_when_not_in_element_list() {
    let results = search_for(180.06339, 5.0, "CNO");
    assert!(results.iter().all(|c| c.composition.iter().all(|&(s, _)| s != "H")));
}

#[test]
fn test_max_results_keeps_best_candidate() {
    let elem_list = elements::parse_element_string("CHNOPS");
    let mut constraints = bounds::default_constraints(&elem_list);
    bounds::compute_bounds(759.57781, 2.0, &mut constraints);
    let params = |max_results| SearchParams {
        target_mass: 759.57781,
        tolerance_ppm: 2.0,
        element_constraints: constraints.clone(),
        adduct: None,
        apply_rules: false,
        strict_rules: false,
        max_results,
        calc_isotopes: false,
    };
    let best_of_all = search::search(&params(100_000))[0].mass_error_da.abs();
    assert_eq!(search::search(&params(1))[0].mass_error_da.abs(), best_of_all);
}

#[test]
fn test_adduct_masses_match_reference() {
    use mfcalc::adduct::find_adduct;
    // Ion masses: atom/molecule mass minus electron (AME2016 via pyteomics)
    for (name, reference) in [("[M+K]+", 38.963_158_10), ("[M+NH4]+", 18.033_825_55), ("[M+Na]+", 22.989_220_70)] {
        let adj = find_adduct(name).unwrap().mass_adjustment;
        assert!((adj - reference).abs() < 1e-6, "{name}: {adj} vs {reference}");
    }
}

#[test]
fn test_isotope_peaks_grouped_by_nucleon_shift() {
    // C34H69ClO2 monoisotopic mass 544.49861 sits at a .5 rounding boundary
    let elem_refs = elements::parse_element_string("CHOCl");
    let pattern = isotope::isotope_distribution(&[("C", 34), ("H", 69), ("Cl", 1), ("O", 2)], &elem_refs);
    let shifts: Vec<f64> = pattern.iter().take(4).map(|p| p.mass - 544.49861).collect();
    assert!(
        shifts.iter().enumerate().all(|(k, s)| (s - k as f64 * 1.0017).abs() < 0.01),
        "shifts {shifts:?}"
    );
}

fn passes_loose_rules(composition: Vec<(&'static str, u32)>) -> bool {
    let candidate = mfcalc::types::FormulaCandidate {
        dbe: mfcalc::types::calculate_dbe(&composition),
        formula_string: mfcalc::types::build_formula_string(&composition),
        monoisotopic_mass: composition
            .iter()
            .map(|&(s, n)| n as f64 * elements::get_element(s).unwrap().monoisotopic_mass)
            .sum(),
        composition,
        mass_error_da: 0.0,
        mass_error_ppm: 0.0,
        isotope_pattern: None,
        adduct_name: None,
        observed_mz: None,
    };
    mfcalc::rules::validate(&candidate, false)
}

// Kind & Fiehn 2007, Table 2 extended range: P/C 0-2, S/C 0-3
#[test]
fn test_rules_accept_sulfur_ratio_within_extended_range() {
    assert!(passes_loose_rules(vec![("C", 2), ("H", 2), ("S", 5)]));
}

#[test]
fn test_rules_reject_phosphorus_ratio_above_extended_range() {
    assert!(!passes_loose_rules(vec![("C", 1), ("H", 3), ("O", 3), ("P", 3)]));
}
