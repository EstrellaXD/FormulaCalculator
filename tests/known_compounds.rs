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
