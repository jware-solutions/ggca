pub mod common;

use crate::common::{
    assert_eq_results, compute_no_truncate, compute_with_top_n,
    get_sorted_by_correlation_abs_desc_gene_gem, get_tuples_from_result, merge_with_adjustment,
    ResultTupleSimple, ResultTupleWithoutAdj,
};
use ggca::{adjustment::AdjustmentMethod, correlation::CorrelationMethod};
use lazy_static::lazy_static;

// Datasets's paths
const DF1_PATH: &str = "tests/small_files/mRNA_for_CNA.csv"; // mRNA = 200 rows
const DF2_PATH: &str = "tests/small_files/CNA.csv"; // miRNA = 300 rows

// Total number of combinations evaluated in a "all vs all" analysis (200 * 300)
const TOTAL_COMBINATIONS_EVALUATED: usize = 60_000;

lazy_static! {
    /// Expected correlations for Kendall method, sorted by p-value descending
    static ref EXPECTED_KENDALL: ResultTupleWithoutAdj = vec![
        ("TM9SF4".to_string(), "OACYLP".to_string(), -0.45247298818523, 8.046158713408351e-17),
        ("TM9SF4".to_string(), "C18orf20".to_string(), -0.456175756128336, 4.5131298922977604e-17),
        ("TM9SF4".to_string(), "SERPINB2".to_string(), -0.456175756128336, 4.5131298922977604e-17),
        ("PI4KB".to_string(), "THEM4".to_string(), 0.464223633461795, 3.90954707333821e-18),
        ("MKKS".to_string(), "LOC643406".to_string(), 0.468485042732911, 1.53033062609261e-19),
        ("MKKS".to_string(), "PRNP".to_string(), 0.481834936382762, 1.32644852388271e-20),
        ("C20orf30".to_string(), "LOC643406".to_string(), 0.48617299400547, 6.35305382533189e-21),
        ("C20orf30".to_string(), "PRNP".to_string(), 0.490074784643379, 2.9295217919757602e-21),
        ("TM9SF4".to_string(), "NCOA5".to_string(), 0.523127128041391, 2.56368643663311e-23),
        ("C18orf32".to_string(), "NFATC1".to_string(), 0.544185780299279, 1.56555442332436e-23),
        ("C18orf32".to_string(), "OACYLP".to_string(), 0.549285980386605, 4.84760003748e-24),
        ("C18orf32".to_string(), "C18orf20".to_string(), 0.549802069921003, 4.3990071793089995e-24),
        ("C18orf32".to_string(), "SERPINB2".to_string(), 0.549802069921003, 4.3990071793089995e-24),
        ("TM9SF4".to_string(), "SLC13A3".to_string(), 0.534637516627216, 3.64348178546081e-24),
        ("TM9SF4".to_string(), "SLPI".to_string(), 0.532626379173813, 3.37174107305418e-24),
        ("TM9SF4".to_string(), "hsa-mir-1257".to_string(), 0.541463497707136, 1.85264855161345e-24),
        ("TM9SF4".to_string(), "FAM65C".to_string(), 0.54883570813472, 2.81631222187641e-25),
        ("TM9SF4".to_string(), "KIAA1755".to_string(), 0.557344951078564, 3.18464201343994e-26),
        ("TM9SF4".to_string(), "DSN1".to_string(), 0.561757583597623, 9.39749876440102e-27),
        ("TM9SF4".to_string(), "BLCAP".to_string(), 0.567937387581985, 4.07679726357826e-27),
        ("TM9SF4".to_string(), "UQCC".to_string(), 0.570884950344156, 1.4143490488538799e-27),
        ("TM9SF4".to_string(), "RBM12".to_string(), 0.574917972860065, 6.06726053852395e-28)
    ];

    /// Benjamini-Hochberg adjustment for EXPECTED_KENDALL_BH (respecting the order by p-value in an descending order)
    static ref BH_ADJUSTMENT: Vec<f64> = vec![
       2.19440692183864e-13,
        1.28946568351365e-13,
        1.28946568351365e-13,
        1.23459381263312e-14,
        5.10110208697537e-16,
        4.6815830254683904e-17,
        2.38239518449946e-17,
        1.1718087167903e-17,
        1.0987227585570501e-19,
        7.22563579995859e-20,
        2.4238000187400002e-20,
        2.39945846144127e-20,
        2.39945846144127e-20,
        2.39945846144127e-20,
        2.39945846144127e-20,
        1.58798447281153e-20,
        2.81631222187641e-21,
        3.8215704161279297e-22,
        1.40962481466015e-22,
        8.15359452715652e-23,
        4.24304714656164e-23,
        3.6403563231143704e-23
    ];

    /// Bonferroni adjustment for EXPECTED_KENDALL_BH (respecting the order by p-value in an descending order)
    static ref BONFERRONI_ADJUSTMENT: Vec<f64> = vec![
        4.8276952280450095e-12,
        2.7078779353786598e-12,
        2.7078779353786598e-12,
        2.34572824400293e-13,
        9.181983756555661e-15,
        7.95869114329626e-16,
        3.81183229519913e-16,
        1.75771307518546e-16,
        1.5382118619798698e-18,
        9.39332653994616e-19,
        2.908560022488e-19,
        2.6394043075853996e-19,
        2.6394043075853996e-19,
        2.18608907127649e-19,
        2.02304464383251e-19,
        1.11158913096807e-19,
        1.68978733312585e-20,
        1.91078520806396e-21,
        5.638499258640609e-22,
        2.44607835814696e-22,
        8.48609429312328e-23,
        3.6403563231143704e-23
    ];

    /// Benjamini-Yekutieli adjustment for EXPECTED_KENDALL_BH (respecting the order by p-value in an descending order)
    static ref BY_ADJUSTMENT: Vec<f64> = vec![
        2.54097483834256e-12,
        1.4931140729218299e-12,
        1.4931140729218299e-12,
        1.42957615666181e-13,
        5.906731300300831e-15,
        5.42095659330112e-16,
        2.75865253547844e-16,
        1.35687526095654e-16,
        1.27224666310685e-18,
        8.366797687339781e-19,
        2.80659653389729e-19,
        2.77841065642553e-19,
        2.77841065642553e-19,
        2.77841065642553e-19,
        2.77841065642553e-19,
        1.8387786462648902e-19,
        3.2610991250012e-20,
        4.42512014234769e-21,
        1.63225022210266e-21,
        9.44131114852349e-22,
        4.913161697594701e-22,
        4.21528647562068e-22
    ];

    // Kendall result with Benjamini-Hochberg adjustment
    static ref EXPECTED_KENDALL_BH: ResultTupleSimple = merge_with_adjustment(&EXPECTED_KENDALL, &BH_ADJUSTMENT);

    // Kendall result with Bonferroni adjustment
    static ref EXPECTED_KENDALL_BONFERRONI: ResultTupleSimple = merge_with_adjustment(&EXPECTED_KENDALL, &BONFERRONI_ADJUSTMENT);

    // Kendall result with Benjamini-Yekutieli adjustment
    static ref EXPECTED_KENDALL_BY: ResultTupleSimple = merge_with_adjustment(&EXPECTED_KENDALL, &BY_ADJUSTMENT);
}

#[test]
/// Tests Kendall correlation with Benjamini-Hochberg adjustment. No threshold set (uses all the rows)
fn kendall_and_bh_all() {
    // Some parameters
    let is_all_vs_all = true;
    let collect_gem_dataset = Some(true); // Better performance. Keep GEM file in memory

    let (result, number_of_elements_evaluated) = compute_no_truncate(
        DF1_PATH.to_string(),
        DF2_PATH.to_string(),
        CorrelationMethod::Kendall,
        0.0, // No threshold
        2_000_000,
        AdjustmentMethod::BenjaminiHochberg,
        is_all_vs_all,
        collect_gem_dataset,
    );

    // No threshold, no rows filtered
    assert_eq!(number_of_elements_evaluated, TOTAL_COMBINATIONS_EVALUATED);
    assert_eq!(result.len(), number_of_elements_evaluated); // No element should be filtered
}

#[test]
/// Tests Kendall correlation with Benjamini-Hochberg adjustment. Correlation threshold set to 0.45
fn kendall_and_bh_cor_0_45() {
    // Some parameters
    let is_all_vs_all = true;
    let collect_gem_dataset = Some(true); // Better performance. Keep GEM file in memory

    let (result, number_of_elements_evaluated) = compute_no_truncate(
        DF1_PATH.to_string(),
        DF2_PATH.to_string(),
        CorrelationMethod::Kendall,
        0.45,
        2_000_000,
        AdjustmentMethod::BenjaminiHochberg,
        is_all_vs_all,
        collect_gem_dataset,
    );

    assert_eq!(number_of_elements_evaluated, TOTAL_COMBINATIONS_EVALUATED); // The number of evaluated elements mustn't be modified
    assert_eq!(result.len(), 22);

    let collected_as_tuples = get_tuples_from_result(&result);

    // NOTE: with Kendall we had some repeated correlation values that produced different orders so the tests failed.
    // Sorting by cor (desc), gene (asc) and GEM (asc) we ensure that the order of both vectors coincide
    let collected_as_tuples = get_sorted_by_correlation_abs_desc_gene_gem(&collected_as_tuples);
    let sorted = get_sorted_by_correlation_abs_desc_gene_gem(&EXPECTED_KENDALL_BH);

    assert_eq_results(&collected_as_tuples, &sorted);
}

#[test]
/// Tests Kendall correlation with Benjamini-Hochberg adjustment. No threshold, Keeps top 10 combinations by correlation
fn kendall_and_bh_top_10() {
    // Some parameters
    let is_all_vs_all = true;
    let keep_top_n = Some(10); // Keep top 10
    let collect_gem_dataset = Some(true); // Better performance. Keep GEM file in memory

    let (result, total_combinations_count, number_of_elements_evaluated) = compute_with_top_n(
        DF1_PATH.to_string(),
        DF2_PATH.to_string(),
        CorrelationMethod::Kendall,
        0.0, // No threshold
        2_000_000,
        AdjustmentMethod::BenjaminiHochberg,
        is_all_vs_all,
        collect_gem_dataset,
        keep_top_n,
    );

    // No threshold, no rows filtered
    assert_eq!(number_of_elements_evaluated, TOTAL_COMBINATIONS_EVALUATED); // The number of evaluated elements mustn't be modified
    assert_eq!(total_combinations_count, TOTAL_COMBINATIONS_EVALUATED);
    assert_eq!(result.len(), 10); // Keeps only 10 elements

    let collected_as_tuples = get_tuples_from_result(&result);

    // If top N is defined it should be sorted by correlation descending.
    // NOTE: with Kendall we had some repeated correlation values that produced different orders so the tests failed.
    // Sorting by cor (desc), gene (asc) and GEM (asc) we ensure that the order of both vectors coincide
    let collected_as_tuples = get_sorted_by_correlation_abs_desc_gene_gem(&collected_as_tuples);
    let sorted = get_sorted_by_correlation_abs_desc_gene_gem(&EXPECTED_KENDALL_BH)[0..10].to_vec();

    // Check with a zip to prevent issues with floating point equality
    assert_eq_results(&collected_as_tuples, &sorted);
}

#[test]
/// Tests Kendall correlation with Benjamini-Hochberg adjustment. Correlation threshold set to 0.45 keeping only top 10 combinations by correlation
fn kendall_and_bh_cor_0_45_top_10() {
    // Some parameters
    let is_all_vs_all = true;
    let keep_top_n = Some(10); // Keep top 10
    let collect_gem_dataset = Some(true); // Better performance. Keep GEM file in memory

    let (result, total_combinations_count, number_of_elements_evaluated) = compute_with_top_n(
        DF1_PATH.to_string(),
        DF2_PATH.to_string(),
        CorrelationMethod::Kendall,
        0.45,
        2_000_000,
        AdjustmentMethod::BenjaminiHochberg,
        is_all_vs_all,
        collect_gem_dataset,
        keep_top_n,
    );

    assert_eq!(number_of_elements_evaluated, TOTAL_COMBINATIONS_EVALUATED); // The number of evaluated elements mustn't be modified
    assert_eq!(total_combinations_count, EXPECTED_KENDALL_BH.len()); // If keep_top_n were None this would be the result
    assert_eq!(result.len(), 10);
}

#[test]
/// Tests Kendall correlation with Benjamini-Hochberg adjustment. Correlation threshold set to 1 (no results)
fn kendall_and_bh_cor_1() {
    // Some parameters
    let is_all_vs_all = true;
    let collect_gem_dataset = Some(true); // Better performance. Keep GEM file in memory

    let (result, number_of_elements_evaluated) = compute_no_truncate(
        DF1_PATH.to_string(),
        DF2_PATH.to_string(),
        CorrelationMethod::Kendall,
        1.0,
        2_000_000,
        AdjustmentMethod::BenjaminiHochberg,
        is_all_vs_all,
        collect_gem_dataset,
    );

    assert_eq!(number_of_elements_evaluated, TOTAL_COMBINATIONS_EVALUATED); // The number of evaluated elements mustn't be modified
    assert_eq!(result.len(), 0);
}

#[test]
/// Tests Kendall correlation with Benjamini-Hochberg adjustment. Only matching genes/GEMs
fn kendall_and_bh_only_matching() {
    // Some parameters
    let is_all_vs_all = false; // Keeps only matching genes/GEMs
    let collect_gem_dataset = Some(true); // Better performance. Keep GEM file in memory

    let (result, number_of_elements_evaluated) = compute_no_truncate(
        DF1_PATH.to_string(),
        DF2_PATH.to_string(),
        CorrelationMethod::Kendall,
        1.0,
        2_000_000,
        AdjustmentMethod::BenjaminiHochberg,
        is_all_vs_all,
        collect_gem_dataset,
    );

    // There are three matching Gene/GEM
    assert_eq!(number_of_elements_evaluated, 3);
    assert_eq!(result.len(), 0);
}

#[test]
/// Tests Kendall correlation with Bonferroni adjustment. Correlation threshold set to 0.45
fn kendall_and_bonferroni_cor_0_45() {
    // Some parameters
    let is_all_vs_all = true;
    let collect_gem_dataset = Some(true); // Better performance. Keep GEM file in memory

    let (result, number_of_elements_evaluated) = compute_no_truncate(
        DF1_PATH.to_string(),
        DF2_PATH.to_string(),
        CorrelationMethod::Kendall,
        0.45,
        2_000_000,
        AdjustmentMethod::Bonferroni,
        is_all_vs_all,
        collect_gem_dataset,
    );

    // The adjustment method should not modify the number of resulting combinations
    assert_eq!(number_of_elements_evaluated, TOTAL_COMBINATIONS_EVALUATED);
    assert_eq!(result.len(), 22);

    let collected_as_tuples = get_tuples_from_result(&result);

    // Bonferroni does not sort, so sorts both vec of tuples by correlation descending and compares.
    // NOTE: with Kendall we had some repeated correlation values that produced different orders so the tests failed.
    // Sorting by cor (desc), gene (asc) and GEM (asc) we ensure that the order of both vectors coincide
    let result_sorted = get_sorted_by_correlation_abs_desc_gene_gem(&collected_as_tuples);
    let expected_sorted = get_sorted_by_correlation_abs_desc_gene_gem(&EXPECTED_KENDALL_BONFERRONI);

    assert_eq_results(&result_sorted, &expected_sorted);
}

#[test]
/// Tests Kendall correlation with Benjamini-Yekutieli adjustment. Correlation threshold set to 0.45
fn kendall_and_by_cor_0_45() {
    // Some parameters
    let is_all_vs_all = true;
    let collect_gem_dataset = Some(true); // Better performance. Keep GEM file in memory

    let (result, number_of_elements_evaluated) = compute_no_truncate(
        DF1_PATH.to_string(),
        DF2_PATH.to_string(),
        CorrelationMethod::Kendall,
        0.45,
        2_000_000,
        AdjustmentMethod::BenjaminiYekutieli,
        is_all_vs_all,
        collect_gem_dataset,
    );

    // The adjustment method should not modify the number of resulting combinations
    assert_eq!(number_of_elements_evaluated, TOTAL_COMBINATIONS_EVALUATED);
    assert_eq!(result.len(), 22);

    let collected_as_tuples = get_tuples_from_result(&result);

    // NOTE: with Kendall we had some repeated correlation values that produced different orders so the tests failed.
    // Sorting by cor (desc), gene (asc) and GEM (asc) we ensure that the order of both vectors coincide
    let collected_as_tuples = get_sorted_by_correlation_abs_desc_gene_gem(&collected_as_tuples);
    let sorted = get_sorted_by_correlation_abs_desc_gene_gem(&EXPECTED_KENDALL_BY);

    assert_eq_results(&collected_as_tuples, &sorted);
}
