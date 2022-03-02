pub mod common;

use crate::common::{
    assert_eq_results, compute_no_truncate, compute_with_top_n, get_sorted_by_correlation_abs_desc,
    get_tuples_from_result, merge_with_adjustment, ResultTupleSimple, ResultTupleWithoutAdj,
};
use ggca::{adjustment::AdjustmentMethod, correlation::CorrelationMethod};
use lazy_static::lazy_static;

// Datasets's paths
const DF1_PATH: &str = "tests/small_files/mRNA.csv"; // mRNA = 600 rows
const DF2_PATH: &str = "tests/small_files/miRNA.csv"; // miRNA = 299 rows

// Total number of combinations evaluated in a "all vs all" analysis (600 * 299)
const TOTAL_COMBINATIONS_EVALUATED: usize = 179400;

lazy_static! {
    /// Expected correlations for Pearson method, sorted by p-value descending
    static ref EXPECTED_PEARSON: ResultTupleWithoutAdj = vec![
        ("UBE2L6".to_string(), "hsa-miR-552".to_string(), -0.6013218, 5.827716E-09),
        ("ZFHX4".to_string(), "hsa-miR-194*".to_string(), -0.6044824, 4.619085E-09),
        ("MEIS1".to_string(), "hsa-miR-141".to_string(), -0.6047322, 4.534499E-09),
        ("MEIS1".to_string(), "hsa-miR-194".to_string(), -0.6093747, 3.207092E-09),
        ("CDO1".to_string(), "hsa-miR-194".to_string(), -0.6135806, 2.331999E-09),
        ("ZFHX4".to_string(), "hsa-miR-192*".to_string(), -0.614248, 2.216075E-09),
        ("ZFHX4".to_string(), "hsa-miR-192".to_string(), -0.6150258, 2.087884E-09),
        ("MEIS1".to_string(), "hsa-miR-200c".to_string(), -0.6207857, 1.336221E-09),
        ("FCGR2C".to_string(), "hsa-miR-223".to_string(), 0.6264482, 8.539654E-10),
        ("ZFHX4".to_string(), "hsa-miR-194".to_string(), -0.6345078, 4.444592E-10),
        ("MEIS1".to_string(), "hsa-miR-192*".to_string(), -0.6533207, 8.966458E-11),
        ("CDO1".to_string(), "hsa-miR-194*".to_string(), -0.6589037, 5.455469E-11),
    ];

    /// Benjamini-Hochberg adjustment for EXPECTED_PEARSON_BH (respecting the order by p-value in an descending order)
    static ref BH_ADJUSTMENT: Vec<f64> = vec![
        8.414626E-05,
        7.533308E-05,
        7.533308E-05,
        6.392803E-05,
        5.229507E-05,
        5.229507E-05,
        5.229507E-05,
        4.794362E-05,
        3.830035E-05,
        2.657866E-05,
        8.042913E-06,
        8.042913E-06,
    ];

    /// Bonferroni adjustment for EXPECTED_PEARSON_BH (respecting the order by p-value in an descending order)
    static ref BONFERRONI_ADJUSTMENT: Vec<f64> = vec![
        0.001045492,
        0.0008286639,
        0.0008134892,
        0.0005753522,
        0.0004183606,
        0.0003975639,
        0.0003745664,
        0.0002397181,
        0.0001532014,
        7.973597E-05,
        1.608583E-05,
        9.787111E-06,
    ];

    /// Benjamini-Yekutieli adjustment for EXPECTED_PEARSON_BH (respecting the order by p-value in an descending order)
    static ref BY_ADJUSTMENT: Vec<f64> = vec![
        0.0010665195,
        0.000954816,
        0.000954816,
        0.0008102616,
        0.0006628187,
        0.0006628187,
        0.0006628187,
        0.0006076658,
        0.0004854413,
        0.0003368736,
        0.0001019406,
        0.0001019406,
    ];

    // Pearson result with Benjamini-Hochberg adjustment
    static ref EXPECTED_PEARSON_BH: ResultTupleSimple = merge_with_adjustment(&EXPECTED_PEARSON, &BH_ADJUSTMENT);

    // Pearson result with Bonferroni adjustment
    static ref EXPECTED_PEARSON_BONFERRONI: ResultTupleSimple = merge_with_adjustment(&EXPECTED_PEARSON, &BONFERRONI_ADJUSTMENT);

    // Pearson result with Benjamini-Yekutieli adjustment
    static ref EXPECTED_PEARSON_BY: ResultTupleSimple = merge_with_adjustment(&EXPECTED_PEARSON, &BY_ADJUSTMENT);
}

#[test]
/// Tests Pearson correlation with Benjamini-Hochberg adjustment. No threshold set (uses all the rows)
fn pearson_and_bh_all() {
    // Some parameters
    let is_all_vs_all = true;
    let collect_gem_dataset = Some(true); // Better performance. Keep GEM file in memory

    let (result, number_of_elements_evaluated) = compute_no_truncate(
        DF1_PATH.to_string(),
        DF2_PATH.to_string(),
        CorrelationMethod::Pearson,
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
/// Tests Pearson correlation with Benjamini-Hochberg adjustment. Correlation threshold set to 0.6
fn pearson_and_bh_cor_0_6() {
    // Some parameters
    let is_all_vs_all = true;
    let collect_gem_dataset = Some(true); // Better performance. Keep GEM file in memory

    let (result, number_of_elements_evaluated) = compute_no_truncate(
        DF1_PATH.to_string(),
        DF2_PATH.to_string(),
        CorrelationMethod::Pearson,
        0.6,
        2_000_000,
        AdjustmentMethod::BenjaminiHochberg,
        is_all_vs_all,
        collect_gem_dataset,
    );

    assert_eq!(number_of_elements_evaluated, TOTAL_COMBINATIONS_EVALUATED); // The number of evaluated elements mustn't be modified
    assert_eq!(result.len(), EXPECTED_PEARSON.len());

    let collected_as_tuples = get_tuples_from_result(&result);
    assert_eq_results(&collected_as_tuples, &EXPECTED_PEARSON_BH);
}

#[test]
/// Tests Pearson correlation with Benjamini-Hochberg adjustment. Correlation threshold set to 0.6 keeping only top 10 combinations by correlation
fn pearson_and_bh_cor_0_6_top_10() {
    // Some parameters
    let is_all_vs_all = true;
    let keep_top_n = Some(10); // Keep top 10
    let collect_gem_dataset = Some(true); // Better performance. Keep GEM file in memory

    let (result, total_combinations_count, number_of_elements_evaluated) = compute_with_top_n(
        DF1_PATH.to_string(),
        DF2_PATH.to_string(),
        CorrelationMethod::Pearson,
        0.6,
        2_000_000,
        AdjustmentMethod::BenjaminiHochberg,
        is_all_vs_all,
        collect_gem_dataset,
        keep_top_n,
    );

    assert_eq!(number_of_elements_evaluated, TOTAL_COMBINATIONS_EVALUATED); // The number of evaluated elements mustn't be modified
    assert_eq!(total_combinations_count, EXPECTED_PEARSON.len()); // If keep_top_n were None this would be the result
    assert_eq!(result.len(), 10);
}

#[test]
/// Tests Pearson correlation with Benjamini-Hochberg adjustment. No threshold, Keeps top 10 combinations by correlation
fn pearson_and_bh_top_10() {
    // Some parameters
    let is_all_vs_all = true;
    let keep_top_n = Some(10); // Keep top 10
    let collect_gem_dataset = Some(true); // Better performance. Keep GEM file in memory

    let (result, total_combinations_count, number_of_elements_evaluated) = compute_with_top_n(
        DF1_PATH.to_string(),
        DF2_PATH.to_string(),
        CorrelationMethod::Pearson,
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

    // If top N is defined it should be sorted by correlation descending
    let result_sorted = get_sorted_by_correlation_abs_desc(&collected_as_tuples);
    let expected_sorted = get_sorted_by_correlation_abs_desc(&EXPECTED_PEARSON_BH);

    // Check with a zip to prevent issues with floating point equality
    assert_eq_results(&result_sorted, &expected_sorted);
}

#[test]
/// Tests Pearson correlation with Benjamini-Hochberg adjustment. Correlation threshold set to 1 (no results)
fn pearson_and_bh_cor_1() {
    // Some parameters
    let is_all_vs_all = true;
    let collect_gem_dataset = Some(true); // Better performance. Keep GEM file in memory

    let (result, number_of_elements_evaluated) = compute_no_truncate(
        DF1_PATH.to_string(),
        DF2_PATH.to_string(),
        CorrelationMethod::Pearson,
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
/// Tests Pearson correlation with Benjamini-Hochberg adjustment. Only matching genes/GEMs
fn pearson_and_bh_only_matching() {
    // Some parameters
    let is_all_vs_all = false; // Keeps only matching genes/GEMs
    let collect_gem_dataset = Some(true); // Better performance. Keep GEM file in memory

    let (result, number_of_elements_evaluated) = compute_no_truncate(
        DF1_PATH.to_string(),
        DF2_PATH.to_string(),
        CorrelationMethod::Pearson,
        1.0,
        2_000_000,
        AdjustmentMethod::BenjaminiHochberg,
        is_all_vs_all,
        collect_gem_dataset,
    );

    // As there is no matching genes/GEMs no combinations are evaluated
    assert_eq!(number_of_elements_evaluated, 0);
    assert_eq!(result.len(), 0);
}

#[test]
/// Tests Pearson correlation with Bonferroni adjustment. Correlation threshold set to 0.6
fn pearson_and_bonferroni_cor_0_6() {
    // Some parameters
    let is_all_vs_all = true;
    let collect_gem_dataset = Some(true); // Better performance. Keep GEM file in memory

    let (result, number_of_elements_evaluated) = compute_no_truncate(
        DF1_PATH.to_string(),
        DF2_PATH.to_string(),
        CorrelationMethod::Pearson,
        0.6,
        2_000_000,
        AdjustmentMethod::Bonferroni,
        is_all_vs_all,
        collect_gem_dataset,
    );

    // The adjustment method should not modify the number of resulting combinations
    assert_eq!(number_of_elements_evaluated, TOTAL_COMBINATIONS_EVALUATED);
    assert_eq!(result.len(), EXPECTED_PEARSON.len());

    let collected_as_tuples = get_tuples_from_result(&result);

    // Bonferroni does not sort, so sorts both vec of tuples by correlation descending and compares
    let result_sorted = get_sorted_by_correlation_abs_desc(&collected_as_tuples);
    let expected_sorted = get_sorted_by_correlation_abs_desc(&EXPECTED_PEARSON_BONFERRONI);

    assert_eq_results(&result_sorted, &expected_sorted);
}

#[test]
/// Tests Pearson correlation with Benjamini-Yekutieli adjustment. Correlation threshold set to 0.6
fn pearson_and_by_cor_0_6() {
    // Some parameters
    let is_all_vs_all = true;
    let collect_gem_dataset = Some(true); // Better performance. Keep GEM file in memory

    let (result, number_of_elements_evaluated) = compute_no_truncate(
        DF1_PATH.to_string(),
        DF2_PATH.to_string(),
        CorrelationMethod::Pearson,
        0.6,
        2_000_000,
        AdjustmentMethod::BenjaminiYekutieli,
        is_all_vs_all,
        collect_gem_dataset,
    );

    // The adjustment method should not modify the number of resulting combinations
    assert_eq!(number_of_elements_evaluated, TOTAL_COMBINATIONS_EVALUATED);
    assert_eq!(result.len(), EXPECTED_PEARSON.len());

    let collected_as_tuples = get_tuples_from_result(&result);
    assert_eq_results(&collected_as_tuples, &EXPECTED_PEARSON_BY);
}
