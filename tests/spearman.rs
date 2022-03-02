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
    /// Expected correlations for Spearman method, sorted by p-value descending
    static ref EXPECTED_SPEARMAN: ResultTupleWithoutAdj = vec![
        ("EFS".to_string(), "hsa-miR-125b".to_string(), 0.6015163, 5.745374E-09),
        ("THBS4".to_string(), "hsa-miR-145".to_string(), 0.6015794, 5.718886E-09),
        ("MEIS1".to_string(), "hsa-miR-199a-5p".to_string(), 0.602729, 5.256441E-09),
        ("SLC22A17".to_string(), "hsa-miR-125b".to_string(), 0.6060311, 4.118024E-09),
        ("SLC6A4".to_string(), "hsa-miR-31".to_string(), -0.6060886, 4.100453E-09),
        ("RAB26".to_string(), "hsa-miR-375".to_string(), 0.6106425, 2.91482E-09),
        ("SLAMF8".to_string(), "hsa-miR-146b-3p".to_string(), 0.6191673, 1.516138E-09),
        ("THBS4".to_string(), "hsa-miR-145*".to_string(), 0.6224393, 1.173549E-09),
        ("OLFML1".to_string(), "hsa-miR-199a-5p".to_string(), 0.625226, 9.413221E-10),
        ("CDO1".to_string(), "hsa-miR-199a-5p".to_string(), 0.6303515, 6.238948E-10),
        ("CCDC8".to_string(), "hsa-miR-199a-5p".to_string(), 0.6304362, 6.196302E-10),
        ("MEIS1".to_string(), "hsa-miR-125b".to_string(), 0.6357503, 4.012115E-10),
        ("PTPRC".to_string(), "hsa-miR-142-3p".to_string(), 0.638176, 3.280958E-10),
        ("THY1".to_string(), "hsa-miR-199a-5p".to_string(), 0.650012, 1.197849E-10),
        ("GREM2".to_string(), "hsa-miR-1".to_string(), 0.6914202, 2.421017E-12),
    ];

    /// Benjamini-Hochberg adjustment for EXPECTED_SPEARMAN_BH (respecting the order by p-value in an descending order)
    static ref BH_ADJUSTMENT: Vec<f64> = vec![
        6.871468E-05,
        6.871468E-05,
        6.871468E-05,
        6.156446E-05,
        6.156446E-05,
        5.229188E-05,
        3.022168E-05,
        2.631683E-05,
        2.412474E-05,
        1.865445E-05,
        1.865445E-05,
        1.799434E-05,
        1.799434E-05,
        1.074471E-05,
        4.343305E-07,
    ];

    /// Bonferroni adjustment for EXPECTED_SPEARMAN_BH (respecting the order by p-value in an descending order)
    static ref BONFERRONI_ADJUSTMENT: Vec<f64> = vec![
        0.00103072,
        0.001025968,
        0.0009430055,
        0.0007387735,
        0.0007356214,
        0.0005229188,
        0.0002719952,
        0.0002105346,
        0.0001688732,
        0.0001119267,
        0.0001111617,
        7.197735E-05,
        5.886039E-05,
        2.148941E-05,
        4.343305E-07
    ];

    /// Benjamini-Yekutieli adjustment for EXPECTED_SPEARMAN_BH (respecting the order by p-value in an descending order)
    static ref BY_ADJUSTMENT: Vec<f64> = vec![
        0.0008709305,
        0.0008709305,
        0.0008709305,
        0.0007803044,
        0.0007803044,
        0.0006627782,
        0.0003830475,
        0.000333555,
        0.0003057712,
        0.0002364376,
        0.0002364376,
        0.0002280709,
        0.0002280709,
        0.0001361848,
        5.504962E-06,
    ];

    // Spearman result with Benjamini-Hochberg adjustment
    static ref EXPECTED_SPEARMAN_BH: ResultTupleSimple = merge_with_adjustment(&EXPECTED_SPEARMAN, &BH_ADJUSTMENT);

    // Spearman result with Bonferroni adjustment
    static ref EXPECTED_SPEARMAN_BONFERRONI: ResultTupleSimple = merge_with_adjustment(&EXPECTED_SPEARMAN, &BONFERRONI_ADJUSTMENT);

    // Spearman result with Benjamini-Yekutieli adjustment
    static ref EXPECTED_SPEARMAN_BY: ResultTupleSimple = merge_with_adjustment(&EXPECTED_SPEARMAN, &BY_ADJUSTMENT);
}

#[test]
/// Tests Spearman correlation with Benjamini-Hochberg adjustment. No threshold set (uses all the rows)
fn spearman_and_bh_all() {
    // Some parameters
    let is_all_vs_all = true;
    let collect_gem_dataset = Some(true); // Better performance. Keep GEM file in memory

    let (result, number_of_elements_evaluated) = compute_no_truncate(
        DF1_PATH.to_string(),
        DF2_PATH.to_string(),
        CorrelationMethod::Spearman,
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
/// Tests Spearman correlation with Benjamini-Hochberg adjustment. Correlation threshold set to 0.6
fn spearman_and_bh_cor_0_6() {
    // Some parameters
    let is_all_vs_all = true;
    let collect_gem_dataset = Some(true); // Better performance. Keep GEM file in memory

    let (result, number_of_elements_evaluated) = compute_no_truncate(
        DF1_PATH.to_string(),
        DF2_PATH.to_string(),
        CorrelationMethod::Spearman,
        0.6,
        2_000_000,
        AdjustmentMethod::BenjaminiHochberg,
        is_all_vs_all,
        collect_gem_dataset,
    );

    assert_eq!(number_of_elements_evaluated, TOTAL_COMBINATIONS_EVALUATED); // The number of evaluated elements mustn't be modified
    assert_eq!(result.len(), EXPECTED_SPEARMAN.len());

    let collected_as_tuples = get_tuples_from_result(&result);

    assert_eq_results(&collected_as_tuples, &EXPECTED_SPEARMAN_BH);
}

#[test]
/// Tests Spearman correlation with Benjamini-Hochberg adjustment. Correlation threshold set to 0.6 keeping only top 10 combinations by correlation
fn spearman_and_bh_cor_0_6_top_10() {
    // Some parameters
    let is_all_vs_all = true;
    let keep_top_n = Some(10); // Keep top 10
    let collect_gem_dataset = Some(true); // Better performance. Keep GEM file in memory

    let (result, total_combinations_count, number_of_elements_evaluated) = compute_with_top_n(
        DF1_PATH.to_string(),
        DF2_PATH.to_string(),
        CorrelationMethod::Spearman,
        0.6,
        2_000_000,
        AdjustmentMethod::BenjaminiHochberg,
        is_all_vs_all,
        collect_gem_dataset,
        keep_top_n,
    );

    assert_eq!(number_of_elements_evaluated, TOTAL_COMBINATIONS_EVALUATED); // The number of evaluated elements mustn't be modified
    assert_eq!(total_combinations_count, EXPECTED_SPEARMAN.len()); // If keep_top_n were None this would be the result
    assert_eq!(result.len(), 10);
}

#[test]
/// Tests Spearman correlation with Benjamini-Hochberg adjustment. No threshold, Keeps top 10 combinations by correlation
fn spearman_and_bh_top_10() {
    // Some parameters
    let is_all_vs_all = true;
    let keep_top_n = Some(10); // Keep top 10
    let collect_gem_dataset = Some(true); // Better performance. Keep GEM file in memory

    let (result, total_combinations_count, number_of_elements_evaluated) = compute_with_top_n(
        DF1_PATH.to_string(),
        DF2_PATH.to_string(),
        CorrelationMethod::Spearman,
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
    let sorted = get_sorted_by_correlation_abs_desc(&EXPECTED_SPEARMAN_BH)[0..10].to_vec();

    // Check with a zip to prevent issues with floating point equality
    assert_eq_results(&collected_as_tuples, &sorted);
}

#[test]
/// Tests Spearman correlation with Benjamini-Hochberg adjustment. Correlation threshold set to 1 (no results)
fn spearman_and_bh_cor_1() {
    // Some parameters
    let is_all_vs_all = true;
    let collect_gem_dataset = Some(true); // Better performance. Keep GEM file in memory

    let (result, number_of_elements_evaluated) = compute_no_truncate(
        DF1_PATH.to_string(),
        DF2_PATH.to_string(),
        CorrelationMethod::Spearman,
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
/// Tests Spearman correlation with Benjamini-Hochberg adjustment. Only matching genes/GEMs
fn spearman_and_bh_only_matching() {
    // Some parameters
    let is_all_vs_all = false; // Keeps only matching genes/GEMs
    let collect_gem_dataset = Some(true); // Better performance. Keep GEM file in memory

    let (result, number_of_elements_evaluated) = compute_no_truncate(
        DF1_PATH.to_string(),
        DF2_PATH.to_string(),
        CorrelationMethod::Spearman,
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
/// Tests Spearman correlation with Bonferroni adjustment. Correlation threshold set to 0.6
fn spearman_and_bonferroni_cor_0_6() {
    // Some parameters
    let is_all_vs_all = true;
    let collect_gem_dataset = Some(true); // Better performance. Keep GEM file in memory

    let (result, number_of_elements_evaluated) = compute_no_truncate(
        DF1_PATH.to_string(),
        DF2_PATH.to_string(),
        CorrelationMethod::Spearman,
        0.6,
        2_000_000,
        AdjustmentMethod::Bonferroni,
        is_all_vs_all,
        collect_gem_dataset,
    );

    // The adjustment method should not modify the number of resulting combinations
    assert_eq!(number_of_elements_evaluated, TOTAL_COMBINATIONS_EVALUATED);
    assert_eq!(result.len(), EXPECTED_SPEARMAN.len());

    let collected_as_tuples = get_tuples_from_result(&result);

    // Bonferroni does not sort, so sorts both vec of tuples by correlation descending and compares
    let result_sorted = get_sorted_by_correlation_abs_desc(&collected_as_tuples);
    let expected_sorted = get_sorted_by_correlation_abs_desc(&EXPECTED_SPEARMAN_BONFERRONI);

    assert_eq_results(&result_sorted, &expected_sorted);
}

#[test]
/// Tests Spearman correlation with Benjamini-Yekutieli adjustment. Correlation threshold set to 0.6
fn spearman_and_by_cor_0_6() {
    // Some parameters
    let is_all_vs_all = true;
    let collect_gem_dataset = Some(true); // Better performance. Keep GEM file in memory

    let (result, number_of_elements_evaluated) = compute_no_truncate(
        DF1_PATH.to_string(),
        DF2_PATH.to_string(),
        CorrelationMethod::Spearman,
        0.6,
        2_000_000,
        AdjustmentMethod::BenjaminiYekutieli,
        is_all_vs_all,
        collect_gem_dataset,
    );

    // The adjustment method should not modify the number of resulting combinations
    assert_eq!(number_of_elements_evaluated, TOTAL_COMBINATIONS_EVALUATED);
    assert_eq!(result.len(), EXPECTED_SPEARMAN.len());

    let collected_as_tuples = get_tuples_from_result(&result);
    assert_eq_results(&collected_as_tuples, &EXPECTED_SPEARMAN_BY);
}
