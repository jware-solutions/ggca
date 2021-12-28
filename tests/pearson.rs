mod common;

use crate::common::{
    assert_eq_results, get_sorted_by_correlation, get_tuples_from_result, merge_with_adjustment,
    ResultTupleSimple, ResultTupleWithoutAdj,
};
use ggca::{adjustment::AdjustmentMethod, analysis::Analysis, correlation::CorrelationMethod};
use itertools::Itertools;
use lazy_static::lazy_static;

// Datasets's paths
const DF1_PATH: &str = "tests/small_files/mRNA.csv"; // mRNA = 760 rows
const DF2_PATH: &str = "tests/small_files/miRNA.csv"; // miRNA = 601 rows

// Total number of combinations evaluated in a "all vs all" analysis (760 * 601)
const TOTAL_COMBINATIONS_EVALUATED: usize = 456760;

lazy_static! {
    static ref ANALYSIS: Analysis = {
        let gem_contains_cpg = false;
        Analysis::new_from_files(DF1_PATH.to_string(), DF2_PATH.to_string(), gem_contains_cpg)
    };

    /// Expected correlations for Pearson method, sorted by p-value ascending
    static ref EXPECTED_PEARSON: ResultTupleWithoutAdj = vec![
        ("C6orf155".to_string(), "hsa-mir-30a".to_string(), 0.8581767, 8.379477E-219),
        ("BTLA".to_string(), "hsa-mir-150".to_string(), 0.8479504, 2.192392E-208),
        ("APOBEC4".to_string(), "hsa-mir-449a".to_string(), 0.8009502, 8.516449E-169),
        ("APOBEC4".to_string(), "hsa-mir-449b".to_string(), 0.7973792, 3.147156E-166),
        ("C5orf20".to_string(), "hsa-mir-150".to_string(), 0.7920717, 1.660754E-162),
        ("ATP2B2".to_string(), "hsa-mir-885".to_string(), 0.7912955, 5.696732E-162),
        ("C21orf34".to_string(), "hsa-let-7c".to_string(), 0.7882625, 6.691326E-160),
        ("C16orf54".to_string(), "hsa-mir-150".to_string(), 0.7851163, 8.648716E-158),
        ("C21orf34".to_string(), "hsa-mir-99a".to_string(), 0.774241, 9.323279E-151),
        ("C1orf158".to_string(), "hsa-mir-449a".to_string(), 0.7694639, 8.61831E-148),
        ("C1orf158".to_string(), "hsa-mir-449b".to_string(), 0.759847, 4.940956E-142),
        ("APOBEC4".to_string(), "hsa-mir-449c".to_string(), 0.7531277, 3.617836E-138),
        ("C1orf158".to_string(), "hsa-mir-449c".to_string(), 0.7528368, 5.283894E-138),
        ("C20orf85".to_string(), "hsa-mir-449a".to_string(), 0.751877, 1.836797E-137),
        ("AKAP14".to_string(), "hsa-mir-449c".to_string(), 0.7501377, 1.731085E-136),
        ("C21orf34".to_string(), "hsa-mir-125b-2".to_string(), 0.7411945, 1.328645E-131),
        ("ACAP1".to_string(), "hsa-mir-150".to_string(), 0.7387877, 2.533343E-130),
        ("AKAP14".to_string(), "hsa-mir-449b".to_string(), 0.7381176, 5.722785E-130),
        ("AIM2".to_string(), "hsa-mir-155".to_string(), 0.7357108, 1.04683E-128),
        ("BANK1".to_string(), "hsa-mir-150".to_string(), 0.7353945, 1.530239E-128),
        ("AKAP14".to_string(), "hsa-mir-449a".to_string(), 0.7319602, 9.111215E-127),
        ("AMICA1".to_string(), "hsa-mir-150".to_string(), 0.7306179, 4.423763E-126),
        ("ART3".to_string(), "hsa-mir-934".to_string(), 0.7223071, 6.370956E-122),
        ("CA12".to_string(), "hsa-mir-934".to_string(), -0.7220895, 8.146964E-122),
        ("C1orf194".to_string(), "hsa-mir-449b".to_string(), 0.7216007, 1.414585E-121),
        ("BLK".to_string(), "hsa-mir-150".to_string(), 0.7202692, 6.320095E-121),
        ("C1orf194".to_string(), "hsa-mir-449a".to_string(), 0.7199512, 9.025014E-121),
        ("C20orf85".to_string(), "hsa-mir-449b".to_string(), 0.7195478, 1.417124E-120),
        ("AGR3".to_string(), "hsa-mir-190b".to_string(), 0.7180282, 7.698685E-120),
        ("AGR2".to_string(), "hsa-mir-934".to_string(), -0.712854, 2.255013E-117),
        ("AGR2".to_string(), "hsa-mir-577".to_string(), -0.7014995, 3.791311E-112),
    ];

    /// Benjamini-Hochberg adjustment for EXPECTED_PEARSON_BH (respecting the order by p-value in an ascending order)
    static ref BH_ADJUSTMENT: Vec<f64> = vec![
        3.82741E-213,
        5.006985E-203,
        1.296658E-163,
        3.593738E-161,
        1.517132E-157,
        4.336732E-157,
        4.366186E-155,
        4.937985E-153,
        4.731668E-146,
        3.936499E-143,
        2.051665E-137,
        1.377069E-133,
        1.856516E-133,
        5.99268E-133,
        5.271269E-132,
        3.792948E-127,
        6.806646E-126,
        1.452189E-125,
        2.516578E-124,
        3.49476E-124,
        1.981733E-122,
        9.184536E-122,
        1.265216E-117,
        1.550503E-117,
        2.584503E-117,
        1.110295E-116,
        1.526765E-116,
        2.311734E-116,
        1.21257E-115,
        3.433333E-113,
        5.586191E-108,
    ];

    /// Bonferroni adjustment for EXPECTED_PEARSON_BH (respecting the order by p-value in an ascending order)
    static ref BONFERRONI_ADJUSTMENT: Vec<f64> = vec![
        3.82741E-213,
        1.001397E-202,
        3.889973E-163,
        1.437495E-160,
        7.585658E-157,
        2.602039E-156,
        3.05633E-154,
        3.950388E-152,
        4.258501E-145,
        3.936499E-142,
        2.256831E-136,
        1.652483E-132,
        2.413471E-132,
        8.389752E-132,
        7.906904E-131,
        6.068717E-126,
        1.15713E-124,
        2.613939E-124,
        4.781499E-123,
        6.989521E-123,
        4.161638E-121,
        2.020598E-120,
        2.909998E-116,
        3.721207E-116,
        6.461256E-116,
        2.886766E-115,
        4.122266E-115,
        6.472854E-115,
        3.516452E-114,
        1.03E-111,
        1.731719E-106,
    ];

    /// Benjamini-Yekutieli adjustment for EXPECTED_PEARSON_BH (respecting the order by p-value in an ascending order)
    static ref BY_ADJUSTMENT: Vec<f64> = vec![
        5.208772E-212,
6.814071E-202,
1.764638E-162,
4.890765E-160,
2.064684E-156,
5.901915E-156,
5.941999E-154,
6.720167E-152,
6.439388E-145,
5.357233E-142,
2.792137E-136,
1.874071E-132,
2.526557E-132,
8.155516E-132,
7.173739E-131,
5.161872E-126,
9.263253E-125,
1.976302E-124,
3.424844E-123,
4.756065E-123,
2.696966E-121,
1.249935E-120,
1.72185E-116,
2.1101E-116,
3.517283E-116,
1.511015E-115,
2.077794E-115,
3.146068E-115,
1.650202E-114,
4.672468E-112,
7.60232E-107,
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
fn test_pearson_and_bh_all() {
    // Some parameters
    let is_all_vs_all = true;
    let keep_top_n = None; // Keep all the results
    let collect_gem_dataset = Some(true); // Better performance. Keep GEM file in memory

    let (result, number_of_elements_evaluated) = ANALYSIS
        .compute(
            CorrelationMethod::Pearson,
            0.0, // No threshold
            2_000_000,
            AdjustmentMethod::BenjaminiHochberg,
            is_all_vs_all,
            collect_gem_dataset,
            keep_top_n,
        )
        .unwrap();

    // No threshold, no rows filtered
    assert_eq!(number_of_elements_evaluated, TOTAL_COMBINATIONS_EVALUATED);
    assert_eq!(result.len(), number_of_elements_evaluated); // No element should be filtered
}

#[test]
/// Tests Pearson correlation with Benjamini-Hochberg adjustment. Correlation threshold set to 0.7
fn test_pearson_and_bh_cor_0_7() {
    // Some parameters
    let is_all_vs_all = true;
    let keep_top_n = None; // Keep all the results
    let collect_gem_dataset = Some(true); // Better performance. Keep GEM file in memory

    let (result, number_of_elements_evaluated) = ANALYSIS
        .compute(
            CorrelationMethod::Pearson,
            0.7,
            2_000_000,
            AdjustmentMethod::BenjaminiHochberg,
            is_all_vs_all,
            collect_gem_dataset,
            keep_top_n,
        )
        .unwrap();

    assert_eq!(number_of_elements_evaluated, TOTAL_COMBINATIONS_EVALUATED); // The number of evaluated elements mustn't be modified
    assert_eq!(result.len(), 31);

    let collected_as_tuples = get_tuples_from_result(&result);
    assert_eq_results(&collected_as_tuples, &EXPECTED_PEARSON_BH);
}

#[test]
/// Tests Pearson correlation with Benjamini-Hochberg adjustment. No threshold, Keeps top 10 results by correlation
fn test_pearson_and_bh_top_10() {
    // Some parameters
    let is_all_vs_all = true;
    let keep_top_n = Some(10); // Keep top 10
    let collect_gem_dataset = Some(true); // Better performance. Keep GEM file in memory

    let (result, number_of_elements_evaluated) = ANALYSIS
        .compute(
            CorrelationMethod::Pearson,
            0.0, // No threshold
            2_000_000,
            AdjustmentMethod::BenjaminiHochberg,
            is_all_vs_all,
            collect_gem_dataset,
            keep_top_n,
        )
        .unwrap();

    // No threshold, no rows filtered
    assert_eq!(number_of_elements_evaluated, TOTAL_COMBINATIONS_EVALUATED); // The number of evaluated elements mustn't be modified
    assert_eq!(result.len(), 10); // Keeps only 10 elements

    let collected_as_tuples = get_tuples_from_result(&result);

    // If top N is defined it should be sorted by correlation ascending
    let sorted: ResultTupleSimple = EXPECTED_PEARSON_BH
        .iter()
        .sorted_by(|a, b| b.2.partial_cmp(&a.2).unwrap())
        .cloned()
        .take(10)
        .collect();

    // Check with a zip to prevent issues with floating point equality
    assert_eq_results(&collected_as_tuples, &sorted);
}

#[test]
/// Tests Pearson correlation with Benjamini-Hochberg adjustment. Correlation threshold set to 1 (no results)
fn test_pearson_and_bh_cor_1() {
    // Some parameters
    let is_all_vs_all = true;
    let keep_top_n = None; // Keep all the results
    let collect_gem_dataset = Some(true); // Better performance. Keep GEM file in memory

    let (result, number_of_elements_evaluated) = ANALYSIS
        .compute(
            CorrelationMethod::Pearson,
            1.0,
            2_000_000,
            AdjustmentMethod::BenjaminiHochberg,
            is_all_vs_all,
            collect_gem_dataset,
            keep_top_n,
        )
        .unwrap();

    assert_eq!(number_of_elements_evaluated, TOTAL_COMBINATIONS_EVALUATED); // The number of evaluated elements mustn't be modified
    assert_eq!(result.len(), 0);
}

#[test]
/// Tests Pearson correlation with Benjamini-Hochberg adjustment. Only matching genes/GEMs
fn test_pearson_and_bh_only_matching() {
    // Some parameters
    let is_all_vs_all = false; // Keeps only matching genes/GEMs
    let keep_top_n = None; // Keep all the results
    let collect_gem_dataset = Some(true); // Better performance. Keep GEM file in memory

    let (result, number_of_elements_evaluated) = ANALYSIS
        .compute(
            CorrelationMethod::Pearson,
            1.0,
            2_000_000,
            AdjustmentMethod::BenjaminiHochberg,
            is_all_vs_all,
            collect_gem_dataset,
            keep_top_n,
        )
        .unwrap();

    // As there is no matching genes/GEMs no combinations are evaluated
    assert_eq!(number_of_elements_evaluated, 0);
    assert_eq!(result.len(), 0);
}

#[test]
/// Tests Pearson correlation with Bonferroni adjustment. Correlation threshold set to 0.7
fn test_pearson_and_bonferroni_cor_0_7() {
    // Some parameters
    let is_all_vs_all = true;
    let keep_top_n = None; // Keep all the results
    let collect_gem_dataset = Some(true); // Better performance. Keep GEM file in memory

    let (result, number_of_elements_evaluated) = ANALYSIS
        .compute(
            CorrelationMethod::Pearson,
            0.7,
            2_000_000,
            AdjustmentMethod::Bonferroni,
            is_all_vs_all,
            collect_gem_dataset,
            keep_top_n,
        )
        .unwrap();

    // The adjustment method should not modify the number of resulting combinations
    assert_eq!(number_of_elements_evaluated, TOTAL_COMBINATIONS_EVALUATED);
    assert_eq!(result.len(), 31);

    let collected_as_tuples = get_tuples_from_result(&result);

    // Bonferroni does not sort, so sorts both vec of tuples by correlation descending and compares
    let result_sorted = get_sorted_by_correlation(&collected_as_tuples);
    let expected_sorted = get_sorted_by_correlation(&EXPECTED_PEARSON_BONFERRONI);

    assert_eq_results(&result_sorted, &expected_sorted);
}

#[test]
/// Tests Pearson correlation with Benjamini-Yekutieli adjustment. Correlation threshold set to 0.7
fn test_pearson_and_by_cor_0_7() {
    // Some parameters
    let is_all_vs_all = true;
    let keep_top_n = None; // Keep all the results
    let collect_gem_dataset = Some(true); // Better performance. Keep GEM file in memory

    let (result, number_of_elements_evaluated) = ANALYSIS
        .compute(
            CorrelationMethod::Pearson,
            0.7,
            2_000_000,
            AdjustmentMethod::BenjaminiYekutieli,
            is_all_vs_all,
            collect_gem_dataset,
            keep_top_n,
        )
        .unwrap();

	// The adjustment method should not modify the number of resulting combinations
    assert_eq!(number_of_elements_evaluated, TOTAL_COMBINATIONS_EVALUATED);
    assert_eq!(result.len(), 31);

    let collected_as_tuples = get_tuples_from_result(&result);
    assert_eq_results(&collected_as_tuples, &EXPECTED_PEARSON_BY);
}