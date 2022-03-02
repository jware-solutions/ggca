pub mod common;

use crate::common::compute_no_truncate;
use ggca::{adjustment::AdjustmentMethod, correlation::CorrelationMethod};

// Datasets's paths
const DF1_PATH: &str = "tests/small_files/mRNA_with_NaNs.csv"; // mRNA = 10 rows
const DF2_PATH: &str = "tests/small_files/miRNA_with_NaNs.csv"; // miRNA = 10 rows

// Total number of combinations evaluated in a "all vs all" analysis (10 * 10 - 10 miRNA rows with problematic std)
const TOTAL_COMBINATIONS_EVALUATED: usize = 90;

#[test]
/// Tests that NaNs correlations or p-values are filtered
fn nans_not_panic() {
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

    assert_eq!(number_of_elements_evaluated, TOTAL_COMBINATIONS_EVALUATED);

    // The below miRNA has no std which leads to NaNs values during correlation analysis.
    // So the specified miRNA should be filtered and not included in final result
    assert!(result
        .iter()
        .find(|elem| elem.gem == "hsa-miR-21*")
        .is_none())
}
