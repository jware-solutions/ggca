use approx::assert_relative_eq;
use ggca::adjustment::{get_adjustment_method, AdjustmentMethod};

/// Computes Benjamini-Hochberg and Benjamini-Yekutieli adjustment.
/// To test Bonferroni use `test_adjustment_bh_or_by()` instead as it does not need to sort by p-value
pub fn test_adjustment_bh_or_by(x: Vec<f64>, expected: Vec<f64>, method: &AdjustmentMethod) {
    // To preserve
    let mut ranked = x.iter().enumerate().collect::<Vec<(usize, &f64)>>();

    // Ascending by p-value
    // ranked.sort_by(|result_1, result_2| result_1.1.partial_cmp(result_2.1).unwrap());
    // Descending by p-value
    ranked.sort_by(|result_1, result_2| result_2.1.partial_cmp(result_1.1).unwrap());

    // Ranked
    let ranked = ranked.iter().enumerate();

    // Computes adjustment
    let mut adjustment_struct = get_adjustment_method(method, x.len() as f64);
    let mut adjusted = ranked
        .map(|(rank, (initial_order, p_value))| {
            let q_value = adjustment_struct.adjust(**p_value, rank);

            (*initial_order, **p_value, q_value)
        })
        .collect::<Vec<(usize, f64, f64)>>();

    // Sort by initial order
    adjusted.sort_by(|result_1, result_2| result_1.0.partial_cmp(&result_2.0).unwrap());

    // Asserts
    adjusted
        .iter()
        .zip(expected.iter())
        .for_each(|((_, _, adjusted), expected_adjusted)| {
            assert_relative_eq!(adjusted, expected_adjusted, epsilon = 1e-8);
        });
}
