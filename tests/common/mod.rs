use approx::assert_relative_eq;
use ggca::correlation::CorResult;
use itertools::Itertools;

/// Vec of tuples with Gene, GEM, correlation, p-value and adjusted p-value
pub type ResultTupleSimple = Vec<(String, String, f64, f64, f64)>;

/// Vec of tuples with Gene, GEM, correlation, p-value
pub type ResultTupleWithoutAdj = Vec<(String, String, f64, f64)>;

/// Generates a Vec of tuples with all the data without the CpG Sited ID
pub fn get_tuples_from_result(result: &Vec<CorResult>) -> ResultTupleSimple {
    result
        .iter()
        .map(|elem| {
            (
                elem.gene.clone(),
                elem.gem.clone(),
                elem.correlation.unwrap(),
                elem.p_value.unwrap(),
                elem.adjusted_p_value.unwrap(),
            )
        })
        .collect()
}

/// Merges correlation and p-values result with adjustment. It's responsibility of developers to provide correct order
pub fn merge_with_adjustment(
    result: &ResultTupleWithoutAdj,
    adjustment_result: &Vec<f64>,
) -> ResultTupleSimple {
    result
        .iter()
        .zip(adjustment_result.iter())
        .map(|(elem, adj)| (elem.0.clone(), elem.1.clone(), elem.2, elem.3, *adj))
        .collect()
}

/// Asserts equality between 2 vec of tuples. Checks with a zip to prevent issues with floating point equality
pub fn assert_eq_results(result: &ResultTupleSimple, expected: &ResultTupleSimple) {
    result.iter().zip(expected.iter()).for_each(|(a, b)| {
        assert_eq!(a.0, b.0);
        assert_eq!(a.1, b.1);
        assert_relative_eq!(a.2, b.2, epsilon = 1e-7); // R cor.test only provides a 6/7 digit precision for correlation values
        assert_relative_eq!(a.3, b.3, epsilon = 1e-9); // P-value
        assert_relative_eq!(a.4, b.4, epsilon = 1e-9); // Adjusted p-value
    });
}

/// Sorts a vec of tuples by abs correlation descending
pub fn get_sorted_by_correlation_abs_desc(result: &ResultTupleSimple) -> ResultTupleSimple {
    result
        .iter()
        .sorted_by(|a, b| b.2.abs().partial_cmp(&a.2.abs()).unwrap())
        .cloned()
        .collect()
}

/// Sorts a vec of tuples by abs correlation descending and then by gene and GEM
pub fn get_sorted_by_correlation_abs_desc_gene_gem(
    result: &ResultTupleSimple,
) -> ResultTupleSimple {
    result
        .iter()
        .sorted_by(|a, b| {
            // If it's the same correlation value, sorts by gene and GEM
            if b.2 == a.2 {
                // If it's the same gene, sorts by GEM (ascending)
                if a.0 == b.0 {
                    a.1.partial_cmp(&b.1).unwrap()
                } else {
                    // Sorts by gene (ascending)
                    a.0.partial_cmp(&b.0).unwrap()
                }
            } else {
                // Sorts by correlation (descending)
                b.2.abs().partial_cmp(&a.2.abs()).unwrap()
            }
        })
        .cloned()
        .collect()
}
