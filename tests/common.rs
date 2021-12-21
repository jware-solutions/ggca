use assert_float_eq::{
    afe_abs, afe_absolute_error_msg, afe_is_absolute_eq, assert_float_absolute_eq,
};
use ggca::correlation::CorResult;

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
        assert_float_absolute_eq!(a.2, b.2);
        assert_float_absolute_eq!(a.3, b.3);
        assert_float_absolute_eq!(a.4, b.4);
    });
}
