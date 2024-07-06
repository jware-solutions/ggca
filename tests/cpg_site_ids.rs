pub mod common;

use approx::assert_relative_eq;
use ggca::{
    adjustment::AdjustmentMethod,
    analysis::Analysis,
    correlation::{CorResult, CorrelationMethod},
    types::VecOfResults,
};
use lazy_static::lazy_static;

/// Vec of tuples with Gene, GEM, CpG Site ID, correlation, p-value and adjusted p-value
pub type ResultTupleWithCpGAndAdj = Vec<(String, String, String, f64, f64, f64)>;

/// Vec of tuples with Gene, GEM, CpG Site ID, correlation, p-value and adjusted p-value
pub type ResultTupleWithCpG = Vec<(String, String, String, f64, f64, f64)>;

// Datasets's paths
const DF1_PATH: &str = "tests/medium_files/methylation_gene.csv"; // mRNA = 41 rows
const DF2_PATH: &str = "tests/medium_files/methylation_gem.csv"; // miRNA = 1505 rows

// Total number of combinations evaluated in a "all vs all" analysis (41 * 1505)
const TOTAL_COMBINATIONS_EVALUATED: usize = 61664;

lazy_static! {
    /// Expected correlations for Pearson method, sorted by p-value descending
    static ref EXPECTED_PEARSON_BH: ResultTupleWithCpGAndAdj = vec![
        ("HOXB2".to_string(), "HOXB2".to_string(), "cg19047660".to_string(), -0.7794346419940181, 5.182521109223952e-127, 4.552350166370168e-125),
        ("HOXB2".to_string(), "HOXB2".to_string(), "cg20401567".to_string(), -0.7661066942629469, 3.5745533879330546e-120, 3.135437554957381e-118),
        ("HOXB2".to_string(), "HOXB2".to_string(), "cg25732028".to_string(), -0.7136550757437177, 3.785550752768964e-97, 3.315798318448088e-95),
        ("HOXB2".to_string(), "HOXB2".to_string(), "cg26841048".to_string(), -0.708808761712803, 2.7808040285098246e-95, 2.4322765902699266e-93),
        ("HOXB2".to_string(), "HOXB2".to_string(), "cg23530553".to_string(), -0.6945849268761468, 5.077431929353865e-90, 4.43477000696426e-88),
        ("HOXB2".to_string(), "HOXB2".to_string(), "cg09313705".to_string(), -0.6844456544670825, 1.8707961677761645e-86, 1.631694128567884e-84),
        ("HOXB3".to_string(), "HOXB3".to_string(), "cg23551720".to_string(), -0.6843012206126322, 2.0978664365574456e-86, 1.8271586997722927e-84),
        ("HOXB2".to_string(), "HOXB2".to_string(), "cg22777724".to_string(), -0.682121094493548, 1.1728470402808636e-85, 1.0200626218882817e-83),
        ("HOXC10".to_string(), "HOXC10".to_string(), "cg13615998".to_string(), -0.6770494794524133, 6.068080574243507e-84, 5.2701707116922764e-82),
        ("HOXB2".to_string(), "HOXB2".to_string(), "cg20811755".to_string(), -0.6723087139998514, 2.2603293259084862e-82, 1.9603508797865105e-80),
        ("HOXB13".to_string(), "HOXB13".to_string(), "cg26506288".to_string(), 0.6622638655917419, 0.0, 0.0),
        ("HOXB4".to_string(), "HOXB4".to_string(), "cg21460081".to_string(), -0.6612863267633795, 7.875538776645934e-79, 6.820747515773804e-77),
        ("HOXB2".to_string(), "HOXB2".to_string(), "cg21780393".to_string(), -0.6571818156170428, 1.5037559041499435e-77, 1.3005274063604784e-75),
        ("HOXB3".to_string(), "HOXB2".to_string(), "cg19047660".to_string(), -0.6567202567874807, 2.089104498269013e-77, 1.8042372518383815e-75),
        ("HOXC10".to_string(), "HOXC10".to_string(), "cg20402783".to_string(), -0.6565575038635021, 2.3455657334701182e-77, 2.0228946208209982e-75)
    ];
}

/// Asserts equality between 2 vec of tuples. Checks with a zip to prevent issues with floating point equality
pub fn assert_eq_results(result: &ResultTupleWithCpG, expected: &ResultTupleWithCpG) {
    // Prints both
    result.iter().zip(expected.iter()).for_each(|(a, b)| {
        assert_eq!(a.0, b.0); // mRNA
        assert_eq!(a.1, b.1); // GEM
        assert_eq!(a.2, b.2); // CpG Site ID
        assert_relative_eq!(a.3, b.3, epsilon = 1e-7); // R cor.test only provides a 6/7 digit precision for correlation values
        assert_relative_eq!(a.4, b.4, epsilon = 1e-9); // P-value
        assert_relative_eq!(a.5, b.5, epsilon = 1e-9); // Adjusted p-value
    });
}

/// Generates a Vec of tuples with all the data
pub fn get_tuples_with_cpg_from_result(result: &Vec<CorResult>) -> ResultTupleWithCpG {
    result
        .iter()
        .map(|elem| {
            (
                elem.gene.clone(),
                elem.gem.clone(),
                elem.cpg_site_id.as_ref().unwrap().clone(),
                elem.correlation.unwrap(),
                elem.p_value.unwrap(),
                elem.adjusted_p_value.unwrap(),
            )
        })
        .collect()
}

/// Computes an analysis with specific parameters (with CpG Site IDs)
fn compute_with_cpg(
    gene_file_path: String,
    gem_file_path: String,
    correlation_method: CorrelationMethod,
    correlation_threshold: f64,
    sort_buf_size: usize,
    adjustment_method: AdjustmentMethod,
    is_all_vs_all: bool,
    collect_gem_dataset: Option<bool>,
    keep_top_n: Option<usize>,
) -> (VecOfResults, usize, usize) {
    Analysis {
        gene_file_path,
        gem_file_path,
        gem_contains_cpg: true,
        correlation_method,
        correlation_threshold,
        sort_buf_size,
        adjustment_method,
        is_all_vs_all,
        collect_gem_dataset,
        keep_top_n,
    }
    .compute()
    .unwrap()
}

#[test]
/// Tests Pearson correlation with Benjamini-Hochberg adjustment. No threshold set (uses all the rows)
fn cpg_site_ids_pearson_and_bh_all() {
    // Some parameters
    let is_all_vs_all = true;
    let collect_gem_dataset = Some(true); // Better performance. Keep GEM file in memory
    let keep_top_n = None;

    let (result, _total_row_count, number_of_elements_evaluated) = compute_with_cpg(
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
    assert_eq!(number_of_elements_evaluated, TOTAL_COMBINATIONS_EVALUATED);
    assert_eq!(result.len(), number_of_elements_evaluated); // No element should be filtered
}

#[test]
/// Tests Pearson correlation with Benjamini-Hochberg adjustment. Correlation threshold set to 0.6
fn cpg_site_ids_pearson_and_bh_cor_0_6() {
    // Some parameters
    let is_all_vs_all = true;
    let collect_gem_dataset = Some(true); // Better performance. Keep GEM file in memory
    let keep_top_n = Some(15);

    let (result, _total_row_count, number_of_elements_evaluated) = compute_with_cpg(
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
    assert_eq!(result.len(), EXPECTED_PEARSON_BH.len());

    let collected_as_tuples = get_tuples_with_cpg_from_result(&result);
    assert_eq_results(&collected_as_tuples, &EXPECTED_PEARSON_BH);
}

#[test]
/// Tests Pearson correlation with Benjamini-Hochberg adjustment. Only matching genes/GEMs
fn cpg_site_ids_pearson_and_bh_only_matching() {
    // Some parameters
    let is_all_vs_all = false; // Keeps only matching genes/GEMs
    let collect_gem_dataset = Some(true); // Better performance. Keep GEM file in memory
    let keep_top_n = None;

    let (result, _total_row_count, number_of_elements_evaluated) = compute_with_cpg(
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

    assert_eq!(number_of_elements_evaluated, 1299);
    assert_eq!(result.len(), 27);
}
