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
        ("HOXC10".to_string(), "HOXC4".to_string(), "cg09720701".to_string(), 0.6013793532330399, 6.0402258309228715e-62, 8.097054035609303e-59),
        ("HOXB3".to_string(), "HOXB2".to_string(), "cg25732028".to_string(), -0.6030351863107195, 2.301408578647104e-62, 3.153645746526556e-59),
        ("HOXB4".to_string(), "HOXB4".to_string(), "cg24114154".to_string(), -0.6069504764067456, 2.2976160543131224e-63, 3.2200044630264634e-60),
        ("HOXB13".to_string(), "HOXB13".to_string(), "cg26847748".to_string(), 0.6075361415493422, 1.623254098395665e-63, 2.3278218772900068e-60),
        ("HOXB2".to_string(), "HOXB2".to_string(), "cg13872950".to_string(), -0.6121575462765358, 1.0200352453630917e-64, 1.4976060326207069e-61),
        ("HOXC13".to_string(), "HOXC13".to_string(), "cg25936147".to_string(), -0.6122318077744325, 9.75307370019449e-65, 1.4668622845092513e-61),
        ("HOXB4".to_string(), "HOXB3".to_string(), "cg23551720".to_string(), -0.6126998098290409, 7.3501749589166866e-65, 1.1331029716665963e-61),
        ("HSD11B1".to_string(), "HS6ST1".to_string(), "cg20873136".to_string(), -0.613526262114304, 4.455150985975635e-65, 7.0441648820308095e-62),
        ("HOXB13".to_string(), "MIR3185".to_string(), "cg09000583".to_string(), 0.6165664052668972, 6.973436020846134e-66, 1.1316051547090948e-62),
        ("HOXB13".to_string(), "HOXB13".to_string(), "cg09000583".to_string(), 0.6165664052668972, 6.973436020846134e-66, 1.1316051547090948e-62),
        ("HOXC11".to_string(), "HOXC8".to_string(), "cg05971894".to_string(), 0.6186801775237464, 1.8977516643002807e-66, 3.250637739650348e-63),
        ("HOXC11".to_string(), "HOXC8".to_string(), "cg03218797".to_string(), 0.618838373223325, 1.7209392373992004e-66, 3.031999918142408e-63),
        ("HSD11B1".to_string(), "HOXC12".to_string(), "cg02066277".to_string(), -0.618871284933843, 1.6862658674107686e-66, 3.031999918142408e-63),
        ("HOXB4".to_string(), "HOXB3".to_string(), "cg02311193".to_string(), -0.6195278708279479, 1.1229473603646496e-66, 2.0983462433189622e-63),
        ("HOXC11".to_string(), "HOXC11".to_string(), "cg17273416".to_string(), -0.6197657983158967, 9.688876433477502e-67, 1.8670464887311146e-63),
        ("HOXB13".to_string(), "HOXB13".to_string(), "cg19057986".to_string(), 0.6343046982076374, 9.18554102852667e-71, 1.8271522644615115e-67),
        ("HOXB13".to_string(), "MIR3185".to_string(), "cg19057986".to_string(), 0.6343046982076374, 9.18554102852667e-71, 1.8271522644615115e-67),
        ("HOXB2".to_string(), "HOXB2".to_string(), "cg11955385".to_string(), -0.6348636355653415, 6.369716762741172e-71, 1.354421429164385e-67),
        ("HOXB3".to_string(), "HOXB2".to_string(), "cg22777724".to_string(), -0.6350728746346477, 5.552896795573257e-71, 1.2229065285793904e-67),
        ("HOXB3".to_string(), "HOXB2".to_string(), "cg23530553".to_string(), -0.6389292467567655, 4.3426798515303703e-72, 9.91803742091736e-69),
        ("HOXB3".to_string(), "HOXB2".to_string(), "cg09313705".to_string(), -0.6402376057932351, 1.8140745700856504e-72, 4.302426703452367e-69),
        ("HOXB2".to_string(), "HOXB3".to_string(), "cg23551720".to_string(), -0.6464311235372376, 2.7462450168005155e-74, 6.773778108639479e-71),
        ("HOXC10".to_string(), "HOXC8".to_string(), "cg05971894".to_string(), 0.6476611838539502, 1.180980500479288e-74, 3.0343325658981175e-71),
        ("HOXB13".to_string(), "HOXB13".to_string(), "cg21865150".to_string(), 0.6482956224777938, 7.630406120787938e-75, 2.045745056662032e-71),
        ("HOXC13".to_string(), "HOXC13".to_string(), "cg23922739".to_string(), 0.6499562999977739, 2.420233693011877e-75, 6.783695020267472e-72),
        ("HOXB2".to_string(), "HOXB2".to_string(), "cg21097733".to_string(), -0.6501901157202621, 2.057750042255915e-75, 6.042338028841369e-72),
        ("HOXB3".to_string(), "HOXB2".to_string(), "cg20401567".to_string(), -0.652068511364378, 5.55950809843457e-76, 1.7141075369093467e-72),
        ("HOXC10".to_string(), "HOXC8".to_string(), "cg03218797".to_string(), 0.6535849358109052, 1.9197411565957302e-76, 6.230469404227322e-73),
        ("HOXB2".to_string(), "HOXB2".to_string(), "cg22807449".to_string(), -0.6541660597634334, 1.2751696418241987e-76, 4.368447821858188e-73),
        ("HOXB2".to_string(), "HOXB2".to_string(), "cg06942183".to_string(), -0.6563924505612371, 2.637630191783614e-77, 9.567460479184986e-74),
        ("HOXB3".to_string(), "HOXB2".to_string(), "cg26841048".to_string(), -0.6564837896848054, 2.4718012705255004e-77, 9.526322096605279e-74),
        ("HOXC10".to_string(), "HOXC10".to_string(), "cg20402783".to_string(), -0.6565575038635187, 2.345565733442542e-77, 9.526322096605279e-74),
        ("HOXB3".to_string(), "HOXB2".to_string(), "cg19047660".to_string(), -0.6567202567874764, 2.089104498275663e-77, 9.201609984405034e-74),
        ("HOXB2".to_string(), "HOXB2".to_string(), "cg21780393".to_string(), -0.6571818156170274, 1.5037559041664506e-77, 7.132892621116924e-74),
        ("HOXB4".to_string(), "HOXB4".to_string(), "cg21460081".to_string(), -0.6612863267633756, 7.875538776668605e-79, 4.046976859370774e-75),
        ("HOXB13".to_string(), "HOXB13".to_string(), "cg26506288".to_string(), 0.6622638655917432, 3.874605552328039e-79, 2.17203342526142e-75),
        ("HOXB2".to_string(), "HOXB2".to_string(), "cg20811755".to_string(), -0.6723087139998452, 2.2603293259189796e-82, 1.3938094755346795e-78),
        ("HOXC10".to_string(), "HOXC10".to_string(), "cg13615998".to_string(), -0.6770494794524183, 6.0680805742207e-84, 4.157579116986059e-80),
        ("HOXB2".to_string(), "HOXB2".to_string(), "cg22777724".to_string(), -0.6821210944935389, 1.1728470402892026e-85, 9.040304986549174e-82),
        ("HOXB3".to_string(), "HOXB3".to_string(), "cg23551720".to_string(), -0.6843012206126284, 2.09786643656353e-86, 1.8480405134893358e-82),
        ("HOXB2".to_string(), "HOXB2".to_string(), "cg09313705".to_string(), -0.6844456544670487, 1.8707961678261294e-86, 1.8480405134893358e-82),
        ("HOXB2".to_string(), "HOXB2".to_string(), "cg23530553".to_string(), -0.6945849268761372, 5.0774319293938145e-90, 6.261895249882803e-86),
        ("HOXB2".to_string(), "HOXB2".to_string(), "cg26841048".to_string(), -0.7088087617127972, 2.7808040285235805e-95, 4.286887490371952e-91),
        ("HOXB2".to_string(), "HOXB2".to_string(), "cg25732028".to_string(), -0.7136550757437042, 3.78555075281375e-97, 7.78107338738357e-93),
        ("HOXB2".to_string(), "HOXB2".to_string(), "cg20401567".to_string(), -0.7661066942629341, 3.574553387966863e-120, 1.1021063005779432e-115),
        ("HOXB2".to_string(), "HOXB2".to_string(), "cg19047660".to_string(), -0.7794346419940004, 5.182521109284191e-127, 3.195749816829003e-122),
    ];
}

/// Asserts equality between 2 vec of tuples. Checks with a zip to prevent issues with floating point equality
pub fn assert_eq_results(result: &ResultTupleWithCpG, expected: &ResultTupleWithCpG) {
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
        keep_top_n: None,
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

    let (result, _total_row_count, number_of_elements_evaluated) = compute_with_cpg(
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
fn cpg_site_ids_pearson_and_bh_cor_0_6() {
    // Some parameters
    let is_all_vs_all = true;
    let collect_gem_dataset = Some(true); // Better performance. Keep GEM file in memory

    let (result, _total_row_count, number_of_elements_evaluated) = compute_with_cpg(
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

    let (result, _total_row_count, number_of_elements_evaluated) = compute_with_cpg(
        DF1_PATH.to_string(),
        DF2_PATH.to_string(),
        CorrelationMethod::Pearson,
        0.6,
        2_000_000,
        AdjustmentMethod::BenjaminiHochberg,
        is_all_vs_all,
        collect_gem_dataset,
    );

    assert_eq!(number_of_elements_evaluated, 1299);
    assert_eq!(result.len(), 27);
}
