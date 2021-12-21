use assert_float_eq::{assert_float_absolute_eq, afe_is_absolute_eq, afe_absolute_error_msg, afe_abs};
use ggca::{adjustment::AdjustmentMethod, analysis::Analysis, correlation::CorrelationMethod};
use lazy_static::lazy_static;

// Datasets's paths
const DF1_PATH: &str = "tests/small_files/mRNA.csv"; // mRNA = 760 rows
const DF2_PATH: &str = "tests/small_files/miRNA.csv"; // miRNA = 601 rows

type ResultTuple = Vec<(String, String, f64, f64, f64)>;

lazy_static! {
    static ref ANALYSIS: Analysis = {
        let gem_contains_cpg = false;
        Analysis::new_from_files(DF1_PATH.to_string(), DF2_PATH.to_string(), gem_contains_cpg)
    };
}

#[test]
/// Tests Pearson correlation with Benjamini-Hochberg adjustment. No threshold set (uses all the rows)
fn test_pearson_and_bh() {
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
    assert_eq!(number_of_elements_evaluated, 456760); // 760 * 601
    assert_eq!(result.len(), number_of_elements_evaluated);

	// Expected correlations sorted by p-value ascending
	let expected: ResultTuple = vec![
        ("C6orf155".to_string(), "hsa-mir-30a".to_string(), 0.8581767, 8.379477E-219, 3.82741E-213),
		("BTLA".to_string(), "hsa-mir-150".to_string(), 0.8479504, 2.192392E-208, 5.006985E-203),
		("APOBEC4".to_string(), "hsa-mir-449a".to_string(), 0.8009502, 8.516449E-169, 1.296658E-163),
		("APOBEC4".to_string(), "hsa-mir-449b".to_string(), 0.7973792, 3.147156E-166, 3.593738E-161),
		("C5orf20".to_string(), "hsa-mir-150".to_string(), 0.7920717, 1.660754E-162, 1.517132E-157),
		("ATP2B2".to_string(), "hsa-mir-885".to_string(), 0.7912955, 5.696732E-162, 4.336732E-157),
		("C21orf34".to_string(), "hsa-let-7c".to_string(), 0.7882625, 6.691326E-160, 4.366186E-155),
		("C16orf54".to_string(), "hsa-mir-150".to_string(), 0.7851163, 8.648716E-158, 4.937985E-153),
		("C21orf34".to_string(), "hsa-mir-99a".to_string(), 0.774241, 9.323279E-151, 4.731668E-146),
		("C1orf158".to_string(), "hsa-mir-449a".to_string(), 0.7694639, 8.61831E-148, 3.936499E-143),
		("C1orf158".to_string(), "hsa-mir-449b".to_string(), 0.759847, 4.940956E-142, 2.051665E-137),
		("APOBEC4".to_string(), "hsa-mir-449c".to_string(), 0.7531277, 3.617836E-138, 1.377069E-133),
		("C1orf158".to_string(), "hsa-mir-449c".to_string(), 0.7528368, 5.283894E-138, 1.856516E-133),
		("C20orf85".to_string(), "hsa-mir-449a".to_string(), 0.751877, 1.836797E-137, 5.99268E-133),
		("AKAP14".to_string(), "hsa-mir-449c".to_string(), 0.7501377, 1.731085E-136, 5.271269E-132),
		("C21orf34".to_string(), "hsa-mir-125b-2".to_string(), 0.7411945, 1.328645E-131, 3.792948E-127),
		("ACAP1".to_string(), "hsa-mir-150".to_string(), 0.7387877, 2.533343E-130, 6.806646E-126),
		("AKAP14".to_string(), "hsa-mir-449b".to_string(), 0.7381176, 5.722785E-130, 1.452189E-125),
		("AIM2".to_string(), "hsa-mir-155".to_string(), 0.7357108, 1.04683E-128, 2.516578E-124),
		("BANK1".to_string(), "hsa-mir-150".to_string(), 0.7353945, 1.530239E-128, 3.49476E-124),
		("AKAP14".to_string(), "hsa-mir-449a".to_string(), 0.7319602, 9.111215E-127, 1.981733E-122),
		("AMICA1".to_string(), "hsa-mir-150".to_string(), 0.7306179, 4.423763E-126, 9.184536E-122),
		("ART3".to_string(), "hsa-mir-934".to_string(), 0.7223071, 6.370956E-122, 1.265216E-117),
		("CA12".to_string(), "hsa-mir-934".to_string(), -0.7220895, 8.146964E-122, 1.550503E-117),
		("C1orf194".to_string(), "hsa-mir-449b".to_string(), 0.7216007, 1.414585E-121, 2.584503E-117),
		("BLK".to_string(), "hsa-mir-150".to_string(), 0.7202692, 6.320095E-121, 1.110295E-116),
		("C1orf194".to_string(), "hsa-mir-449a".to_string(), 0.7199512, 9.025014E-121, 1.526765E-116),
		("C20orf85".to_string(), "hsa-mir-449b".to_string(), 0.7195478, 1.417124E-120, 2.311734E-116),
		("AGR3".to_string(), "hsa-mir-190b".to_string(), 0.7180282, 7.698685E-120, 1.21257E-115),
		("AGR2".to_string(), "hsa-mir-934".to_string(), -0.712854, 2.255013E-117, 3.433333E-113),
		("AGR2".to_string(), "hsa-mir-577".to_string(), -0.7014995, 3.791311E-112, 5.586191E-108),
        ];

    let collected_as_tuples: ResultTuple = result
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
        .collect();

	// Check with a zip to prevent issues with floating point equality
	collected_as_tuples.iter().zip(expected.iter()).for_each(|(a, b)| {
		assert_eq!(a.0, b.0);
		assert_eq!(a.1, b.1);
		assert_float_absolute_eq!(a.2, b.2);
		assert_float_absolute_eq!(a.3, b.3);
		assert_float_absolute_eq!(a.4, b.4);
	});
}
