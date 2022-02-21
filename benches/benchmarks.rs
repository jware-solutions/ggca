use criterion::{
    criterion_group, criterion_main, measurement::WallTime, BenchmarkGroup, BenchmarkId, Criterion,
};
use ggca::{adjustment::AdjustmentMethod, analysis::Analysis, correlation::CorrelationMethod};

// Datasets's paths
const GENE_SMALL_FILE_PATH: &str = "tests/small_files/mRNA.csv"; // mRNA = 600 rows
const GEM_SMALL_FILE_PATH: &str = "tests/small_files/miRNA.csv"; // miRNA = 299 rows

const GENE_METHYLATION_FILE_PATH: &str = "tests/medium_files/methylation_gene.csv"; // mRNA = 41 rows
const METHYLATION_FILE_PATH: &str = "tests/medium_files/methylation_gem.csv"; // miRNA = 1505 rows

/// Tries all the combinations of correlation methods, threshold, GEM collection and p-value adjustment methods
/// # Args
/// * `gene_file_path`: Genes file path
/// * `gem_file_path`: GEM file path
/// * `gem_contains_cpg`: True if the GEM dataset contains a CpG Site ID column
/// * `datasets_description`: Description to set as ID of the test group
/// * `use_both_collect_strategies`: True if want to test collecting GEM dataset in RAM and disk. False to test only disk
/// * `group`: Group instance
fn bench_group(
    gene_file_path: &str,
    gem_file_path: &str,
    gem_contains_cpg: bool,
    datasets_description: &str,
    use_both_collect_strategies: bool,
    group: &mut BenchmarkGroup<WallTime>,
) {
    // For every correlation method
    for cor_method in [
        CorrelationMethod::Pearson,
        CorrelationMethod::Spearman,
        CorrelationMethod::Kendall,
    ] {
        // Tests some important thresholds
        for threshold in [0.0, 0.6, 1.0] {
            // Tests all the adjustment methods
            for adjustment in [
                AdjustmentMethod::BenjaminiHochberg,
                AdjustmentMethod::BenjaminiYekutieli,
                AdjustmentMethod::Bonferroni,
            ] {
                // Test with GEM file in disk and RAM
                for collect_gem_in_ram in [false, true] {
                    if !collect_gem_in_ram {
                        continue;
                    }; // TODO: remove!!!
                    let collect_desc = if collect_gem_in_ram { "RAM" } else { "disk" };

                    group.bench_with_input(
						BenchmarkId::new(
							format!("{datasets_description}"),
							format!(
								"{} (threshold = {:.1}) using {} adjustment. GEM in {}",
								&cor_method, threshold, &adjustment, collect_desc
							),
						),
						&(&cor_method, threshold, &adjustment, collect_gem_in_ram),
						|b, (correlation_method, threshold, adjustment_method, collect_gem_in_ram)| {
							b.iter(|| Analysis {
								gene_file_path: gene_file_path.to_string(),
								gem_file_path: gem_file_path.to_string(),
								gem_contains_cpg,
								correlation_method: (*correlation_method).clone(),
								correlation_threshold: *threshold,
								sort_buf_size: 2_000_000,
								adjustment_method: (*adjustment_method).clone(),
								is_all_vs_all: true,
								collect_gem_dataset: Some(*collect_gem_in_ram),
								keep_top_n: None,
							})
						},
					);
                }
            }
        }
    }
}

pub fn ggca_benchmarks(c: &mut Criterion) {
    let group_desc = "Small mRNA/miRNA files";
    let mut group_small = c.benchmark_group(group_desc);
    bench_group(
        GENE_SMALL_FILE_PATH,
        GEM_SMALL_FILE_PATH,
        false,
        group_desc,
        true,
        &mut group_small,
    );
    group_small.finish();

    let group_desc = "Small mRNA/Methylation files";
    let mut group_methylation = c.benchmark_group(group_desc);
    bench_group(
        GENE_METHYLATION_FILE_PATH,
        METHYLATION_FILE_PATH,
        true,
        group_desc,
        true,
        &mut group_methylation,
    );
    group_methylation.finish();

    // c.bench_function("small files", |b| {
    //     b.iter(|| {
    //         Analysis {
    //             gene_file_path: DF1_PATH.to_string(),
    //             gem_file_path: DF2_PATH.to_string(),
    //             gem_contains_cpg: false,
    //             correlation_method: CorrelationMethod::Pearson,
    //             correlation_threshold: 0.0, // No threshold
    //             sort_buf_size: 2_000_000,
    //             adjustment_method: AdjustmentMethod::BenjaminiHochberg,
    //             is_all_vs_all: true,
    //             collect_gem_dataset: Some(true),
    //             keep_top_n: None,
    //         }
    //     })
    // });
}

criterion_group!(benches, ggca_benchmarks);
criterion_main!(benches);
