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
/// * `group`: Group instance
fn bench_group(
    gene_file_path: &str,
    gem_file_path: &str,
    gem_contains_cpg: bool,
    group: &mut BenchmarkGroup<WallTime>,
) {
    // Correlation task are heavy, so
    group.sample_size(10);

    // For every correlation method
    for cor_method in [
        CorrelationMethod::Pearson,
        CorrelationMethod::Spearman,
        CorrelationMethod::Kendall,
    ] {
        // Tests some important thresholds
        for threshold in [0.0, 0.6] {
            // Tests all the adjustment methods
            for adjustment in [
                AdjustmentMethod::BenjaminiHochberg,
                AdjustmentMethod::BenjaminiYekutieli,
                AdjustmentMethod::Bonferroni,
            ] {
                // Tests performance of correlation sorting to keep best combinations
                for keep_top_cors in [None, Some(10)] {
                    // Tests with GEM file in disk and RAM
                    for collect_gem_in_ram in [false, true] {
                        let collect_desc = if collect_gem_in_ram { "RAM" } else { "disk" };
                        let keep_top_desc = if let Some(top) = keep_top_cors {
                            format!("Top {}", top)
                        } else {
                            "Keep all".to_string()
                        };
                        group.bench_with_input(
                            BenchmarkId::from_parameter(format!(
                                // "{} (threshold = {:.1} | {}) using {} adjustment. GEM in {}",
                                "{} | {} (t = {:.1} | {}) on {}",
                                &cor_method, &adjustment, threshold, keep_top_desc, collect_desc
                            )),
                            &(
                                &cor_method,
                                threshold,
                                &adjustment,
                                collect_gem_in_ram,
                                keep_top_cors,
                            ),
                            |b,
                             (
                                correlation_method,
                                threshold,
                                adjustment_method,
                                collect_gem_in_ram,
                                keep_top_n,
                            )| {
                                b.iter(|| {
                                    Analysis {
                                        gene_file_path: gene_file_path.to_string(),
                                        gem_file_path: gem_file_path.to_string(),
                                        gem_contains_cpg,
                                        correlation_method: (*correlation_method).clone(),
                                        correlation_threshold: *threshold,
                                        sort_buf_size: 2_000_000,
                                        adjustment_method: (*adjustment_method).clone(),
                                        is_all_vs_all: true,
                                        collect_gem_dataset: Some(*collect_gem_in_ram),
                                        keep_top_n: *keep_top_n,
                                    }
                                    .compute()
                                    .unwrap();
                                })
                            },
                        );
                    }
                }
            }
        }
    }
}

pub fn ggca_benchmarks(c: &mut Criterion) {
    // Test small files without CpG Site IDs
    let mut group_small = c.benchmark_group("Small mRNA/miRNA files");
    bench_group(
        GENE_SMALL_FILE_PATH,
        GEM_SMALL_FILE_PATH,
        false,
        &mut group_small,
    );
    group_small.finish();

    // Test small files with CpG Site IDs
    let mut group_methylation = c.benchmark_group("Small mRNA/Methy. (CpG) files");
    bench_group(
        GENE_METHYLATION_FILE_PATH,
        METHYLATION_FILE_PATH,
        true,
        &mut group_methylation,
    );
    group_methylation.finish();
}

criterion_group!(benches, ggca_benchmarks);
criterion_main!(benches);
