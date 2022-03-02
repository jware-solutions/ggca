use ggca::adjustment::AdjustmentMethod;
use ggca::analysis::Analysis;
use ggca::correlation::CorrelationMethod;
use pyo3::PyResult;
use std::time::Instant;

fn main() -> PyResult<()> {
    // Datasets's paths
    let gene_file_path = "mrna.csv".to_string();
    let gem_file_path = "methylation_with_cpgs.csv".to_string();

    // Some parameters
    let gem_contains_cpg = true; // Second column in df2 contains CpG Site IDs
    let is_all_vs_all = false; // Only matching genes
    let keep_top_n = Some(10); // Keeps the top 10 of correlation (sorting by abs values)
    let collect_gem_dataset = None;

    let now = Instant::now();

    // Creates and run an analysis
    let analysis = Analysis {
        gene_file_path,
        gem_file_path,
        gem_contains_cpg,
        correlation_method: CorrelationMethod::Pearson,
        correlation_threshold: 0.8,
        sort_buf_size: 2_000_000,
        adjustment_method: AdjustmentMethod::Bonferroni,
        is_all_vs_all,
        collect_gem_dataset,
        keep_top_n,
    };

    let (result, _total_combinations_count, number_of_elements_evaluated) = analysis.compute()?;

    let seconds = now.elapsed().as_secs();

    for cor_p_value in result.iter() {
        println!("{}", cor_p_value);
    }

    println!("Finished in -> {} seconds", seconds);
    println!(
        "Number of elements -> {} of {} combinations evaluated",
        result.len(),
        number_of_elements_evaluated
    );

    Ok(())
}
