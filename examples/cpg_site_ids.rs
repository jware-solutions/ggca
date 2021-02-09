use ggca::adjustment::AdjustmentMethod;
use ggca::analysis::Analysis;
use ggca::correlation::CorrelationMethod;
use pyo3::PyResult;
use std::time::Instant;

fn main() -> PyResult<()> {
    // Datasets's paths
    let df1_path = "mrna.csv";
    let df2_path = "methylation_with_cpgs.csv";

    // Some parameters
    let gem_contains_cpg = true; // Second column in df2 contains CpG Site IDs
    let is_all_vs_all = false; // Only matching genes
    let keep_top_n = Some(10); // Keeps the top 10 of correlation (sorting by abs values)
    let collect_gem_dataset = None;

    let analysis =
        Analysis::new_from_files(df1_path.to_string(), df2_path.to_string(), gem_contains_cpg);

    let now = Instant::now();
    let (result, number_of_elements_evaluated) = analysis.compute(
        CorrelationMethod::Pearson,
        0.8,
        2_000_000,
        AdjustmentMethod::Bonferroni,
        is_all_vs_all,
        collect_gem_dataset,
        keep_top_n,
    )?;

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
