use std::time::Instant;
use ggca::experiment::{new_from_files, Computation};
use ggca::adjustment::AdjustmentMethod;
use ggca::correlation::CorrelationMethod;
use pyo3::PyResult;


fn main() -> PyResult<()> {
    // File's paths
    let df1_path = "mrna.csv";
    let df2_path = "mirna.csv";

    let now = Instant::now();
    
    let experiment = new_from_files(df1_path.to_string(), df2_path.to_string());
	let result = experiment.compute(CorrelationMethod::Pearson, 0.7, 2_000_000, AdjustmentMethod::BenjaminiYekutieli)?;
	
    println!("Finished in -> {} seconds", now.elapsed().as_secs());
    println!("Number of elements -> {}", result.len());

    for cor_p_value in result.iter() {
        println!("{}", cor_p_value);
    }

    Ok(())
}
