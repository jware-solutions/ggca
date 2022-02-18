//! # Gene GEM Correlation Analysis (GGCA)
//! 
//! Computes efficiently the correlation (Pearson, Spearman or Kendall) and the p-value (two-sided) between all the pairs from two datasets. It also supports [CpG Site IDs][cpg-site].
//! 
//! ## Installation
//! 
//! 1. Install [GSL][gsl] >= 2.6 in your system. 
//! 1. Add crate to `Cargo.toml`: `ggca = "0.4.0"`
//! 
//! 
//! ## Usage
//! 
//! **Basic example**:
//! 
//! ```
//! use ggca::adjustment::AdjustmentMethod;
//! use ggca::analysis::Analysis;
//! use ggca::correlation::CorrelationMethod;
//! 
//! // File's paths
//! let df1_path = "mrna.csv";
//! let df2_path = "mirna.csv";
//! 
//! // Some parameters
//! let gem_contains_cpg = false;
//! let is_all_vs_all = true;
//! let keep_top_n = Some(10); // Keeps the top 10 of correlation (sorting by abs values)
//! let collect_gem_dataset = None; // Better performance. Keep small GEM files in memory
//! 
//! let analysis = Analysis::new_from_files(df1_path.to_string(), df2_path.to_string(), false);
//! let (result, number_of_elements_evaluated) = analysis.compute(
//! 	CorrelationMethod::Pearson,
//! 	0.7,
//! 	2_000_000,
//! 	AdjustmentMethod::BenjaminiHochberg,
//! 	is_all_vs_all,
//! 	collect_gem_dataset,
//! 	keep_top_n,
//! ).unwrap();
//! 
//! println!("Number of elements -> {} of {} combinations evaluated", result.len(), number_of_elements_evaluated);
//! 
//! for cor_p_value in result.iter() {
//! 	println!("{}", cor_p_value);
//! }
//! ```
//! 
//! **With CpG Site IDs**:
//! 
//! ```
//! use ggca::adjustment::AdjustmentMethod;
//! use ggca::analysis::Analysis;
//! use ggca::correlation::CorrelationMethod;
//! 
//! // Datasets's paths
//! let df1_path = "mrna.csv";
//! let df2_path = "methylation_with_cpgs.csv";
//! 
//! // Some parameters
//! let gem_contains_cpg = true; // Second column in df2 contains CpG Site IDs
//! let is_all_vs_all = false; // Only matching genes
//! let keep_top_n = Some(10); // Keeps the top 10 of correlation (sorting by abs values)
//! let collect_gem_dataset = None;
//! 
//! let analysis =
//! 	Analysis::new_from_files(df1_path.to_string(), df2_path.to_string(), gem_contains_cpg);
//! 
//! let (result, number_of_elements_evaluated) = analysis.compute(
//! 	CorrelationMethod::Pearson,
//! 	0.8,
//! 	2_000_000,
//! 	AdjustmentMethod::Bonferroni,
//! 	is_all_vs_all,
//! 	collect_gem_dataset,
//! 	keep_top_n,
//! ).unwrap();
//! 
//! for cor_p_value in result.iter() {
//! 	println!("{}", cor_p_value);
//! }
//! 
//! println!(
//! 	"Number of elements -> {} of {} combinations evaluated",
//! 	result.len(),
//! 	number_of_elements_evaluated
//! );
//! ```
//! 
//! 
//! ## More examples
//! 
//! You can check the [examples][examples-folder] folder for more types of analysis!
//! 
//! 
//! [gsl]: https://www.gnu.org/software/gsl/
//! [cpg-site]: https://en.wikipedia.org/wiki/CpG_site
//! [examples-folder]: https://github.com/jware-solutions/ggca/tree/master/examples


pub mod analysis;
pub mod adjustment;
pub mod correlation;
pub mod dataset;
pub mod types;

use adjustment::AdjustmentMethod;
use analysis::Analysis;
use correlation::{CorResult, CorrelationMethod};
use dataset::GGCAError;
use pyo3::wrap_pyfunction;
use pyo3::{create_exception, prelude::*};
use types::VecOfResults;

create_exception!(ggca, GGCADiffSamplesLength, pyo3::exceptions::PyException);
create_exception!(ggca, GGCADiffSamples, pyo3::exceptions::PyException);

#[pyfunction]
fn correlate(
    py: Python,
    gene_file_path: String,
    gem_file_path: String,
    correlation_method: i32,
    correlation_threshold: f64,
    sort_buf_size: usize,
    adjustment_method: i32,
    all_vs_all: bool,
    gem_contains_cpg: bool,
    collect_gem_dataset: Option<bool>,
    keep_top_n: Option<usize>,
) -> PyResult<(VecOfResults, usize)> {
    py.allow_threads(|| {
        let experiment = Analysis::new_from_files(gene_file_path, gem_file_path, gem_contains_cpg);
        let correlation_method = match correlation_method {
            1 => CorrelationMethod::Spearman,
            2 => CorrelationMethod::Kendall,
            _ => CorrelationMethod::Pearson,
        };

        let adjustment_method = match adjustment_method {
            1 => AdjustmentMethod::BenjaminiHochberg,
            2 => AdjustmentMethod::BenjaminiYekutieli,
            _ => AdjustmentMethod::Bonferroni,
        };

        experiment.compute(
            correlation_method,
            correlation_threshold,
            sort_buf_size,
            adjustment_method,
            all_vs_all,
            collect_gem_dataset,
            keep_top_n,
        )
    })
}

/// A Python module implemented in Rust.
#[pymodule]
fn ggca(py: Python, m: &PyModule) -> PyResult<()> {
    m.add_function(wrap_pyfunction!(correlate, m)?)?;
    m.add_class::<CorResult>()?;
    m.add("GGCAError", py.get_type::<GGCAError>())?;
    m.add(
        "GGCADiffSamplesLength",
        py.get_type::<GGCADiffSamplesLength>(),
    )?;
    m.add("GGCADiffSamples", py.get_type::<GGCADiffSamples>())?;
    Ok(())
}