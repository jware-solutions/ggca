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
//! let gene_file_path = "mrna.csv".to_string();
//! let gem_file_path = "mirna.csv".to_string();
//!
//! // Some parameters
//! let gem_contains_cpg = false;
//! let is_all_vs_all = true;
//! let keep_top_n = Some(10); // Keeps the top 10 of correlation (sorting by abs values)
//! let collect_gem_dataset = None; // Better performance. Keep small GEM files in memory
//!
//! // Creates and run an analysis
//! let analysis = Analysis {
//!     gene_file_path,
//!     gem_file_path,
//!     gem_contains_cpg: false,
//!     correlation_method: CorrelationMethod::Pearson,
//!     correlation_threshold: 0.7,
//!     sort_buf_size: 2_000_000,
//!     adjustment_method: AdjustmentMethod::BenjaminiHochberg,
//!     is_all_vs_all,
//!     collect_gem_dataset,
//!     keep_top_n,
//! };
//! 
//! let (result, _total_combinations_count, number_of_elements_evaluated) = analysis.compute().unwrap();
//!
//! println!(
//!     "Number of elements -> {} of {} combinations evaluated",
//!     result.len(),
//!     number_of_elements_evaluated
//! );
//!
//! for cor_p_value in result.iter() {
//!     println!("{}", cor_p_value);
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
//! let gene_file_path = "mrna.csv".to_string();
//! let gem_file_path = "methylation_with_cpgs.csv".to_string();
//!
//! // Some parameters
//! let gem_contains_cpg = true; // Second column in df2 contains CpG Site IDs
//! let is_all_vs_all = false; // Only matching genes
//! let keep_top_n = Some(10); // Keeps the top 10 of correlation (sorting by abs values)
//! let collect_gem_dataset = None;
//!
//! let analysis = Analysis {
//!     gene_file_path,
//!     gem_file_path,
//!     gem_contains_cpg,
//!     correlation_method: CorrelationMethod::Pearson,
//!     correlation_threshold: 0.8,
//!     sort_buf_size: 2_000_000,
//!     adjustment_method: AdjustmentMethod::Bonferroni,
//!     is_all_vs_all,
//!     collect_gem_dataset,
//!     keep_top_n,
//! 
//! };
//!
//! let (result, _total_combinations_count, number_of_elements_evaluated) = analysis.compute().unwrap();
//! 
//! println!(
//!     "Number of elements -> {} of {} combinations evaluated",
//!     result.len(),
//!     number_of_elements_evaluated
//! );
//!
//! for cor_p_value in result.iter() {
//!     println!("{}", cor_p_value);
//! }
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

pub mod adjustment;
pub mod analysis;
pub mod correlation;
pub mod dataset;
pub mod types;

use adjustment::AdjustmentMethod;
use analysis::{Analysis, GGCADiffSamples, GGCADiffSamplesLength};
use correlation::{CorResult, CorrelationMethod};
use dataset::GGCAError;
use pyo3::wrap_pyfunction;
use pyo3::{create_exception, prelude::*};
use types::VecOfResults;

// Errors
create_exception!(
    ggca,
    InvalidCorrelationMethod,
    pyo3::exceptions::PyException
);
create_exception!(ggca, InvalidAdjustmentMethod, pyo3::exceptions::PyException);

// NOTE: Python has named arguments, so this linting warning can be disabled without sacrificing maintainability
#[pyfunction]
#[allow(clippy::too_many_arguments)] fn correlate(
    py: Python,
    gene_file_path: String,
    gem_file_path: String,
    correlation_method: i32,
    correlation_threshold: f64,
    sort_buf_size: usize,
    adjustment_method: i32,
    is_all_vs_all: bool,
    gem_contains_cpg: bool,
    collect_gem_dataset: Option<bool>,
    keep_top_n: Option<usize>,
) -> PyResult<(VecOfResults, usize, usize)> {
    py.allow_threads(|| {
        let correlation_method = match correlation_method {
            1 => Ok(CorrelationMethod::Spearman),
            2 => Ok(CorrelationMethod::Kendall),
            3 => Ok(CorrelationMethod::Pearson),
            selected => Err(
                InvalidCorrelationMethod::new_err(
                    format!("Wrong correlation method ({selected}). Only values 1 (Spearman), 2 (Kendall) or 3 (Pearson) are valid")
                )
            )
        }?;

        let adjustment_method = match adjustment_method {
            1 => Ok(AdjustmentMethod::BenjaminiHochberg),
            2 => Ok(AdjustmentMethod::BenjaminiYekutieli),
            3 => Ok(AdjustmentMethod::Bonferroni),
            selected => Err(
                InvalidAdjustmentMethod::new_err(
                    format!("Wrong adjustment method ({selected}). Only values 1 (Benjamini-Hochberg), 2 (Benjamini-Yekutieli) or 3 (Bonferroni) are valid")
                )
            )
        }?;

        // Creates analysis and run
        let analysis = Analysis {
            gene_file_path, 
            gem_file_path, 
            gem_contains_cpg,
            correlation_method,
            correlation_threshold,
            sort_buf_size,
            adjustment_method,
            is_all_vs_all,
            collect_gem_dataset,
            keep_top_n,
        };

        analysis.compute()
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
    m.add(
        "InvalidCorrelationMethod",
        py.get_type::<InvalidCorrelationMethod>(),
    )?;
    m.add(
        "InvalidAdjustmentMethod",
        py.get_type::<InvalidAdjustmentMethod>(),
    )?;
    Ok(())
}
