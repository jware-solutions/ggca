//! # Gene GEM Correlation Analysis (GGCA)
//!
//! Computes efficiently the correlation (Pearson, Spearman or Kendall) and the p-value (two-sided) between all the pairs from two datasets. It also supports [CpG Site IDs][cpg-site].
//!
//! ## Installation
//!
//! 1. Add crate to `Cargo.toml`: `ggca = "1.0.0"`
//!
//!
//! ## Usage
//!
//! **Basic example**:
//!
//! ```ignore
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
//! ```ignore
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
#[pyo3(signature = (gene_file_path, gem_file_path, correlation_method, correlation_threshold, sort_buf_size, adjustment_method, is_all_vs_all, gem_contains_cpg, collect_gem_dataset=None, keep_top_n=None))]
#[allow(clippy::too_many_arguments)]
fn correlate(
    py: Python,
    gene_file_path: String,
    gem_file_path: String,
    correlation_method: CorrelationMethod,
    correlation_threshold: f64,
    sort_buf_size: usize,
    adjustment_method: AdjustmentMethod,
    is_all_vs_all: bool,
    gem_contains_cpg: bool,
    collect_gem_dataset: Option<bool>,
    keep_top_n: Option<usize>,
) -> PyResult<(VecOfResults, usize, usize)> {
    py.allow_threads(|| {
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
fn ggca(py: Python<'_>, m: &Bound<'_, PyModule>) -> PyResult<()> {
    // Functions
    m.add_function(wrap_pyfunction!(correlate, m)?)?;

    // Enums
    m.add_class::<CorrelationMethod>()?;
    m.add_class::<AdjustmentMethod>()?;

    // Structs
    m.add_class::<CorResult>()?;

    // Errors
    m.add("GGCAError", py.get_type_bound::<GGCAError>())?;
    m.add(
        "GGCADiffSamplesLength",
        py.get_type_bound::<GGCADiffSamplesLength>(),
    )?;
    m.add("GGCADiffSamples", py.get_type_bound::<GGCADiffSamples>())?;
    m.add(
        "InvalidCorrelationMethod",
        py.get_type_bound::<InvalidCorrelationMethod>(),
    )?;
    m.add(
        "InvalidAdjustmentMethod",
        py.get_type_bound::<InvalidAdjustmentMethod>(),
    )?;
    Ok(())
}
