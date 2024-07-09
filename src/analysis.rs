extern crate extsort;
use crate::dataset::{Dataset, GGCAError};
use crate::types::{TupleExpressionValues, VecOfResults};
use crate::{
    adjustment::{get_adjustment_method, AdjustmentMethod},
    correlation::{get_correlation_method, CorResult, Correlation, CorrelationMethod},
};
use extsort::ExternalSorter;
use itertools::Itertools;
use log::warn;
// Do not remove, it's used for tee()
use pyo3::{create_exception, prelude::*};
use rayon::iter::{IntoParallelIterator, ParallelIterator};
use std::fs;
use std::sync::Mutex;

create_exception!(ggca, GGCADiffSamplesLength, pyo3::exceptions::PyException);
create_exception!(ggca, GGCADiffSamples, pyo3::exceptions::PyException);

fn get_correlation_result(
    tuple_1: TupleExpressionValues,
    tuple_2: TupleExpressionValues,
    correlation_method_struct: &dyn Correlation,
) -> CorResult {
    // Gene and GEM
    let gene = tuple_1.0;
    let gem = tuple_2.0;
    let cpg_site_id = tuple_2.1;

    let (correlation, p_value) =
        correlation_method_struct.correlate(tuple_1.2.as_slice(), tuple_2.2.as_slice());

    CorResult {
        gene,
        gem,
        cpg_site_id,
        correlation: Some(correlation),
        p_value: Some(p_value),
        adjusted_p_value: None,
    }
}

/// Generates a cartesian product without checking for equal genes
fn cartesian_all_vs_all(
    tuple_1: TupleExpressionValues,
    tuple_2: TupleExpressionValues,
    correlation_method_struct: &dyn Correlation,
) -> CorResult {
    get_correlation_result(tuple_1, tuple_2, correlation_method_struct)
}

/// Generates a cartesian product checking for equal genes saving computation time
fn cartesian_equal_genes(
    tuple_1: TupleExpressionValues,
    tuple_2: TupleExpressionValues,
    correlation_method_struct: &dyn Correlation,
) -> CorResult {
    // Gene and GEM
    let gene = tuple_1.0.clone();
    let gem = tuple_2.0.clone();
    let cpg_site_id = tuple_2.1.clone();

    if gene == gem {
        get_correlation_result(tuple_1, tuple_2, correlation_method_struct)
    } else {
        CorResult {
            gene,
            gem,
            cpg_site_id,
            correlation: None,
            p_value: None,
            adjusted_p_value: None,
        }
    }
}

/// Independent struct to log warnings in case some combinations are filtered
struct ConstantInputError {
    number_of_cor_filtered: Mutex<usize>,
}

impl ConstantInputError {
    fn new() -> Self {
        let _ = env_logger::try_init(); // To log warnings
        ConstantInputError {
            number_of_cor_filtered: Mutex::new(0),
        }
    }

    /// Checks if p-value is NaN
    fn p_value_is_nan(&self, cor_result: &CorResult) -> bool {
        if cor_result.p_value.unwrap().is_nan() {
            *self.number_of_cor_filtered.lock().unwrap() += 1;
            true
        } else {
            false
        }
    }

    /// Logs a warning indicating (if needed) that some correlations were filtered
    fn warn_if_needed(&mut self) {
        let number_of_filtered = *self.number_of_cor_filtered.lock().unwrap();
        if number_of_filtered > 0 {
            warn!(
                "{} combinations produced NaNs values as an input array is constant. The correlation coefficient is not defined so that/those combination/s were be filtered.",
                number_of_filtered
            );
        }
    }
}

/// Represents a correlation analysis
#[derive(Clone, Debug)]
pub struct Analysis {
    /// mRNA file path
    pub gene_file_path: String,
    /// Gene Expression Modulator (GEM) file path
    pub gem_file_path: String,
    /// Indicates if the second column of GEM dataset contains CpG Site IDs
    pub gem_contains_cpg: bool,
    /// Correlation method (Pearson, Spearman or Kendall)
    pub correlation_method: CorrelationMethod,
    /// The threshold to discard all results whose correlation statistic values are below this value
    pub correlation_threshold: f64,
    /// Number of elements to hold in memory during external sorting. Bigger values gives better performance but uses more memory
    pub sort_buf_size: usize,
    /// P-values adjustment method (Benjamini-Hochberg, Benjamini-Yekutieli, Bonferroni or None)
    pub adjustment_method: AdjustmentMethod,
    /// True if all Genes must be evaluated with all GEMs. Otherwise only matching Genes/GEM will be evaluated (useful for CNA or Methylation analysis)
    pub is_all_vs_all: bool,
    /// True to make the GEM dataset available in memory. This has a HUGE impact in analysis performance. Specify a boolean value to force or use None to allocate in memory automatically when GEM dataset size is small (<= 100MB)
    pub collect_gem_dataset: Option<bool>,
    /// Specify a number of results to keep or None to return all the resulting combinations
    pub keep_top_n: Option<usize>,
}

impl Analysis {
    fn run_analysis(
        &self,
        dataset_1: Dataset,
        dataset_2: Dataset,
        number_of_samples: usize,
        _should_collect_gem_dataset: bool,
    ) -> PyResult<(VecOfResults, usize, usize)> {
        // Cartesian product computing correlation and p-value (two-sided)
        let correlation_method_struct =
            get_correlation_method(&self.correlation_method, number_of_samples);
        let correlation_function = if self.is_all_vs_all {
            cartesian_all_vs_all
        } else {
            cartesian_equal_genes
        };

        // Right part of iproduct must implement Clone. For more info read:
        // https://users.rust-lang.org/t/iterators-over-csv-files-with-iproduct/51947
        let d1_vec = dataset_1.lazy_matrix.collect_vec();
        let d2_vec = dataset_2.lazy_matrix.collect_vec();

        // Dataset 1 is consumed (into_iter) due it's the external iterator
        // Dataset 2 is referenced (iter) due each thread needs to iterate over its elements
        let correlations_and_p_values = d1_vec
            .into_par_iter()
            .map(|x| {
                d2_vec
                    .iter()
                    .map(|y| {
                        correlation_function(x.clone(), y.clone(), &*correlation_method_struct)
                    })
                    .collect_vec()
            })
            .flatten();

        // Filtering by equal genes (if needed) and NaN values
        let mut nan_errors = ConstantInputError::new();
        let filtered: Box<dyn Iterator<Item = CorResult>> = if self.is_all_vs_all {
            let filtered_nan = correlations_and_p_values
                .filter(|cor_result| !nan_errors.p_value_is_nan(cor_result));
            Box::new(filtered_nan.collect::<Vec<_>>().into_iter())
        } else {
            let filtered_nan = correlations_and_p_values.filter(|cor_result| {
                cor_result.gene == cor_result.gem && !nan_errors.p_value_is_nan(cor_result)
            });
            Box::new(filtered_nan.collect::<Vec<_>>().into_iter())
        };

        // Counts element for future p-value adjustment
        let (filtered, filtered_aux) = filtered.tee();
        let number_of_evaluated_combinations = filtered_aux.count();

        // Sorting (for future adjustment). Note: consumes iterator
        let sorted: Box<dyn Iterator<Item = CorResult>> = match self.adjustment_method {
            AdjustmentMethod::Bonferroni => Box::new(filtered),
            _ => {
                // Benjamini-Hochberg and Benjamini-Yekutieli needs sort by p-value (descending order) to
                // make the adjustment
                let sorter = ExternalSorter::new()
                    .with_parallel_sort()
                    .with_segment_size(self.sort_buf_size);
                Box::new(sorter.sort_by(filtered, |result_1, result_2| {
                    result_2
                        .p_value
                        .unwrap()
                        .partial_cmp(&result_1.p_value.unwrap())
                        .unwrap()
                })?)
            }
        };

        // Ranking
        let ranked = sorted.enumerate();

        // Adjustment
        let mut adjustment_struct = get_adjustment_method(
            &self.adjustment_method,
            number_of_evaluated_combinations as f64,
        );
        let adjusted = ranked.map(|(rank, mut cor_and_p_value)| {
            // Unwrap is safe as None values were filtered before (they are None when is_all_vs_all is false)
            let p_value = cor_and_p_value.p_value.unwrap();
            let q_value = adjustment_struct.adjust(p_value, rank);
            cor_and_p_value.adjusted_p_value = Some(q_value);

            cor_and_p_value
        });

        // Filtering. Improves performance adjusting only the final combinations, note that
        // ranking is preserved
        let filtered = adjusted.filter(|cor_and_p_value| {
            // Unwrap is safe as None correlation values were filtered before (they are None
            // when is_all_vs_all == false)
            cor_and_p_value.abs_correlation() >= self.correlation_threshold
        });

        // Keep top N if needed
        let (limited, total_combinations_count): (Box<dyn Iterator<Item = CorResult>>, usize) =
            match self.keep_top_n {
                Some(top_n) => {
                    let (filtered, filtered_aux_count) = filtered.tee();

                    // Sorts by correlation in descending order
                    let sorter = ExternalSorter::new()
                        .with_parallel_sort()
                        .with_segment_size(self.sort_buf_size);
                    let sorted_cor_desc =
                        sorter.sort_by(filtered, |combination_1, combination_2| {
                            // Unwrap is safe as correlation values are all valid in this stage of algorithm (None
                            // were discarded in all vs all checking stage)
                            combination_2
                                .abs_correlation()
                                .partial_cmp(&combination_1.abs_correlation())
                                .unwrap()
                        })?;
                    (
                        Box::new(sorted_cor_desc.take(top_n)),
                        filtered_aux_count.count(),
                    )
                }
                None => (Box::new(filtered), number_of_evaluated_combinations),
            };

        let result = limited.collect::<VecOfResults>();

        // Generates warnings if needed
        nan_errors.warn_if_needed();

        Ok((
            result,
            total_combinations_count,
            number_of_evaluated_combinations,
        ))
    }

    /// Gets DataFrame and the total number of samples
    fn datasets_and_shapes(
        &self,
        df1_path: &str,
        df2_path: &str,
        gem_contains_cpg: bool,
    ) -> PyResult<(Dataset, Dataset, usize)> {
        let dataset_1 = Dataset::new(df1_path, false)?;
        let mut dataset_1_headers = dataset_1.headers.clone();
        let dataset_1_aux = Dataset::new(df1_path, false)?;
        let mut d1_lazy_matrix_aux = dataset_1_aux.lazy_matrix;
        let one_row = d1_lazy_matrix_aux.next();

        match one_row {
            None => Err(GGCAError::new_err(
                "Error reading first row of mRNA file, is it empty?",
            )),
            Some(row) => {
                let number_of_samples = row.2.len();

                let matrix_2 = Dataset::new(df2_path, gem_contains_cpg)?;
                let m2_headers = matrix_2.headers.clone();

                // In case of an extra column for CpG Site IDs it shouldn't threw an error
                // Index column will be remove after
                let mut m2_headers_to_compare = m2_headers;
                let mut advice_message = "";
                if gem_contains_cpg {
                    m2_headers_to_compare.remove(1); // Removes CpG Site IDs column
                    advice_message = "(without taking the CpG Site IDs column into account)";
                }

                if dataset_1_headers.len() != m2_headers_to_compare.len() {
                    Err(GGCADiffSamplesLength::new_err(format!(
                        "Length of samples in datasets are different:
							\tSamples in mRNA dataset -> {}
							\tSamples in GEM dataset -> {} {}",
                        dataset_1_headers.len(),
                        m2_headers_to_compare.len(),
                        advice_message
                    )))
                } else {
                    // Needs to remove index column as these are not samples and could have different names
                    dataset_1_headers.remove(0);
                    m2_headers_to_compare.remove(0);

                    if dataset_1_headers != m2_headers_to_compare {
                        Err(GGCADiffSamples::new_err(
							"Samples in both datasets are different but they have the same length, maybe they are in different order"
						))
                    } else {
                        Ok((dataset_1, matrix_2, number_of_samples))
                    }
                }
            }
        }
    }

    /// Computes analysis and returns a vec of CorResult, the number of combinations before truncating by 'keep_top_n' parameter and the number of combinations evaluated
    pub fn compute(&self) -> PyResult<(VecOfResults, usize, usize)> {
        let (m1, m2, number_of_samples) = self.datasets_and_shapes(
            self.gene_file_path.as_str(),
            self.gem_file_path.as_str(),
            self.gem_contains_cpg,
        )?;

        let should_collect_gem_dataset = match self.collect_gem_dataset {
            Some(collect_gem_dataset_value) => collect_gem_dataset_value,
            None => {
                // If None, it's defined automatically from the size of the dataset:
                // It'll be collected if the GEM dataset size is less than 100 MB
                let metadata = fs::metadata(self.gem_file_path.as_str())?;
                let size_in_mb = metadata.len() / 1_048_576;
                size_in_mb <= 100
            }
        };

        self.run_analysis(m1, m2, number_of_samples, should_collect_gem_dataset)
    }
}
