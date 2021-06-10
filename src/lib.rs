pub mod adjustment;
pub mod correlation;
pub mod dataset;
pub mod types;

pub mod analysis {
    extern crate extsort;
    use crate::dataset::{Dataset, GGCAError};
    use crate::types::{CollectedMatrix, TupleExpressionValues, VecOfResults};
    use crate::{
        adjustment::{get_adjustment_method, AdjustmentMethod},
        correlation::{get_correlation_method, CorResult, Correlation, CorrelationMethod},
    };
    use extsort::ExternalSorter;
    use itertools::Itertools; // Do not remove, it's used for tee()
    use pyo3::wrap_pyfunction;
    use pyo3::{create_exception, prelude::*};
    use std::fs;

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

    pub struct Analysis {
        file_1_path: String,
        file_2_path: String,
        gem_contains_cpg: bool,
    }

    impl Analysis {
        /// Creates an instance of Analysis from both datasets path
        /// # Args
        /// * `file_1_path`: Path of mRNA dataset
        /// * `file_2_path`: Path of GEM (miRNA, CNA, Methylation) dataset 
        /// * `gem_contains_cpg`: True if second column of GEM dataset must be considered as a CpG Site ID. False for normal datasets
        pub fn new_from_files(
            file_1_path: String,
            file_2_path: String,
            gem_contains_cpg: bool,
        ) -> Analysis {
            Analysis {
                file_1_path,
                file_2_path,
                gem_contains_cpg,
            }
        }

        /// Generic cartesian product for Iterator and collected LazyMatrix
        /// See [this question](https://users.rust-lang.org/t/use-iterator-or-collected-vec/55324/2)
        fn cartesian_product<X: Clone, Y, J: IntoIterator<Item = Y>>(
            &self,
            i: impl Iterator<Item = X>,
            j: J,
        ) -> impl Iterator<Item = (X, Y)>
        where
            J::IntoIter: Clone,
        {
            i.cartesian_product(j.into_iter())
        }

        fn run_analysis(
            &self,
            dataset_1: Dataset,
            len_m1: u64,
            dataset_2: Dataset,
            len_m2: u64,
            number_of_samples: usize,
            correlation_method: CorrelationMethod,
            correlation_threshold: f64,
            sort_buf_size: u64,
            adjustment_method: AdjustmentMethod,
            is_all_vs_all: bool,
            collect_gem_dataset: bool,
            keep_top_n: Option<usize>,
        ) -> PyResult<(VecOfResults, u64)> {
            // Cartesian product computing correlation and p-value (two-sided)
            let correlation_method_struct =
                get_correlation_method(correlation_method, number_of_samples);
            let correlation_function = if is_all_vs_all {
                cartesian_all_vs_all
            } else {
                cartesian_equal_genes
            };

            // Right part of iproduct must implement Clone. For more info read:
            // https://users.rust-lang.org/t/iterators-over-csv-files-with-iproduct/51947
            let cross_product: Box<
                dyn Iterator<Item = (TupleExpressionValues, TupleExpressionValues)>,
            > = if collect_gem_dataset {
                Box::new(self.cartesian_product(
                    dataset_1.lazy_matrix,
                    dataset_2.lazy_matrix.collect::<CollectedMatrix>(),
                ))
            } else {
                Box::new(self.cartesian_product(dataset_1.lazy_matrix, dataset_2.lazy_matrix))
            };

            let correlations_and_p_values = cross_product.map(|(tuple_1, tuple_2)| {
                correlation_function(tuple_1, tuple_2, &*correlation_method_struct)
            });

            // Filtering by equal genes (if needed)
            // Counts element for future p-value adjustment
            let (filtered, number_of_evaluated_combinations): (
                Box<dyn Iterator<Item = CorResult>>,
                u64,
            ) = if is_all_vs_all {
                (Box::new(correlations_and_p_values), len_m1 * len_m2)
            } else {
                let (res, count_aux) = correlations_and_p_values
                    .filter(|cor_result| cor_result.gene == cor_result.gem)
                    .tee();
                (Box::new(res), count_aux.count() as u64)
            };

            // Sorting (for future adjustment). Note: consumes iterator
            let sorted: Box<dyn Iterator<Item = CorResult>> = match adjustment_method {
                AdjustmentMethod::Bonferroni => Box::new(filtered),
                _ => {
                    // Benjamini-Hochberg and Benjamini-Yekutieli needs sort by p-value (ascending order) to
                    // make the adjustment
                    let sorter = ExternalSorter::new().with_segment_size(sort_buf_size as usize);
                    Box::new(sorter.sort_by(filtered, |result_1, result_2| {
                        result_1.p_value.partial_cmp(&result_2.p_value).unwrap()
                    })?)
                }
            };

            // Ranking
            let ranked = sorted.enumerate();

            // Filtering. Improves performance adjusting only the final combinations, note that
            // ranking is preserved
            let filtered = ranked.filter(|(_, cor_and_p_value)| {
                // Unwrap is safe as None correlation values will be filtered before (they are None
                // when is_all_vs_all == false)
                cor_and_p_value.correlation.unwrap().abs() >= correlation_threshold
            });

            // Adjustment
            let mut adjustment_struct =
                get_adjustment_method(adjustment_method, number_of_evaluated_combinations as f64);
            let adjusted = filtered.map(|(rank, mut cor_and_p_value)| {
                let p_value = cor_and_p_value.p_value;
                // Unwrap is safe as None p-values will be filtered before (they are None
                // when is_all_vs_all == false)
                let q_value = adjustment_struct.adjust(p_value.unwrap(), rank);

                cor_and_p_value.adjusted_p_value = Some(q_value);
                cor_and_p_value
            });

            // Keep top N if needed
            let limited: Box<dyn Iterator<Item = CorResult>> = match keep_top_n {
                Some(top_n) => {
                    // Sorts by correlation in descending order
                    let sorter = ExternalSorter::new().with_segment_size(sort_buf_size as usize);
                    let sorted_cor_desc = sorter.sort_by(adjusted, |combination_1, combination_2| {
                        // Unwrap is safe as Correlation values are all valid in this stage of algorithm (None
                        // were discarded in all vs all checking stage)
                        combination_2.correlation.unwrap().abs().partial_cmp(
                            &combination_1.correlation.unwrap().abs()
                        ).unwrap()
                    })?;
                    Box::new(sorted_cor_desc.take(top_n))
                },
                None => Box::new(adjusted),
            };

            Ok((
                limited.collect::<VecOfResults>(),
                number_of_evaluated_combinations,
            ))
        }

        /// Gets DataFrame and its shape (to compute, for example, the number of evaluated combinations)
        fn datasets_and_shapes(
            &self,
            df1_path: &str,
            df2_path: &str,
            gem_contains_cpg: bool,
        ) -> PyResult<(Dataset, u64, Dataset, u64, usize)> {
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

                    let len_m1 = d1_lazy_matrix_aux.count() + 1; // Plus discarded element by next()

                    let matrix_2 = Dataset::new(df2_path, gem_contains_cpg)?;
                    let m2_headers = matrix_2.headers.clone();
                    let m2_aux = Dataset::new(df2_path, gem_contains_cpg)?;

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
                            let len_m2 = m2_aux.lazy_matrix.count();

                            Ok((
                                dataset_1,
                                len_m1 as u64,
                                matrix_2,
                                len_m2 as u64,
                                number_of_samples,
                            ))
                        }
                    }
                }
            }
        }

        /// Computes analysis and returns a vec of CorResult
        /// # Args
        /// * `correlation_method`: Correlation method (Pearson, Spearman or Kendall)
        /// * `correlation_threshold`: Threshold to discard all results which are below this value
        /// * `sort_buf_size`: Number of elements to hold in memory during external sorting. Bigger values gives better performance but uses more memory
        /// * `adjustment_method`: P-values adjustment method (Benjamini-Hochberg, Benjamini-Yekutieli, Bonferroni or None)
        /// * `is_all_vs_all`: True if all Genes must be evaluated with all GEMs. Otherwise only matching Genes/GEM will be evaluated (useful for CNA or Methylation analysis)
        /// * `collect_gem_dataset`: True to make the GEM dataset available in memory. This has a HUGE impact in analysis performance. Specify a boolean value to force or use None to allocate in memory automatically when GEM dataset size is small (<= 100MB)
        /// * `keep_top_n`: Specify a number of results to keep or None to return all the resulting combinations
        pub fn compute(
            &self,
            correlation_method: CorrelationMethod,
            correlation_threshold: f64,
            sort_buf_size: u64,
            adjustment_method: AdjustmentMethod,
            is_all_vs_all: bool,
            collect_gem_dataset: Option<bool>,
            keep_top_n: Option<usize>,
        ) -> PyResult<(VecOfResults, u64)> {
            let (m1, len_m1, m2, len_m2, number_of_samples) = self.datasets_and_shapes(
                self.file_1_path.as_str(),
                self.file_2_path.as_str(),
                self.gem_contains_cpg,
            )?;

            let should_collect_gem_dataset = match collect_gem_dataset {
                Some(collect_gem_dataset_value) => collect_gem_dataset_value,
                None => {
                    // If None, it's defined automatically from the size of the dataset:
                    // It'll be collected if the GEM dataset size is less than 100 MB
                    let metadata = fs::metadata(self.file_2_path.as_str())?;
                    let size_in_mb = metadata.len() / 1_048_576;
                    size_in_mb <= 100
                }
            };

            self.run_analysis(
                m1,
                len_m1,
                m2,
                len_m2,
                number_of_samples,
                correlation_method,
                correlation_threshold,
                sort_buf_size,
                adjustment_method,
                is_all_vs_all,
                should_collect_gem_dataset,
                keep_top_n,
            )
        }
    }

    /// A Python module implemented in Rust.
    #[pymodule]
    fn ggca(py: Python, m: &PyModule) -> PyResult<()> {
        m.add_function(wrap_pyfunction!(correlate, m)?)?;
        m.add("GGCAError", py.get_type::<GGCAError>())?;
        m.add(
            "GGCADiffSamplesLength",
            py.get_type::<GGCADiffSamplesLength>(),
        )?;
        m.add("GGCADiffSamples", py.get_type::<GGCADiffSamples>())?;
        Ok(())
    }

    #[pyfunction]
    fn correlate(
        py: Python,
        file_1_path: String,
        file_2_path: String,
        correlation_method: i32,
        correlation_threshold: f64,
        sort_buf_size: u64,
        adjustment_method: i32,
        all_vs_all: bool,
        gem_contains_cpg: bool,
        collect_gem_dataset: Option<bool>,
        keep_top_n: Option<usize>,
    ) -> PyResult<(VecOfResults, u64)> {
        py.allow_threads(|| {
            let experiment = Analysis::new_from_files(file_1_path, file_2_path, gem_contains_cpg);
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
}
