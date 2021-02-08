pub mod types;
pub mod adjustment;
pub mod correlation;
pub mod dataset;

pub mod experiment {
    extern crate extsort;
    use crate::dataset::{Dataset, GGCAError};
    use crate::types::{VecOfResults, TupleExpressionValues, LazyMatrixInner, CollectedMatrix};
    use crate::{
        adjustment::{get_adjustment_method, AdjustmentMethod},
        correlation::{CorResult, CorrelationMethod, get_correlation_method, Correlation},
    };
    use extsort::ExternalSorter;
    use itertools::iproduct;
    use itertools::Itertools;
    use pyo3::wrap_pyfunction;
    use pyo3::{create_exception, prelude::*};

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

    pub trait Computation {
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
        ) -> PyResult<(VecOfResults, u64)> {
            // Cartesian product computing correlation and p-value
            println!("Seleccionando metodo");
            let correlation_method_struct =
                get_correlation_method(correlation_method, number_of_samples);
            let correlation_function = if is_all_vs_all {
                cartesian_all_vs_all
            } else {
                cartesian_equal_genes
            };

            println!("Aplicando cartesiano");
            // Right part of iproduct must implement Clone. For more info read:
            // https://users.rust-lang.org/t/iterators-over-csv-files-with-iproduct/51947
            // let matrix_2_collected = dataset_2.lazy_matrix.collect::<CollectedMatrix>();
            let correlations_and_p_values =
                // iproduct!(dataset_1.lazy_matrix, matrix_2_collected).map(|(tuple_1, tuple_2)| {
                iproduct!(dataset_1.lazy_matrix, dataset_2.lazy_matrix).map(|(tuple_1, tuple_2)| {
                    // println!("Ejecutando iproduct");
                    correlation_function(tuple_1, tuple_2, &*correlation_method_struct)
                });

            // Filtering by equal genes (if needed)
            // Counts element for future p-value adjustment
            let (filtered, number_of_evaluated_combinations): (
                Box<dyn Iterator<Item = CorResult>>,
                u64,
            ) = if is_all_vs_all {
                println!("Se seleccionó todos vs todos");
                (Box::new(correlations_and_p_values), len_m1 * len_m2)
            } else {
                println!("Se seleccionó solo coincidentes");
                let (res, count_aux) = correlations_and_p_values
                    .filter(|cor_result| cor_result.gene == cor_result.gem)
                    .tee();
                (Box::new(res), count_aux.count() as u64)
            };

            // println!("Numero de elementos -> {}",filtered.count());

            // return Ok((
            //     // adjusted.collect::<VecOfResults>(),
            //     Vec::new(), // TODO: remove
            //     number_of_evaluated_combinations,
            // ));

            // Sorting (for future adjustment). Note: consumes iterator
            let sorted: Box<dyn Iterator<Item = CorResult>> = match adjustment_method {
                AdjustmentMethod::Bonferroni => Box::new(filtered),
                _ => {
                    // Benjamini-Hochberg and Benjamini-Yekutieli needs sort by p-value to
                    // make the adjustment
                    let sorter = ExternalSorter::new().with_segment_size(sort_buf_size as usize);
                    Box::new(sorter.sort(filtered)?)
                }
            };

            // Ranking
            println!("Aplicando ranking");
            let ranked = sorted.enumerate();

            // Filtering. Improves performance adjusting only the final combinations, note that
            // ranking is preserved
            println!("Aplicando filtro");
            let filtered = ranked.filter(|(_, cor_and_p_value)| {
                // Unwrap is safe as None correlation values will be filtered before (they are None
                // when is_all_vs_all == false)
                cor_and_p_value.correlation.unwrap().abs() >= correlation_threshold
            });

            // Adjustment
            println!("Aplicando ajuste");
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

            // println!("Returning {} elements", adjusted.count());

            Ok((
                adjusted.collect::<VecOfResults>(),
                number_of_evaluated_combinations,
            ))
        }

        /// Gets DataFrame and its shape (to compute, for example, the number of evaluated combinations)
        fn df_and_shape(
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
                    println!("Hace un remove de headers to compar");
                    if gem_contains_cpg {
                        m2_headers_to_compare.remove(1); // Removes CpG Site IDs column
                        advice_message = "(without taking the CpG Site IDs column into account)";
                    }

                    println!("Compara longitudes");
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
                        println!("Hace un remove de los headers");
                        // Needs to remove index column as these are not samples and could have different names
                        dataset_1_headers.remove(0);
                        m2_headers_to_compare.remove(0);
                        
                        println!("Compara los samples");
                        if dataset_1_headers != m2_headers_to_compare {
                            println!("Hubo un problema");
                            Err(GGCADiffSamples::new_err(
                                "Samples in both datasets are different but they have the same length, maybe they are in different order"
                            ))
                        } else {
                            println!("Hace un count de m2_aux.1");
                            let len_m2 = m2_aux.lazy_matrix.count();
                            
                            println!("Sale de df_and_shape");
                            Ok((dataset_1, len_m1 as u64, matrix_2, len_m2 as u64, number_of_samples))
                        }
                    }
                }
            }
        }

        fn compute(
            &self,
            correlation_method: CorrelationMethod,
            correlation_threshold: f64,
            sort_buf_size: u64,
            adjustment_method: AdjustmentMethod,
            is_all_vs_all: bool,
        ) -> PyResult<(VecOfResults, u64)>;
    }

    pub struct ExperimentFromFiles {
        file_1_path: String,
        file_2_path: String,
        gem_contains_cpg: bool,
    }

    impl Computation for ExperimentFromFiles {
        fn compute(
            &self,
            correlation_method: CorrelationMethod,
            correlation_threshold: f64,
            sort_buf_size: u64,
            adjustment_method: AdjustmentMethod,
            is_all_vs_all: bool,
        ) -> PyResult<(VecOfResults, u64)> {
            let (m1, len_m1, m2, len_m2, number_of_samples) = self.df_and_shape(
                self.file_1_path.as_str(),
                self.file_2_path.as_str(),
                self.gem_contains_cpg,
            )?;

            println!("Comienza el analysis!");
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
            )
        }
    }

    pub fn new_from_files(
        file_1_path: String,
        file_2_path: String,
        gem_contains_cpg: bool,
    ) -> ExperimentFromFiles {
        ExperimentFromFiles {
            file_1_path,
            file_2_path,
            gem_contains_cpg,
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
    pub fn correlate(
        py: Python,
        file_1_path: String,
        file_2_path: String,
        correlation_method: i32,
        correlation_threshold: f64,
        sort_buf_size: u64,
        adjustment_method: i32,
        all_vs_all: bool,
        gem_contains_cpg: bool,
    ) -> PyResult<(VecOfResults, u64)> {
        py.allow_threads(|| {
            let experiment = new_from_files(file_1_path, file_2_path, gem_contains_cpg);
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
            )
        })
    }
}
