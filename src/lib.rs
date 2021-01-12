pub mod adjustment;
pub mod correlation;

pub mod experiment {
    extern crate extsort;
    use crate::correlation::{get_correlation_method, CorrelationMethod};
    use crate::{
        adjustment::{get_adjustment_method, AdjustmentMethod},
        correlation::Correlation,
    };
    use csv::ReaderBuilder;
    use extsort::*;
    use itertools::iproduct;
    use itertools::Itertools;
    use pyo3::wrap_pyfunction;
    use pyo3::{create_exception, prelude::*};
    use serde_derive::{Deserialize, Serialize};
    use std::{cmp::Ordering, iter::Enumerate};
    use std::{
        fmt::Debug,
        io::{Read, Write},
    };

    type VecOfResults = Vec<CorResult>;
    type TupleExpressionValues = (String, Option<String>, Vec<f64>);
    type LazyMatrix = Box<dyn Iterator<Item = TupleExpressionValues>>;
    type CollectedMatrix = Vec<TupleExpressionValues>;

    create_exception!(ggca, GGCAError, pyo3::exceptions::PyException);
    create_exception!(ggca, GGCADiffSamplesLength, pyo3::exceptions::PyException);
    create_exception!(ggca, GGCADiffSamples, pyo3::exceptions::PyException);

    #[pyclass]
    #[derive(Clone, PartialEq, Serialize, Deserialize, Debug)]
    pub struct CorResult {
        #[pyo3(get, set)]
        gene: String,
        #[pyo3(get, set)]
        gem: String,
        #[pyo3(get, set)]
        cpg_site_id: Option<String>,
        #[pyo3(get, set)]
        correlation: Option<f64>,
        #[pyo3(get, set)]
        p_value: Option<f64>,
        #[pyo3(get, set)]
        adjusted_p_value: Option<f64>,
    }

    impl std::fmt::Display for CorResult {
        // This trait requires `fmt` with this exact signature.
        fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
            // Write strictly the first element into the supplied output
            // stream: `f`. Returns `fmt::Result` which indicates whether the
            // operation succeeded or failed. Note that `write!` uses syntax which
            // is very similar to `println!`.
            let cpg_site_id = match &self.cpg_site_id {
                Some(cpg) => cpg.clone(),
                None => String::from("-"),
            };

            write!(
                f,
                "Gene -> {} | GEM -> {} | CpG Site ID -> {}
                \tCor -> {}
                \tP-value -> {:+e}
                \tAdjusted p-value -> {:+e}",
                self.gene,
                self.gem,
                cpg_site_id,
                self.correlation.unwrap_or(0.0),
                self.p_value.unwrap_or(0.0),
                self.adjusted_p_value.unwrap_or(0.0)
            )
        }
    }

    impl Eq for CorResult {}

    impl Ord for CorResult {
        fn cmp(&self, other: &Self) -> Ordering {
            self.partial_cmp(&other).unwrap()
        }
    }

    impl PartialOrd for CorResult {
        fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
            self.p_value.partial_cmp(&other.p_value)
        }
    }

    impl Sortable for CorResult {
        fn encode<W: Write>(&self, writer: &mut W) {
            let serialized = bincode::serialize(self).unwrap();
            writer.write_all(&serialized[..]).unwrap();
        }

        fn decode<R: Read>(reader: &mut R) -> Option<Self> {
            bincode::deserialize_from(reader).ok()
        }
    }

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
            matrix_1: LazyMatrix,
            len_m1: u64,
            matrix_2: LazyMatrix,
            len_m2: u64,
            number_of_samples: usize,
            correlation_method: CorrelationMethod,
            correlation_threshold: f64,
            sort_buf_size: u64,
            adjustment_method: AdjustmentMethod,
            is_all_vs_all: bool,
        ) -> PyResult<(VecOfResults, u64)> {
            // Right part of iproduct must implement Clone. For more info read:
            // https://users.rust-lang.org/t/iterators-over-csv-files-with-iproduct/51947
            let matrix_2_collected = matrix_2.collect::<CollectedMatrix>();

            // Cartesian product computing correlation and p-value
            let correlation_method_struct =
                get_correlation_method(correlation_method, number_of_samples);
            let correlation_function = if is_all_vs_all {
                cartesian_all_vs_all
            } else {
                cartesian_equal_genes
            };
            let correlations_and_p_values =
                iproduct!(matrix_1, matrix_2_collected).map(|(tuple_1, tuple_2)| {
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

            // Sorting (for future adjustment)
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

            Ok((
                adjusted.collect::<VecOfResults>(),
                number_of_evaluated_combinations,
            ))
        }

        fn csv_builder(
            &self,
            path: &str,
        ) -> PyResult<(
            Vec<String>,
            Enumerate<csv::StringRecordsIntoIter<std::fs::File>>,
        )> {
            let reader_builder = ReaderBuilder::new().delimiter(b'\t').from_path(path);

            match reader_builder {
                Err(er) => Err(GGCAError::new_err(format!(
                    "The dataset '{}' has thrown an error: {}",
                    path, er
                ))),
                Ok(mut reader) => {
                    let headers: csv::StringRecord = reader.headers().unwrap().to_owned();
                    let headers = headers
                        .into_iter()
                        .map(|header| header.to_string())
                        .collect::<Vec<String>>();
                    let reader_iterator = reader.into_records().enumerate();
                    Ok((headers, reader_iterator))
                }
            }
        }

        fn get_df(&self, path: &str) -> PyResult<(Vec<String>, LazyMatrix)> {
            // Build the CSV reader and iterate over each record.
            let (headers, reader) = self.csv_builder(path)?;
            let dataframe_parsed = reader.map(|(row_idx, record_result)| {
                let record = record_result.unwrap();
                let mut it = record.into_iter();
                let gene_or_gem = it.next().unwrap().to_string();
                let values = it
                    .enumerate()
                    .map(|(column_idx, cell)| {
                        // It's added 1 as we're not considering headers row or index column
                        cell.parse::<f64>().expect(&format!(
                            "Row {} column {} (0 indexed) has an invalid value -> {}.
                                \nFirst column must be the Gene/GEM and the rest the samples",
                            row_idx + 1,
                            column_idx + 1,
                            cell
                        ))
                    })
                    .collect::<Vec<f64>>();

                (gene_or_gem, None, values)
            });

            Ok((headers, Box::new(dataframe_parsed)))
        }

        fn get_df_with_cpg(&self, path: &str) -> PyResult<(Vec<String>, LazyMatrix)> {
            // Build the CSV reader and iterate over each record.
            let (headers, reader) = self.csv_builder(path)?;
            let dataframe_parsed = reader
                .map(|(row_idx, record_result)| {
                    let record = record_result.unwrap();
                    let mut it = record.into_iter();
                    let gene_or_gem = it.next().unwrap().to_string();
                    let cpg_site_id = it.next().unwrap().to_string();

                    let values = it
                        .enumerate()
                        .map(|(column_idx, cell)| {
                            cell.parse::<f64>().expect(&format!(
                                "Row {} column {} (0 indexed) has an invalid value -> {}
                                \nFirst column must be the Gene, the second the CpG Site ID and the rest the samples",
                                row_idx + 1, // + 1 for the header row
                                column_idx + 2, // +2 for Gene and CpG columns
                                cell
                            ))
                        })
                        .collect::<Vec<f64>>();

                    (gene_or_gem, Some(cpg_site_id), values)
                });

            Ok((headers, Box::new(dataframe_parsed)))
        }

        /// Gets DataFrame and its shape (to compute, for example, the number of evaluated combinations)
        fn df_and_shape(
            &self,
            df1_path: &str,
            df2_path: &str,
            gem_contains_cpg: bool,
        ) -> PyResult<(LazyMatrix, u64, LazyMatrix, u64, usize)> {
            let (mut m1_headers, m1) = self.get_df(df1_path)?;
            let (_m1_aux_headers, mut m1_aux) = self.get_df(df1_path)?;
            let one_row = m1_aux.next();

            match one_row {
                None => Err(GGCAError::new_err(
                    "Error reading first row of mRNA file, is it empty?",
                )),
                Some(row) => {
                    let number_of_samples = row.2.len();

                    let len_m1 = m1_aux.count() + 1; // Plus discarded element by next()

                    let ((m2_headers, m2), m2_aux) = if !gem_contains_cpg {
                        (self.get_df(df2_path)?, self.get_df(df2_path)?)
                    } else {
                        (
                            self.get_df_with_cpg(df2_path)?,
                            self.get_df_with_cpg(df2_path)?,
                        )
                    };

                    // In case of an extra column for CpG Site IDs it shouldn't threw an error
                    // Index column will be remove after
                    let mut m2_headers_to_compare = m2_headers;
                    let mut advice_message = "";
                    if gem_contains_cpg {
                        m2_headers_to_compare.remove(1); // Removes CpG Site IDs column
                        advice_message = "(without taking the CpG Site IDs column into account)";
                    }

                    if m1_headers.len() != m2_headers_to_compare.len() {
                        Err(GGCADiffSamplesLength::new_err(format!(
                            "Length of samples in datasets are different:
                                \tSamples in mRNA dataset -> {}
                                \tSamples in GEM dataset -> {} {}",
                            m1_headers.len(),
                            m2_headers_to_compare.len(),
                            advice_message
                        )))
                    } else {
                        // Needs to remove index column as these are not samples and could have different names
                        m1_headers.remove(0);
                        m2_headers_to_compare.remove(0);

                        if m1_headers != m2_headers_to_compare {
                            Err(GGCADiffSamples::new_err(
                                "Samples in both datasets are different but they have the same length, maybe they are in different order"
                            ))
                        } else {
                            let len_m2 = m2_aux.1.count();

                            Ok((m1, len_m1 as u64, m2, len_m2 as u64, number_of_samples))
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
