pub mod adjustment;
pub mod correlation;

pub mod experiment {
    use pyo3::prelude::*;
    use pyo3::wrap_pyfunction;
    extern crate external_sort;
    use crate::adjustment::{get_adjustment_method, AdjustmentMethod};
    use crate::correlation::{get_correlation_method, CorrelationMethod};
    use csv::ReaderBuilder;
    use external_sort::{ExternalSorter, ExternallySortable};
    use itertools::iproduct;
    use serde::{Deserialize, Serialize};
    use std::cmp::Ordering;

    type TupleExpressionValues = (String, Vec<f64>);
    pub type Batch = Vec<TupleExpressionValues>;
    type LazyMatrix = Box<dyn Iterator<Item = TupleExpressionValues>>;

    #[derive(Clone, PartialEq, Serialize, Deserialize, Debug)]
    struct CorResult {
        gene: String,
        gem: String,
        r: f32,
        p_value: f64,
        p_value_adjusted: Option<f64>,
    }

    impl Eq for CorResult {}

    impl Ord for CorResult {
        // Sorts in descending order
        fn cmp(&self, other: &Self) -> Ordering {
            self.partial_cmp(&other).unwrap()
        }
    }

    impl PartialOrd for CorResult {
        fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
            self.p_value.partial_cmp(&other.p_value)
        }
    }

    impl ExternallySortable for CorResult {
        fn get_size(&self) -> u64 {
            1
        }
    }

    pub trait Computation {
        fn all_vs_all(
            &self,
            m1: LazyMatrix,
            len_m1: u64,
            m3: LazyMatrix,
            len_m3: u64,
            number_of_columns: usize,
            correlation_method: CorrelationMethod,
            correlation_threhold: f32,
            sort_chunk_size: u64,
            adjustment_method: AdjustmentMethod,
        ) {
            let total_number_of_elements: u64 = len_m1 * len_m3;

            // We need a collected object for right-side of the iproduct macro
            // TODO: make dinamic checking for m1 and m3 length
            let m3_collected = m3.collect::<Vec<TupleExpressionValues>>();

            let correlation_struct = get_correlation_method(correlation_method, number_of_columns);
            let correlations_and_p_values = iproduct!(m1, m3_collected).map(|(tuple1, tuple3)| {
                // Gene and GEM
                let gene = tuple1.0;
                let gem = tuple3.0;

                let (r, p_value) =
                    correlation_struct.correlate(tuple1.1.as_slice(), tuple3.1.as_slice());

                CorResult {
                    gene,
                    gem,
                    r: r as f32,
                    p_value,
                    p_value_adjusted: None,
                }
            });

            // Sorting
            let external_sorter: ExternalSorter<CorResult> =
                ExternalSorter::new(sort_chunk_size, None);
            let sorted = external_sorter.sort(correlations_and_p_values).unwrap();

            // Ranking
            let ranked = sorted.enumerate();

            // Filtering
            let filtered = ranked.filter(|(_, cor_and_p_value)| {
                cor_and_p_value.as_ref().unwrap().r.abs() >= correlation_threhold
            });

            // Adjustment
            let mut adjustment_struct =
                get_adjustment_method(adjustment_method, total_number_of_elements as f64);
            let adjusted = filtered.map(|(rank, mut cor_and_p_value)| {
                let p_value = cor_and_p_value.as_ref().unwrap().p_value;
                let q_value = adjustment_struct.adjust(p_value, rank);

                cor_and_p_value.as_mut().unwrap().p_value_adjusted = Some(q_value);
                cor_and_p_value
            });

            let mut number_of_result_elements = 0;
            for elem in adjusted {
                let valid_elem = elem.unwrap();
                println!(
                    "{} x {} -> Cor: {} | p-value: {:+e} | adjusted p-value {:+e}",
                    valid_elem.gene,
                    valid_elem.gem,
                    valid_elem.r,
                    valid_elem.p_value,
                    valid_elem.p_value_adjusted.unwrap()
                );
                number_of_result_elements += 1;
            }

            println!("Cantidad final de datos -> {}", number_of_result_elements);
        }

        fn get_df(&self, path: &str) -> LazyMatrix {
            // Build the CSV reader and iterate over each record.
            let rdr = ReaderBuilder::new()
                .delimiter(b'\t')
                .from_path(path)
                .unwrap();

            let dataframe_parsed =
                rdr.into_records()
                    .enumerate()
                    .map(|(row_idx, record_result)| {
                        let record = record_result.unwrap();
                        let mut it = record.into_iter();
                        let gene_or_gem = it.next().unwrap().to_string();
                        let values = it
                            .enumerate()
                            .map(|(column_idx, cell)| {
                                // It's added 2 as we're not considering headers row or index column
                                cell.parse::<f64>().expect(&format!(
                                    "Row {} column {} has an invalid value",
                                    row_idx + 2,
                                    column_idx + 2
                                ))
                            })
                            .collect::<Vec<f64>>();

                        (gene_or_gem, values)
                    });

            Box::new(dataframe_parsed)
        }

        fn get_both_df_and_shape(
            &self,
            df1_path: &str,
            df2_path: &str,
        ) -> (LazyMatrix, u64, LazyMatrix, u64, usize) {
            let m1 = self.get_df(df1_path);
            let mut m1_aux = self.get_df(df1_path);
            let number_of_columns = m1_aux.next().unwrap().1.len();
            let len_m1 = m1_aux.count() + 1; // Plus discarded element by next()

            let m3 = self.get_df(df2_path);
            let m3_aux = self.get_df(df2_path);
            let len_m3 = m3_aux.count();

            (m1, len_m1 as u64, m3, len_m3 as u64, number_of_columns)
        }

        fn compute(
            &self,
            correlation_method: CorrelationMethod,
            correlation_threhold: f32,
            sort_chunk_size: u64,
            adjustment_method: AdjustmentMethod,
        );
    }

    pub struct ExperimentFromFiles {
        file_1_path: String,
        file_2_path: String,
    }

    impl Computation for ExperimentFromFiles {
        fn compute(
            &self,
            correlation_method: CorrelationMethod,
            correlation_threhold: f32,
            sort_chunk_size: u64,
            adjustment_method: AdjustmentMethod,
        ) {
            let (m1, len_m1, m3, len_m3, number_of_columns) =
                self.get_both_df_and_shape(self.file_1_path.as_str(), self.file_2_path.as_str());

            self.all_vs_all(
                m1,
                len_m1,
                m3,
                len_m3,
                number_of_columns,
                correlation_method,
                correlation_threhold,
                sort_chunk_size,
                adjustment_method,
            )
        }
    }

    pub fn new_from_files(file_1_path: String, file_2_path: String) -> ExperimentFromFiles {
        ExperimentFromFiles {
            file_1_path,
            file_2_path,
        }
    }

    /// A Python module implemented in Rust.
    #[pymodule]
    fn pruebas_correlation(_py: Python, m: &PyModule) -> PyResult<()> {
        m.add_function(wrap_pyfunction!(prueba, m)?)?;

        Ok(())
    }

    #[pyfunction]
    pub fn prueba(
        file_1_path: String,
        file_2_path: String,
        correlation_method: i32,
        correlation_threhold: f32,
        sort_chunk_size: u64,
        adjustment_method: i32,
    ) -> PyResult<()> {
        let experiment = new_from_files(file_1_path, file_2_path);
        let correlation_method = match correlation_method {
            1 => CorrelationMethod::Spearman,
            2 => CorrelationMethod::Kendall,
            _ => CorrelationMethod::Pearson
        };

        let adjustment_method = match adjustment_method {
            1 => AdjustmentMethod::BenjaminiHochberg,
            2 => AdjustmentMethod::BenjaminiYekutieli,
            _ => AdjustmentMethod::Bonferroni
        };

        experiment.compute(correlation_method, correlation_threhold, sort_chunk_size, adjustment_method);
        Ok(())
    }
}