use crate::types::{LazyMatrixInner, TupleExpressionValues};
use csv::{Reader, ReaderBuilder};
use pyo3::create_exception;
use pyo3::exceptions::PyException;
use pyo3::PyResult;
use rayon::iter::ParallelBridge;
use rayon::iter::ParallelIterator;
use std::fs::File;

create_exception!(ggca, GGCAError, PyException);

/// Creates CSV Reader from a path with the Tab delimiter
/// # Args
/// * `path`: Path of file to create the Reader
fn reader_from_path(path: &str) -> PyResult<Reader<File>> {
    let reader_builder = ReaderBuilder::new()
        .buffer_capacity(16_384) // 16 KB
        .delimiter(b'\t')
        .from_path(path);

    match reader_builder {
        Err(er) => Err(GGCAError::new_err(format!(
            "The dataset '{}' has thrown an error: {}",
            path, er
        ))),
        Ok(reader) => Ok(reader),
    }
}

pub struct LazyMatrix {
    path: String,
    gem_contains_cpg: bool,
    inner: LazyMatrixInner,
}

impl LazyMatrix {
    /// Creates an instance of LazyMatrix from a path
    /// # Args
    /// * `path`: Path of file to create the LazyMatrix
    /// * `gem_contains_cpg`: True if second column of GEM dataset must be considered as a CpG Site ID. False for normal datasets
    fn new(path: &str, gem_contains_cpg: bool) -> PyResult<Self> {
        let lazy_matrix = Self::lazy_matrix(path, gem_contains_cpg)?;

        Ok(LazyMatrix {
            path: path.to_string(),
            gem_contains_cpg,
            inner: lazy_matrix,
        })
    }

    /// Parses the datasets cells (strings) to float using efficient implementation
    /// from [fast-float](https://github.com/aldanor/fast-float-rust) crate.
    /// Returns an Map iterator with the values of Gene/GEM, CpG Site ID (optional) and a vec with expression values
    /// In case of error, panics informing the line and column with invalid format.
    /// # Args
    /// * `path`: Path of the dataset
    /// * `gem_contains_cpg`: True if second column of GEM dataset must be considered as a CpG Site ID. False for normal datasets
    fn lazy_matrix(path: &str, with_cpg_site_id: bool) -> PyResult<LazyMatrixInner> {
        // Build the CSV reader and iterate over each record.
        let reader = reader_from_path(path)?;
        let res_lazy_matrix = reader
            .into_records()
            .par_bridge()
            .map(move |record_result| {
                // Gets current record and its line
                let record = record_result.unwrap();
                let line = record.position().unwrap().line();
                let mut cells_it = record.into_iter();

                // Gets Gene/GEM and, if needed, CpG Site ID
                let gene_or_gem = cells_it.next().unwrap().to_string();
                let cpg_site_id = if with_cpg_site_id {
                    Some(cells_it.next().unwrap().to_string())
                } else {
                    None
                };

                // Casts all cells to float
                let values = cells_it
                    .enumerate()
                    .map(|(column_idx, cell)| {
                        // + 1 or +2 as it doesn't take index (and CpG Site ID) column/s into account
                        fast_float::parse(cell)
                            // Avoids using expect() to prevent calling format()
                            .unwrap_or_else(|_| {
                                panic!(
                                    "Line {} column {} has an invalid value -> '{}'.
                        \nFirst column must be the Gene/GEM and the rest the samples",
                                    line,
                                    column_idx + (if with_cpg_site_id { 2 } else { 1 }),
                                    cell
                                )
                            })
                    })
                    .collect::<Vec<f64>>();

                (gene_or_gem, cpg_site_id, values)
            });

        let aux = res_lazy_matrix.collect::<Vec<_>>();
        Ok(Box::new(aux.into_iter()))
    }
}

impl Iterator for LazyMatrix {
    type Item = TupleExpressionValues;
    fn next(&mut self) -> Option<Self::Item> {
        self.inner.next()
    }
}

impl Clone for LazyMatrix {
    fn clone(&self) -> Self {
        Self::new(self.path.as_str(), self.gem_contains_cpg).unwrap()
    }
}

pub struct Dataset {
    pub headers: Vec<String>,
    pub lazy_matrix: LazyMatrix,
}

impl Dataset {
    /// Creates an instance of Dataset from a path
    /// # Args
    /// * `path`: Path of file to create the Dataset
    /// * `gem_contains_cpg`: True if second column of GEM dataset must be considered as a CpG Site ID. False for normal datasets
    pub fn new(path: &str, gem_contains_cpg: bool) -> PyResult<Self> {
        let mut reader = reader_from_path(path)?;
        let headers = Self::headers_from_reader(&mut reader);

        Ok(Dataset {
            headers,
            lazy_matrix: LazyMatrix::new(path, gem_contains_cpg)?,
        })
    }

    /// Gets Dataset headers from its reader
    /// # Args
    /// * `reader`: CSV Reader to extract headers
    fn headers_from_reader(reader: &mut Reader<File>) -> Vec<String> {
        let headers: csv::StringRecord = reader.headers().unwrap().to_owned();
        headers
            .into_iter()
            .map(|header| header.to_string())
            .collect::<Vec<String>>()
    }
}
