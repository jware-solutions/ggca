use std::fs::File;
use csv::{Reader, ReaderBuilder};
use pyo3::PyResult;
use crate::types::{LazyMatrixInner, TupleExpressionValues};
use pyo3::create_exception;

create_exception!(ggca, GGCAError, pyo3::exceptions::PyException);


fn reader_from_path(path: &str) -> PyResult<Reader<File>> {
    let reader_builder = ReaderBuilder::new().delimiter(b'\t').from_path(path);

    match reader_builder {
        Err(er) => Err(GGCAError::new_err(format!(
            "The dataset '{}' has thrown an error: {}",
            path, er
        ))),
        Ok(reader) => Ok(reader)
    }
}


fn headers_from_reader(reader: &mut Reader<File>) -> Vec<String> {
    let headers: csv::StringRecord = reader.headers().unwrap().to_owned();
    headers
        .into_iter()
        .map(|header| header.to_string())
        .collect::<Vec<String>>()

    // let string_record_wrapper = StringRecordsWrapper {
    //     filename: path.to_string(),
    //     inner: reader.into_records()
    // };
    // let string_record_iterator = string_record_wrapper.enumerate();

    // (headers, reader.into_records().enumerate())
}

pub struct LazyMatrix {
    path: String,
    gem_contains_cpg: bool,
    inner: LazyMatrixInner
}

impl LazyMatrix {
    pub fn new(path: &str, gem_contains_cpg: bool) -> PyResult<Self> {
        let lazy_matrix = if gem_contains_cpg {
            Self::get_df_with_cpg(path)?
        } else {
            Self::get_df(path)?
        };

        Ok(LazyMatrix {
            path: path.to_string(),
            gem_contains_cpg,
            inner: lazy_matrix
        })
    }


    fn get_df(path: &str) -> PyResult<LazyMatrixInner> {
            // Build the CSV reader and iterate over each record.
            // println!("Entra en get_df");
            let reader = reader_from_path(path)?;
            // let (headers, reader) = Self::headers_from_reader(path)?;
            let dataframe_parsed = reader.into_records().enumerate().map(|(row_idx, record_result)| {
                let record = record_result.unwrap();
                let mut it = record.into_iter();
                let gene_or_gem = it.next().unwrap().to_string();
                let lazy_matrix = it
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

                (gene_or_gem, None, lazy_matrix)
            });

            // println!("Sale de get_df");
            Ok(Box::new(dataframe_parsed))
        }

        fn get_df_with_cpg(path: &str) -> PyResult<LazyMatrixInner> {
            // Build the CSV reader and iterate over each record.
            // println!("Entra en get_df_with_cpg");
            let reader = reader_from_path(path)?;
            // let (headers, reader) = Self::headers_from_reader(path)?;
            let dataframe_parsed = reader.into_records().enumerate()
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
            // println!("Sale de get_df_with_cpg");
            Ok(Box::new(dataframe_parsed))
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
        // let mut new_inner: LazyMatrixInner = /*  */;
        // LazyMatrix {
        //     inner: new_inner
        // }
    }
}

// struct StringRecordsWrapper {
//     filename: String,
//     inner: csv::StringRecordsIntoIter<std::fs::File>,
// }

// impl Iterator for StringRecordsWrapper {
//     type Item = Result<csv::StringRecord, csv::Error>;
//     fn next(&mut self) -> Option<Self::Item> {
//         self.inner.next()
//     }
// }

// impl Clone for StringRecordsWrapper {
//     fn clone(&self) -> Self {
//         let new_reader = Dataset::builder(self.filename.as_str()).unwrap();
//         let mut new_inner: csv::StringRecordsIntoIter<std::fs::File> = new_reader.into_records();
//         new_inner.reader_mut().seek(*self.inner.reader().position());
//         StringRecordsWrapper {
//             filename: self.filename.clone(),
//             inner: new_inner
//         }
//     }
// }

pub struct Dataset {
    pub headers: Vec<String>,
    // file: csv::StringRecordsIntoIter<std::fs::File>,
    pub lazy_matrix: LazyMatrix
}

impl Dataset {
    pub fn new(path: &str, gem_contains_cpg: bool) -> PyResult<Self> {
        let mut reader = reader_from_path(path)?;
        let headers = headers_from_reader(&mut reader);

        Ok(Dataset {
            headers,
            lazy_matrix: LazyMatrix::new(path, gem_contains_cpg)?
        })
    }
}