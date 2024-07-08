use bincode::{deserialize, serialize};
use extsort::Sortable;
use pyo3::prelude::*;
use pyo3::types::PyBytes;
use pyo3::types::PyDict;
use pyo3::types::PyTuple;
use serde_derive::{Deserialize, Serialize};
use statrs::distribution::ContinuousCDF;
use statrs::distribution::Normal;
use statrs::distribution::StudentsT;
use std::cmp::Ordering;
use std::{
    fmt::Debug,
    io::{Read, Write},
};

/// Calculates the Pearson correlation coefficient between two arrays of f64 values.
/// Returns the correlation coefficient.
/// # Arguments
/// * `x` - Array of f64 values
/// * `y` - Array of f64 values
/// # Example
/// ```
/// use ggca::correlation::pearson_correlation;
///
/// let x = vec![1.0, 2.0, 3.0, 4.0, 5.0];
/// let y = vec![1.0, 2.0, 3.0, 4.0, 5.0];
/// let correlation = pearson_correlation(&x, &y);
///
/// assert_eq!(correlation, 1.0);
/// ```
pub fn pearson_correlation(x: &[f64], y: &[f64]) -> f64 {
    let n = x.len() as f64;
    let mut sum_x = 0.0;
    let mut sum_y = 0.0;
    let mut sum_xy = 0.0;
    let mut sum_x2 = 0.0;
    let mut sum_y2 = 0.0;

    for i in 0..x.len() {
        sum_x += x[i];
        sum_y += y[i];
        sum_xy += x[i] * y[i];
        sum_x2 += x[i] * x[i];
        sum_y2 += y[i] * y[i];
    }

    let numerator = sum_xy - (sum_x * sum_y / n);
    let denominator = ((sum_x2 - (sum_x * sum_x / n)) * (sum_y2 - (sum_y * sum_y / n))).sqrt();

    numerator / denominator
}

/// Calculates the Spearman correlation coefficient between two arrays of f64 values.
/// Returns the correlation coefficient.
///
/// # Arguments
/// * `x` - Array of f64 values
/// * `y` - Array of f64 values
///
/// # Example
/// ```
/// use ggca::correlation::spearman_correlation;
///
/// let x = vec![1.0, 2.0, 3.0, 4.0, 5.0];
/// let y = vec![1.0, 2.0, 3.0, 4.0, 5.0];
/// let correlation = spearman_correlation(&x, &y);
///
/// assert_eq!(correlation, 1.0);
/// ```
pub fn spearman_correlation(x: &[f64], y: &[f64]) -> f64 {
    // Calculate ranks for x and y
    let rank_x = rank_vector_avg(x);
    let rank_y = rank_vector_avg(y);

    // Calculate Pearson correlation between ranks https://en.wikipedia.org/wiki/Spearman%27s_rank_correlation_coefficient#Definition_and_calculation
    pearson_correlation(&rank_x, &rank_y)
}

/// Rank vector elements. Ties are assigned the average rank.
fn rank_vector_avg(v: &[f64]) -> Vec<f64> {
    let n = v.len();
    let mut indexed: Vec<(usize, &f64)> = v.iter().enumerate().collect();

    // Sort by value
    indexed.sort_by(|a, b| a.1.partial_cmp(b.1).unwrap_or(std::cmp::Ordering::Equal));

    let mut ranks = vec![0.0; n];
    let mut i = 0;
    while i < n {
        let mut j = i + 1;
        while j < n && (indexed[j].1 - indexed[i].1).abs() < f64::EPSILON {
            j += 1;
        }

        // Calculate average rank for tied elements
        let rank_sum: f64 = (i + 1..=j).sum::<usize>() as f64;
        let avg_rank = rank_sum / (j - i) as f64;

        // Assign average rank to all tied elements
        for k in i..j {
            ranks[indexed[k].0] = avg_rank;
        }

        i = j;
    }

    ranks
}

/// Represents an correlation analysis result. Includes Gene, GEM, CpG Site ID (if specified) correlation statistic,
/// p-value and adjusted p-value.
#[pyclass(module = "ggca")]
#[derive(Clone, PartialEq, Serialize, Deserialize, Debug)]
pub struct CorResult {
    /// Gene name
    #[pyo3(get, set)]
    pub gene: String,
    /// Gene Expression Modulator (GEM) name
    #[pyo3(get, set)]
    pub gem: String,
    #[pyo3(get, set)]
    /// CpG Site ID
    pub cpg_site_id: Option<String>,
    /// Correlation statistic (Pearson, Spearman or Kendall, as selected)
    #[pyo3(get, set)]
    pub correlation: Option<f64>,
    /// P-value
    #[pyo3(get, set)]
    pub p_value: Option<f64>,
    /// Adjusted p-value (Benjamini-Hochberg, Benjamini-Yekutieli or Bonferroni, as selected)
    #[pyo3(get, set)]
    pub adjusted_p_value: Option<f64>,
}

#[pymethods]
impl CorResult {
    #[new]
    #[pyo3(signature = (gene, gem, cpg_site_id=None, correlation=None, p_value=None, adjusted_p_value=None))]
    fn new(
        gene: String,
        gem: String,
        cpg_site_id: Option<String>,
        correlation: Option<f64>,
        p_value: Option<f64>,
        adjusted_p_value: Option<f64>,
    ) -> Self {
        CorResult {
            gene,
            gem,
            cpg_site_id,
            correlation,
            p_value,
            adjusted_p_value,
        }
    }

    // Adds support for pickle
    pub fn __setstate__(&mut self, py: Python, state: PyObject) -> PyResult<()> {
        match state.extract::<Bound<'_, PyTuple>>(py) {
            Ok(args) => {
                let gene_bytes = args
                    .get_item(0)
                    .unwrap()
                    .extract::<Bound<'_, PyBytes>>()
                    .unwrap();
                self.gene = deserialize(gene_bytes.as_bytes()).unwrap();

                let gem_bytes = args
                    .get_item(1)
                    .unwrap()
                    .extract::<Bound<'_, PyBytes>>()
                    .unwrap();
                self.gem = deserialize(gem_bytes.as_bytes()).unwrap();

                let cpg_site_id_bytes = args
                    .get_item(2)
                    .unwrap()
                    .extract::<Bound<'_, PyBytes>>()
                    .unwrap();
                self.cpg_site_id = deserialize(cpg_site_id_bytes.as_bytes()).unwrap();

                let correlation_bytes = args
                    .get_item(3)
                    .unwrap()
                    .extract::<Bound<'_, PyBytes>>()
                    .unwrap();
                self.correlation = deserialize(correlation_bytes.as_bytes()).unwrap();

                let p_value_bytes = args
                    .get_item(4)
                    .unwrap()
                    .extract::<Bound<'_, PyBytes>>()
                    .unwrap();
                self.p_value = deserialize(p_value_bytes.as_bytes()).unwrap();

                let adjusted_p_value_bytes = args
                    .get_item(5)
                    .unwrap()
                    .extract::<Bound<'_, PyBytes>>()
                    .unwrap();
                self.adjusted_p_value = deserialize(adjusted_p_value_bytes.as_bytes()).unwrap();

                Ok(())
            }
            Err(e) => Err(e),
        }
    }

    // Adds support for pickle
    pub fn __getstate__(&self, py: Python) -> PyResult<PyObject> {
        let obj = (
            PyBytes::new_bound(py, &serialize(&self.gene).unwrap()),
            PyBytes::new_bound(py, &serialize(&self.gem).unwrap()),
            PyBytes::new_bound(py, &serialize(&self.cpg_site_id).unwrap()),
            PyBytes::new_bound(py, &serialize(&self.correlation).unwrap()),
            PyBytes::new_bound(py, &serialize(&self.p_value).unwrap()),
            PyBytes::new_bound(py, &serialize(&self.adjusted_p_value).unwrap()),
        )
            .to_object(py);
        Ok(obj)
    }

    /// Returns the CpG Site ID description. Empty String if it's None
    fn cpg_site_id_description(&self) -> &str {
        match &self.cpg_site_id {
            Some(cpg) => cpg,
            None => "",
        }
    }

    /// Gets the absolute correlation. Panics if the correlation value is None
    pub fn abs_correlation(&self) -> f64 {
        self.correlation.unwrap().abs()
    }

    // Will be auto-generated by PyO3
    pub fn __str__(&self) -> String {
        format!("{}", self)
    }

    // Will be auto-generated by PyO3
    pub fn __repr__(&self) -> String {
        format!(
            r#"CorResult("{}", "{}", "{}", {}, {:+e}, {:+e})"#,
            self.gene,
            self.gem,
            self.cpg_site_id_description(),
            self.correlation.unwrap_or(0.0),
            self.p_value.unwrap_or(0.0),
            self.adjusted_p_value.unwrap_or(0.0)
        )
    }

    /// Returns a dictionary with the CorResult fields
    pub fn __dict__(&self, py: Python) -> PyResult<PyObject> {
        let dict = PyDict::new_bound(py);
        dict.set_item("gene", self.gene.clone())?;
        dict.set_item("gem", self.gem.clone())?;
        dict.set_item("cpg_site_id", self.cpg_site_id_description())?;
        dict.set_item("correlation", self.correlation.unwrap_or(0.0))?;
        dict.set_item("p_value", self.p_value.unwrap_or(0.0))?;
        dict.set_item("adjusted_p_value", self.adjusted_p_value.unwrap_or(0.0))?;
        Ok(dict.to_object(py))
    }
}

impl std::fmt::Display for CorResult {
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        write!(
            f,
            r#"Gene: "{}" | GEM: "{}" | CpG Site ID: "{}"
    Cor: {}
    P-value: {:+e}
    Adjusted p-value: {:+e}"#,
            self.gene,
            self.gem,
            self.cpg_site_id_description(),
            self.correlation.unwrap_or(0.0),
            self.p_value.unwrap_or(0.0),
            self.adjusted_p_value.unwrap_or(0.0)
        )
    }
}

impl Eq for CorResult {}

impl Sortable for CorResult {
    fn encode<W: Write>(&self, writer: &mut W) {
        let serialized = bincode::serialize(self).unwrap();
        writer.write_all(&serialized[..]).unwrap();
    }

    fn decode<R: Read>(reader: &mut R) -> Option<Self> {
        bincode::deserialize_from(reader).ok()
    }
}

pub trait Correlation {
    fn correlate(&self, x: &[f64], y: &[f64]) -> (f64, f64);
}

pub struct Pearson {
    degrees_of_freedom: f64,
}

impl Pearson {
    fn new(n: usize) -> Self {
        Pearson {
            degrees_of_freedom: (n - 2) as f64,
        }
    }
}

impl Correlation for Pearson {
    fn correlate(&self, x: &[f64], y: &[f64]) -> (f64, f64) {
        let r = pearson_correlation(x, y);

        // P-value (two-sided)
        // Based on R's cor.test method (https://github.com/SurajGupta/r-source/blob/a28e609e72ed7c47f6ddfbb86c85279a0750f0b7/src/library/stats/R/cor.test.R#L21)
        let statistic = self.degrees_of_freedom.sqrt() * r / (1.0 - r.powi(2)).sqrt();

        // In case of NaN, the p-value is also NaN to filter it out later
        if statistic.is_nan() {
            (r, f64::NAN)
        } else {
            let cdf = StudentsT::new(0.0, 1.0, self.degrees_of_freedom) // TODO: add instantiation of Normal struct in Pearson struct
                .unwrap()
                .cdf(statistic);
            let ccdf = 1.0 - cdf;
            let p_value = 2.0 * cdf.min(ccdf);

            (r, p_value)
        }
    }
}

struct Spearman {
    degrees_of_freedom: f64,
}

impl Spearman {
    fn new(n: usize) -> Self {
        Spearman {
            degrees_of_freedom: (n - 2) as f64,
        }
    }
}

impl Correlation for Spearman {
    fn correlate(&self, x: &[f64], y: &[f64]) -> (f64, f64) {
        let rs = spearman_correlation(x, y);

        // P-value (two-sided)
        // Same behavior as Python Scipy's spearmanr method
        let t = rs * (self.degrees_of_freedom / ((rs + 1.0) * (1.0 - rs))).sqrt();

        // In case of NaN, the p-value is also NaN to filter it out later
        if t.is_nan() {
            (rs, f64::NAN)
        } else {
            let cdf = StudentsT::new(0.0, 1.0, self.degrees_of_freedom) // TODO: add instantiation of Normal struct in Spearman struct
                .unwrap()
                .cdf(t.abs());
            let ccdf = 1.0 - cdf;
            let p_value = 2.0 * ccdf;

            (rs, p_value)
        }
    }
}

struct Kendall {}

impl Kendall {
    fn new(_n: usize) -> Self {
        Kendall {}
    }
}

impl Correlation for Kendall {
    fn correlate(&self, x: &[f64], y: &[f64]) -> (f64, f64) {
        // TODO: replace with new implementation
        let (tau, significance) = kendalls::tau_b_with_comparator(x, y, |a: &f64, b: &f64| {
            a.partial_cmp(b).unwrap_or(Ordering::Greater)
        })
        .unwrap();

        // P-value (two-sided)
        let cdf = Normal::new(0.0, 1.0).unwrap().cdf(-significance.abs()); // TODO: add instantiation of Normal struct in Kendall struct
        let p_value = 2.0 * cdf;

        (tau, p_value)
    }
}

#[pyclass(eq, eq_int)]
#[derive(Clone, Debug, PartialEq)]
pub enum CorrelationMethod {
    Spearman = 1,
    Kendall = 2,
    Pearson = 3,
}

impl std::fmt::Display for CorrelationMethod {
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        let description = match &self {
            CorrelationMethod::Spearman => "Spearman",
            CorrelationMethod::Kendall => "Kendall",
            CorrelationMethod::Pearson => "Pearson",
        };

        write!(f, "{description}")
    }
}

pub fn get_correlation_method(
    correlation_method: &CorrelationMethod,
    number_of_samples: usize,
) -> Box<dyn Correlation> {
    match correlation_method {
        CorrelationMethod::Pearson => Box::new(Pearson::new(number_of_samples)),
        CorrelationMethod::Spearman => Box::new(Spearman::new(number_of_samples)),
        CorrelationMethod::Kendall => Box::new(Kendall::new(number_of_samples)),
    }
}
