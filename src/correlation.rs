use bincode::{deserialize, serialize};
use distrs::StudentsT;
use extsort::Sortable;
use ordered_float::OrderedFloat;
use pyo3::prelude::*;
use pyo3::types::{PyBytes, PyTuple};
use pyo3::ToPyObject;
use serde_derive::{Deserialize, Serialize};
use std::cmp::Ordering;
use std::collections::BTreeMap;
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
fn pearson_correlation(x: &[f64], y: &[f64]) -> f64 {
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
fn spearman_correlation(x: &[f64], y: &[f64]) -> f64 {
    let n = x.len() as f64;
    let mut rank_x: Vec<f64> = x.iter().map(|v| *v).collect();
    let mut rank_y: Vec<f64> = y.iter().map(|v| *v).collect();

    rank_x.sort_by(|a, b| a.partial_cmp(b).unwrap());
    rank_y.sort_by(|a, b| a.partial_cmp(b).unwrap());

    let mut rank_map_x: BTreeMap<OrderedFloat<f64>, f64> = BTreeMap::new();
    let mut rank_map_y: BTreeMap<OrderedFloat<f64>, f64> = BTreeMap::new();

    for i in 0..n as usize {
        rank_map_x.insert(OrderedFloat(rank_x[i]), (i + 1) as f64);
        rank_map_y.insert(OrderedFloat(rank_y[i]), (i + 1) as f64);
    }

    let mut sum_d2 = 0.0;
    for i in 0..n as usize {
        let d = rank_map_x[&OrderedFloat(x[i])] - rank_map_y[&OrderedFloat(y[i])];
        sum_d2 += d * d;
    }

    1.0 - (6.0 * sum_d2) / (n * (n * n - 1.0))
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

use std::f64::consts::SQRT_2;

fn gaussian_p(x: f64, sigma: f64) -> f64 {
    let mu = 0.0;
    if sigma == 0.0 {
        return if x < mu { 0.0 } else { 1.0 };
    }

    let t = (x - mu) / (sigma * SQRT_2);
    let y = erfc(-t);
    y / 2.0
}

fn erfc(x: f64) -> f64 {
    let mut sum = 0.0;
    let mut term = 1.0;
    let mut n = 1;

    loop {
        sum += term;
        n += 2;
        term *= -x * x / (n - 1) as f64;

        if term.abs() < 1e-15 * sum.abs() {
            break;
        }
    }

    2.0 * std::f64::consts::E.powf(-x * x) / SQRT_2 / std::f64::consts::PI * sum
}

#[pymethods]
impl CorResult {
    #[new]
    #[pyo3(signature = (*args))]
    fn new(args: &PyTuple) -> Self {
        if args.len() >= 2 {
            CorResult {
                gene: args.get_item(0).unwrap().extract::<String>().unwrap(),
                gem: args.get_item(1).unwrap().extract::<String>().unwrap(),
                cpg_site_id: args
                    .get_item(2)
                    .unwrap()
                    .extract::<Option<String>>()
                    .unwrap(),
                correlation: args.get_item(3).unwrap().extract::<Option<f64>>().unwrap(),
                p_value: args.get_item(4).unwrap().extract::<Option<f64>>().unwrap(),
                adjusted_p_value: args.get_item(5).unwrap().extract::<Option<f64>>().unwrap(),
            }
        } else {
            CorResult {
                gene: String::from(""),
                gem: String::from(""),
                cpg_site_id: None,
                correlation: None,
                p_value: None,
                adjusted_p_value: None,
            }
        }
    }

    // Adds support for pickle
    pub fn __setstate__(&mut self, py: Python, state: PyObject) -> PyResult<()> {
        match state.extract::<&PyTuple>(py) {
            Ok(args) => {
                let gene_bytes = args.get_item(0).unwrap().extract::<&PyBytes>().unwrap();
                self.gene = deserialize(gene_bytes.as_bytes()).unwrap();
                let gem_bytes = args.get_item(1).unwrap().extract::<&PyBytes>().unwrap();
                self.gem = deserialize(gem_bytes.as_bytes()).unwrap();
                let cpg_site_id_bytes = args.get_item(2).unwrap().extract::<&PyBytes>().unwrap();
                self.cpg_site_id = deserialize(cpg_site_id_bytes.as_bytes()).unwrap();
                let correlation_bytes = args.get_item(3).unwrap().extract::<&PyBytes>().unwrap();
                self.correlation = deserialize(correlation_bytes.as_bytes()).unwrap();
                let p_value_bytes = args.get_item(4).unwrap().extract::<&PyBytes>().unwrap();
                self.p_value = deserialize(p_value_bytes.as_bytes()).unwrap();
                let adjusted_p_value_bytes =
                    args.get_item(5).unwrap().extract::<&PyBytes>().unwrap();
                self.adjusted_p_value = deserialize(adjusted_p_value_bytes.as_bytes()).unwrap();
                Ok(())
            }
            Err(e) => Err(e),
        }
    }

    // Adds support for pickle
    pub fn __getstate__(&self, py: Python) -> PyResult<PyObject> {
        let obj = (
            PyBytes::new(py, &serialize(&self.gene).unwrap()),
            PyBytes::new(py, &serialize(&self.gem).unwrap()),
            PyBytes::new(py, &serialize(&self.cpg_site_id).unwrap()),
            PyBytes::new(py, &serialize(&self.correlation).unwrap()),
            PyBytes::new(py, &serialize(&self.p_value).unwrap()),
            PyBytes::new(py, &serialize(&self.adjusted_p_value).unwrap()),
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
    n: usize,
    degrees_of_freedom: f64,
}

impl Pearson {
    fn new(n: usize) -> Self {
        Pearson {
            n,
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
        
        let cdf = StudentsT::cdf(statistic, self.degrees_of_freedom);
        let ccdf = 1.0 - cdf;
        let p_value = 2.0 * cdf.min(ccdf);

        (r, p_value)
    }
}

struct Spearman {
    n: usize,
    degrees_of_freedom: f64,
}

impl Spearman {
    fn new(n: usize) -> Self {
        Spearman {
            n,
            degrees_of_freedom: (n - 2) as f64,
        }
    }
}

impl Correlation for Spearman {
    fn correlate(&self, x: &[f64], y: &[f64]) -> (f64, f64) {
        let rs = spearman_correlation(x, y);

        // P-value (two-sided)
        // Same behavior as Python Scipy's spearmanr method
        // let t = rs * (self.degrees_of_freedom / ((rs + 1.0) * (1.0 - rs))).sqrt();
        // let ccdf = tdist_Q(t.abs(), self.degrees_of_freedom);
        // let p_value = 2.0 * ccdf;

        let t = rs * (self.degrees_of_freedom / ((rs + 1.0) * (1.0 - rs))).sqrt();
        let ccdf = 1.0 - StudentsT::cdf(t.abs(), self.degrees_of_freedom);
        let p_value = 2.0 * ccdf;

        (rs, p_value)
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
        let (tau, significance) = kendalls::tau_b_with_comparator(x, y, |a: &f64, b: &f64| {
            a.partial_cmp(b).unwrap_or(Ordering::Greater)
        })
        .unwrap();

        // P-value (two-sided)
        let cdf = gaussian_p(-significance.abs(), 1.0);
        let p_value = 2.0 * cdf;

        (tau, p_value)
    }
}

#[derive(Clone, Debug)]
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
