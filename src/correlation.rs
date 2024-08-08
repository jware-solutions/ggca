use bincode::{deserialize, serialize};
use extsort::Sortable;
use itertools::Itertools;
use pyo3::prelude::*;
use pyo3::types::PyBytes;
use pyo3::types::PyDict;
use pyo3::types::PyTuple;
use rayon::slice::ParallelSliceMut;
use serde_derive::{Deserialize, Serialize};
use statrs::distribution::ContinuousCDF;
use statrs::distribution::Normal;
use statrs::distribution::StudentsT;
use std::cmp::Ordering;
use std::{
    fmt::Debug,
    io::{Read, Write},
};

fn pairs_comparator(a: &f64, b: &f64) -> Ordering {
    a.partial_cmp(b).unwrap_or(Ordering::Greater)
}

/// It's the same as [Kendalls crate](https://github.com/zolkko/kendalls/) Tau-b method. But here we apply a
/// parallel sort to improve performance.
#[allow(clippy::many_single_char_names)]
pub fn tau_b_with_comparator(x: &[f64], y: &[f64]) -> (f64, f64) {
    let n = x.len();

    let mut pairs = x.iter().cloned().zip(y.iter().cloned()).collect_vec();

    pairs.par_sort_by(|pair1, pair2| {
        let res = pairs_comparator(&pair1.0, &pair2.0);
        if res == Ordering::Equal {
            pairs_comparator(&pair1.1, &pair2.1)
        } else {
            res
        }
    });

    let mut v1_part_1 = 0usize;
    let mut v2_part_1 = 0isize;

    let mut tied_x_pairs = 0usize;
    let mut tied_xy_pairs = 0usize;
    let mut vt = 0usize;
    let mut consecutive_x_ties = 1usize;
    let mut consecutive_xy_ties = 1usize;

    for i in 1..n {
        let prev = &pairs[i - 1];
        let curr = &pairs[i];
        if curr.0 == prev.0 {
            consecutive_x_ties += 1;
            if curr.1 == prev.1 {
                consecutive_xy_ties += 1;
            } else {
                tied_xy_pairs += sum(consecutive_xy_ties - 1);
                consecutive_xy_ties = 1;
            }
        } else {
            update_x_group(
                &mut vt,
                &mut tied_x_pairs,
                &mut tied_xy_pairs,
                &mut v1_part_1,
                &mut v2_part_1,
                consecutive_x_ties,
                consecutive_xy_ties,
            );
            consecutive_x_ties = 1;
            consecutive_xy_ties = 1;
        }
    }

    update_x_group(
        &mut vt,
        &mut tied_x_pairs,
        &mut tied_xy_pairs,
        &mut v1_part_1,
        &mut v2_part_1,
        consecutive_x_ties,
        consecutive_xy_ties,
    );

    let mut swaps = 0usize;
    let mut pairs_dest: Vec<(f64, f64)> = vec![(Default::default(), Default::default()); n];

    let mut segment_size = 1usize;
    while segment_size < n {
        for offset in (0..n).step_by(2 * segment_size) {
            let mut i = offset;
            let i_end = n.min(i + segment_size);
            let mut j = i_end;
            let j_end = n.min(j + segment_size);
            let mut copy_location = offset;

            while i < i_end && j < j_end {
                let a = &pairs[i].1;
                let b = &pairs[j].1;

                if pairs_comparator(a, b) == Ordering::Greater {
                    pairs_dest[copy_location] = pairs[j];
                    j += 1;
                    swaps += i_end - i;
                } else {
                    pairs_dest[copy_location] = pairs[i];
                    i += 1;
                }

                copy_location += 1;
            }

            while i < i_end {
                pairs_dest[copy_location] = pairs[i];
                i += 1;
                copy_location += 1
            }

            while j < j_end {
                pairs_dest[copy_location] = pairs[j];
                j += 1;
                copy_location += 1
            }
        }

        std::mem::swap(&mut pairs, &mut pairs_dest);

        segment_size <<= 1;
    }

    let mut v1_part_2 = 0usize;
    let mut v2_part_2 = 0isize;
    let mut tied_y_pairs = 0usize;
    let mut consecutive_y_ties = 1usize;
    let mut vu = 0usize;

    for j in 1..n {
        let prev = &pairs[j - 1];
        let curr = &pairs[j];
        if curr.1 == prev.1 {
            consecutive_y_ties += 1;
        } else {
            update_y_group(
                &mut vu,
                &mut tied_y_pairs,
                &mut v1_part_2,
                &mut v2_part_2,
                consecutive_y_ties,
            );
            consecutive_y_ties = 1;
        }
    }

    update_y_group(
        &mut vu,
        &mut tied_y_pairs,
        &mut v1_part_2,
        &mut v2_part_2,
        consecutive_y_ties,
    );

    // Generates T1 and T2 for significance
    let v1 = (v1_part_1 * v1_part_2) as f64;
    let v2 = (v2_part_1 * v2_part_2) as f64;

    // Prevents overflow on subtraction
    let num_pairs_f: f64 = ((n * (n - 1)) as f64) / 2.0; // sum(n - 1).as_();
    let tied_x_pairs_f: f64 = tied_x_pairs as f64;
    let tied_y_pairs_f: f64 = tied_y_pairs as f64;
    let tied_xy_pairs_f: f64 = tied_xy_pairs as f64;
    let swaps_f: f64 = (2 * swaps) as f64;

    // Note that tot = con + dis + (xtie - ntie) + (ytie - ntie) + ntie
    //               = con + dis + xtie + ytie - ntie
    //
    //           C-D = tot - xtie - ytie + ntie - 2 * dis
    let concordant_minus_discordant =
        num_pairs_f - tied_x_pairs_f - tied_y_pairs_f + tied_xy_pairs_f - swaps_f;

    // non_tied_pairs_multiplied = ((n0 - n1) * (n0 - n2)).sqrt()
    let non_tied_pairs_multiplied = (num_pairs_f - tied_x_pairs_f) * (num_pairs_f - tied_y_pairs_f);

    let tau_b = concordant_minus_discordant / non_tied_pairs_multiplied.sqrt();

    // Significance
    let v0 = (n * (n - 1)) * (2 * n + 5);
    let n_f = n as f64;

    let v0_isize = v0 as isize;
    let vt_isize = vt as isize;
    let vu_isize = vu as isize;
    let var_s = (v0_isize - vt_isize - vu_isize) as f64 / 18.0
        + v1 / (2.0 * n_f * (n_f - 1.0))
        + v2 / (9.0 * n_f * (n_f - 1.0) * (n_f - 2.0));

    let s = tau_b * non_tied_pairs_multiplied.sqrt();
    let z = s / var_s.sqrt();

    // Limit range to fix computational errors
    (tau_b.clamp(-1.0, 1.0), z)
}

#[inline]
fn sum(n: usize) -> usize {
    n * (n + 1_usize) / 2_usize
}

/// Updated vt, v1_part_1, v2_part_1, tied_x_pairs, tied_xy_pairs variables with current tied group in X
fn update_x_group(
    vt: &mut usize,
    tied_x_pairs: &mut usize,
    tied_xy_pairs: &mut usize,
    v1_part_1: &mut usize,
    v2_part_1: &mut isize,
    consecutive_x_ties: usize,
    consecutive_xy_ties: usize,
) {
    *vt += consecutive_x_ties * (consecutive_x_ties - 1) * (2 * consecutive_x_ties + 5);
    *v1_part_1 += consecutive_x_ties * (consecutive_x_ties - 1);

    let consecutive_x_ties_i = consecutive_x_ties as isize;
    *v2_part_1 += consecutive_x_ties_i * (consecutive_x_ties_i - 1) * (consecutive_x_ties_i - 2);

    *tied_x_pairs += sum(consecutive_x_ties - 1);
    *tied_xy_pairs += sum(consecutive_xy_ties - 1);
}

/// Updated vu, tied_y_pairs, v1_part_2 and v2_part_2 variables with current tied group in Y
fn update_y_group(
    vu: &mut usize,
    tied_y_pairs: &mut usize,
    v1_part_2: &mut usize,
    v2_part_2: &mut isize,
    consecutive_y_ties: usize,
) {
    *vu += consecutive_y_ties * (consecutive_y_ties - 1) * (2 * consecutive_y_ties + 5);
    *v1_part_2 += consecutive_y_ties * (consecutive_y_ties - 1);

    let consecutive_y_ties_i = consecutive_y_ties as isize;
    *v2_part_2 += consecutive_y_ties_i * (consecutive_y_ties_i - 1) * (consecutive_y_ties_i - 2);

    *tied_y_pairs += sum(consecutive_y_ties - 1);
}

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

pub trait Correlation: Sync {
    fn correlate(&self, x: &[f64], y: &[f64]) -> (f64, f64);
}

pub struct Pearson {
    degrees_of_freedom: f64,
    /// Student's T distribution already instantiated for p-value calculation
    student_distribution: StudentsT,
}

impl Pearson {
    fn new(n: usize) -> Self {
        let degrees_of_freedom = (n - 2) as f64;
        Pearson {
            degrees_of_freedom,
            student_distribution: StudentsT::new(0.0, 1.0, degrees_of_freedom).unwrap(),
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
            let cdf = self.student_distribution.cdf(statistic);
            let ccdf = 1.0 - cdf;
            let p_value = 2.0 * cdf.min(ccdf);

            (r, p_value)
        }
    }
}

struct Spearman {
    degrees_of_freedom: f64,
    /// Student's T distribution already instantiated for p-value calculation
    student_distribution: StudentsT,
}

impl Spearman {
    fn new(n: usize) -> Self {
        let degrees_of_freedom = (n - 2) as f64;
        Spearman {
            degrees_of_freedom,
            student_distribution: StudentsT::new(0.0, 1.0, degrees_of_freedom).unwrap(),
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
            let cdf = self.student_distribution.cdf(t.abs());
            let ccdf = 1.0 - cdf;
            let p_value = 2.0 * ccdf;

            (rs, p_value)
        }
    }
}

// FIXME: this is not passing tests when --release is set
struct Kendall {
    /// Normal distribution already instantiated for p-value calculation
    normal_distribution: Normal,
}

impl Kendall {
    fn new(_n: usize) -> Self {
        Kendall {
            normal_distribution: Normal::new(0.0, 1.0).unwrap(),
        }
    }
}

impl Correlation for Kendall {
    fn correlate(&self, x: &[f64], y: &[f64]) -> (f64, f64) {
        let (tau, significance) = tau_b_with_comparator(x, y);

        // P-value (two-sided)
        let cdf = self.normal_distribution.cdf(-significance.abs());
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
