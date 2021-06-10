use pyo3::prelude::*;
use rgsl::{randist::t_distribution::{tdist_P, tdist_Q}, statistics::{correlation, spearman}};
use serde_derive::{Deserialize, Serialize};
use std::{
    fmt::Debug,
    io::{Read, Write},
};
use extsort::Sortable;


#[pyclass]
#[derive(Clone, PartialEq, Serialize, Deserialize, Debug)]
pub struct CorResult {
    #[pyo3(get, set)]
    pub gene: String,
    #[pyo3(get, set)]
    pub gem: String,
    #[pyo3(get, set)]
    pub cpg_site_id: Option<String>,
    #[pyo3(get, set)]
    pub correlation: Option<f64>,
    #[pyo3(get, set)]
    pub p_value: Option<f64>,
    #[pyo3(get, set)]
    pub adjusted_p_value: Option<f64>,
}

impl std::fmt::Display for CorResult {
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
        let r = correlation(x, 1, y, 1, self.n);

        // P-value (two-sided)
        // Based on R's cor.test method (https://github.com/SurajGupta/r-source/blob/a28e609e72ed7c47f6ddfbb86c85279a0750f0b7/src/library/stats/R/cor.test.R#L21)
        let statistic = self.degrees_of_freedom.sqrt() * r / (1.0 - r.powi(2)).sqrt();
        let p_value = 2.0 * tdist_P(statistic, self.degrees_of_freedom).min(tdist_Q(statistic, self.degrees_of_freedom));

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
        let mut vec = Vec::with_capacity(2 * self.n);
        let workspace: &mut [f64] = vec.as_mut_slice();
        let rs = spearman(x, 1, y, 1, self.n, workspace);

        // P-value (two-sided)
        // Same behavior as Python Scipy's spearmanr method
        // let t = r * np.sqrt((dof/((rs+1.0)*(1.0-rs))).clip(0))
        let t = rs * (self.degrees_of_freedom / ((rs + 1.0) * (1.0 - rs))).sqrt();
        let ccdf = tdist_Q(t.abs(), self.degrees_of_freedom);
        let p_value = 2.0 * ccdf;

        (rs, p_value)
    }
}

struct Kendall {
    // n: usize
}

impl Kendall {
    fn new(_n: usize) -> Self {
        // Kendall { n }
        Kendall {}
    }
}

impl Correlation for Kendall {
    fn correlate(&self, x: &[f64], y: &[f64]) -> (f64, f64) {
        // Hotfix until https://github.com/zolkko/kendalls/issues/2 is solved
        let gil = Python::acquire_gil();
        let py = gil.python();
        let scipy = PyModule::import(py, "scipy.stats").expect("Scipy was not found");
        scipy.call1("kendalltau", (x.to_vec(), y.to_vec(),)).unwrap().extract().unwrap()

        // let tau = kendalls::tau_b_with_comparator(x, y, |a: &f64, b: &f64| {
        //     a.partial_cmp(&b).unwrap_or(Ordering::Greater)
        // }).unwrap();

        // P-value (two-sided)
        
        // FIXME: significance value seems to be wrong. Issue: https://github.com/zolkko/kendalls/issues/2
        // let significance = kendalls::significance(tau, x.len()); // If this line is correct, use self.n instead of x.len()
        // let significance: f64 = -2.09764910069;
        // let cdf = gaussian_P(-significance.abs(), 1.0);
        // let p_value = 2.0 * cdf;

        // println!("Tau -> {} | Z -> {} P-value -> {}", tau, significance, p_value);

        // (tau, p_value)
    }
}

#[derive(Clone)]
pub enum CorrelationMethod {
    Spearman = 1,
    Kendall = 2,
    Pearson = 3,
}

pub fn get_correlation_method(
    correlation_method: CorrelationMethod,
    number_of_samples: usize,
) -> Box<dyn Correlation> {
    match correlation_method {
        CorrelationMethod::Pearson => Box::new(Pearson::new(number_of_samples)),
        CorrelationMethod::Spearman => Box::new(Spearman::new(number_of_samples)),
        CorrelationMethod::Kendall => Box::new(Kendall::new(number_of_samples)),
    }
}
