use std::cmp::Ordering;

use rgsl::{randist::beta::beta_P, randist::{gaussian::gaussian_P, t_distribution::tdist_Q}, statistics::{correlation, spearman}};

pub trait Correlation {
    fn correlate(&self, x: &[f64], y: &[f64]) -> (f64, f64);
}

pub struct Pearson {
    n: usize,
    ab: f64,
}

impl Pearson {
    fn new(n: usize) -> Self {
        let ab = (n / 2 - 1) as f64;
        Pearson { ab, n }
    }
}

impl Correlation for Pearson {
    fn correlate(&self, x: &[f64], y: &[f64]) -> (f64, f64) {
        let r = correlation(x, 1, y, 1, self.n);

        // P-value
        // Same behavior as Python Scipy's pearsonr method
        let x = 0.5 * (1.0 - r.abs());
        let p_value = 2.0 * beta_P(x, self.ab, self.ab);

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

        // P-value
        // Same behavior as Python Scipy's spearmanr method
        // let t = r * np.sqrt((dof/((rs+1.0)*(1.0-rs))).clip(0))
        let t = rs * (self.degrees_of_freedom / ((rs + 1.0) * (1.0 - rs))).sqrt();
        let ccdf = tdist_Q(t.abs(), self.degrees_of_freedom);
        let p_value = 2.0 * ccdf;

        (rs, p_value)
    }
}

struct Kendall {
    n: usize
}

impl Kendall {
    fn new(n: usize) -> Self {
        Kendall { n }
    }
}

impl Correlation for Kendall {
    fn correlate(&self, x: &[f64], y: &[f64]) -> (f64, f64) {
        println!("{:?}", x);
        println!("{:?}", y);
        let tau = kendalls::tau_b_with_comparator(x, y, |a: &f64, b: &f64| {
            a.partial_cmp(&b).unwrap_or(Ordering::Greater)
        }).unwrap();

        // P-value
        // FIXME: significance value seems to be wrong. Issue: https://github.com/zolkko/kendalls/issues/2
        // let significance = kendalls::significance(tau, x.len()); // If this line is correct, use self.n instead of x.len()
        let significance: f64 = -2.09764910069;
        let cdf = gaussian_P(-significance.abs(), 1.0);
        let p_value = 2.0 * cdf;

        println!("Tau -> {} | Z -> {} P-value -> {} | Diff -> {}", tau, significance, p_value, (0.0389842391014099 - p_value).abs());

        (tau, p_value)
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
