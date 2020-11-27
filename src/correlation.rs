use rgsl::{randist::beta::beta_P,  statistics::{correlation, spearman}, randist::t_distribution::tdist_Q};

pub trait Correlation {
	fn correlate(&self, x: &[f64], y: &[f64]) -> (f32, f64);
}

pub struct Pearson {
	n: usize,
	ab: f64
}

impl Pearson {
	fn new(n: usize) -> Self {
		let ab = (n / 2 - 1) as f64;
		Pearson { ab, n }
	}
}

impl Correlation for Pearson {
	fn correlate(&self, x: &[f64], y: &[f64]) -> (f32, f64) {
		let r = correlation(x, 1, y, 1, self.n);

        // P-value
        // Same behaviour as Python Scipy's pearsonr method
        let x = 0.5 * (1.0 - r.abs());
		let p_value = 2.0 * beta_P(x, self.ab, self.ab);

		(r as f32, p_value)
	}
}

struct Spearman {
	n: usize,
	degrees_of_freedom: f64
}

impl Spearman {
	fn new(n: usize) -> Self {
		Spearman { n, degrees_of_freedom: (n - 2) as f64 }
	}
}

impl Correlation for Spearman {
	fn correlate(&self, x: &[f64], y: &[f64]) -> (f32, f64) {
		let mut vec = Vec::with_capacity(2 * self.n);
		let workspace: &mut[f64] = vec.as_mut_slice();
		let rs = spearman(x, 1, y, 1, self.n, workspace);

        // P-value
        // Same behaviour as Python Scipy's spearmanr method
		// let t = r * np.sqrt((dof/((rs+1.0)*(1.0-rs))).clip(0))
		let t = rs * (self.degrees_of_freedom / ((rs + 1.0) * (1.0 - rs))).sqrt();
		let ccdf = tdist_Q(t.abs(), self.degrees_of_freedom);
		let p_value = 2.0 * ccdf;

		if p_value == 0.0 {
			println!("r -> {} | t -> {} | degrees_of_freedom -> {} | ccdf -> {} | p_value -> {}", rs, t, self.degrees_of_freedom, ccdf, p_value);
			std::process::exit(0)
		}

		(rs as f32, p_value)
	}
}


pub enum CorrelationMethod {
	Pearson,
	Spearman,
	Kendall
}

pub fn get_correlation_method(correlation_method: CorrelationMethod, n: usize) -> Box<dyn Correlation> {
	match correlation_method {
		CorrelationMethod::Pearson => Box::new(Pearson::new(n)),
		CorrelationMethod::Spearman => Box::new(Spearman::new(n)),
		_ => Box::new(Pearson::new(n))
	}
}