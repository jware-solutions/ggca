pub trait Adjustment {
	// fn new(total_number_of_elements: u64) -> Self;
	fn adjust(&mut self, p_value: f64, rank: usize) -> f64;
}

pub struct BenjaminiHochberg {
	total_number_of_elements: f64, // To prevent casting on every adjustment
	previous_p_value: f64,
}

impl BenjaminiHochberg {
	fn new(total_number_of_elements: f64) -> Self {
		BenjaminiHochberg {
			total_number_of_elements,
			previous_p_value: -1.0
		}
	}
}

impl Adjustment for BenjaminiHochberg {
	fn adjust(&mut self, p_value: f64, rank: usize) -> f64 {
		let valid_rank = rank + 1;
		let q_value = p_value * (self.total_number_of_elements / valid_rank as f64);
		let q_value = q_value.min(1.0).max(self.previous_p_value);
		self.previous_p_value = q_value;
		q_value
	}
}

pub struct Bonferroni {
	total_number_of_elements: f64,
}

impl Bonferroni {
	fn new(total_number_of_elements: f64) -> Self {
		Bonferroni { total_number_of_elements }
	}
}

impl Adjustment for Bonferroni {
	fn adjust(&mut self, p_value: f64, _: usize) -> f64 {
		p_value * self.total_number_of_elements
	}
}

pub enum AdjustmentMethod {
	Bonferroni,
	BenjaminiHochberg,
	BenjaminiYekutieli
}

pub fn get_adjustment_method(adjustment_method: AdjustmentMethod, total_number_of_elements: f64) -> Box<dyn Adjustment> {
	match adjustment_method {
		AdjustmentMethod::Bonferroni => Box::new(Bonferroni::new(total_number_of_elements)),
		_ => Box::new(BenjaminiHochberg::new(total_number_of_elements))
	}
}