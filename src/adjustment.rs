use pyo3::prelude::*;

pub trait Adjustment {
    fn adjust(&mut self, p_value: f64, rank: usize) -> f64;
}

pub struct Bonferroni {
    total_number_of_elements: f64,
}

impl Bonferroni {
    fn new(total_number_of_elements: f64) -> Self {
        Bonferroni {
            total_number_of_elements,
        }
    }
}

impl Adjustment for Bonferroni {
    fn adjust(&mut self, p_value: f64, _: usize) -> f64 {
        (p_value * self.total_number_of_elements).min(1.0)
    }
}

pub struct BenjaminiHochberg {
    total_number_of_elements: f64, // To prevent casting on every adjustment
    previous_max_p_value: f64,
}

impl BenjaminiHochberg {
    fn new(total_number_of_elements: f64) -> Self {
        BenjaminiHochberg {
            total_number_of_elements,
            previous_max_p_value: f64::MAX,
        }
    }
}

impl Adjustment for BenjaminiHochberg {
    fn adjust(&mut self, p_value: f64, rank: usize) -> f64 {
        let valid_rank = self.total_number_of_elements - rank as f64;
        let q_value = p_value * (self.total_number_of_elements / valid_rank);
        let q_value = q_value.min(self.previous_max_p_value).min(1.0);
        self.previous_max_p_value = q_value;
        q_value
    }
}

pub struct BenjaminiYekutieli {
    total_number_of_elements: f64, // To prevent casting on every adjustment
    previous_max_p_value: f64,
    accumulator: f64,
}

impl BenjaminiYekutieli {
    fn new(total_number_of_elements: f64) -> Self {
        let accumulator =
            (1..total_number_of_elements as u64 + 1).fold(0.0, |acc, x| acc + (1.0 / x as f64));

        BenjaminiYekutieli {
            total_number_of_elements,
            previous_max_p_value: f64::MAX,
            accumulator,
        }
    }
}

impl Adjustment for BenjaminiYekutieli {
    fn adjust(&mut self, p_value: f64, rank: usize) -> f64 {
        let valid_rank = self.total_number_of_elements - rank as f64;
        let adjustment_factor = self.accumulator * self.total_number_of_elements / valid_rank;
        let q_value = p_value * adjustment_factor;
        let q_value = q_value.min(self.previous_max_p_value).min(1.0);
        self.previous_max_p_value = q_value;
        q_value
    }
}

#[pyclass(eq, eq_int)]
#[derive(Clone, Debug, PartialEq)]
pub enum AdjustmentMethod {
    BenjaminiHochberg = 1,
    BenjaminiYekutieli = 2,
    Bonferroni = 3,
}

impl std::fmt::Display for AdjustmentMethod {
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        let description = match &self {
            AdjustmentMethod::BenjaminiHochberg => "Benjamini-Hochberg",
            AdjustmentMethod::BenjaminiYekutieli => "Benjamini-Yekutieli",
            AdjustmentMethod::Bonferroni => "Bonferroni",
        };

        write!(f, "{description}")
    }
}

pub fn get_adjustment_method(
    adjustment_method: &AdjustmentMethod,
    total_number_of_elements: f64,
) -> Box<dyn Adjustment> {
    match adjustment_method {
        AdjustmentMethod::Bonferroni => Box::new(Bonferroni::new(total_number_of_elements)),
        AdjustmentMethod::BenjaminiHochberg => {
            Box::new(BenjaminiHochberg::new(total_number_of_elements))
        }
        AdjustmentMethod::BenjaminiYekutieli => {
            Box::new(BenjaminiYekutieli::new(total_number_of_elements))
        }
    }
}
