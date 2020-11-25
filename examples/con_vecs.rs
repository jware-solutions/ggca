use std::cmp::Ordering;
use csv::ReaderBuilder;
use itertools::iproduct;
use rgsl::{randist::beta::beta_P, statistics::correlation};
use serde::{Serialize, Deserialize};
extern crate external_sort;
use external_sort::{ExternalSorter, ExternallySortable};
use std::time::Instant;

type Matrix = Vec<Vec<f64>>;

#[derive(Clone, PartialEq, Serialize, Deserialize, Debug)]
struct CorResult {
    r: f64,
    p_value: f64,
    p_value_adjusted: Option<f64>
}

impl Eq for CorResult { }

impl Ord for CorResult {
    // Sorts in descending order
    fn cmp(&self, other: &Self) -> Ordering {
        self.partial_cmp(&other).unwrap()
    }
}

impl PartialOrd for CorResult {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        self.p_value.partial_cmp(&other.p_value)
        // other.p_value.partial_cmp(&self.p_value)
    }
}

impl ExternallySortable for CorResult {
    fn get_size(&self) -> u64 { 1 }
}

fn all_vs_all(m1: Matrix, m3: Matrix) {
    let n = m1[0].len();
    let stride = 1;
    let ab = (n / 2 - 1) as f64;
    
    let total_number_of_elements = (m1.len() * m3.len()) as u64;

    let correlations_and_p_values = iproduct!(m1, m3).map(|(tuple1, tuple3)| {        
        // Correlation
        let r = correlation(tuple1.as_slice(), stride, tuple3.as_slice(), stride, n);

        // P-value
        // Same behaviour as Python Scipy's pearsonr method
        let x = 0.5 * (1.0 - r.abs());
        let p_value = 2.0 * beta_P(x, ab, ab);
        CorResult{r, p_value, p_value_adjusted: None}
    });

    // for elem in correlations_and_p_values.clone() {
    //     println!("{:?}", elem);
    // }
    // println!("{}", correlations_and_p_values.clone().count());

    // Sorting
    let external_sorter: ExternalSorter<CorResult> = ExternalSorter::new(total_number_of_elements, None);
    let sorted = external_sorter.sort(correlations_and_p_values).unwrap();

    // for elem in sorted {
    //     println!("{:?}", elem);
    // }

    // Ranking 
    let ranked = sorted.enumerate();

    // Filtering
    let correlation_threhold = 0.7;
    let filtered = ranked.filter(|(_, cor_and_p_value)| cor_and_p_value.as_ref().unwrap().r.abs() >= correlation_threhold);
    
    // println!("Cantidad de tuplas a ajustar -> {}", filtered.count());
    
    // Adjustment
    let mut previous_p_value = 999999.0;
    let adjusted = filtered.map(|(rank, mut cor_and_p_value)| {
        let valid_rank = rank + 1;
        let p_value = cor_and_p_value.as_ref().unwrap().p_value;
        // println!("p_value -> {}", p_value);
        // println!("valid_rank -> {}", valid_rank);
        let q_value = p_value * (total_number_of_elements / valid_rank as u64) as f64;
        // println!("q_value -> {}", q_value);
        let q_value = q_value.min(1.0).min(previous_p_value);
        // println!("q_value despues del min -> {}", q_value);
        previous_p_value = q_value;
        
        cor_and_p_value.as_mut().unwrap().p_value_adjusted = Some(q_value);
        cor_and_p_value
    });

    for elem in adjusted {
        let valid_elem = elem.unwrap();
        println!("Cor: {} | p-value: {:+e} | adjusted p-value {:+e}", valid_elem.r, valid_elem.p_value, valid_elem.p_value_adjusted.unwrap());
    }

    // println!("Cantidad final de datos -> {}", adjusted.count());
}

fn get_df(path: &str) -> Matrix {
    // Build the CSV reader and iterate over each record.
    let mut rdr =  ReaderBuilder::new()
        .delimiter(b'\t')
        .from_path(path).unwrap();
    
    rdr.records().map(|result| {
        let record = result.unwrap();
        record.iter().map(|x| x.parse::<f64>().expect(x)).collect::<Vec<f64>>()
    }).collect()
}

fn main() {
    let m1 = get_df("/home/genaro/Descargas/ParaRust/mrna_rust.csv");
    println!("Dimensiones -> {} x {}", m1.len(), m1[0].len());
    let m3 = get_df("/home/genaro/Descargas/ParaRust/mirna_rust.csv");
    println!("Dimensiones -> {} x {}", m3.len(), m3[0].len());
    
    let now = Instant::now();
    all_vs_all(m1, m3);
    // all_vs_all_with_iterators(m1, m3);
    println!("Tiempo del experimento -> {} segundos", now.elapsed().as_secs());
}
