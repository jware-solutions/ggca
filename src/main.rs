use std::{cmp::Ordering, process::exit};
use csv::ReaderBuilder;
use rand::prelude::*;
use itertools::iproduct;
use rgsl::{statistics::correlation, randist::t_distribution::tdist_P};
use serde::{Serialize, Deserialize};
extern crate external_sort;
use external_sort::{ExternalSorter, ExternallySortable};

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
        other.p_value.partial_cmp(&self.p_value).unwrap()
    }
}

impl PartialOrd for CorResult {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        other.p_value.partial_cmp(&self.p_value)
    }
}

impl ExternallySortable for CorResult {
    fn get_size(&self) -> u64 { 1 }
}

fn all_vs_all(m1: Matrix, m3: Matrix) {
    let n_float = m1[0].len() as f64;
    let stride = 1;
    let sqrt = (n_float - 2.0).sqrt();
    
    let total_number_of_elements = (m1.len() * m3.len()) as u64;

    let correlations_and_p_values = iproduct!(m1, m3).map(|(tuple1, tuple3)| {        
        // Correlation
        let r = correlation(tuple1.as_slice(), stride, tuple3.as_slice(), stride, tuple1.len());
        
        // P-value
        let t = (r * sqrt) / (1.0 - r.powi(2)).sqrt();
        let ccdf = 1.0 - tdist_P(t, n_float - 2.0);
        let p_value = 2.0 * ccdf;

        CorResult{r, p_value, p_value_adjusted: None}
    });

    // for elem in correlations_and_p_values.clone() {
    //     println!("{:?}", elem);
    // }
    // println!("{}", correlations_and_p_values.clone().count());
    // let first = correlations_and_p_values.next().unwrap();
    // println!("{:#?}", first);

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
    let filtered = ranked.filter(|(_, cor_and_p_value)| cor_and_p_value.as_ref().unwrap().r >= correlation_threhold);
    
    // println!("Cantidad de tuplas a ajustar -> {}", filtered.count());
    
    // Adjustment
    let mut previous_p_value = 999999.0;
    let adjusted = filtered.map(|(rank, mut cor_and_p_value)| {
        let p_value = cor_and_p_value.as_ref().unwrap().p_value;
        let q_value = p_value * (total_number_of_elements / rank as u64) as f64;
        // let q_value = if q_value < 1.0 { q_value } else { 1.0 };
        let q_value = q_value.min(1.0).min(previous_p_value);
		// let q_value = if q_value < previous_max_p_value { q_value } else { previous_max_p_value };
		// let q_value = min(q_value, previous_max_p_value);
        previous_p_value = q_value;
        
        cor_and_p_value.as_mut().unwrap().p_value_adjusted = Some(q_value);
        cor_and_p_value
    });

    // for elem in adjusted {
        // println!("{:#?}", elem);
    // }

    println!("Cantidad final de datos -> {}", adjusted.count());
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
    // for result in rdr.records() {
    //     // The iterator yields Result<StringRecord, Error>, so we check the
    //     // error here.
    //     let record = result.unwrap();
    //     let mut iter = record.iter();
    //     iter.next();
    //     let casted: Vec<f64> = iter
    //     // for elem in casted {
    //     //     println!("{:?}", elem);
    //     // }
    //     println!("{:?}", casted);
    //     break;
    // }
}

fn main() {
    let m1 = get_df("/home/genaro/Descargas/ParaRust/mrna_rust.csv");
    // println!("Dimensiones -> {} x {}", m1.len(), m1[0].len());
    let m3 = get_df("/home/genaro/Descargas/ParaRust/mirna_rust.csv");
    // println!("Dimensiones -> {} x {}", m3.len(), m3[0].len());
    // exit(0);
    // let m1 = (1..50)
    // .map(|_| {
    //     vec![70.0, 65.0, 50.0, 25.0, 30.0, 62.0, 49.0, 30.0, 42.0, 19.0, 22.0]
    // })
    // .collect::<Matrix>();

    // let m3 = (1..50)
    // .map(|j| {
    //     let mut rng = rand::thread_rng();
    //     let i = j as f64;
    //     let mut ar = vec![28.0 * i, 26.0 * i, 27.0 * i, 17.0 * i, 24.0 * i, 27.0 * i, 22.0 * i, 23.0 * i, 18.0 * i, 17.0 * i, 21.0 * i];
    //     ar.shuffle(&mut rng);
    //     ar
    // })
    // .collect::<Matrix>();

    all_vs_all(m1, m3);
}
