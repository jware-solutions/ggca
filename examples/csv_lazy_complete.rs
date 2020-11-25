// ESTE ARCHIVO CONTIENE UN INTENTO DE MANEJAR UN BOX DE BOXES AL LEVANTAR EL ARCHIVO 
// ES SOLO EXPERIMENTAL PARA VER SI APORTABA ALGUNA MEJORA AL CONSUMO DE MEMORIA, PERO EN VISTA DE QUE FUNCIONA
// BASTANTE BIEN LA VERSION CON ITERADORES DE VECTORES NO VALE MUCHO LA PENA SEGUIR INVIRTIENDO TIEMPO EN ESTE ARCHIVO

use std::cmp::Ordering;
use csv::{ReaderBuilder, StringRecord, Error};
use itertools::iproduct;
use rgsl::{randist::beta::beta_P, statistics::correlation};
use serde::{Serialize, Deserialize};
extern crate external_sort;
use external_sort::{ExternalSorter, ExternallySortable};
use std::time::Instant;

type Matrix = Box<dyn Iterator<Item = Box<dyn Iterator<Item = f64>>>>;

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

fn all_vs_all(m1: Matrix, len_m1: u64, m3: Matrix, len_m3: u64, number_of_columns: usize) {
    let stride = 1;
    let ab = (number_of_columns / 2 - 1) as f64;

	let total_number_of_elements: u64 = len_m1 * len_m3;
	println!("Numero total de combinaciones a evaluar -> {}", total_number_of_elements);

	// We need a collected object for right-side of the iproduct macro
	let m3_collected = m3.collect::<Vec<Vec<f64>>>();

    let correlations_and_p_values = iproduct!(m1, m3_collected).map(|(tuple1, tuple3)| {        
        // Correlation
        let r = correlation(tuple1.as_slice(), stride, tuple3.as_slice(), stride, number_of_columns);

        // P-value
        // Same behaviour as Python Scipy's pearsonr method
        let x = 0.5 * (1.0 - r.abs());
		let p_value = 2.0 * beta_P(x, ab, ab);

        CorResult{r, p_value, p_value_adjusted: None}
    });

    // println!("correlations_and_p_values.count -> {}", correlations_and_p_values.count());

    // for elem in correlations_and_p_values {
    //     println!("{:?}", elem);
    // }
    
    // Sorting
    let chunk_size = 1_000_000;
    // let external_sorter: ExternalSorter<CorResult> = ExternalSorter::new(total_number_of_elements.clone(), None);
    let external_sorter: ExternalSorter<CorResult> = ExternalSorter::new(chunk_size, None);
    let sorted = external_sorter.sort(correlations_and_p_values).unwrap();
    
    // println!("sorted.count -> {}", sorted.count());
    
    // // for elem in sorted {
    // //     println!("{:?}", elem);
    // // }

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
	
	let mut number_of_result_elements = 0;
    for elem in adjusted {
        let valid_elem = elem.unwrap();
		println!("Cor: {} | p-value: {:+e} | adjusted p-value {:+e}", valid_elem.r, valid_elem.p_value, valid_elem.p_value_adjusted.unwrap());
		number_of_result_elements += 1;
    }

    println!("Cantidad final de datos -> {}", number_of_result_elements);
}

fn get_df(path: &str) -> Matrix {
    // Build the CSV reader and iterate over each record.
    let rdr = ReaderBuilder::new()
        .delimiter(b'\t')
		.from_path(path).unwrap();
	
	let parse_record = |record_result: Result<StringRecord, Error>| {
        Box::new(record_result.unwrap().into_iter().map(|cell| {
            cell.parse::<f64>().unwrap()
        })) as Box<dyn Iterator<Item = f64>>
    };

	// Box<Iterator<Item = Box<Iterator<Item = f64>>>>
    Box::new(rdr.into_records().map(parse_record))
}

fn main() {
    // Chicos
    // let m1_path = "/home/genaro/Descargas/ParaRust/mrna_rust.csv";
    // let m3_path = "/home/genaro/Descargas/ParaRust/mirna_rust.csv";

    // Medianos
    // let m1_path = "/home/genaro/Descargas/ParaRust/mrna_rust_mediano.csv";
    // let m3_path = "/home/genaro/Descargas/ParaRust/mirna_rust_mediano.csv";

    // Grandes
    // let m1_path = "/home/genaro/Descargas/ParaRust/mrna_rust_gigante.csv";
    // let m3_path = "/home/genaro/Descargas/ParaRust/mirna_rust_gigante.csv";
    
    // Masivos
    let m1_path = "/home/genaro/Descargas/ParaRust/mrna_rust_gigante.csv";
    let m3_path = "/home/genaro/Descargas/ParaRust/cna_rust_gigante.csv";


    let m1 = get_df(m1_path);
    let mut m1_aux = get_df(m1_path);
    let number_of_columns = m1_aux.next().unwrap().count();
    let len_m1 = m1_aux.count() + 1; // Plus discarded element by next()

	let m3 = get_df(m3_path);
    let m3_aux = get_df(m3_path);
    let len_m3 = m3_aux.count();
    
    println!("Dimensiones de m1 -> {} x {}", len_m1, number_of_columns);
    println!("Dimensiones de m3 -> {} x {}", len_m3, number_of_columns);
    
    let now = Instant::now();
    

	// Grandes
	all_vs_all(m1, len_m1 as u64, m3, len_m3 as u64, number_of_columns);
	
    println!("Tiempo del experimento -> {} segundos", now.elapsed().as_secs());
}
