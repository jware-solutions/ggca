use std::time::Instant;
use pruebas_correlation::{new_from_files, Computation};
use pruebas_correlation::adjustment::AdjustmentMethod;
use pruebas_correlation::correlation::CorrelationMethod;


fn main() {
    // Chicos (reducido)
    let m1_path = "/home/genaro/Descargas/ParaRust/mrna_rust.csv";
    let m3_path = "/home/genaro/Descargas/ParaRust/mirna_rust.csv";

    // Medianos
    // let m1_path = "/home/genaro/Descargas/ParaRust/mrna_rust_mediano.csv";
    // let m3_path = "/home/genaro/Descargas/ParaRust/mirna_rust_mediano.csv";

    // Grandes
    // let m1_path = "/home/genaro/Descargas/ParaRust/mrna_rust_gigante.csv";
    // let m3_path = "/home/genaro/Descargas/ParaRust/mirna_rust_gigante.csv";
    
    // Masivos
    // let m1_path = "/home/genaro/Descargas/ParaRust/mrna_rust_gigante.csv";
    // let m3_path = "/home/genaro/Descargas/ParaRust/cna_rust_gigante.csv";


    let experiment = new_from_files(m1_path.to_string(), m3_path.to_string());
    let now = Instant::now();
    
	// all_vs_all(m1, len_m1 as u64, m3, len_m3 as u64, number_of_columns, CorrelationMethod::Kendall, 0.7, 2_000_000, AdjustmentMethod::BenjaminiHochberg);
	experiment.compute(CorrelationMethod::Pearson, 0.7, 2_000_000, AdjustmentMethod::BenjaminiYekutieli);
	
    println!("Tiempo del experimento -> {} segundos", now.elapsed().as_secs());
}
