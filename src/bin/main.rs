use std::time::Instant;
use ggca::experiment::{new_from_files, Computation};
use ggca::adjustment::AdjustmentMethod;
use ggca::correlation::CorrelationMethod;


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


    let now = Instant::now();
    
    let experiment = new_from_files(m1_path.to_string(), m3_path.to_string());
	experiment.compute(CorrelationMethod::Pearson, 0.7, 2_000_000, AdjustmentMethod::BenjaminiYekutieli);
	
    println!("Tiempo del experimento -> {} segundos", now.elapsed().as_secs());
}
