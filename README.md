# Gene GEM Correlation Analysis (GGCA)

[![CI](https://github.com/jware-solutions/ggca/actions/workflows/ci.yml/badge.svg)](https://github.com/jware-solutions/ggca/actions/workflows/ci.yml)

Computes efficiently the correlation (Pearson, Spearman or Kendall) and the p-value (two-sided) between all the pairs from two datasets. It also supports [CpG Site IDs][cpg-site].

**IMPORTANT**: GGCA is the heart of a platform called Multiomix. On the official website you will be able to use this library in a fast and agile way through a friendly graphical interface (along with many extra features!). Go to https://multiomix.org/ to get started now!

[Python PyPi][pypi-site] | [Rust Crate][crate-site]


## Index

- [Requirements](#requirements)
- [Usage](#usage)
	- [Python](#python)
	- [Rust](#rust)
- [Contributing](#contributing)
- [Considerations](#considerations)


## Requirements

You need to install [GSL][gsl] >= 2.6 in your system to use this library.


## Usage

There are a few examples in `examples` folder for both languages.


### Python

1. Install: `pip install ggca`
1. Configure and call the `correlate` method:

```python
import ggca

# Possible Correlation methods
SPEARMAN = 1
KENDALL = 2
PEARSON = 3

# Possible P-values adjustment methods
BENJAMINI_HOCHBERG = 1
BENJAMINI_YEKUTIELI = 2
BONFERRONI = 3

mrna_file_path = "mrna.csv"
gem_file_path = "mirna.csv"

try:
	(result_combinations, evaluated_combinations) = ggca.correlate(
		mrna_file_path,
		gem_file_path,
		correlation_method=PEARSON,
		correlation_threshold=0.5,
		sort_buf_size=2_000_000,
		adjustment_method=BENJAMINI_HOCHBERG,
		all_vs_all=True,
		gem_contains_cpg=False,
		collect_gem_dataset=None,
		keep_top_n=2  # Keeps only top 2 elements
	)

	print(f'Number of resulting combinations: {len(result_combinations)} of {evaluated_combinations} evaluated combinations')
	for combination in result_combinations:
		print(
			combination.gene,
			combination.gem,
			combination.correlation,
			combination.p_value,
			combination.adjusted_p_value
		)
except ggca.GGCADiffSamplesLength as ex:
	print('Raised GGCADiffSamplesLength:', ex)
except ggca.GGCADiffSamples as ex:
	print('Raised GGCADiffSamples:', ex)
except ggca.InvalidCorrelationMethod as ex:
	print('Raised InvalidCorrelationMethod:', ex)
except ggca.InvalidAdjustmentMethod as ex:
	print('Raised InvalidAdjustmentMethod:', ex)
except ggca.GGCAError as ex:
	print('Raised GGCAError:', ex)
```


### Rust

1. Add crate to `Cargo.toml`: `ggca = { version = "0.4.1", default-features = false  }`
1. Create an analysis and run it:

```rust
use ggca::adjustment::AdjustmentMethod;
use ggca::analysis::Analysis;
use ggca::correlation::CorrelationMethod;

// File's paths
let df1_path = "mrna.csv";
let df2_path = "mirna.csv";

// Some parameters
let gem_contains_cpg = false;
let is_all_vs_all = true;
let keep_top_n = Some(10); // Keeps the top 10 of correlation (sorting by abs values)
let collect_gem_dataset = None; // Better performance. Keep small GEM files in memory

let analysis = Analysis::new_from_files(df1_path.to_string(), df2_path.to_string(), false);
let (result, number_of_elements_evaluated) = analysis.compute(
	CorrelationMethod::Pearson,
	0.7,
	2_000_000,
	AdjustmentMethod::BenjaminiHochberg,
	is_all_vs_all,
	collect_gem_dataset,
	keep_top_n,
)?;

println!("Number of elements -> {} of {} combinations evaluated", result.len(), number_of_elements_evaluated);

for cor_p_value in result.iter() {
	println!("{}", cor_p_value);
}
```

Note that [env_logger][env-logger] crate is used to provide some warning in case some mRNA/GEM combinations produce NaN values (for instance, because the input array has 0 std). In that case, you can add RUST_LOG=warn to your command to produce warnings in the stderr. E.g:

`RUST_LOG=warn cargo test --no-default-features --tests`

or 

`RUST_LOG=warn cargo run --example basic --no-default-features`


## Contributing

All kind of help is welcome! Feel free o submit an issue or a PR. There are some TODOs which are listed below:

- [ ] Parallelize iterators to improve performance
- [ ] Make Rust enums accessible from Python
- [ ] Add support for Windows OS
- [X] Add Rust documentation
- [X] Add tests
- [X] Add MyPy support


### Developing

We provide a Docker image to execute all the commands listed below:

- Build for rust: cargo build [--release]
- Build for Python (uses Maturin):
	1. First of all, remove the `target/wheels` folder: `rm -rf ./target/wheels`
	1. Build wheels: `docker run --rm -v $(pwd):/io jwaresolutions/ggca-build:0.2.3 maturin build --release --skip-auditwheel --manylinux=2014`
	1. Repair wheels to include some missing `.so` files: `docker run --rm -v $(pwd):/io jwaresolutions/ggca-build:0.2.3 ./repair-wheels.sh`
- Only for development:
	- To run the examples you must run with this command due to [an issue][pyo3-issue] with Pyo3: `cargo run [--release] --no-default-features --example <example>`
	- In Python you can test with Maturin:
		- Maturin develop (installs the package in your current virtualenv): `maturin develop`
		- Maturin build (just builds the wheel for your current Python version): `maturin build --manylinux=1-unchecked`
- To publish in PyPi:
	1. Install twine: `pip install twine`
	1. Upload: `twine upload ./target/wheels/wheelhouse/*`


### Tests

All the correlation, p-values and adjusted p-values were taken from [cor.test][r-cor-test] and [p.adjust][r-p-adjust] functions from the R programming language and [statsmodels][statsmodels] package for Python language.

Data in `small_files` folder was retrieved with random sampling from the *Colorectal Adenocarcinoma (TCGA, Nature 2012)* dataset. This dataset can be downloaded from [cBioPortal datasets page][cbioportal-datasets-page] or [this direct link][colorectal-dataset].

All the correlations results were compared directly with R-Multiomics output (old version of [multiomix.org][multiomix] only available for R lang).


### Performance

We use [criterion.rs][criterion] to perform benchmarks. In case you have made a contribution you can check that no regression was added to the project. Just run `cargo bench --no-default-features` before and after your changes to perform a statistical analysis of performance.


## Considerations

If you use any part of our code, or the tool itself is useful for your research, please consider citing:

```
@article{camele2022multiomix,
  title={Multiomix: a cloud-based platform to infer cancer genomic and epigenomic events associated with gene expression modulation},
  author={Camele, Genaro and Menazzi, Sebastian and Chanfreau, Hern{\'a}n and Marraco, Agustin and Hasperu{\'e}, Waldo and Butti, Matias D and Abba, Martin C},
  journal={Bioinformatics},
  volume={38},
  number={3},
  pages={866--868},
  year={2022},
  publisher={Oxford University Press}
}
```


[pypi-site]: https://pypi.org/project/ggca/
[crate-site]: https://crates.io/crates/ggca
[cpg-site]: https://en.wikipedia.org/wiki/CpG_site
[gsl]: https://www.gnu.org/software/gsl/
[pyo3-issue]: https://github.com/PyO3/pyo3/issues/1084
[r-cor-test]: https://www.rdocumentation.org/packages/stats/versions/3.6.2/topics/cor.test
[r-p-adjust]: https://www.rdocumentation.org/packages/stats/versions/3.6.2/topics/p.adjust
[statsmodels]: https://www.statsmodels.org/dev/generated/statsmodels.stats.multitest.multipletests.html
[cbioportal-datasets-page]: https://www.cbioportal.org/datasets
[colorectal-dataset]: https://cbioportal-datahub.s3.amazonaws.com/coadread_tcga_pub.tar.gz
[multiomix]: https://www.multiomix.org
[env-logger]: https://docs.rs/env_logger/latest/env_logger/
[criterion]: https://github.com/bheisler/criterion.rs
