# Gene GEM Correlation Analysis (GGCA)


## Requirements

1. [GSL][gsl] >= 2.6
1. If Kendall is going to be used it's needed to make some extra stuff until [Kendall crate issue][kendalls-issue] is fixed:
	- Install [Scipy][scipy] in your system. Or...
	- Use a virtual environment:
		1. Create a new Python virtualenv: `python3 -m venv venv`
		1. Activate virtualenv: `source venv/bin/activate`
		1. Install Scipy: `pip install scipy`
		1. Use this library with the virtualenv activated

## Usage

There are a few examples in `examples` folder for both languages.


### Python

1. Install: `pip install ggca`
1. Calls correlate method:

```python
"""Possible Correlation methods"""
SPEARMAN = 1
KENDALL = 2
PEARSON = 3

"""Possible P-values adjustment methods"""
BENJAMINI_HOCHBERG = 1
BENJAMINI_YEKUTIELI = 2
BONFERRONI = 3

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
		combination.cpg_site_id,
		combination.correlation,
		combination.p_value,
		combination.adjusted_p_value
	)
```

### Rust

1. Add crate to `Cargo.toml`: `ggca = "0.1.9"`
1. Create an experiment and run it:

```rust
// File's paths
let df1_path = "mrna.csv";
let df2_path = "mirna.csv";

// Some parameters
let gem_contains_cpg = false;
let is_all_vs_all = true;
let keep_top_n = Some(10); // Keeps the top 10 of correlation (sorting by abs values)
let collect_gem_dataset = None; // Better performance. Keep small GEM files in memory

let experiment = new_from_files(df1_path.to_string(), df2_path.to_string(), false);
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


## Contributing

All kind of help is welcome! Feel free o submit an issue or a PR. There are some TODOs which are listed below:

- [ ] Fix Kendall p-value (waiting [#2][kendalls-issue] to be fixed)
- [ ] Add Rust documentation
- [ ] Make Rust enums accessible from Python
- [ ] Add tests
- [ ] Add MyPy support

### Developing

We provide a Docker image to execute all the commands listed below:

- Build for rust: cargo build [--release]
- Build for Python (uses Maturin):
	1. First of all, remove the `target/wheels` folder: `rm -rf ./target/wheels`
	1. Build wheels: `docker run --rm -v $(pwd):/io jwaresolutions/ggca-build maturin build --release --skip-auditwheel --manylinux=2014`
	1. Repair wheels to include some missing `.so` files: `docker run --rm -v $(pwd):/io jwaresolutions/ggca-build ./repair-wheels.sh`
- Only for development:
	- To run the examples you must run with this command due to [an issue][pyo3-issue] with Pyo3: `cargo run [--release] --no-default-features --example <example>`
	- In Python you can test with Maturin:
		- Maturin develop (installs the package in your current virtualenv): `maturin develop`
		- Maturin build (just builds the wheel for your current Python version): `maturin build --manylinux=1-unchecked`
- To publish in PyPi:
	1. Install twine: `pip install twine`
	1. Upload: `twine upload -u <username> -p <password> ./target/wheels/wheelhouse/*`

[gsl]: https://www.gnu.org/software/gsl/
[pyo3-issue]: https://github.com/PyO3/pyo3/issues/1084
[kendalls-issue]: https://github.com/zolkko/kendalls/issues/2
[scipy]: https://www.scipy.org/