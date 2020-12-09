## With only Lib

1. Build: cargo build [--release]
1. Maturin develop: maturin develop
1. Maturin build: maturin build --manylinux=1-unchecked
1. Maturin publish: maturin publish --manylinux=1-unchecked --username [username] --password [password]
1. Maturin build manylinux: docker run --rm -v $(pwd):/io ggca-build maturin build --skip-auditwheel
1. Maturin publish manylinux: docker run --rm -v $(pwd):/io ggca-build maturin publish --skip-auditwheel --username [username] --password [password]

## With Bin/Lib

See issue https://github.com/PyO3/pyo3/issues/1084

1. Build: cargo build [--release] --no-default-features
1. Maturin develop: maturin develop --cargo-extra-args="--no-default-features"
1. Maturin build/publish: maturin build --cargo-extra-args="--no-default-features" --manylinux 1-unchecked

Don't check compliant libraries as GSL is not in the ManyLinux list
1. Maturin build manylinux: docker run --rm -v $(pwd):/io ggca-build:5 build --skip-auditwheel
1. Maturin publish manylinux: docker run --rm -v $(pwd):/io ggca-build:5 publish --skip-auditwheel --username [username] --password [password]


## Pasos para solucionar el problema

1. Primero hice un chequeo de funcionalidad ejecutando python3 setup.py bdist_wheel
1. Importe en un proyecto aparte el wheel generado en el punto anterior: pip install /path/to/whl
1. Hacerlo con Setuptools_rust siguiendo los pasos de: https://github.com/PyO3/setuptools-rust
1. Hacer imagen docker con las dependencias necesarias
1. Hacer un build: docker run --rm -v $(pwd):/io ggca-build:3 /io/build-wheels.sh
1. Instalar twine
1. Hacer upload de los wheels con twine: twine upload dist/*manylinux2014_x86_64.whl