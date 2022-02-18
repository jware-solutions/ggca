FROM konstin2/maturin:v0.12.9

# Installs GSL 2.6
WORKDIR /
COPY gsl-2.6.tar.gz /gsl-2.6.tar.gz
RUN tar xvzf gsl-2.6.tar.gz
WORKDIR /gsl-2.6
RUN ./configure && make && make install && make clean

WORKDIR /io

ENTRYPOINT [""]
