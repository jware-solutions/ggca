FROM konstin2/maturin

# Solves CentOS 6 EOL
RUN echo "https://vault.centos.org/6.10/os/x86_64/" > /var/cache/yum/x86_64/6/base/mirrorlist.txt
RUN echo "http://vault.centos.org/6.10/extras/x86_64/" > /var/cache/yum/x86_64/6/extras/mirrorlist.txt
RUN echo "http://vault.centos.org/6.10/updates/x86_64/" > /var/cache/yum/x86_64/6/updates/mirrorlist.txt
RUN echo "http://vault.centos.org/6.10/sclo/x86_64/rh" > /var/cache/yum/x86_64/6/centos-sclo-rh/mirrorlist.txt
RUN echo "http://vault.centos.org/6.10/sclo/x86_64/sclo" > /var/cache/yum/x86_64/6/centos-sclo-sclo/mirrorlist.txt
RUN yum install -y gsl-devel

RUN echo $LD_LIBRARY_PATH

ENTRYPOINT [""]
