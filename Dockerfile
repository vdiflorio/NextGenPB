FROM rockylinux:9

LABEL maintainer="vincenzo.diflorio@iit.it"

RUN echo "Installing NextGenPB 0.0"

#ENV CFLAGS="-Ofast -mtune=native -march=native"
ENV CFLAGS="-O2 -mtune=generic"
ENV CXXFLAGS="$CFLAGS"
ENV FCFLAGS="$CFLAGS"

# Install dependencies and clean up
RUN dnf upgrade -y && \
    dnf --enablerepo=crb install -y epel-release && \
    dnf --enablerepo=crb install -y \
        gcc gcc-c++ gcc-gfortran make python3 texinfo nano gawk procps wget \
        openssh-clients p11-kit diffutils which git rsync zip unzip bzip2 \
        glibc-static patch xz perl-locale perl-Unicode-Normalize \
        infiniband-diags libibverbs libibverbs-utils rdma-core boost boost-devel \
        zlib-devel json-c jansson jansson-devel gmp gmp-devel mpfr mpfr-devel \
        openssl-devel tar glpk-utils bison openblas-devel autoconf automake \
        libtool scotch scotch-devel redhat-rpm-config libasan bzip2 bzip2-devel \
        MUMPS MUMPS-devel qhull-devel arpack-devel octave octave-devel cmake tbb \
        tbb-devel CGAL-devel && \
    dnf clean all


###OPENMPI (if you want a specific version of openmpi to match local verion of openmpi)
RUN wget https://download.open-mpi.org/release/open-mpi/v4.0/openmpi-4.0.0.tar.bz2 && \
    tar xf openmpi-4.0.0.tar.bz2 && \
    rm -rf openmpi-4.0.0.tar.bz2 && \
    cd openmpi-4.0.0 && \
    mkdir build && cd build && \
    ../configure --prefix="/opt/openmpi" \
                 --disable-silent-rules \
                 --enable-mpi-fortran=all \
                 --disable-mpi-cxx \
                 --enable-ipv6 --without-x \
                 --without-libltdl --disable-dlopen \
                 --enable-pretty-print-stacktrace \
                 --enable-openib-rdmacm --disable-sphinx && \
    make -j$(nproc) && \
    make install && \
    cd /opt && rm -rf openmpi-4.0.0

ENV PATH=/opt/openmpi/bin:$PATH
ENV LD_LIBRARY_PATH=/opt/openmpi/lib:$LD_LIBRARY_PATH

### If you don't need to match the local verion of openmpi

#RUN dnf --enablerepo crb install -y openmpi openmpi-devel mpi
#ENV PATH=/usr/lib64/openmpi/bin:$PATH
#ENV LD_LIBRARY_PATH=/usr/lib64/openmpi/lib:$LD_LIBRARY_PATH

### NanoShaper ###
WORKDIR /opt
RUN git clone https://gitlab.iit.it/SDecherchi/nanoshaper.git && \
    cd nanoshaper && \
    cp CMakeLists_so.txt CMakeLists.txt && \
    cd build_lib && \
    cmake .. -DCMAKE_BUILD_TYPE="Release" && \
    make -j$(nproc)

### LIS ###
WORKDIR /opt
ENV LIS_VERSION="2.1.6"
RUN wget https://www.ssisc.org/lis/dl/lis-${LIS_VERSION}.zip && \
    unzip lis-${LIS_VERSION}.zip && cd lis-${LIS_VERSION} && \
    chmod 755 configure config/install-sh && \
    ./configure --prefix="/opt/lis" --enable-mpi --enable-shared --disable-static && \
    make -j$(nproc) && make install && \
    cd /opt && rm -rf lis-${LIS_VERSION} lis-${LIS_VERSION}.zip


### P4EST ###
ENV P4EST_VERSION="2.8.6"
RUN wget https://p4est.github.io/release/p4est-${P4EST_VERSION}.tar.gz && \
    tar -xzf p4est-${P4EST_VERSION}.tar.gz && \
    rm -rf p4est-${P4EST_VERSION}.tar.gz && \
    cd /opt/p4est-${P4EST_VERSION} && \
    mkdir build && \
    cd build && \
    CC=mpicc CXX=mpic++ ../configure --prefix="/opt/p4est" \
                --enable-mpi \
                --enable-shared \
                --disable-static \
                --disable-vtk-binary \
                --without-blas \
                CPPFLAGS="-I/opt/openmpi/include" && \
    make -j$(nproc) && \
    make install && \
    cd /opt && rm -rf p4est-${P4EST_VERSION}



### Octave-File-IO ###
WORKDIR /opt
ENV OCTAVE_FILE_IO_VERSION="1.0.91"
RUN wget https://github.com/carlodefalco/octave_file_io/archive/refs/tags/v${OCTAVE_FILE_IO_VERSION}.tar.gz && \
    tar -xf v${OCTAVE_FILE_IO_VERSION}.tar.gz && \
    rm -rf v${OCTAVE_FILE_IO_VERSION}.tar.gz && \
    cd octave_file_io-${OCTAVE_FILE_IO_VERSION} && \
    ./autogen.sh && mkdir build && cd build && \
    ../configure --prefix="/opt/octave_file_io" \
                 --with-octave-home=/usr/bin \
                 CC=mpicc CXX=mpicxx && \
    make -j$(nproc) && make install && \
    cd /opt && rm -rf octave_file_io-${OCTAVE_FILE_IO_VERSION}


###BIMPP###
WORKDIR /opt
ENV pkgname=bimpp
ENV pkgver=patch
RUN git clone --branch nextgenPB https://github.com/vdiflorio/bimpp.git $pkgname-$pkgver
WORKDIR /opt/$pkgname-$pkgver
RUN mkdir build
RUN ./autogen.sh
WORKDIR /opt/$pkgname-$pkgver/build
RUN ../configure --prefix=/opt/bimpp/ \
                CPPFLAGS="-I/opt/openmpi/include -I/usr/include/ -I/usr/include/MUMPS/ -I/opt/p4est/include -DOMPI_SKIP_MPICXX -DHAVE_OCTAVE_44 -DBIM_TIMING" \
                LDFLAGS="-L/usr/lib64 -L/opt/openmpi/lib -L/lib64" \
                 --with-blas-lapack="-lopenblas"  \
                 --with-octave_file_io-home=/opt/octave_file_io \
                 --with-octave-home=/usr/bin \
                 --with-p4est-home=/opt/p4est \
                 --with-lis-home=/opt/lis \
                 --with-mumps-home=/usr/lib64 \
F77=mpif90 CXX=mpicxx MPICC=mpicc CC=mpicc \
--with-mumps-extra-libs="-L/usr/lib64 -L/opt/openmpi/lib -L/usr/lib64/lib -lscotch -lmpi \
 -lmpi_usempif08 -lmpi_usempi_ignore_tkr -lmpi_mpifh -lmpi -lopenblas  -lgfortran"
RUN make -j$(nproc)
RUN make install
WORKDIR /opt

### NextGenPB ###
WORKDIR /usr/local/nextgenPB
RUN git clone --branch main https://github.com/vdiflorio/nextgenPB.git . && \
    cp /usr/local/nextgenPB/local_setting/local_settings_rocky.mk /usr/local/nextgenPB/src/local_settings.mk && \
    cd src && make clean all

# Clean up unnecessary files
RUN rm -rf /opt/nanoshaper/{example,src_client,test}

WORKDIR /opt

# Create a volume for passing input data (potfile, pqrfile, etc.)
ENV PATH=/usr/local/nextgenPB/src:$PATH

VOLUME ["/App"]
WORKDIR /App
