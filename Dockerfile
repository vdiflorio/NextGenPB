FROM rockylinux:9

LABEL maintainer="vincenzo.diflorio@iit.it"

RUN echo "Installing NextGenPB 0.0"

#ENV CFLAGS="-Ofast -mtune=native -march=native"
#ENV CFLAGS="-O3 -mtune=generic -march=generic"
ENV CFLAGS="-O2 -mtune=generic"
ENV CXXFLAGS="$CFLAGS"
ENV FCFLAGS="$CFLAGS"

# Install dependencies and clean up
RUN dnf upgrade -y && \
    dnf --enablerepo crb install -y \
    gcc gcc-c++ gcc-gfortran make python3 texinfo nano\
    gawk procps wget openssh-clients p11-kit diffutils which \
    git rsync zip unzip bzip2 glibc-static patch xz perl-locale \
    perl-Unicode-Normalize infiniband-diags libibverbs libibverbs-utils rdma-core \
    boost boost-devel \
    zlib-devel json-c  jansson jansson-devel \
    gmp gmp-devel mpfr mpfr-devel \
    openssl-devel tar && \
    dnf clean all
RUN dnf install -y epel-release
RUN dnf --enablerepo crb install -y glpk-utils bison \
         openblas-devel autoconf automake libtool scotch scotch-devel \
         redhat-rpm-config libasan bzip2 bzip2-devel MUMPS MUMPS-devel \
         qhull-devel arpack-devel octave octave-devel \
         cmake tbb tbb-devel CGAL-devel

#RUN dnf --enablerepo crb install -y openmpi openmpi-devel mpi
#ENV PATH=/usr/lib64/openmpi/bin:$PATH
#ENV LD_LIBRARY_PATH=/usr/lib64/openmpi/lib:$LD_LIBRARY_PATH

###OPENMPI 
WORKDIR /opt
ENV pkgname=openmpi
ENV pkgver=5.0.0
ENV mver=v5.0
RUN wget https://download.open-mpi.org/release/open-mpi/$mver/$pkgname-$pkgver.tar.bz2 && \
    tar xf $pkgname-$pkgver.tar.bz2 && \
    rm -rf $pkgname-$pkgver.tar.bz2
   
WORKDIR /opt/$pkgname-$pkgver
RUN mkdir build
WORKDIR /opt/$pkgname-$pkgver/build/
RUN ../configure --prefix="/opt/openmpi" \
                 --disable-silent-rules \
                 --enable-mpi-fortran=all \
                 --disable-mpi-cxx \
                 --enable-ipv6 --without-x   \
                 --without-libltdl --disable-dlopen   \
                 --enable-pretty-print-stacktrace  \
                 --enable-openib-rdmacm  \
                 --disable-sphinx
RUN make -j10
RUN make install
WORKDIR /opt
RUN rm -rf /opt/$pkgname-$pkgver 
ENV PATH=/opt/openmpi/bin:$PATH
ENV LD_LIBRARY_PATH=/opt/openmpi/lib:$LD_LIBRARY_PATH

###NanoShaper###
WORKDIR /opt/
RUN git clone https://gitlab.iit.it/SDecherchi/nanoshaper.git
WORKDIR /opt/nanoshaper/
RUN cp  CMakeLists_so.txt CMakeLists.txt 
WORKDIR /opt/nanoshaper/build_lib
RUN cmake .. -DCMAKE_BUILD_TYPE="Release"
RUN make -j10


##LIS###
WORKDIR /opt/
ENV pkgname=lis
ENV pkgver=2.1.6
ENV archive=$pkgname-$pkgver.zip
RUN wget https://www.ssisc.org/lis/dl/$archive
RUN unzip $archive
WORKDIR /opt/$pkgname-$pkgver
RUN chmod 755 configure config/install-sh
RUN ./configure --prefix="/opt/$pkgname" \
                --enable-mpi \
                --enable-shared \
                --disable-static
RUN make -j10
RUN make -j4 install
WORKDIR /opt/
RUN rm -rf $pkgname-$pkgver $archive


###P4EST###
ENV pkgname=p4est
ENV pkgver=2.8.6
ENV archive=$pkgname-$pkgver.tar.gz
RUN wget https://p4est.github.io/release/$archive && \
    tar -xzf $archive
RUN    rm -rf $archive
WORKDIR /opt/$pkgname-$pkgver
ENV CC=mpicc
ENV CXX=mpic++
RUN mkdir /opt/$pkgname-$pkgver/build
WORKDIR /opt/$pkgname-$pkgver/build
RUN ../configure --prefix="/opt/$pkgname" \
                --enable-mpi \
                --enable-shared \
                --disable-static \
                --disable-vtk-binary \
                --without-blas \
                CPPFLAGS="-I/opt/openmpi/include"
RUN make -j10
RUN make install
WORKDIR /opt/
RUN rm -rf $pkgname-$pkgver

###OCTAVE-FILE-IO###
ENV pkgname=octave_file_io
ENV pkgver=1.0.91
ENV archive=v$pkgver.tar.gz
RUN wget https://github.com/carlodefalco/octave_file_io/archive/refs/tags/$archive && \
    tar -xf $archive && \
    rm -rf $archive
WORKDIR /opt/$pkgname-$pkgver
RUN ./autogen.sh
RUN mkdir build
WORKDIR /opt/$pkgname-$pkgver/build
RUN ../configure --prefix="/opt/$pkgname" \
                --with-octave-home=/usr/bin \
                CC=mpicc CXX=mpicxx
RUN make -j10
RUN make install
WORKDIR /opt
RUN rm -rf /opt/$pkgname-$pkgver

RUN mkdir /usr/local/nextgenPB
RUN git config --global \
  url."https://vdiflorio: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@github.com/".insteadOf \
  "https://github.com/"
RUN git clone --branch main --single-branch https://github.com/vdiflorio/nextgenPB.git /usr/local/nextgenPB/

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
RUN make -j10
RUN make install
WORKDIR /opt

RUN rm -rf /opt/$pkgname-$pkgver
RUN rm -rf /opt/nanoshaper/example && \
    rm -rf /opt/nanoshaper/src_client && \
    rm -rf /opt/nanoshaper/test && \
    rm -rf /usr/local/nextgenPB/d* && \
    rm -rf /usr/local/nextgenPB/*.bz2

###NGPB###
WORKDIR  /usr/local/nextgenPB
RUN cp local_settings_docker_generic.mk /usr/local/nextgenPB/src/local_settings.mk
WORKDIR  /usr/local/nextgenPB/src
RUN make clean all

WORKDIR /opt

# Create a volume for passing input data (potfile, pqrfile, etc.)
ENV PATH=/usr/local/nextgenPB/src:$PATH

VOLUME ["/App"]
WORKDIR /App
