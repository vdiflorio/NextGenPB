FROM rockylinux:9

LABEL maintainer="vincenzo.diflorio@iit.it"

RUN echo "Installing NextGenPB 0.0"

# Install dependencies and clean up
RUN dnf upgrade -y && \
    dnf --enablerepo crb install -y \
    gcc gcc-c++ gcc-gfortran make python3 texinfo nano\
    gawk procps wget openssh-clients p11-kit diffutils which \
    git rsync zip unzip bzip2 glibc-static patch xz perl-locale \
    perl-Unicode-Normalize \
    boost boost-devel \
    mpi openmpi zlib-devel json-c  jansson jansson-devel\
    gmp gmp-devel mpfr mpfr-devel \
    openssl-devel tar && \
    dnf clean all


###MK###
WORKDIR /opt/
RUN git clone https://github.com/carlodefalco/mk.git
WORKDIR /opt/mk
RUN git checkout nextgenPB

ENV CFLAGS="-O3 -mtune=native"
ENV CXXFLAGS=$CFLAGS
ENV FCFLAGS=$CFLAGS

###CMAKE (at least 3.1 for tbb)###
WORKDIR /opt/
RUN wget https://cmake.org/files/v3.19/cmake-3.19.0-Linux-x86_64.tar.gz && \
    tar -xzf cmake-3.19.0-Linux-x86_64.tar.gz && \
    rm cmake-3.19.0-Linux-x86_64.tar.gz


###TBB###
WORKDIR /opt/
RUN wget https://github.com/oneapi-src/oneTBB/archive/v2021.4.0.tar.gz && \
    tar -xvzf v2021.4.0.tar.gz && \
    mv oneTBB-2021.4.0 tbb
WORKDIR /opt/tbb
RUN /opt/cmake-3.19.0-Linux-x86_64/bin/cmake \
        -DCMAKE_BUILD_TYPE=Release \
        -DCMAKE_INSTALL_PREFIX=/opt/tbb/lib \
        -DTBB_TEST=OFF .
RUN make -j4
RUN make install
WORKDIR /opt/
RUN rm v2021.4.0.tar.gz


###NanoShaper###
WORKDIR /opt/
RUN git clone https://gitlab.iit.it/SDecherchi/nanoshaper.git


###CGAL-5.6.1###
ENV CGALVER="5.6.1"
ENV CGALFILE="CGAL-5.6.1.tar.xz"
ENV CGALDIR="CGAL-5.6.1"
WORKDIR /opt/
RUN wget https://github.com/CGAL/cgal/releases/download/v$CGALVER/$CGALFILE && \
   tar -xvf ${CGALFILE} && \
   rm ${CGALFILE}
###Patching###
RUN cp /opt/mk/pkgs/nanoshaper/wp* /opt/${CGALDIR}
WORKDIR /opt/${CGALDIR}
RUN [ -r wp2.diff ] && patch -p0 -i wp2.diff && \
    [ -r wp3.diff ] && patch -p0 -i wp3.diff


###Install CGAL###
RUN  mkdir build
WORKDIR /opt/${CGALDIR}/build
RUN /opt/cmake-3.19.0-Linux-x86_64/bin/cmake .. \
        -DCMAKE_INSTALL_PREFIX=/opt/${CGALDIR}/ \
        -DWITH_examples=false \
        -DWITH_CGAL_ImageIO=false \
        -DCMAKE_BUILD_TYPE=Release
RUN make -j4 install


###NS Install NS###
WORKDIR /opt/nanoshaper/
RUN git checkout Nanoshaper_devel_TBB
RUN cp  CMakeLists_so.txt CMakeLists.txt 
WORKDIR /opt/nanoshaper/build_lib
RUN /opt/cmake-3.19.0-Linux-x86_64/bin/cmake .. \
      -DCGAL_DIR=/opt/${CGALDIR}/lib64/cmake/CGAL \
      -DTBB_DIR=/opt/tbb/lib/lib64/cmake/TBB  \
      -DCMAKE_BUILD_TYPE="Release"
RUN make clean 
RUN make -j4

###BIMPP LIBRARY###
RUN echo "export PATH=/usr/lib64/openmpi/bin:$PATH" >> ~/.bashrc
RUN echo "export LD_LIBRARY_PATH=/usr/lib64/openmpi/lib:$LD_LIBRARY_PATH" >> ~/.bashrc
ENV PATH=/usr/lib64/openmpi/bin:$PATH
ENV LD_LIBRARY_PATH=/usr/lib64/openmpi/lib:$LD_LIBRARY_PATH

RUN dnf install -y epel-release
RUN dnf --enablerepo crb install -y glpk-utils bison \
         openblas-devel autoconf automake libtool scotch scotch-devel \
         redhat-rpm-config libasan bzip2 bzip2-devel
###MUMPS###
RUN dnf install -y MUMPS MUMPS-devel

###OCTAVE###
RUN dnf --enablerepo crb install -y qhull-devel arpack-devel octave octave-devel

##LIS###
WORKDIR /opt/
ENV pkgname=lis
ENV pkgver=2.1.3
ENV archive=$pkgname-$pkgver.zip
RUN wget https://www.ssisc.org/lis/dl/$archive
RUN unzip $archive
WORKDIR /opt/$pkgname-$pkgver
RUN chmod 755 configure config/install-sh
RUN ./configure --prefix="/opt/$pkgname" \
                --enable-mpi \
                --enable-shared \
                --disable-static
RUN make -j4
RUN make -j4 install
WORKDIR /opt/
RUN rm -rf $pkgname-$pkgver $archive


###P4EST###
ENV pkgname=p4est
ENV pkgver=2.8.6
ENV archive=$pkgname-$pkgver.tar.gz
RUN wget https://p4est.github.io/release/$archive && \
    tar -xzf $archive && \
    rm -rf $archive
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
                CPPFLAGS="-I/usr/include/openmpi-x86_64/"
RUN make -j4
RUN make install
WORKDIR /opt/
RUN rm -rf $pkgname-$pkgver


###SCALAPACK###
ENV pkgname=scalapack
ENV pkgver=2.1.0
ENV archive=$pkgname-$pkgver.tgz
RUN wget https://www.netlib.org/$pkgname/$archive && \
    tar -xf $archive && \
    rm -rf $archive
WORKDIR /opt/$pkgname-$pkgver
RUN mkdir build
WORKDIR /opt/$pkgname-$pkgver/build
ENV FFLAGS="-std=legacy"
RUN /opt/cmake-3.19.0-Linux-x86_64/bin/cmake \
          -D CMAKE_INSTALL_PREFIX="/opt/$pkgname" \
          -D CMAKE_SKIP_INSTALL_RPATH=ON \
          -D CMAKE_SKIP_RPATH=ON \
          -D BUILD_TESTING=OFF \
          -D BUILD_SHARED_LIBS=ON \
          -D LAPACK_LIBRARIES=openblas \
          -D BLAS_LIBRARIES=openblas \
          -D MPI_BASE_DIR=/usr/lib64/openmpi \
          ..
RUN make -j4 VERBOSE=1
RUN make install
WORKDIR /opt/
RUN rm -rf /opt/$pkgname-$pkgver


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
RUN make -j4
RUN make install
WORKDIR /opt
RUN rm -rf /opt/$pkgname-$pkgver


###NextGenPB###
RUN mkdir /usr/local/nextgenPB
RUN git config --global \
  url."https://vdiflorio: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@github.com/".insteadOf \
  "https://github.com/"
RUN git clone https://github.com/vdiflorio/nextgenPB.git /usr/local/nextgenPB/


###BIMPP###
WORKDIR /opt
ENV pkgname=bimpp
ENV pkgver=patch
ENV archive=$pkgname-$pkgver.tar.bz2
#RUN cp /usr/local/nextgenPB/$archive /opt && mkdir $pkgname-patch && \
#    tar -xf $archive -C $pkgname-patch --strip-components=1 && \
#    rm -rf $archive
RUN cp /usr/local/nextgenPB/$archive /opt && \
    tar -xf $archive && \
    rm -rf $archive
WORKDIR /opt/$pkgname-patch
RUN mkdir build
RUN ./autogen.sh
WORKDIR /opt/$pkgname-patch/build
RUN ../configure --prefix=/opt/bimpp/ \
                CPPFLAGS="-I/usr/include/openmpi-x86_64 -I/usr/include/ -I/usr/include/MUMPS/ -I/opt/p4est/include -DOMPI_SKIP_MPICXX -DHAVE_OCTAVE_44 -DBIM_TIMING" \
                LDFLAGS="-L/usr/lib64 -L/usr/lib64/openmpi/lib -L/lib64" \
                 --with-blas-lapack="-lopenblas"  \
                 --with-octave_file_io-home=/opt/octave_file_io \
                 --with-octave-home=/usr/bin \
                 --with-p4est-home=/opt/p4est \
                 --with-lis-home=/opt/lis \
                 --with-mumps-home=/usr/lib64 \
F77=mpif90 CXX=mpicxx MPICC=mpicc CC=mpicc \
CXXFLAGS="-O3" \
--with-mumps-extra-libs="-L/usr/lib64 -L/lib64 -L/usr/lib64/lib -lscotch -lmpi \
 -lmpi_usempif08 -lmpi_usempi_ignore_tkr -lmpi_mpifh -lmpi -lopenblas  -lgfortran"
RUN make -j4
RUN make install
WORKDIR /opt


###NGPB###
WORKDIR  /usr/local/nextgenPB
RUN git checkout fast_density
RUN cp local_settings_docker.mk /usr/local/nextgenPB/src/local_settings.mk
WORKDIR  /usr/local/nextgenPB/src
RUN make clean all

