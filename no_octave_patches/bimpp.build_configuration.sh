#!/bin/sh
# Build bimpp locally, WITHOUT Octave, against the no-octave octave_file_io shim.
#
# REPO is the repository root (this script lives in <REPO>/bimpp/build). The
# path may contain '&', so every value is double-quoted.
REPO="$(cd "$(dirname "$0")/../.." && pwd)"

# NB: prefix is 'install_local', not 'install': on a case-insensitive macOS
# filesystem 'install' collides with bimpp's existing INSTALL doc file.
../configure --prefix="$REPO/bimpp/install_local" --disable-3d-interpolation  --disable-option-checking \
LDFLAGS="-L/opt/local/lib -Wl,-rpath,/opt/local/lib/libgcc -Wl,-rpath,/opt/local/lib/gcc14" \
CPPFLAGS="-I/opt/local/include/ -I/opt/local/include/gcc14 -DOMPI_SKIP_MPICXX -DOFIO_NO_OCTAVE -DBIM_TIMING" \
 --with-blas-lapack="-lopenblas"  \
 --with-octave_file_io-home="$REPO/octave_file_io/install_no_octave" \
 --with-octave-home=/nonexistent \
 --with-p4est-home=/opt/p4est/2.8.7 \
 --with-lis-home=/opt/lis/2.1.6 \
 --with-mumps-home=/opt/local \
 F77=mpif90 CXX=mpicxx MPICC=mpicc CC=mpicc \
CXXFLAGS="-std=c++17 -O3 -mtune=native -march=native" \
--with-mumps-extra-libs="-L/opt/local/lib -lptscotch -lscotch -lmpi -Wl,-flat_namespace -Wl,-commons,use_dylibs \
-L/opt/local/lib/openmpi-mp -lmpi_usempif08 -lmpi_usempi_ignore_tkr -lmpi_mpifh -lopenblas -L/opt/local/lib/gcc14 -lgfortran"
