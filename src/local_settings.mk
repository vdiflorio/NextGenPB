CPPFLAGS=\
-I/opt/octave_file_io/1.0.91/include/ \
-I/opt/p4est/2.8.5/include/ \
-I/opt/nanoshaper/src \
-I../include -I../addons \
-I/opt/bimpp1/include \
-I/opt/local/include \
-I/opt/lis/2.1.3/include \
-DHAVE_OCTAVE_44 -DOMPI_SKIP_MPICXX -DBIM_TIMING -DUSE_MPI

CXXFLAGS= -std=c++17 -ggdb  -O3 -fsanitize=address
#CXXFLAGS= -Ofast -mtune=native

LDFLAGS=-L/opt/octave_file_io/1.0.91/lib \
-L/opt/p4est/2.8.5/lib \
-L/opt/lis/2.1.3/lib \
-L/opt/nanoshaper/build \
-L/opt/bimpp1/lib \
-L/opt/p4est/2.8.5/lib \
-Wl,-rpath,/opt/nanoshaper/build \
-Wl,-rpath,/opt/bimpp1/lib  \
-Wl,-rpath,/opt/octave_file_io/1.0.91/lib \
-Wl,-rpath,/opt/lis/2.1.3/lib \
-Wl,-rpath,/opt/p4est/2.8.5/lib

LIBS=-ld_classic -lNanoShaper -lbim -lbimmumps -lbimlis -lbimp4est -lbimlinalg -llis -ldmumps -lmumps_common \
-lscotcherr -lbz2 -lmpi_usempif08 -lmpi_usempi_ignore_tkr -lmpi_mpifh -lmpi -lp4est -lsc
