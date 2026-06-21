# local_settings for NextGenPB built WITHOUT Octave (Ubuntu container).
#
# Same as local_settings_ubuntu.mk but for the no-Octave chain:
#   - octave_file_io is the text-backend shim built with -DOFIO_NO_OCTAVE
#     (liboctave_file_io.so has NO dependency on real Octave).
#   - bimpp is built against that shim (configure.ac no-Octave branch).
#   - no mkoctfile, no octave includes, no -loctave/-loctinterp, no HAVE_OCTAVE_44.
#
# The 'x86_64-linux-gnu' token is rewritten to the real DEB_HOST_MULTIARCH triplet
# by the Dockerfile (so it links on arm64 too).
CPPFLAGS=\
-I/usr/include/ \
-I/opt/octave_file_io/include/ \
-I/opt/nanoshaper/src \
-I/usr/local/nextgenPB/include -I/usr/local/nextgenPB/addons \
-I/opt/bimpp/include \
-I/opt/lis/include \
-DOFIO_NO_OCTAVE -DOMPI_SKIP_MPICXX -DBIM_TIMING -DUSE_MPI


CXXFLAGS= -O2 -mtune=generic
#CXXFLAGS= -Ofast -mtune=native -march=native
CXX=mpicxx

LDFLAGS=-L/opt/octave_file_io/lib \
-L/usr/lib/x86_64-linux-gnu \
-L/opt/lis/lib \
-L/opt/nanoshaper/build_lib \
-L/opt/bimpp/lib \
-Wl,-rpath,/opt/nanoshaper/build_lib \
-Wl,-rpath,/opt/bimpp/lib  \
-Wl,-rpath,/opt/octave_file_io/lib \
-Wl,-rpath,/opt/lis/lib

LIBS=-lNanoShaper -lbim -lbimmumps -lbimlis -lbimp4est \
     -lbimlinalg -loctave_file_io \
     -llis -ldmumps -lmumps_common \
     -lscotcherr -lbz2 -lmpi_usempif08 \
     -lmpi_usempi_ignore_tkr -lmpi_mpifh -lmpi -lp4est -lsc
