CXX=mpic++

CXXFLAGS= -std=c++17 -ggdb -O0

CPPFLAGS=-I$(mkP4estInc) -I$(mkMumpsInc) -I$(mkLisInc) -I/u/software/octave_file_io/1.0.91/include/ \
-I/u/software/octave/6.4.0/include/octave-6.4.0 \
-I/u/software/nanoshaper/nanoshaper/src -I../include -I../addons \
-I/u/software/bimpp/dev/include -DHAVE_OCTAVE_44 -DOMPI_SKIP_MPICXX -DBIM_TIMING -DUSE_MPI

LDFLAGS=-L/u/software/octave_file_io/1.0.91/lib \
-L/u/software/nanoshaper/nanoshaper/build_so \
-L/u/software/bimpp/dev/lib -L$(mkLisLib) -L$(mkMumpsLib) -L$(mkScotchLib)

LIBS=-lNanoShaper -lbim -lbimmumps -lbimlis -lbimp4est -lbimlinalg -llis -ldmumps -lmumps_common \
-lscotcherr -lbz2 -lmpi_usempif08 -lmpi_usempi_ignore_tkr -lmpi_mpifh -lmpi

include local_settings.mk

all : poisson_boltzmann

%.o : %.cpp
	$(CXX) -c $(CXXFLAGS) $(CPPFLAGS) $<

poisson_boltzmann : poisson_boltzmann.o $(patsubst %.cpp, %.o, $(wildcard *.cpp))
	$(CXX) $(CXXFLAGS) $(LDFLAGS) $? -o $@ $(LIBS)

.PHONY : clean distclean

clean :
	$(RM) *.o

distclean : clean
	$(RM) poisson_boltzmann


