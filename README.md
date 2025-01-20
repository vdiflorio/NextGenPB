# NextGenPB 

Linearized Poisson-Boltzmann solver on Octree mesh.

Solves the linearized PBE 

$$
-\mathrm{div} \left( \varepsilon_0 \varepsilon_r \nabla \varphi \right) + \kappa^2 \varphi = \rho^f
$$

on a rectangular domain.

---

# Table of Contents
1. [Installation](#installation)
   1.1 [Dependencies](#dependencies)
   1.2 [Setup](#setup)
2. [Usage](#usage)
3. [Contributing](#contributing)

---
# Installation
The following dependencies are required to build and run the program:  

### Core Libraries

- **`lis`**
- **`p4est`**
- **`bim++`**  
- **`NanoShaper`**  

### Post-Processing Tools

- **Paraview**  
- **GNU Octave**


## Installation on macOS
### Prerequisites

Use MacPorts (or Homebrew) to install the required packages: 

```bash
sudo port install openmpi cmake onetbb boost cgal5 nlohmann-json jansson octave
sudo port install mumps +openmpi -mpich
sudo port install lis +openmpi -mpich
sudo port install p4est +openmpi -mpich
```

### Installing Additional Dependencies
**`NanoShaper`**
Clone and build the NanoShaper library:
```bash
git clone https://gitlab.iit.it/SDecherchi/nanoshaper.git
cd nanoshaper
cp CMakeLists_so.txt CMakeLists.txt
cd build_lib
cmake .. -DCGAL_DIR=/opt/local/ -DCMAKE_BUILD_TYPE=Release
make
```

**`octave_file_io`**
Clone and build the octave_file_io library:
```bash
git clone https://github.com/carlodefalco/octave_file_io.git
cd octave_file_io
./autogen.sh
mkdir build
cd build
../configure CXX=mpicxx --prefix=/opt/octave_file_io/1.0.91 --with-octave-home=/opt/local/bin 'LDFLAGS=-Wl,-rpath -Wl,/opt/local/lib/libgcc -Wl,-rpath -Wl,/opt/local/lib/gcc13 -ld_classic'
make
sudo make install
```

**`bim++`**
Clone and build the bim++ library:
```bash
git clone https://github.com/carlodefalco/bimpp.git
cd bimpp
git checkout nextgenPB
./autogen.sh
mkdir build
cd build
../configure --prefix=/opt/bimpp --disable-option-checking \
LDFLAGS="-L/opt/local/lib -Wl,-rpath,/opt/local/lib/libgcc -Wl,-rpath,/opt/local/lib/gcc14" \
CPPFLAGS="-I/opt/local/include/ -I/opt/local/include/gcc14 -DOMPI_SKIP_MPICXX -DHAVE_OCTAVE_44 -DBIM_TIMING" \
 --with-blas-lapack="-lopenblas"  \
 --with-octave_file_io-home=/opt/octave_file_io/1.0.91 \
 --with-octave-home=/opt/local/bin \
 --with-p4est-home=/opt/local \
 --with-lis-home=/opt/local \
 --with-mumps-home=/opt/local \
 F77=mpif90 CXX=mpicxx MPICC=mpicc CC=mpicc \
CXXFLAGS="-std=c++17 -O3 -mtune=native" \
--with-mumps-extra-libs="-L/opt/local/lib -lptscotch -lscotch -lmpi -Wl,-flat_namespace -Wl,-commons,use_dylibs \
-L/opt/local/lib/openmpi-mp -lmpi_usempif08 -lmpi_usempi_ignore_tkr -lmpi_mpifh -lopenblas -L/opt/local/lib/gcc14 -lgfortran"

make
sudo make install
```

## Installation on Rocky Linux and Ubuntu



To compile the program, create a `local_settings.mk` file to add or modify the compilation options specified in the `Makefile`. Then, run:  


```bash
make
```

### command line options

To run the code in serial, using the default `options.pot` file and the `1CCM.pqr` file, 
the necessary commands are

```
mkdir results
cd results
../poisson_boltzmann 
```
 
To run in parallel with MPI, change the last line with

```
mpirun -np <number of processes> ../poisson_boltzmann
```

with `<number of processes>` equal to the desired number of processors. 

To select a different option file add the following option

```
--potfile /your/folder/<file name>.pot
```

While to change the input molecule add the following option

```
--pqrfile /your/folder/<file name>.pqr
```
Files .pqr must not contain the optional feature `CHAIN_ID`, because the code cannot read that column. For example, .pqr files created with PDB2PQR and the –whitespace option are guaranteed to conform to the correct format:

```
pdb2pqr --ff=charmm --whitespace 4ake.pdb 4ake.pqr
```

without the option `--keep-chain`. Otherwise, if the user utilizes online website to generate the file .pqr, it has to be sure to un-tick the option `add/keep chain IDs in the PQR file`.

### Post-processing of Output Files

The user can use the `GNU Octave` script called `export.m` and 
distributed with `bimpp` (it is located in the  folder `script/m`)
for the post-processing. When the user adopts the default names for the files, it 
can write the following command line for the creation of a .vtu file that can be 
opened using Paraview. 

```
   export_tmesh_data ('poisson_boltzmann_surface_%1.1d_%4.4d',
                       {'poisson_boltzmann_surface_%1.1d_%4.4d', 'phi_%1.1d_%4.4d', 
                       'rhs_%1.1d_%4.4d', 'rho_%1.1d_%4.4d'}, {'surface', 'phi', 'rhs', 'rho'},
                       {'poisson_boltzmann_marker_%1.1d_%4.4d', 'epsilon_%1.1d_%4.4d',
                       'reaction_%1.1d_%4.4d'}, {'marker', 'epsilon', 'reaction'},
                       'output_name', 0, processors_number)
```

where `output_name` can be substituted with an arbitrary name for the `.vtu ` file and `processors_number` 
must be replaced by the number of processors used for the simulation. 

The user can also decide to convert only some among the outputs according to his necessities. 

### Projects

Below an (unordered) list of improvement/fixes/new features that we expect to implement in the near future

1. Separate the effect of each residue on the energy
2. Account for the *Stern Layer* in the model
3. Singularity removal
4. Nonlinear ionic density
5. Coulombi boundary conditions + other choices for boundary conditions
6. Map of Electrostatic potential on triangulated molecular surface (*.offc* file format)
7. Allow more flexible format for *.pqr* files (in particular with or without ChainID column)
