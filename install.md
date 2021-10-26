# Installation

Instuctions on how to setup, and check the setup of, `EXCEED-DM` using CMake. 

## Build

To build `EXCEED-DM` to `/your/path/EXCEED-DM` with CMake,

    > tar -xvzf EXCEED-DM-vX.Y.Z.tar.gz -C /your/path/EXCEED-DM
    > cd /your/path/EXCEED-DM
    > mkdir build
    > cd build
    > cmake ..
    > make

### Build Options

CMake allows the user to change many of the build parameters at the command line. Simply append these options to the `cmake ..` command in the previous section.

- To specify a different Fortran compiler,

        -DCMAKE_Fortran_COMPILER=<your compiler>

- To specify a build mode, 

        -DCMAKE_BUILD_TYPE=<build-type>

    Options are: `Release`, `Debug`

### Cluster/Supercomputer

Large computing systems should have a way to load the pre-requisite software, e.g. `module load openmpi`. Once all of the prerequisites are loaded cmake should automatically find the correct file locations. 

Notes:
- You might have to specify the compiler directly with a build option.
- Your version of HDF5 should have been compiled with the same compiler you're using. There seems to be some wiggle room with different versions of the same compiler, but expect errors if you compiled HDF5 with GNU Fortran and want to compile the main program with an Intel compiler.

### Check Build

To check that the program was built run
    
    mpirun -np 2 ./build/exdm

from the main folder.

### Clean Build

To remove the current build of the program, run

    cmake --build . --target clean

from the `build/` folder.

## Pre-requisites

- Fortran compiler
- OpenMPI
- FFTW3
- HDF5
- LAPACK
- BLAS
- CMake (v>=3.16)

## Fortran90 compiler

### Linux

#### Install

    sudo apt install gfortran

#### Check Install

    gfortran --version

### Mac

#### Install

#### Check Install

    gfortran --version
    
- Note: attempt compiling a seperate test program, if an error about **-lSystem** shows up try
    
    brew reinstall gcc 

## OpenMPI

### Linux

#### Install

    sudo apt install libopenmpi-dev

#### Check Install

    mpirun --version 

### Mac

#### Install

Follow the instructions at [https://www.open-mpi.org/faq/?category=building#easy-build](https://www.open-mpi.org/faq/?category=building#easy-build).

#### Check Install

    mpirun --version
    
## FFTW3

### Linux

#### Install

    sudo apt install libfftw3-dev

#### Check Install

    sudo find / -name libfftw3*

Should return a list of libraries, but could take a while. Try checking likely folders first like **/usr**.

### Mac

#### Install

    brew install fftw

#### Check Install

    sudo find / -name "libfftw3*"

## HDF5

### Linux

#### Install 

    sudo apt install libhdf5-dev

#### Check Install

    h5fc -showconfig

### Mac

#### Install

Download the latest release of HDF5 from [https://www.hdfgroup.org/downloads/hdf5/](https://www.hdfgroup.org/downloads/hdf5/). Unzip the file with 

    gunzip < hdf5-X.Y.Z.tar.gz | tar xf -       # replace X.Y.Z with the release number 

Then, to install HDF5 at location /usr/local/hdf5, run

    cd hdf5-X.Y.Z           # replace X.Y.Z with the release number
    ./configure --prefix=/usr/local/hdf5 --enable-fortran
    make
    make check              # optional
    make install 
    make check-install      # optional

#### Check Install

    h5fc -showconfig
    
## CMake

### Linux

#### Install

    sudo apt install cmake

#### Check Install

    cmake --version

