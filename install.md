# Installation

Instuctions on how to setup, and check the setup of, `dm-electron` on a local machine. For setup on a cluster ignore the **install** instructions, but make sure to check that the cluster has all the packages installed/modules loaded by running the commands in **check install**

## Pre-requisites

- Fortran90 compiler (gfortran, ifort)
- OpenMPI
- FFTW3
- HDF5
- FoBiS.py

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

Follow the instructions at [https://www.open-mpi.org/faq/?category=building#easy-build](https://www.open-mpi.org/faq/?category=building#easy-build).

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
    
## FoBiS.py

### Linux/Mac

#### Install

    pip install FoBiS.py

#### Check Install

    FoBiS.py --version

## Build

To build the program we use **FoBiS.py** which is a simple build systerm for Fortran code. **FoBiS** requires a **fobos** file which specifies the build options. There is an example in the main folder. The user must specify the location of the hdf5 libraries and fftw libraries.

### Step 0: Download dm-electron

    git clone https://github.com/tanner-trickle/dm-electron
    cd dm-electron

### Step 1: Specify file locations

Inside the **fobos** file there are four variables which need to be specified:
    
    $hdf5_inc_dir
    $hdf5_lib_dir
    $fftw3_inc_dir
    $fftw3_lib_dir

(For running on a cluster simply make new copies of these variables as done for **ubuntu** and **caltech_hpc** in the **fobos** file.)

The directories can be found by finding the location of four files (assuming these weren't previously moved upon installation),

- hdf5 include (inc) directory -> directory of `hdf5.mod`
- hdf5 library (lib) directory -> directory of `libhdf5.a`
- fftw3 include directory -> directory of `fftw3.f03` 
- fftw3 library directory -> directory of `libfftw3.a`

### Step 2: Compile

Finally, run

    FoBiS.py build -mode local-gnu

assuming your Fortran compiler is **gfortran**.

### Check Build

To check that the program was built run
    
    mpirun -np 2 ./build/dme 

### Clean Build

To remove the current build of the program, run

    FoBiS.py clean
