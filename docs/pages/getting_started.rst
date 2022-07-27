===============
Getting Started
===============

-----------
Quick Start
-----------

After downloading :code:`EXCEED-DMvX.Y.Z` from 

- `<https://github.com/tanner-trickle/EXCEED-DM/releases>`_

Run

.. code-block:: none

    tar -xvzf EXCEED-DMvX.Y.Z.tar.gz -C EXCEED-DMvX.Y.Z 
    cd EXCEED-DMvX.Y.Z
    mkdir build && cd build
    cmake ..
    make

--------
Download
--------

Release
=======

To download release version :code:`vX.Y.Z` directly, simply click on :code:`EXCEED-DMvX.Y.Z.tar.gz` in the Releases section of the Github page 

- `<https://github.com/tanner-trickle/EXCEED-DM/releases>`_.

Git Branches
============

To download the :code:`<branch>` branch of the source code from Github,

.. code-block:: none
    
    git clone -b <branch> https://github.com/tanner-trickle/EXCEED-DM.git ./EXCEED-DM

The source code will be in the :code:`EXCEED-DM` folder. Available branches are :code:`main` and :code:`develop`.

-------
Install
-------

Before you begin installation, make sure that you have all the necessary pre-requisite software installed. Commands are shown for installation on a Linux machine via :code:`apt`. For more information follow the links.

Pre-requisites
==============

* Fortran compiler - `gfortran <https://gcc.gnu.org/wiki/GFortran>`_, ifort, etc.

.. code-block:: none

   sudo apt install gfortran

**Libraries**

* `OpenMPI <https://www.open-mpi.org/>`_

.. code-block:: none

   sudo apt install libopenmpi-dev

* `FFTW3 <https://www.fftw.org/>`_

.. code-block:: none

   sudo apt install libfftw3-dev
* `HDF5 (serial, with Fortran support) <https://www.hdfgroup.org/downloads/hdf5/>`_

.. code-block:: none

   sudo apt install libhdf5-serial-dev

* `LAPACK <https://netlib.org/lapack/>`_
.. code-block:: none

   sudo apt install liblapack-dev
* `BLAS <https://netlib.org/blas/>`_
.. code-block:: none

   sudo apt install libblas-dev
* `CMake <https://cmake.org/>`_
.. code-block:: none

   sudo apt install cmake

Build
=====

After the pre-requisite sofware has been installed installed go to the folder where :code:`EXCEED-DM` was downloaded to (containing the :code:`README.md` file). From this folder, run

.. code-block:: none

    mkdir build && cd build
    cmake ..
    make

To delete the build simply delete the contents of the build folder,

.. code-block:: none

   rm -r build/*

Build Options
-------------

Build options are specified with flags when running the :code:`cmake` command, e.g.,

.. code-block:: none

    cmake .. -DCMAKE_BUILD_TYPE=DEBUG

The currently configured flags are:

* :code:`-DCMAKE_BUILD_TYPE`
    * :code:`RELEASE`: All optimizations turned on. 
    * :code:`DEBUG`: No optimizations, maximize the number of errors caught at runtime. 

Test
====

To check that :code:`EXCEED-DM` was installed correctly, from the :code:`EXCEED-DM` folder run

.. code-block:: none

   ./build/exdm

The output should read something like

.. code-block:: none

     --------------------------------------------------------------------------------

         EXCEED-DM - v1.0.0

         Running on 1 processors
         Compiled with GCC version 11.1.0

         Started at 16:52:19.532 7/26/2022

     --------------------------------------------------------------------------------

     No input file specified, aborting.
