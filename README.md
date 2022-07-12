<p align="center">
    <img src="https://github.com/tanner-trickle/EXCEED-DM/blob/develop/docs/media/exdm-prelim-logo.png?raw=true" alt= "EXCEED-DM-logo"/>
</p>

# EXCEED-DM (EXtended Calculation of Electronic Excitations for Direct detection of Dark Matter)

[![License: GPL v3](https://img.shields.io/badge/License-GPLv3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)
[![DOI](https://zenodo.org/badge/354900532.svg)](https://zenodo.org/badge/latestdoi/354900532)

**NOTICE: This program is in active development and the code, as well as required inputs and outputs, are subject to change. Star and watch the repository to stay up to date!**

`EXCEED-DM` is an OpenMPI Fortran program that computes Dark Matter (DM) induced electronic transition rates. 

[**Home**](https://exceed-dm.caltech.edu)

[**Documentation**](https://tanner-trickle.github.io/EXCEED-DM) (undergoing revision/conversion to Sphinx)

## Features

`EXCEED-DM` provides a complete framework for computing DM-electron interaction rates. Given an electronic configuration, `EXCEED-DM` computes the relevant electronic matrix elements, then particle physics specific rates from these matrix elements. This allows for separation between approximations regarding the electronic state configuration, and the specific calculation being performed. 

- **Variety of Electronic State Approximations**
    - `EXCEED-DM` provides a variety of different bases to specify the initial, and final, electronic states to accurately characterize the target. For example, plane wave (PW) (spin-dependent) bases and Slate type orbital (STO) bases.

- **Calculations**

    - **Scattering**: Given a DM model (DM masses, interaction potential) `EXCEED-DM` computes the expected number of interactions per kg-year binned in energy and momentum deposition.
        - With the electronic state approximations separated from the scattering rate calculation, new transition form factors can be easily added as functions of the electronic matrix elements.
    - **Absorption**: Given a DM model (e.g., scalar, pseudoscalar, vector DM) `EXCEED-DM` computes the expected number of interactions per kg-year.
    - **Dielectric**: For some processes the dielectric will screen the interaction rate. The complex dielectric can be computed and subsequently used in scattering rate calculations.

    - **Check out the examples folder to see specific calculations in action!**

`EXCEED-DM` is :

- **fast** - Being parallelized with OpenMPI means that `EXCEED-DM` can take full advantage of large computing clusters. 
- **DFT calculator independent** - `EXCEED-DM` depends only on the output of DFT calculations, which means that any DFT calculator can be used to compute the targets electronic properties and then used as input. 

## Getting Started

Follow these instructions to compile and run `EXCEED-DM` on a fresh Linux distribution. For installation on other systems, or if something goes wrong, see the `Getting Started` section of the documentation.

1) Install preliminary software

- Fortran compiler (`sudo apt install gfortran`)
- OpenMPI (`sudo apt install libopenmpi-dev`)
- fftw3 (`sudo apt install libfftw3-dev`)
- hdf5 (`sudo apt install libhdf5-serial-dev`)
- LAPACK (`sudo apt install liblapack-dev`)
- BLAS (`sudo apt install libblas-dev`)
- CMake (`sudo apt install cmake`)

Note : It's recommended to run `sudo apt update` before, and on a completely fresh Linux distribution.

2) Download, then extact the latest release with
    
        > tar -xvzf EXCEED-DM-vX.Y.Z.tar.gz -C /your/specific/path/

3) Compile the main program, `exdm`,

        > cd /your/specific/path
        > mkdir build
        > cd build
        > cmake ..
        > make

4) Test the installation from `/your/specific/path',

        > mpirun -np 2 ./build/exdm ./examples/Si/inputs/scatter_vc_test_input.txt

If installed correctly you should see something similar to,

        ----------------------------------------------------------------------

            EXCEED-DM - v1.0.0

            Running on 2 processors
            Compiled with GCC version 9.3.0

        ----------------------------------------------------------------------

## Papers using `EXCEED-DM`

- CDEX Collaboration, [Constraints on Sub-GeV Dark Matter-Electron Scattering from the CDEX-10 Experiment]
- H-Y. Chen, A. Mitridate, T. Trickle, Z. Zhang, M. Bernardi, K. M. Zurek, [Dark Matter Direct Detection in Materials with Spin-Orbit Coupling]
- A. Mitridate, T. Trickle, Z. Zhang, K. M. Zurek, [Dark Matter Absorption via Electronic Excitations]
- S. M. Griffin, K. Inzani, T. Trickle, Z. Zhang and K. M. Zurek, [Extended Calculation of Dark Matter-Electron Scattering in Crystal Targets]
- T. Trickle, Z. Zhang, K. M. Zurek, K. Inzani and S. Griffin, [Multi-Channel Direct Detection of Light Dark Matter: Theoretical Framework]
- S. M. Griffin, K. Inzani, T. Trickle, Z. Zhang and K. M. Zurek, [Multichannel direct detection of light dark matter: Target comparison]

[Constraints on Sub-GeV Dark Matter-Electron Scattering from the CDEX-10 Experiment]: https://arxiv.org/abs/2206.04128 
[Dark Matter Direct Detection in Materials with Spin-Orbit Coupling]: https://arxiv.org/abs/2202.11716
[Dark Matter Absorption via Electronic Excitations]: https://link.springer.com/article/10.1007/JHEP09(2021)123 
[Extended Calculation of Dark Matter-Electron Scattering in Crystal Targets]: https://journals.aps.org/prd/abstract/10.1103/PhysRevD.104.095015
[Multi-Channel Direct Detection of Light Dark Matter: Theoretical Framework]: https://doi.org/10.1007/JHEP03(2020)036
[Multichannel direct detection of light dark matter: Target comparison]: https://doi.org/10.1103/PhysRevD.101.055004

## Attribution

If you use `EXCEED-DM` in your work, please cite,

[![DOI](https://zenodo.org/badge/354900532.svg)](https://zenodo.org/badge/latestdoi/354900532)

along with,

    @article{Griffin:2021znd,
        author = "Griffin, Sin\'ead M. and Inzani, Katherine and Trickle, Tanner and Zhang, Zhengkang and Zurek, Kathryn M.",
        title = "{Extended Calculation of Dark Matter-Electron Scattering in Crystal Targets}",
        eprint = "2105.05253",
        archivePrefix = "arXiv",
        primaryClass = "hep-ph",
        month = "5",
        year = "2021"
    }
