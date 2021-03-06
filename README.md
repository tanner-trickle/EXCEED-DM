<p align="center">
    <img src="https://github.com/tanner-trickle/EXCEED-DM/blob/develop/docs/media/exdm-prelim-logo.png?raw=true" alt= "EXCEED-DM-logo"/>
</p>

# EXCEED-DM (EXtended Calculation of Electronic Excitations for Direct detection of Dark Matter)

[![License: GPL v3](https://img.shields.io/badge/License-GPLv3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)
[![DOI](https://zenodo.org/badge/354900532.svg)](https://zenodo.org/badge/latestdoi/354900532)

**NOTICE: This program is in active development and the code, as well as required inputs and outputs, are subject to change. Star and watch the repository to stay up to date!**

`EXCEED-DM` is an OpenMPI Fortran program that computes Dark Matter (DM) induced electronic transition rates for any gapped crystal target. 

[**Home**](https://exceed-dm.caltech.edu)

[**Documentation**](https://tanner-trickle.github.io/EXCEED-DM)

## Features

`EXCEED-DM` provides a complete framework for computing DM-electron interaction rates. 

- **Scattering**: Given a range of DM masses, mediator form factors, and times of day, `EXCEED-DM` computes the scattering rate per kg-year binned in energy and momentum deposition.
    - **All kinematically allowed transitions are included.** In addition to transitions near the Fermi surface, where density functional theory (DFT) calculations are a necessary component of the calculation, electronic states further below, and above, are modeled semi-analytically and included in the scattering rate calculation.
    - **Daily and annual modulation signals**: No assumptions about the isotropy of the target are made, allowing one to study the daily modulation signal, and the DM velocity distribution parameters can be changed very easily to compute annual modulation. 
    - **Spin-dependent wave functions**: Some targets, such as those with spin-orbit coupling, will have electronic states which are not eigenstates of the spin operator. This means the wave functions have two components instead of one, and particle physics couplings to the spin operator are no longer trivial to compute. `EXCEED-DM` can perform these spin-dependent scattering rate calculations.
      - Note: currently (v >= 0.2.0) only valence to conduction transitions are supported.
    - **Calculate the dielectric**: For some processes the dielectric will screen the interaction rate. The complex dielectric/screening factor can now be computed and used in scattering rate calculations.

- **Absorption**: Given a range of DM masses and times of day, `EXCEED-DM` computes the absorption rate, electronic self energies, and generalized crystal form factors, needed to compute the DM absorption rate.
    - **Scalar, pseudoscalar, and vector DM**: compute the absorption rate from these bosonic DM candidates.
    - Note: currently (v >= 0.2.0) only valence to conduction transitions are supported.

`EXCEED-DM` is :

- **fast** - Being parallelized with OpenMPI means that `EXCEED-DM` can take full advantage of large computing clusters. 
- **DFT calculator independent** - `EXCEED-DM` depends only on the output of DFT calculations, which means that any DFT calculator can be used to compute the targets electronic properties and then used as input. An example converter for `VASP` output can be found in the `utilities` folder. This also allows for modular work flows where particle physicists don't have to worry about DFT details, and material scientists don't have to worry about particle physics details.

and has more features on the way!

## Extras

Inside the `utilities` folder are other programs meant to aid in using `EXCEED-DM`:

- `vasp_converter/` - python program to convert the output files from VASP calculations to the input files needed for `EXCEED-DM`. Performs the all electron reconstruction with `pawpyseed`.
- `binned_wfc/` - Fortran program to compute the square magnitudes of the Bloch wave functions, binned in momentum transfer.  
- `post_analysis/` - python scripts and example notebooks for post-processing the output of `EXCEED-DM`.
- `core_elec_config/` - Calculations involving core electrons require a 'core electron configuration' file that depends on the target material. These files can be generated with the notebook in `core_elec_config/` by simply passing the [Materials Project](https://materialsproject.org/) ID.

## Getting Started

Follow these instructions to compile and run `EXCEED-DM` on a fresh Linux distribution. For installation on other systems, or if something goes wrong, see `install-cmake.md` for more detailed instructions.

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

            EXCEED-DM - v0.3.0

            Running on 2 processors
            Compiled with GCC version 9.3.0

        ----------------------------------------------------------------------

        Loading control parameters...

        ----------------------------------------------------------------------
            -------
            Control
            -------

## Support 

- Installation instructions can be found in `install-cmake.md`.
- More detailed usage intstructions can be found in the user manual (**in preparation**).
- Input files needed for example calculations can be found in `examples/`.
- Output files of example calculations can be found in `examples/<material name>/outputs/` and here:
    - Si/Ge - [![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.4737654.svg)](https://doi.org/10.5281/zenodo.4737654) 
- Documentation can be found [here](https://tanner-trickle.github.io/EXCEED-DM), as well as the `docs/` folder.
- Larger DFT input files can be found here: 
    - Si/Ge - [![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.4735777.svg)](https://doi.org/10.5281/zenodo.4735777)
    - Si/Ge (compatible with v>=0.2.8) - [![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.6015637.svg)](https://doi.org/10.5281/zenodo.6015637)


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

