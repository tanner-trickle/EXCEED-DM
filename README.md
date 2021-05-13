# EXCEED-DM (EXtended Calculation of Electronic Excitations for Direct detection of Dark Matter)

[![License: GPL v3](https://img.shields.io/badge/License-GPLv3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)
[![DOI](https://zenodo.org/badge/354900532.svg)](https://zenodo.org/badge/latestdoi/354900532)

**NOTICE: This program is in active development and the code, as well as required inputs and outputs, are subject to change. Star and watch the repository to stay up to date!**

`EXCEED-DM` is an OpenMPI Fortran program that computes Dark Matter (DM) induced electronic transition rates for any gapped crystal target. 

## Features

`EXCEED-DM` provides a complete framework for computing DM-electron interaction rates. Given a range of DM masses, mediator form factors, and times of day, `EXCEED-DM` computes the scattering rate per kg-year binned in energy and momentum deposition.

- **All kinematically allowed transitions are included.** In addition to transitions near the Fermi surface, where density functional theory (DFT) calculations are a necessary component of the calculation, electronic states further below, and above, are modeled semi-analytically and included in the scattering rate calculation.
- **Compute daily and annual modulation signals.** No assumptions about the isotropy of the target are made, allowing one to study the daily modulation signal, and the DM velocity distribution parameters can be changed very easily to compute annual modulation. 

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

Follow these instructions to compile and run `EXCEED-DM` on a fresh Ubuntu (18.04, 20.04) distribution. For installation on other systems, or if something goes wrong, see `install.md` for more detailed instructions.

1) Install preliminary software
    - Fortran compiler (`sudo apt install gfortran`) 
    - OpenMPI (`sudo apt install libopenmpi-dev`)
    - fftw3 (`sudo apt install libfftw3-dev`)
    - hdf5 (`sudo apt install libhdf5-serial-dev`)
    - FoBiS.py (`sudo pip3 install FoBiS.py`)

Note : It's recommended to run `sudo apt update` before, and on a completely fresh Ubuntu installation pip will need to be installed (`sudo apt install python3-pip`).

2) Download, then extact the latest release with
    
        > tar -xvzf EXCEED-DM-vX.Y.Z.tar.gz -C /your/specific/path/

3) From `/your/specific/path`, compile the main program, `exdm`,

        > FoBiS.py build -mode ubuntu-gnu

To compile `exdm` on a cluster a few file paths specific to the cluster will have to be put in to the `fobos` file. See that file and `install.md` for more instructions.

4) Test the installation

        > mpirun -np 2 ./build/exdm

If installed correctly you should see something similar to,

         --------------------

            EXCEED-DM - v0.1.0

         --------------------

         Running on            2 processors
         Compiled with GCC version 9.3.0

         ----------

         !! ERROR !!

            Input file for control parameters :  does NOT exist.

         !!!!!!!!!!!

5) Perform example calculation. Open `examples/Si/inputs/vc_test_input.txt` and change the file paths as necessary (DFT input file is stored in `./examples/Si/dft/Si/Si_2x2x2_AE.hdf5`). Then run

        > mpirun -np 2 ./build/exdm ./examples/Si/inputs/vc_test_input.txt

## Support 

- Installation instructions can be found in `install.md`.
- More detailed usage intstructions can be found in the user manual (**in preparation**).
- Input files needed for example calculations can be found in `examples/`.
- Output files of example calculations can be found in `examples/<material name>/outputs/` and here:
    - Si/Ge - [![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.4737654.svg)](https://doi.org/10.5281/zenodo.4737654) 
- Documentation can be found [here](https://tanner-trickle.github.io/EXCEED-DM), as well as the `docs/` folder.
- Larger DFT input files can be found here: 
    - Si/Ge - [![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.4735777.svg)](https://doi.org/10.5281/zenodo.4735777)

## Papers using `EXCEED-DM`

- T. Trickle, Z. Zhang, K. M. Zurek, K. Inzani and S. Griffin, [Multi-Channel Direct Detection of Light Dark Matter: Theoretical Framework]
- S. M. Griffin, K. Inzani, T. Trickle, Z. Zhang and K. M. Zurek, [Multichannel direct detection of light dark matter: Target comparison]

[Multi-Channel Direct Detection of Light Dark Matter: Theoretical Framework]: https://doi.org/10.1007/JHEP03(2020)036
[Multichannel direct detection of light dark matter: Target comparison]: https://doi.org/10.1103/PhysRevD.101.055004

## Attribution

If you use `EXCEED-DM` in your work, please cite,

[![DOI](https://zenodo.org/badge/354900532.svg)](https://zenodo.org/badge/latestdoi/354900532)

along with,

    @article{griffin2021extended,
      title={Extended Calculation of Dark Matter-Electron Scattering in Crystal Targets}, 
      author={Sinéad M. Griffin and Katherine Inzani and Tanner Trickle and Zhengkang Zhang and Kathryn M. Zurek},
      year={2021},
      eprint={2105.05253},
      archivePrefix={arXiv},
      primaryClass={hep-ph}
    }

