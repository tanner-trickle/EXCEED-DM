<p align="center">
    <img src="https://github.com/tanner-trickle/EXCEED-DM/blob/develop/docs/media/exdm-prelim-logo.png?raw=true" alt= "EXCEED-DM-logo"/>
</p>

# EXCEED-DM (EXtended Calculation of Electronic Excitations for Direct detection of Dark Matter)

[![DOI](https://zenodo.org/badge/354900532.svg)](https://zenodo.org/badge/latestdoi/354900532)

`EXCEED-DM` is an OpenMPI Fortran program that computes Dark Matter (DM) induced electronic transition rates. 

## Overview

[**Documentation**](https://tanner-trickle.github.io/EXCEED-DM)

[**User Manual**](https://arxiv.org/abs/2210.14917)

`EXCEED-DM` provides a complete framework for computing DM-electron interaction rates. Given an electronic configuration, `EXCEED-DM` computes the relevant electronic matrix elements, then particle physics specific rates from these matrix elements. This allows for separation between approximations regarding the electronic state configuration, and the specific calculation being performed. 

- To lean how to install and run `EXCEED-DM`, see the documentation [documentation](https://tanner-trickle.github.io/EXCEED-DM).

- For more detailed information about what `EXCEED-DM` is computing, and to see it used for new physics results, see the [user manual](https://arxiv.org/abs/2210.14917).

    - Results from the [user manual](https://arxiv.org/abs/2210.14917) can be found in a [Zenodo repository](https://zenodo.org/record/7250090#.Y1gwxLaSnZc) 

    - For the large Si and Ge electronic configuration files, used in the [user manual](https://arxiv.org/abs/2210.14917), see this [Zenodo repository](https://zenodo.org/record/7246141#.Y1cKIbaSnZc).

## Features

- **Variety of Electronic State Approximations**
    - `EXCEED-DM` provides a variety of different bases to specify the initial, and final, electronic states to accurately characterize the target. For example, plane wave (PW) (spin-dependent) bases and Slate type orbital (STO) bases.

- **Calculations**

    - **Scattering**: Given a DM model (DM masses, interaction potential, etc.) `EXCEED-DM` computes the expected number of interactions per kg-year binned in energy and momentum deposition.
        - With the electronic state approximations separated from the scattering rate calculation, new transition form factors can be easily added as functions of the electronic matrix elements.
    - **Absorption**: Given a DM model (e.g., scalar, pseudoscalar, vector DM) `EXCEED-DM` computes the expected number of interactions per kg-year.
    - **Dielectric**: For some processes the dielectric will screen the interaction rate. The complex dielectric can be computed and subsequently used in scattering rate calculations.

    - **Check out the examples folder to see specific calculations in action!**

`EXCEED-DM` is :

- **fast** - Being parallelized with OpenMPI means that `EXCEED-DM` can take full advantage of large computing clusters. 
- **DFT calculator independent** - `EXCEED-DM` depends only on the output of DFT calculations, which means that any DFT calculator can be used to compute the targets electronic properties and then used as input. 

## Papers using `EXCEED-DM`

- T. Trickle [EXCEED-DM: Extended Calculation of Electronic Excitations for Direct Detection of Dark Matter]
- B. Dutta, S. Ghosh, T. Li, A. Thompson, A. Verma [Non-standard neutrino interactions in light mediator models at reactor experiments]
- CDEX Collaboration, [Constraints on Sub-GeV Dark Matter-Electron Scattering from the CDEX-10 Experiment]
- H-Y. Chen, A. Mitridate, T. Trickle, Z. Zhang, M. Bernardi, K. M. Zurek, [Dark Matter Direct Detection in Materials with Spin-Orbit Coupling]
- A. Mitridate, T. Trickle, Z. Zhang, K. M. Zurek, [Dark Matter Absorption via Electronic Excitations]
- S. M. Griffin, K. Inzani, T. Trickle, Z. Zhang and K. M. Zurek, [Extended Calculation of Dark Matter-Electron Scattering in Crystal Targets]
- T. Trickle, Z. Zhang, K. M. Zurek, K. Inzani and S. Griffin, [Multi-Channel Direct Detection of Light Dark Matter: Theoretical Framework]
- S. M. Griffin, K. Inzani, T. Trickle, Z. Zhang and K. M. Zurek, [Multichannel direct detection of light dark matter: Target comparison]

[EXCEED-DM: Extended Calculation of Electronic Excitations for Direct Detection of Dark Matter]: https://arxiv.org/abs/2210.14917
[Non-standard neutrino interactions in light mediator models at reactor experiments]: https://arxiv.org/abs/2209.13566 
[Constraints on Sub-GeV Dark Matter-Electron Scattering from the CDEX-10 Experiment]: https://arxiv.org/abs/2206.04128 
[Dark Matter Direct Detection in Materials with Spin-Orbit Coupling]: https://arxiv.org/abs/2202.11716
[Dark Matter Absorption via Electronic Excitations]: https://link.springer.com/article/10.1007/JHEP09(2021)123 
[Extended Calculation of Dark Matter-Electron Scattering in Crystal Targets]: https://journals.aps.org/prd/abstract/10.1103/PhysRevD.104.095015
[Multi-Channel Direct Detection of Light Dark Matter: Theoretical Framework]: https://doi.org/10.1007/JHEP03(2020)036
[Multichannel direct detection of light dark matter: Target comparison]: https://doi.org/10.1103/PhysRevD.101.055004

## Attribution

If you use `EXCEED-DM` in your work, please cite,

    @misc{exdmv1,
        doi = {10.5281/ZENODO.7250321},
        url = {https://zenodo.org/record/7250321},
        author = {Trickle,  Tanner and {Kinzani}},
        title = {tanner-trickle/EXCEED-DM: EXCEED-DMv1.0.0},
        publisher = {Zenodo},
        year = {2022},
        copyright = {Open Access}
    }

    @article{Trickle2022,
        author = "Trickle, Tanner",
        title = "{EXCEED-DM: Extended Calculation of Electronic Excitations for Direct Detection of Dark Matter}",
        eprint = "https://arxiv.org/abs/2210.14917",
        archivePrefix = "arXiv",
        primaryClass = "hep-ph",
        month = "10",
        year = "2022"
    }

    @article{Griffin:2021znd,
        author = "Griffin, Sin\'ead M. and Inzani, Katherine and Trickle, Tanner and Zhang, Zhengkang and Zurek, Kathryn M.",
        title = "{Extended Calculation of Dark Matter-Electron Scattering in Crystal Targets}",
        eprint = "2105.05253",
        archivePrefix = "arXiv",
        primaryClass = "hep-ph",
        month = "5",
        year = "2021"
    }
