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

### Experiment


- DAMIC-M Collaboration (2025) [Probing Benchmark Models of Hidden-Sector Dark Matter with DAMIC-M]
- SENSEI Collaboration [SENSEI: First Direct-Detection Results on sub-GeV Dark Matter from SENSEI at SNOLAB]
- DAMIC-M Collaboration (2023) [First Constraints from DAMIC-M on Sub-GeV Dark-Matter Particles Interacting with Electrons]
- CDEX Collaboration [Constraints on Sub-GeV Dark Matter-Electron Scattering from the CDEX-10 Experiment]

### Theory

- R. Essig, Y. Hochberg, Y. Shoji, A. Singal, G. Suczewski [Low-Energy Compton Scattering in Materials]
- C. E. Dreyer, R. Essig, M. Fernandez-Serra, A. Singal, C. Zhen [Fully ab-initio all-electron calculation of dark matter--electron scattering in crystals with evaluation of systematic uncertainties] 
- G. Krnjaic, T. Trickle [Absorption of Vector Dark Matter Beyond Kinetic Mixing]
- T. Trickle [EXCEED-DM: Extended Calculation of Electronic Excitations for Direct Detection of Dark Matter]
- B. Dutta, S. Ghosh, T. Li, A. Thompson, A. Verma [Non-standard neutrino interactions in light mediator models at reactor experiments]
- H-Y. Chen, A. Mitridate, T. Trickle, Z. Zhang, M. Bernardi, K. M. Zurek, [Dark Matter Direct Detection in Materials with Spin-Orbit Coupling]
- A. Mitridate, T. Trickle, Z. Zhang, K. M. Zurek [Dark Matter Absorption via Electronic Excitations]
- S. M. Griffin, K. Inzani, T. Trickle, Z. Zhang and K. M. Zurek [Extended Calculation of Dark Matter-Electron Scattering in Crystal Targets]
- T. Trickle, Z. Zhang, K. M. Zurek, K. Inzani and S. Griffin [Multi-Channel Direct Detection of Light Dark Matter: Theoretical Framework]
- S. M. Griffin, K. Inzani, T. Trickle, Z. Zhang and K. M. Zurek [Multichannel direct detection of light dark matter: Target comparison]

[Probing Benchmark Models of Hidden-Sector Dark Matter with DAMIC-M]: https://arxiv.org/abs/2503.14617
[SENSEI: First Direct-Detection Results on sub-GeV Dark Matter from SENSEI at SNOLAB]: https://journals.aps.org/prl/abstract/10.1103/PhysRevLett.134.011804
[Fully ab-initio all-electron calculation of dark matter--electron scattering in crystals with evaluation of systematic uncertainties]: https://journals.aps.org/prd/abstract/10.1103/PhysRevD.109.115008
[Low-Energy Compton Scattering in Materials]: https://journals.aps.org/prd/abstract/10.1103/PhysRevD.109.116011
[First Constraints from DAMIC-M on Sub-GeV Dark-Matter Particles Interacting with Electrons]: https://journals.aps.org/prl/abstract/10.1103/PhysRevLett.130.171003
[Absorption of Vector Dark Matter Beyond Kinetic Mixing]: https://journals.aps.org/prd/abstract/10.1103/PhysRevD.108.015024
[EXCEED-DM: Extended Calculation of Electronic Excitations for Direct Detection of Dark Matter]: https://journals.aps.org/prd/abstract/10.1103/PhysRevD.107.035035
[Non-standard neutrino interactions in light mediator models at reactor experiments]: https://link.springer.com/article/10.1007/JHEP03(2023)163
[Constraints on Sub-GeV Dark Matter-Electron Scattering from the CDEX-10 Experiment]: https://journals.aps.org/prl/abstract/10.1103/PhysRevLett.129.221301
[Dark Matter Direct Detection in Materials with Spin-Orbit Coupling]: https://journals.aps.org/prd/abstract/10.1103/PhysRevD.106.015024
[Dark Matter Absorption via Electronic Excitations]: https://link.springer.com/article/10.1007/JHEP09(2021)123 
[Extended Calculation of Dark Matter-Electron Scattering in Crystal Targets]: https://journals.aps.org/prd/abstract/10.1103/PhysRevD.104.095015
[Multi-Channel Direct Detection of Light Dark Matter: Theoretical Framework]: https://doi.org/10.1007/JHEP03(2020)036
[Multichannel direct detection of light dark matter: Target comparison]: https://doi.org/10.1103/PhysRevD.101.055004

## Attribution

If you use `EXCEED-DM` in your work please cite,

    @article{Trickle:2022fwt,
        author = "Trickle, Tanner",
        title = "{Extended calculation of electronic excitations for direct detection of dark matter}",
        eprint = "2210.14917",
        archivePrefix = "arXiv",
        primaryClass = "hep-ph",
        reportNumber = "FERMILAB-PUB-22-767-T",
        doi = "10.1103/PhysRevD.107.035035",
        journal = "Phys. Rev. D",
        volume = "107",
        number = "3",
        pages = "035035",
        year = "2023"
    }

    @article{Griffin:2021znd,
        author = "Griffin, Sin\'ead M. and Inzani, Katherine and Trickle, Tanner and Zhang, Zhengkang and Zurek, Kathryn M.",
        title = "{Extended calculation of dark matter-electron scattering in crystal targets}",
        eprint = "2105.05253",
        archivePrefix = "arXiv",
        primaryClass = "hep-ph",
        doi = "10.1103/PhysRevD.104.095015",
        journal = "Phys. Rev. D",
        volume = "104",
        number = "9",
        pages = "095015",
        year = "2021"
    }
