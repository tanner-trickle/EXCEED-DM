.. EXCEED-DM documentation master file, created by
   sphinx-quickstart on Thu May  5 10:51:11 2022.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Home
====

.. .. image:: ./media/exdm-prelim-logo-modified.png
   .. :align: center

.. toctree::
    :maxdepth: 2

    Getting Started<pages/getting_started.rst>
    Examples<pages/examples_link.md>
    File Specifications<pages/file_spec/index.rst>
    Data<pages/data.rst>
    Papers<pages/papers.rst>
.. Source Documentation<source/index.rst>

----

.. .. image:: ./media/exdm-prelim-logo-modified.png
   .. :align: center

.. image:: https://zenodo.org/badge/DOI/10.5281/zenodo.6097642.svg
   :target: https://doi.org/10.5281/zenodo.6097642

Welcome to the :code:`EXCEED-DM` documentation!

:code:`EXCEED-DM` is an OpenMPI Fortran program that computes Dark Matter (DM) induced electronic transition rates. 

`User Manual <https://arxiv.org/abs/2210.14917>`_

========
Overview
========

:code:`EXCEED-DM` provides a complete framework for computing DM-electron interaction rates. Given an electronic configuration, :code:`EXCEED-DM` computes the relevant electronic matrix elements, then particle physics specific rates from these matrix elements. This allows for separation between approximations regarding the electronic state configuration, and the specific calculation being performed. 

- For more detailed information about what :code:`EXCEED-DM` is computing, and to see it used for new physics results, see the `user manual <https://arxiv.org/abs/2210.14917>`_.

- For the large electronic configuration files for Si and Ge targets, see the `Zenodo repository <https://zenodo.org/record/7246141#.Y1cKIbaSnZc>`_.

========
Features
========

- **Variety of Electronic State Approximations**
   - :code:`EXCEED-DM` provides a variety of different bases to specify the initial, and final, electronic states to accurately characterize the target. For example, plane wave (PW) (spin-dependent) bases and Slate type orbital (STO) bases.

- **Calculations**
    - **Scattering**: Given a DM model (DM masses, interaction potential, etc.) :code:`EXCEED-DM` computes the expected number of interactions per kg-year binned in energy and momentum deposition.
        - With the electronic state approximations separated from the scattering rate calculation, new transition form factors can be easily added as functions of the electronic matrix elements.
    - **Absorption**: Given a DM model (e.g., scalar, pseudoscalar, vector DM) :code:`EXCEED-DM` computes the expected number of interactions per kg-year.
    - **Dielectric**: For some processes the dielectric will screen the interaction rate. The complex dielectric can be computed and subsequently used in scattering rate calculations.

    - **Check out the examples folder to see specific calculations in action!**

:code:`EXCEED-DM` is:

- **fast** - Being parallelized with OpenMPI means that :code:`EXCEED-DM` can take full advantage of large computing clusters. 
- **DFT calculator independent** - :code:`EXCEED-DM` depends only on the output of DFT calculations, which means that any DFT calculator can be used to compute the targets electronic properties and then used as input. 

===========
Attribution
===========

If you use :code:`EXCEED-DM` in your work, please cite,

.. code-block:: none

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
