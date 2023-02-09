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

