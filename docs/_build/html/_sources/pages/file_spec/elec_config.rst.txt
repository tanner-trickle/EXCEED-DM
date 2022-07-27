Electronic Configuration File
=============================

.. contents:: Table of Contents
   :depth: 3
   :local:
   :backlinks: none

Introduction
------------

The electronic configuration file is an :code:`HDF5` file containing all the information about the electronic states in the target material. 

.. warning::

    :code:`Fortran` is a **column-major** language, which affects reading/writing of multi-dimensional datasets with HDF5. For example, for an :math:`n \times m` matrix to be read in to :code:`Fortran`, it must be saved as an :math:`m \times n` matrix if the language creating the HDF5 file is **row-major** (e.g., :code:`python`).


Electronic State Approximations
-------------------------------


************
Bloch States
************

Plane Wave (PW) basis
*********************

.. include:: elec_state_approx/bloch/PW_basis.rst
