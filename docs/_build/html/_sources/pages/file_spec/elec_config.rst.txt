Electronic Configuration File
=============================

Introduction
------------

The electronic configuration file is an :code:`HDF5` file containing all the information about the electronic states in the target material. Available electronic states are shown below. These data must be kept inside the HDF5 folder specified, e.g., to add initial electronic states in the plane wave basis, put the data specified below in the, :code:`init/bloch/PW_basis/` folder inside the HDF5 file. 

.. warning::

    :code:`Fortran` is a **column-major** language, which affects reading/writing of multi-dimensional datasets with HDF5. For example, for an :math:`n \times m` matrix to be read in to :code:`Fortran`, it must be saved as an :math:`m \times n` matrix if the language creating the HDF5 file is **row-major** (e.g., :code:`python`).


Electronic State Approximations
-------------------------------

************
Bloch States
************

:code:`bloch/PW_basis`: Plane Wave (PW) basis
*********************************************

**Equations**

.. math::

   u_{i, \mathbf{k}, s}(\mathbf{x}) = \sum_{\mathbf{G}} \, e^{i \mathbf{G} \cdot \mathbf{x}} \widetilde{u}_{i, \mathbf{k}, \mathbf{G}, s}

.. include:: elec_state_approx/bloch/PW_basis.rst

:code:`bloch/STO_basis`: Slater Type Orbital (STO) basis
********************************************************

**Equations**

.. math::

   u_{\kappa, n, \ell, m, \mathbf{k}}(\mathbf{x}) & = \sqrt{\Omega} \sum_{\mathbf{r}} e^{- i \mathbf{k} \cdot         \mathbf{y}_{\mathbf{r}, \kappa}} \sum_j C_{j, \ell, n, \kappa} R_\text{STO}(y_{\mathbf{r}; \kappa}, Z_{j, \ell,        \kappa}, n_{j, \ell, \kappa}) Y_l^m(\hat{\mathbf{y}}_{\mathbf{r}, \kappa})\, , \\
   R_\text{STO}(r; Z, n) & = a_0^{-3/2} \frac{(2Z)^{n + \frac{1}{2}}}{\sqrt{(2n)!}} \left( \frac{r}{a_0}              \right)^{n - 1} e^{-Zr / a_0}

.. include:: elec_state_approx/bloch/STO_basis.rst

:code:`bloch/single_PW`: Single PW
**********************************

**Equations**

.. math::

   u_{\mathbf{p}}(\mathbf{x}) = e^{i \mathbf{G}(\mathbf{p}) \cdot \mathbf{x}} 

.. include:: elec_state_approx/bloch/single_PW.rst
