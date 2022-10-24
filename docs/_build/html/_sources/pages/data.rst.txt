====
Data
====

Electronic Configuration Files
==============================

While the :code:`examples` folder contains many example electronic configuration files useful for testing and learning what :code:`EXCEED-DM` can do, they are not suitable for publication quality results. The datasets are simply too small; they do no adequately represent the collection of electronic states in the target. The publication quality electronic configuration files are too big to host on the Github repository, typically being :math:`\mathcal{O}(\text{GB})` in size. There are two ways around this: 

1) Create your own electronic configuration file to the specifications given in :doc:`File Specifications - Electronic Configuration File </pages/file_spec/elec_config>`. This is the best option if you are studying a novel target.
2) Use a publicly available electronic configuration file. The files for Si and Ge are hosted on `Zenodo <https://zenodo.org/>`_, with links given below. It is **highly recommeded** that you make your electronic configuration files public so that results can be easily reproduced.

.. note:: Remember to cite the associated Zenodo repository if you use these publically available electronic configuration files.


**
Si
**

- Link: (add Zenodo link when available)
    - File: :code:`Si/scatter/Si_scatter_elec_config.hdf5`
    - Description: 
        - Electronic states assumed to be spin-degenerate, i.e., one-component wave functions. Used for binned scattering rate and dielectric calculations.
        - Initial States: 
            - Modelled with a combination of **STO basis** states and **PW basis** states. **STO basis** is used for the low energy, "core" states, while **PW basis** is used for the valence states.
            - STO basis
                - 10 states (1 :math:`\mathbf{k}` point (:math:`\mathbf{k} = 0`), 10 bands). States are the electrons in the :math:`1s \rightarrow 2p` orbitals, for both Si in the unit cell. The sum over lattice vectors, :math:`\mathbf{r}` extends to cells :math:`\pm 1` away from the center, i.e., includes 27 cells in total. Each state is sampled on a :math:`128 \times 128 \times 128` uniform grid in the unit cell. STO basis coefficients, e.g., :math:`C_{j, l, n, \kappa}` are taken from the tabulated values `here <https://linkinghub.elsevier.com/retrieve/pii/S0092640X8371003X>`_. Energy of each state is taken from the Materials Project database, material ID mp-149.
                - Notes: The reason for the small number of :math:`\mathbf{k}` points is due to runtime considerations, one has to choose between a larger sampling grid, i.e., sample the high momentum contributions in the unit cell, or more :math:`\mathbf{k}` points. Since the main features of these states are at high momentum, this is prioritized.
            - **PW basis**
                - 4000 states (:math:`10 \times 10 \times 10 \, \mathbf{k}` points, 4 bands). Computed with DFT (VASP), see Ref. :cite:`Griffin:2021znd` for more details. Uniform sampling in the 1BZ. Each state was expanded to an :math:`E_\text{cut} = \text{keV}` and then all-electron reconstructed with :code:`pawpyseed` to :math:`E_\text{cut} = 2 \, \text{keV}`. Lowest energy state at :math:`-11.814 \, \text{eV}`.


        - Final States:
            - Modelled with a combination of **PW basis** and **single PW basis** states. **PW basis** is used for the lower energy "conduction" bands, **single PW basis** is used for higher energy states being approximated as "free".
            - **PW basis**
                - 60000 states (:math:`10 \times 10 \times 10 \, \mathbf{k}` points, 60 bands). Computed with DFT (VASP), see Ref. :cite:`Griffin:2021znd` for more details. Uniform sampling in the 1BZ. Each state was expanded to an :math:`E_\text{cut} = \text{keV}` and then all-electron reconstructed with :code:`pawpyseed` to :math:`E_\text{cut} = 2 \, \text{keV}`. Lowest energy state at :math:`1.11 \, \text{eV}`. All bands which have an :math:`E_{i \mathbf{k}} < 60 \, \text{eV}`, for any :math:`\mathbf{k}`, are included.
            - **single PW basis**
                - 40000 states (:math:`10 \times 10 \times 400` grid in :math:`(\theta, \phi, p)` space). Uniformly sampled on the sphere in :math:`(\theta, \phi)`, logarithmically sampled in :math:`E = p^2/2m_e` between :math:`E_\text{min} = 60 \, \text{eV}` and :math:`E_\text{max} = 400 \, \text{eV}`.


- Link: (add Zenodo link when available)
    - File: :code:`Si/abs/Si_abs_elec_config.hdf5`
    - Description: 
        - Electronic states assumed to be spin-degenerate, i.e., one-component wave functions. Used for absorption rate calculations.
        - Initial States: 
            - Modelled with a combination of **STO basis** states and **PW basis** states. **STO basis** is used for the low energy, "core" states, while **PW basis** is used for the valence states.
            - STO basis
                - 10000 states (:math:`10 \times 10 \times 10 \, \mathbf{k}` points, 10 bands). States are the electrons in the :math:`1s \rightarrow 2p` orbitals, for both Si in the unit cell. :math:`\mathbf{k}` grid is uniformly sampled over the 1BZ. The sum over lattice vectors, :math:`\mathbf{r}` extends to cells :math:`\pm 1` away from the center, i.e., includes 27 cells in total. Each state is sampled on a :math:`128 \times 128 \times 128` uniform grid in the unit cell. STO basis coefficients, e.g., :math:`C_{j, l, n, \kappa}` are taken from the tabulated values `here <https://linkinghub.elsevier.com/retrieve/pii/S0092640X8371003X>`_. Energy of each state is taken from the Materials Project database, material ID mp-149.
                - Notes: A larger number of :math:`\mathbf{k}` vectors can, and must be, used here is because transitions must be vertical. This limits the number of transitions, relative to a scattering rate calculation.
            - **PW basis**
                - 4000 states (:math:`10 \times 10 \times 10 \, \mathbf{k}` points, 4 bands). Computed with DFT (VASP), see Ref. :cite:`Griffin:2021znd` for more details. Uniform sampling in the 1BZ. Each state was expanded to an :math:`E_\text{cut} = \text{keV}` and then all-electron reconstructed with :code:`pawpyseed` to :math:`E_\text{cut} = 2 \, \text{keV}`. Lowest energy state at :math:`-11.814 \, \text{eV}`.


        - Final States:
            - Modelled with a combination of **PW basis** and **single PW basis** states. **PW basis** is used for the lower energy "conduction" bands, **single PW basis** is used for higher energy states being approximated as "free".
            - **PW basis**
                - 60000 states (:math:`10 \times 10 \times 10 \, \mathbf{k}` points, 60 bands). Computed with DFT (VASP), see Ref. :cite:`Griffin:2021znd` for more details. Uniform sampling in the 1BZ. Each state was expanded to an :math:`E_\text{cut} = \text{keV}` and then all-electron reconstructed with :code:`pawpyseed` to :math:`E_\text{cut} = 2 \, \text{keV}`. Lowest energy state at :math:`1.11 \, \text{eV}`. All bands which have an :math:`E_{i \mathbf{k}} < 60 \, \text{eV}`, for any :math:`\mathbf{k}`, are included.
            - **single PW basis**
                - 2152000 states (:math:`10 \times 10 \times 10 \, \mathbf{k}` grid). :math:`\mathbf{k}` points are sampled uniformly in the 1BZ. For each :math:`\mathbf{k}`, all :math:`\mathbf{G}` were included such that :math:`60 \, \text{eV} < |\mathbf{k} + \mathbf{G}|^2/2m_e < \text{keV}`. Different :math:`\mathbf{G}` correspond to different bands when the parabolic dispersion relation gets folded in to the 1BZ.

**
Ge
**

- Link: (add Zenodo link when available)
    - File: :code:`Ge/scatter/Ge_scatter_elec_config.hdf5`
    - Description: 
        - Electronic states assumed to be spin-degenerate, i.e., one-component wave functions. Used for binned scattering rate and dielectric calculations.
        - Initial States: 
            - Modelled with a combination of **STO basis** states and **PW basis** states. **STO basis** is used for the low energy, "core" states, while **PW basis** is used for the valence states.
            - STO basis
                - 28 states (1 :math:`\mathbf{k}` point (:math:`\mathbf{k} = 0`), 28 bands). States are the electrons in the :math:`1s \rightarrow 3d` orbitals, for both Ge in the unit cell. The sum over lattice vectors, :math:`\mathbf{r}` extends to cells :math:`\pm 1` away from the center, i.e., includes 27 cells in total. Each state is sampled on a :math:`128 \times 128 \times 128` uniform grid in the unit cell. STO basis coefficients, e.g., :math:`C_{j, l, n, \kappa}` are taken from the tabulated values `here <https://linkinghub.elsevier.com/retrieve/pii/S0092640X8371003X>`_. Energy of each state is taken from the Materials Project database, material ID mp-32.
                - Notes: The reason for the small number of :math:`\mathbf{k}` points is due to runtime considerations, one has to choose between a larger sampling grid, i.e., sample the high momentum contributions in the unit cell, or more :math:`\mathbf{k}` points. Since the main features of these states are at high momentum, this is prioritized.
            - **PW basis**
                - 4000 states (:math:`10 \times 10 \times 10 \, \mathbf{k}` points, 4 bands). Computed with DFT (VASP), see Ref. :cite:`Griffin:2021znd` for more details. Uniform sampling in the 1BZ. Each state was expanded to an :math:`E_\text{cut} = \text{keV}` and then all-electron reconstructed with :code:`pawpyseed` to :math:`E_\text{cut} = 2 \, \text{keV}`. Lowest energy state at :math:`-11.814 \, \text{eV}`.


        - Final States:
            - Modelled with a combination of **PW basis** and **single PW basis** states. **PW basis** is used for the lower energy "conduction" bands, **single PW basis** is used for higher energy states being approximated as "free".
            - **PW basis**
                - 82000 states (:math:`10 \times 10 \times 10 \, \mathbf{k}` points, 82 bands). Computed with DFT (VASP), see Ref. :cite:`Griffin:2021znd` for more details. Uniform sampling in the 1BZ. Each state was expanded to an :math:`E_\text{cut} = \text{keV}` and then all-electron reconstructed with :code:`pawpyseed` to :math:`E_\text{cut} = 2 \, \text{keV}`. Lowest energy state at :math:`0.67 \, \text{eV}`. All bands which have an :math:`E_{i \mathbf{k}} < 60 \, \text{eV}`, for any :math:`\mathbf{k}`, are included.
            - **single PW basis**
                - 40000 states (:math:`10 \times 10 \times 400` grid in :math:`(\theta, \phi, p)` space). Uniformly sampled on the sphere in :math:`(\theta, \phi)`, logarithmically sampled in :math:`E = p^2/2m_e` between :math:`E_\text{min} = 60 \, \text{eV}` and :math:`E_\text{max} = 400 \, \text{eV}`.

- Link: (add Zenodo link when available)
    - File: :code:`Ge/abs/Ge_abs_elec_config.hdf5`
    - Description: 
        - Electronic states assumed to be spin-degenerate, i.e., one-component wave functions. Used for absorption rate calculations.
        - Initial States: 
            - Modelled with a combination of **STO basis** states and **PW basis** states. **STO basis** is used for the low energy, "core" states, while **PW basis** is used for the valence states.
            - STO basis
                - 28000 states (:math:`10 \times 10 \times 10 \, \mathbf{k}` points, 28 bands). States are the electrons in the :math:`1s \rightarrow 3d` orbitals, for both Ge in the unit cell. :math:`\mathbf{k}` grid is uniformly sampled over the 1BZ. The sum over lattice vectors, :math:`\mathbf{r}` extends to cells :math:`\pm 1` away from the center, i.e., includes 27 cells in total. Each state is sampled on a :math:`128 \times 128 \times 128` uniform grid in the unit cell. STO basis coefficients, e.g., :math:`C_{j, l, n, \kappa}` are taken from the tabulated values `here <https://linkinghub.elsevier.com/retrieve/pii/S0092640X8371003X>`_. Energy of each state is taken from the Materials Project database, material ID mp-32.
                - Notes: A larger number of :math:`\mathbf{k}` vectors can, and must be, used here is because transitions must be vertical. This limits the number of transitions, relative to a scattering rate calculation.
            - **PW basis**
                - 4000 states (:math:`10 \times 10 \times 10 \, \mathbf{k}` points, 4 bands). Computed with DFT (VASP), see Ref. :cite:`Griffin:2021znd` for more details. Uniform sampling in the 1BZ. Each state was expanded to an :math:`E_\text{cut} = \text{keV}` and then all-electron reconstructed with :code:`pawpyseed` to :math:`E_\text{cut} = 2 \, \text{keV}`. Lowest energy state at :math:`-11.814 \, \text{eV}`.

        - Final States:
            - Modelled with a combination of **PW basis** and **single PW basis** states. **PW basis** is used for the lower energy "conduction" bands, **single PW basis** is used for higher energy states being approximated as "free".
            - **PW basis**
                - 82000 states (:math:`10 \times 10 \times 10 \, \mathbf{k}` points, 82 bands). Computed with DFT (VASP), see Ref. :cite:`Griffin:2021znd` for more details. Uniform sampling in the 1BZ. Each state was expanded to an :math:`E_\text{cut} = \text{keV}` and then all-electron reconstructed with :code:`pawpyseed` to :math:`E_\text{cut} = 2 \, \text{keV}`. Lowest energy state at :math:`1.11 \, \text{eV}`. All bands which have an :math:`E_{i \mathbf{k}} < 60 \, \text{eV}`, for any :math:`\mathbf{k}`, are included.
            - **single PW basis**
                - 2586000 states (:math:`10 \times 10 \times 10 \, \mathbf{k}` grid). :math:`\mathbf{k}` points are sampled uniformly in the 1BZ. For each :math:`\mathbf{k}`, all :math:`\mathbf{G}` were included such that :math:`60 \, \text{eV} < |\mathbf{k} + \mathbf{G}|^2/2m_e < \text{keV}`. Different :math:`\mathbf{G}` correspond to different bands when the parabolic dispersion relation gets folded in to the 1BZ.


:code:`EXCEED-DM` Results
=========================

Below are links to datasets used for previously published results. See the refereces for more details.

***************************************************************************************************
EXCEED-DMv1.0.0: Extended Calculation of Electronic Excitations for Direct Detection of Dark Matter 
***************************************************************************************************

- Ref: (add arXiv link when available)
    - Link: (add Zenodo link when available)
    - Description:
        Output of all the calculations performed in the :code:`v1.0.0` user manual. Specifically, for Si and Ge targets,
            - Numerically computed dielectric (to be used to screen the scattering rate calculation).
            - Binned scatter rate of DM fermion in kinetically mixed dark photon model with different screenings: no screening, an analytic model of screening, and screened with the aforementioned numerically computed dielectric.
            - Binned scatter rate of DM fermion in a model where the scattering potential depends on the electron velocity (light mediator).
            - Extended absorption rate calculation for scalar, pseudoscalar, and vector DM.
            - Annual modulation of binned scatter rate of DM fermion in kinetically mixed dark photon model.

.. bibliography:: ../bibliography.bib
