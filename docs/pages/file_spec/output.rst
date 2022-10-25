Output File
===========

Introduction
------------

The output file is and :code:`HDF5` file and contains the output of a calculation. Overlap with variables in the input file is minimized, i.e., if a variable is present in the input file it may not be written to the output file.

Below we discuss all of the variables in each output group hierarchically (groups are in **bold**, datasets are in :code:`code`). Similar to the input file, if the **Dim** of a variable is **[ : ]**,  this means the array has variable dimension, and may change from calculation to calculation.

.. note:: 

   To further understand the output file, go through this document with the output file open in :code:`HDFView`, a graphical user interface to :code:`HDF5` files: `download <https://www.hdfgroup.org/downloads/hdfview/>`_

.. warning::

    :code:`Fortran` is a **column-major** language, which affects reading/writing of multi-dimensional datasets with HDF5. For example, an :math:`n \times m` matrix, written to an HDF5 file with :code:`Fortran`, will be read in to a **row-major** language (e.g., :code:`python`) as an :math:`m \times n` matrix. Dimensions listed below are the variables dimensions when written with :code:`EXCEED-DM`

General Groups
--------------

* :code:`exdm_version` - Version number of :code:`EXCEED-DM`.

* **dm_model**

  * :code:`FIF_id` - String ID of the scattering form factor used.

  * :code:`mX` - Dark matter masses, :math:`m_\chi`.

    * **Units**: :math:`\text{eV}`

    * **Dim**: [ : ]

  * :code:`med_FF` - Mediator form factor powers, :math:`\beta`.

    * **Formula**: :math:`\mathcal{F}_\text{med} = \left( \frac{q_0}{q} \right)^\beta` 

    * **Dim**: [ : ]

  * :code:`particle_type` - Type of DM particle

  * :code:`rho_X` - DM density, :math:`\rho_\chi`

    * **Units**: :math:`\text{GeV}/\text{cm}^3`

* **experiment**

  * :code:`M` - Mass of the target.

    * **Units**: :math:`\text{kg}`

  * :code:`T` - Exposure time of the target.

    * **Units**: :math:`\text{year}`

* **material**

  * :code:`band_gap` - Band gap of the target material. 

    * **Units**: :math:`\text{eV}`

  * :code:`name` - Name of the target material.

  * :code:`pc_vol` - Primitive cell/unit cell volume.

    * **Units**: :math:`\text{Å}^3`

  * :code:`rho_T` - Target density.

    * **Units**: :math:`\text{g}/\text{cm}^3`

* **timing**

  * :code:`dt_compute` - Time taken to perform sums in the calculation.  

    * **Units**: :math:`\text{s}`

  * :code:`dt_total` - Total run time of the calculation.  

    * **Units**: :math:`\text{s}`

  * :code:`start_date` - Start year/month/day of the calculation.

    * **Dim**: [3]

  * :code:`end_date` - End year/month/day of the calculation.

    * **Dim**: [3]


Calculation Specific Groups
---------------------------

Binned Scatter Rate
*******************

* **binned_scatter_rate**

  * **model_<n>/v_e_<v>/mass_<m>** - Folder structure for a **binned_scatter_rate** calculation. Note that if only a single DM model, or Earth velocity vector, is specified these folders will be absent. For example, if only a light mediator with one Earth velocity vector is chosen the folder structure will be **binned_scatter_rate/mass_<m>** for each mass. However if a heavy and light mediator are calculated the folder structure will be **binned_scatter_rate/model_<n>/mass_<m>**. Similarly for mulitple Earth velocity vectors.

    * :code:`total_binned_scatter_rate` - Total number of events in each :math:`q, \omega` bin, assuming :math:`\overline{\sigma}_e = 1` and an exposure of :code:`M_kg x T_year (kg-yr)`, where :code:`M_kg, T_year` are specified in the **input file**.

      .. note:: 

        Note the absence of units on :math:`\overline{\sigma}_e` above. To compute the total number of events simply multiply the data by :math:`\overline{\sigma}_e` in units of :math:`\text{cm}^2`. For example to compute the number of events with :math:`\overline{\sigma}_e = 10^{-40} \, \text{cm}^2`, mulitply the data by :math:`10^{-40}`. 

      The event rate, in units of events/kg-year, is then calculated by dividing the entries by the exposure in kg-year. For example, if **M_kg = 1, T_year = 1** (the default) in the **experiment** group in the **input file**, the entries are both, a) the total number of events in a kg-year and b) the number of events/kg-year. If **M_kg = 10, T_year = 1**, the entries are the total number of events in 10 kg-year, and the event rate is easily calculated (simply divide by 10 in this example). 

      * **Dim**: [ :, : ], [:math:`N_q`, :math:`N_\omega`]

      * **Units**: :math:`\text{cm}^{-2}`

    * :code:`i_<i>/binned_scatter_rate` - Total number of events in each :math:`q, \omega` bin from each initial state group labelled by :code:`i`, e.g., the band number for Bloch states. The sum of the entries in each :code:`i` sum to the :code:`total_binned_scatter_rate`.

      * **Dim**: [ :, : ], [:math:`N_q`, :math:`N_\omega`]

      * **Units**: :math:`\text{cm}^{-2}`

* **numerics_binned_scatter_rate**

  * :code:`E_bin_width` - Width of the bins in :math:`\omega` space.

    * **Units**: :math:`\text{eV}`

  * :code:`q_bin_width` - Width of the bins in :math:`q` space.

    * **Units**: :math:`\text{keV}`

Absorption Rate
***************

* **absorption_rate**

  * **width_<i>/mass_<m>** - Folder structure for an **absorption_rate** calculation. Note that if only a single width is specified these folders will be absent. For example, if only a single width is chosen the folder structure will be **absorption_rate/mass_<m>** for each mass. However if multiple widths are calculated the folder structure will be **binned_scatter_rate/width_<i>/mass_<m>**.

    * :code:`absorption_rate` - Total number of events assuming :math:`g_e = 1`, and an exposure of :code:`M_kg x T_year (kg-yr)`. The event rate, in units of events/kg-year, is then calculated by dividing the entries by the exposure in kg-year. For example, if **M_kg = 1, T_year = 1** (the default) in the **experiment** group in the **input file**, the entries are both, a) the total number of events in a kg-year and b) the number of events/kg-year. If **M_kg = 10, T_year = 1**, the entries are the total number of events in 10 kg-year, and the event rate is easily calculated (simply divide by 10 in this example). 

* **numerics_absorption_rate**

  * :code:`widths` - List of width parameterizations computed for, [:math:`a`, :math:`b`, :math:`c`].

    * **Formula**: :math:`\delta = \text{min}(c, a + b \omega)`

    * **Units**: [ : , { :math:`\text{eV}`, - , :math:`\text{eV}` } ]

    * **Dim**: [ : , 3]

  * :code:`smear_type` - Defines broadening behavior for the imaginary part of the Greens function.

Dielectric
**********

* **dielectric** 

  * **width_<i>** - Folder containing the dielectric for width parameterization :code:`i`

    * :code:`dielectric_r` - Real part of the averaged dielectric, :math:`\overline{\varepsilon}(\mathbf{q}, \omega)`.

      * **Dim**: [ : , : , : , : ], [:math:`N_q`, :math:`N_\theta`, :math:`N_\phi`, :math:`N_\omega`]

    * :code:`dielectric_c` - Imaginary part of the averaged dielectric, :math:`\overline{\varepsilon}(\mathbf{q}, \omega)`.

      * **Dim**: [ : , : , : , : ], [:math:`N_q`, :math:`N_\theta`, :math:`N_\phi`, :math:`N_\omega`]

* **numerics_dielectric**

  * :code:`E_bin_width` - Width of the bins in :math:`\omega` space.

    * **Units**: :math:`\text{eV}`

  * :code:`q_bin_width` - Width of the bins in :math:`q` space.

    * **Units**: :math:`\text{keV}`

  * :code:`widths` - List of width parameterizations computed for, [:math:`a`, :math:`b`, :math:`c`].

    * **Formula**: :math:`\delta = \text{min}(c, a + b \omega)`

    * **Units**: [ : , { :math:`\text{eV}`, - , :math:`\text{eV}` } ]

    * **Dim**: [ : , 3]

  * :code:`smear_type` - Defines broadening behavior for the imaginary part of the Greens function.
