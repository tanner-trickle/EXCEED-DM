
## astroph_model

* Dark matter SHM velocity distribution parameter, \( v_0 \) <br />Units : km/s

        v_0_km_per_sec = 2.200000E+02

* List of Earth velocity vectors in the galactic frame, \( \mathbf{v}_e \). Every three elements defines a \( \mathbf{v}_e \) vector<br />Units : km/sDim : [:, 3]

        v_e_km_per_sec = 0.000000E+00 0.000000E+00 2.300000E+02

* Dark matter SHM velocity distribution parameter, \( v_\mathrm{esc} \) <br />Units : km/s

        v_esc_km_per_sec = 6.000000E+02

* Specify the velocity distribution to use in the calculation.

        vel_distribution_name = 'SHM'

## control

* Which calculation to perform

        calculation = ''

* Filename to store default input parameters in markdown format to

        default_input_markdown_filename = './default_inputs.md'

* Filename to store input parameters in markdown format to

        input_markdown_filename = './inputs.md'

* Folder to store the ouput data

        out_folder = './'

* Small description of calculation which will be appended to 'EXDM_out_' to set the output filename

        run_description = ''

* Toggle to save the default input parameters to a markdown file

        save_default_inputs_markdown = F

* Toggle to save the input parameters to a markdown file

        save_inputs_markdown = F

* Toggle output printing to the console

        verbose = T

## dm_model

* Form factor ID to compute for

        FIF_id = 'SI'

* Dark matter masses, $m_\chi$<br />Units : eV<br />Dim : [:]

        mX = 0.000000E+00

* Add linearly spaced dark matter masses. [ N, mX_min, mX_max] (eV)

        mX_linspace = 0.000000E+00 1.000000E+00 1.000000E+00

* Add logarithmically spaced dark matter masses. [ N, mX_min, mX_max] (eV)

        mX_logspace = 0.000000E+00 1.000000E+00 1.000000E+00

* Mediator form factor powers, \( \beta \)<br />Formula: \( \mathcal{F}_\text{med} = \left( \frac{\alpha m_e}{q} \right)^\beta \)<br />Dim : [:]

        med_FF = 0.000000E+00

* Type of incoming dark matter particle

        particle_type = 'fermion'

* Local dark matter density (GeV/cm^3)

        rho_X_GeV_per_cm3 = 4.000000E-01

## elec_config_input

* File containing the electronic configuration

        filename = ''

## experiment

* Mass of the experimental target <br />Units: kg

        M_kg = 1.000000E+00

* Exposure time of the experimental target <br />Units: yr

        T_year = 1.000000E+00

## material

* Lattice vectors of the target material, \( \mathbf{a}_i \). Each row is a different lattice vector<br />Units : \( \mathrm{\AA} \)<br />Dim : [3, 3]

        a_vecs_Ang = 1.000000E+00 0.000000E+00 0.000000E+00 0.000000E+00 1.000000E+00 0.000000E+00 0.000000E+00 0.000000E+00 1.000000E+00

* Band gap of the target material, \( E_g \)<br />Units: \( \mathrm{g}/\mathrm{cm}^3 \)

        band_gap = 0.000000E+00

* Materials Project ID of the target material

        materials_project_ID = ''

* Name of the target material

        name = ''

* Mass density of the target material, \( \rho_T \)<br />Units: \( \mathrm{g}/\mathrm{cm}^3 \)

        rho_T_g_per_cm3 = 1.000000E+00

## numerics_absorption_rate

* Defines broadening behavior for the imaginary part of the Greens function.

        smear_type = 'lorentz'

* List of widths to compute for<br />Dim : [:, 3]

        widths = 0.000000E+00 1.000000E-01 1.000000E+02

## numerics_binned_scatter_rate

* Width of bins in \( \omega \) space<br />Units : eV

        E_bin_width = 1.000000E+00

* Number of bins in \( \omega \) space

        n_E_bins = 1

* Number of bins in \( q \) space

        n_q_bins = 1

* Width of bins in \( q \) space<br />Units : keV

        q_bin_width = 1.000000E+00

## numerics_dielectric

* Width of bins in $\omega$<br />Units: eV

        E_bin_width = 1.000000E+00

* Number of bins in $\omega$

        n_E_bins = 1

* Number of bins in $q$

        n_q_bins = 1

* Number of bins in $\phi_\mathbf{q}$

        n_q_phi = 1

* Number of bins in $\theta_\mathbf{q}$

        n_q_theta = 1

* Width of bins in $q$<br />Units: keV

        q_bin_width = 1.000000E+00

* Defines broadening behavior for the imaginary part of the Greens function.

        smear_type = 'lorentz'

* List of widths to compute for<br />Dim. - [:, 3]

        widths = 0.000000E+00 1.000000E-01 1.000000E+02

## screening

* Width of the bins in $\omega$. Only used for the numeric screening.<br />Units : eV

        E_bin_width = 1.000000E+00

* Shape parameter for analytic screening, $\alpha$ in Eq. (6) of [https://journals.aps.org/prb/pdf/10.1103/PhysRevB.47.9892](https://journals.aps.org/prb/pdf/10.1103/PhysRevB.47.9892)<br />Units : None

        alpha = 1.000000E+00

* Location of the numerically computed dielectric to use as the screening factor.

        dielectric_filename = ''

* Static dielectric parameter for analytic screening, $\epsilon(0)$ in Eq. (6) of[https://journals.aps.org/prb/pdf/10.1103/PhysRevB.47.9892](https://journals.aps.org/prb/pdf/10.1103/PhysRevB.47.9892)<br

        e0 = 1.000000E+00

* Plasma frequency parameter for analytic screening, $\omega_p$ in Eq. (6) of [https://journals.aps.org/prb/pdf/10.1103/PhysRevB.47.9892](https://journals.aps.org/prb/pdf/10.1103/PhysRevB.47.9892)<br />

        omega_p = 1.000000E+00

* Width of the bins in $q$. Only used for the numeric screening.<br />Units : keV

        q_bin_width = 1.000000E+00

* Thomas-Fermi momentum parameter for analytic screening, $q_\text{TF}$ in Eq. (6) of[https://journals.aps.org/prb/pdf/10.1103/PhysRevB.47.9892](https://journals.aps.org/prb/pdf/10.1103/PhysRevB.47.9892

        q_tf = 1.000000E+00

* Model to use for the screening factor. Default is a screening factor of 1.

        type = ''

* Width ID number to use for the numeric dielectric.

        width_id = 1

