
### astroph_model

* Dark matter SHM velocity distribution parameter, $v_0$.<br /><ul><li><b>Units</b>: $\text{km}/\text{s}$</li></ul>

        v_0_km_per_sec = 2.200000E+02

* List of Earth velocity vectors, $\mathbf{v}_e$.<br /><ul><li><b>Units</b>: $\text{km}/\text{s}$</li><li><b>Dim</b>: [ : , 3]</li></ul>

        v_e_km_per_sec = 0.000000E+00 0.000000E+00 2.300000E+02

* Dark matter SHM velocity distribution parameter, $v_\mathrm{esc}$.<br /><ul><li><b>Units</b>: $\text{km}/\text{s}$</li></ul>

        v_esc_km_per_sec = 6.000000E+02

* Specify the velocity distribution to use in the calculation.

        vel_distribution_name = 'SHM'

### control

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

### dm_model

* Form factor ID to compute for

        FIF_id = 'SI'

* Dark matter masses, $m_\chi$<br /><ul><li><b>Units</b>: $\text{eV}$</li><li><b>Dim</b>: [ : ]</li></ul>

        mX = 0.000000E+00

* Add $N$ linearly spaced dark matter masses between $m_\text{min}$ and $m_\text{max}$: [$N$, $m_\text{min}$, $m_\text{max}$]<br /><ul><li><b>Units</b>: [-, $\text{eV}$, $\text{eV}$]</li><li><b>Dim</b>: [3]</li></ul>

        mX_linspace = 0.000000E+00 1.000000E+00 1.000000E+00

* Add $N$ logarithmically spaced dark matter masses between $m_\text{min}$and $m_\text{max}$: [$N$, $m_\text{min}$, $m_\text{max}$]<br /><ul><li><b>Units</b>: [-, $\text{eV}$, $\text{eV}$]</li><li><b>Dim</b>: [3]</li></ul>

        mX_logspace = 0.000000E+00 1.000000E+00 1.000000E+00

* Mediator form factor powers, $\beta$<br /><ul><li><b>Formula</b>: $\mathcal{F}_\text{med} = \left( \frac{\alpha m_e}{q} \right)^\beta$</li><li><b>Dim</b>: [ : ]</li></ul>

        med_FF = 0.000000E+00

* Type of incoming dark matter particle

        particle_type = 'fermion'

* Local dark matter density<br /><ul><li><b>Units</b>: $\text{GeV}/\text{cm}^3$</li></ul>

        rho_X_GeV_per_cm3 = 4.000000E-01

### elec_config_input

* File containing the electronic configuration

        filename = ''

### experiment

* Mass of the experimental target<br /><ul><li><b>Units</b>: $\text{kg}$</li></ul>

        M_kg = 1.000000E+00

* Exposure time of the experimental target<br /><ul><li><b>Units</b>: $\text{yr}$</li></ul>

        T_year = 1.000000E+00

### material

* Lattice vectors of the target material, $\mathbf{a}_i$<br /><ul><li><b>Units</b>: $\text{Å}$</li><li><b>Dim</b>: [3, 3]</li></ul>

        a_vecs_Ang = 1.000000E+00 0.000000E+00 0.000000E+00 0.000000E+00 1.000000E+00 0.000000E+00 0.000000E+00 0.000000E+00 1.000000E+00

* Band gap of the target material, $E_g$<br /><ul><li><b>Units</b>: $\text{eV}$</li></ul>

        band_gap = 0.000000E+00

* Materials Project ID of the target material

        materials_project_ID = ''

* Name of the target material

        name = ''

* Mass density of the target material, $\rho_T$<br /><ul><li><b>Units</b>: $\text{g}/\text{cm}^3$</li></ul>

        rho_T_g_per_cm3 = 1.000000E+00

### numerics_absorption_rate

* Defines broadening behavior for the imaginary part of the Greens function

        smear_type = 'lorentz'

* List of widths, $\delta \, [\text{eV}]$, to compute for, parameterized as [$a$, $b$, $c$]<br /><ul><li><b>Formula</b>: $\delta = \text{min}(a + b \omega, c)$</li><li><b>Units</b>: [$\text{eV}$, -, $\text{eV}$]</li><li><b>Dim</b>: [ : , 3]</li></ul>

        widths = 0.000000E+00 1.000000E-01 1.000000E+02

### numerics_binned_scatter_rate

* Width of bins in $\omega$ space<br /><ul><li><b>Units</b>: $\text{eV}$</li></ul>

        E_bin_width = 1.000000E+00

* Number of bins in $\omega$ space

        n_E_bins = 1

* Number of bins in $q$ space

        n_q_bins = 1

* Width of bins in $q$ space<br /><ul><li><b>Units</b>: $\text{keV}$</li></ul>

        q_bin_width = 1.000000E+00

### numerics_dielectric

* Width of bins in $\omega$ space<br /><ul><li><b>Units</b>: $\text{eV}$</li></ul>

        E_bin_width = 1.000000E+00

* Number of bins in $\omega$ space

        n_E_bins = 1

* Number of bins in $q$ space

        n_q_bins = 1

* Number of bins in $\phi_\mathbf{q}$

        n_q_phi = 1

* Number of bins in $\theta_\mathbf{q}$

        n_q_theta = 1

* Width of bins in $q$ space<br /><ul><li><b>Units</b>: $\text{keV}$</li></ul>

        q_bin_width = 1.000000E+00

* Defines broadening behavior for the imaginary part of the Greens function.

        smear_type = 'lorentz'

* List of widths, $\delta \, [\text{eV}]$, to compute for, parameterized as [$a$, $b$, $c$]<br /><ul><li><b>Formula</b>: $\delta = \text{min}(a + b \omega, c)$</li><li><b>Units</b>: [$\text{eV}$, -, $\text{eV}$]</li><li><b>Dim</b>: [ : , 3]</li></ul>

        widths = 0.000000E+00 1.000000E-01 1.000000E+02

### screening

* Width of bins in $\omega$ space<br /><ul><li><b>Units</b>: $\text{eV}$</li></ul>

        E_bin_width = 1.000000E+00

* Shape parameter for analytic screening<ul><li><b>Formula</b>: $\alpha$, Eq. (6), <https://journals.aps.org/prb/pdf/10.1103/PhysRevB.47.9892></li></ul>

        alpha = 1.000000E+00

* Location of the numerically computed dielectric to use as the screening factor.

        dielectric_filename = ''

* Static dielectric parameter for analytic screening<ul><li><b>Formula</b>: $\epsilon(0)$, Eq. (6), <https://journals.aps.org/prb/pdf/10.1103/PhysRevB.47.9892></li></ul>

        e0 = 1.000000E+00

* Plasma frequency parameter for analytic screening<ul><li><b>Formula</b>: $\omega_p$, Eq. (6), <https://journals.aps.org/prb/pdf/10.1103/PhysRevB.47.9892></li><li><b>Units</b>: $\text{eV}$</li></ul>

        omega_p = 1.000000E+00

* Width of bins in $q$ space<br /><ul><li><b>Units</b>: $\text{keV}$</li></ul>

        q_bin_width = 1.000000E+00

* Thomas-Fermi momentum parameter for analytic screening<ul><li><b>Formula</b>: $q_\text{TF}$, Eq. (6), <https://journals.aps.org/prb/pdf/10.1103/PhysRevB.47.9892></li><li><b>Units</b>: $\text{keV}$</li></ul>

        q_tf = 1.000000E+00

* Model to use for the screening factor. Default is a screening factor of $1$ (no screening).

        type = ''

* Width ID of the dielectric to use. Dielectric will be taken from 'dielectric/width_<width_id>'inside 'dielectric_filename'

        width_id = 1

