
### astroph_model

* `v_0_km_per_sec`: Dark matter SHM velocity distribution parameter, $v_0$.<br /><ul><li><b>Units</b>: $\text{km}/\text{s}$</li></ul>

      v_0_km_per_sec = 2.300000E+02

* `v_e_km_per_sec`: List of Earth velocity vectors, $\mathbf{v}_e$.<br /><ul><li><b>Units</b>: $\text{km}/\text{s}$</li><li><b>Dim</b>: [ : , 3]</li></ul>

      v_e_km_per_sec = 0.000000E+00 0.000000E+00 2.400000E+02

* `v_esc_km_per_sec`: Dark matter SHM velocity distribution parameter, $v_\mathrm{esc}$.<br /><ul><li><b>Units</b>: $\text{km}/\text{s}$</li></ul>

      v_esc_km_per_sec = 6.000000E+02

* `vel_distribution_name`: Specify the velocity distribution to use in the calculation.

      vel_distribution_name = 'SHM'

### control

* `calculation`: Which calculation to perform

      calculation = ''

* `default_input_markdown_filename`: Filename to store default input parameters in markdown format to

      default_input_markdown_filename = './default_inputs.md'

* `input_markdown_filename`: Filename to store input parameters in markdown format to

      input_markdown_filename = './inputs.md'

* `out_folder`: Folder to store the ouput data

      out_folder = './'

* `run_description`: Small description of calculation which will be appended to 'EXDM_out_' to set the output filename

      run_description = ''

* `save_default_inputs_markdown`: Toggle to save the default input parameters to a markdown file

      save_default_inputs_markdown = F

* `save_inputs_markdown`: Toggle to save the input parameters to a markdown file

      save_inputs_markdown = F

* `verbose`: Toggle output printing to the console

      verbose = T

### dm_model

* `FIF_id`: Form factor ID to compute for

      FIF_id = 'SI'

* `mX`: Dark matter masses, $m_\chi$<br /><ul><li><b>Units</b>: $\text{eV}$</li><li><b>Dim</b>: [ : ]</li></ul>

      mX = 0.000000E+00

* `mX_linspace`: Add $N$ linearly spaced dark matter masses between $m_\text{min}$ and $m_\text{max}$: [$N$, $m_\text{min}$, $m_\text{max}$]<br /><ul><li><b>Units</b>: [-, $\text{eV}$, $\text{eV}$]</li><li><b>Dim</b>: [3]</li></ul>

      mX_linspace = 0.000000E+00 1.000000E+00 1.000000E+00

* `mX_logspace`: Add $N$ logarithmically spaced dark matter masses between $m_\text{min}$and $m_\text{max}$: [$N$, $m_\text{min}$, $m_\text{max}$]<br /><ul><li><b>Units</b>: [-, $\text{eV}$, $\text{eV}$]</li><li><b>Dim</b>: [3]</li></ul>

      mX_logspace = 0.000000E+00 1.000000E+00 1.000000E+00

* `med_FF`: Mediator form factor powers, $\beta$<br /><ul><li><b>Formula</b>: $\mathcal{F}_\text{med} = \left( \frac{\alpha m_e}{q} \right)^\beta$</li><li><b>Dim</b>: [ : ]</li></ul>

      med_FF = 0.000000E+00

* `particle_type`: Type of incoming dark matter particle

      particle_type = 'fermion'

* `rho_X_GeV_per_cm3`: Local dark matter density<br /><ul><li><b>Units</b>: $\text{GeV}/\text{cm}^3$</li></ul>

      rho_X_GeV_per_cm3 = 4.000000E-01

### elec_config_input

* `filename`: File containing the electronic configuration

      filename = ''

### experiment

* `M_kg`: Mass of the experimental target<br /><ul><li><b>Units</b>: $\text{kg}$</li></ul>

      M_kg = 1.000000E+00

* `T_year`: Exposure time of the experimental target<br /><ul><li><b>Units</b>: $\text{yr}$</li></ul>

      T_year = 1.000000E+00

### material

* `a_vecs_Ang`: Lattice vectors of the target material, $\mathbf{a}_i$<br /><ul><li><b>Units</b>: $\text{Å}$</li><li><b>Dim</b>: [3, 3]</li></ul>

      a_vecs_Ang = 1.000000E+00 0.000000E+00 0.000000E+00 0.000000E+00 1.000000E+00 0.000000E+00 0.000000E+00 0.000000E+00 1.000000E+00

* `band_gap`: Band gap of the target material, $E_g$<br /><ul><li><b>Units</b>: $\text{eV}$</li></ul>

      band_gap = 0.000000E+00

* `materials_project_ID`: Materials Project ID of the target material

      materials_project_ID = ''

* `n_T_g_per_cm3_per_AMU`: Number density of the target material, $n_T$<br /><ul><li><b>Units</b>: $\text{g}/\text{cm}^3/\text{AMU}$</li></ul>

      n_T_g_per_cm3_per_AMU = 1.000000E+00

* `name`: Name of the target material

      name = ''

* `rho_T_g_per_cm3`: Mass density of the target material, $\rho_T$<br /><ul><li><b>Units</b>: $\text{g}/\text{cm}^3$</li></ul>

      rho_T_g_per_cm3 = 1.000000E+00

### numerics_TIF_calculator_atomic

* `integration_scheme`: Specific method of sampling the radial direction.

      integration_scheme = 'log'

* `n_r`: Number of radial points to integrate with.

      n_r = 1

* `r_max_a0`: Maximum radius to use in the integration, in units of the Bohr radius, $a_0$.<ul><li><b>[$a_0$]</b></li></ul>

      r_max_a0 = 1.000000E+03

* `r_min_a0`: Minimum radius to use in the integration, in units of the Bohr radius, $a_0$.<ul><li><b>[$a_0$]</b></li></ul>

      r_min_a0 = 1.000000E-03

### numerics_absorption_rate

* `smear_type`: Defines broadening behavior for the imaginary part of the Greens function

      smear_type = 'lorentz'

* `widths`: List of widths, $\delta \, [\text{eV}]$, to compute for, parameterized as [$a$, $b$, $c$]<br /><ul><li><b>Formula</b>: $\delta = \text{min}(a + b \omega, c)$</li><li><b>Units</b>: [$\text{eV}$, -, $\text{eV}$]</li><li><b>Dim</b>: [ : , 3]</li></ul>

      widths = 0.000000E+00 1.000000E-01 1.000000E+02

### numerics_binned_scatter_rate

* `E_bin_width`: Width of bins in $\omega$ space<br /><ul><li><b>Units</b>: $\text{eV}$</li></ul>

      E_bin_width = 1.000000E+00

* `n_E_bins`: Number of bins in $\omega$ space

      n_E_bins = 1

* `n_q_bins`: Number of bins in $q$ space

      n_q_bins = 1

* `q_bin_width`: Width of bins in $q$ space<br /><ul><li><b>Units</b>: $\text{keV}$</li></ul>

      q_bin_width = 1.000000E+00

### numerics_dielectric

* `E_bin_width`: Width of bins in $\omega$ space<br /><ul><li><b>Units</b>: $\text{eV}$</li></ul>

      E_bin_width = 1.000000E+00

* `eta`: Parameter controlling the small $q$ approximation for the dielectric.

      eta = 1.000000E-01

* `n_E_bins`: Number of bins in $\omega$ space

      n_E_bins = 1

* `n_q_bins`: Number of bins in $q$ space

      n_q_bins = 1

* `n_q_phi`: Number of bins in $\phi_\mathbf{q}$

      n_q_phi = 1

* `n_q_theta`: Number of bins in $\theta_\mathbf{q}$

      n_q_theta = 1

* `q_bin_width`: Width of bins in $q$ space<br /><ul><li><b>Units</b>: $\text{keV}$</li></ul>

      q_bin_width = 1.000000E+00

* `smear_type`: Defines broadening behavior for the imaginary part of the Greens function.

      smear_type = 'lorentz'

* `widths`: List of widths, $\delta \, [\text{eV}]$, to compute for, parameterized as [$a$, $b$, $c$]<br /><ul><li><b>Formula</b>: $\delta = \text{min}(a + b \omega, c)$</li><li><b>Units</b>: [$\text{eV}$, -, $\text{eV}$]</li><li><b>Dim</b>: [ : , 3]</li></ul>

      widths = 0.000000E+00 1.000000E-01 1.000000E+02

### screening

* `E_bin_width`: Width of bins in $\omega$ space<br /><ul><li><b>Units</b>: $\text{eV}$</li></ul>

      E_bin_width = 1.000000E+00

* `alpha`: Shape parameter for analytic screening<ul><li><b>Formula</b>: $\alpha$, Eq. (6), <https://journals.aps.org/prb/pdf/10.1103/PhysRevB.47.9892></li></ul>

      alpha = 1.000000E+00

* `dielectric_filename`: Location of the numerically computed dielectric to use as the screening factor.

      dielectric_filename = ''

* `e0`: Static dielectric parameter for analytic screening<ul><li><b>Formula</b>: $\epsilon(0)$, Eq. (6), <https://journals.aps.org/prb/pdf/10.1103/PhysRevB.47.9892></li></ul>

      e0 = 1.000000E+00

* `omega_p`: Plasma frequency parameter for analytic screening<ul><li><b>Formula</b>: $\omega_p$, Eq. (6), <https://journals.aps.org/prb/pdf/10.1103/PhysRevB.47.9892></li><li><b>Units</b>: $\text{eV}$</li></ul>

      omega_p = 1.000000E+00

* `q_bin_width`: Width of bins in $q$ space<br /><ul><li><b>Units</b>: $\text{keV}$</li></ul>

      q_bin_width = 1.000000E+00

* `q_tf`: Thomas-Fermi momentum parameter for analytic screening<ul><li><b>Formula</b>: $q_\text{TF}$, Eq. (6), <https://journals.aps.org/prb/pdf/10.1103/PhysRevB.47.9892></li><li><b>Units</b>: $\text{keV}$</li></ul>

      q_tf = 1.000000E+00

* `type`: Model to use for the screening factor. Default is a screening factor of $1$ (no screening).

      type = ''

* `width_id`: Width ID of the dielectric to use. Dielectric will be taken from 'dielectric/width_<width_id>'inside 'dielectric_filename'

      width_id = 1

