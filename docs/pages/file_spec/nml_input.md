title: namelist Input File
author: Tanner Trickle
date: v0.2.5

\(\require{mediawiki-texvc}\)

**Filetype**: text

**Extension**: `.txt`

**Example**: `./examples/Si/inputs/scatter_vc_test_input.txt`

A collection of `Fortran` namelist entries which specify the run input parameters. Each namelist has a title and some variable entries which set variables within the program. For example,

    &material
        name = 'Si'
    /

sets the name of the target material to 'Si'.

# Namelists

## control

- `process` - str
    - Physics process to compute for. Options : 'scatter', 'absorption', 'dielectric'.
- `calc_mode` - str
    - Calculation mode. Options : 'vc', 'cc', 'cf', 'vf', 'LO'.
- `overwrite_output` - bool
    - Flag to overwrite the output file. If `.FALSE.` an output filename will be generated.
- `timer` - bool
    - Flag to time the calculation.
- `quiet` - bool
    - Flag to shut off the printed output.

## io

- `out_folder` - str
    - Output folder.
- `run_description` - str
    - Short description of the run, will be appended to 'EXDMout_' to create the output filename.
- `out_filename` - str
    - Specify the output filename directly. Overwrites the `out_folder`/`run_description` combo.
- `PW_data_filename` - str
    - Filename of data holding \( u_{i, \mathbf{k}, \mathbf{G}} \) (with an \( s \) index for spin dependent wave functions ) data and corresponding metadata.
- `core_elec_config_filename` - str
    - Filename of core electron configuration data.
- `dielectric_filename` - str
    - Filename where dielectric data is stored to or read from.
- `nml_input_filename` - str
    - Filename where namelist entries are.
- `sto_data_filename` - str
    - Filename of STO coefficients and parameters for RHF wave functions.

## material

- `name` - str 
    - Name of the target material.
- `band_gap` - real
    - Band gap of the target material.
    - Units : \( \text{eV} \)
- `pc_vol_A` - real
    - Volume of the primitive cell
    - Units : \( \AA^3 \)
- `rho_T_g_per_cm3` - real
    - Target density.
    - Units : \( \text{g} \; \text{cm}^{-3} \)

## experiment

- `E_threshold` - real
    - Experimental energy deposited threshold.
    - Units : \( \text{eV} \)
- `n_time` - int
    - Number of times of day data is taken.
- `theta_E` - real
    - Angle between the Earth's rotation axis and the DM wind.
    - Units : \( \text{rad} \)
- `m_T_kg` - real
    - Mass of the target.
    - Units : \( \text{kg} \)
- `exposure_yr` - real
    - Exposure time of the target.
    - Units : \( \text{yr} \)

##dm_model

- `particle_type` - str
    - Dark matter particle type. Options : 'scalar', 'ps', 'vector', 'fermion'.
- `n_mX` - int
    - Number of dark matter masses.
- `log_mX_min` - real
    - Minimum \( \text{log10}(m_\chi / \text{eV}) \).
- `log_mX_max` - real
    - Maximum \( \text{log10}(m_\chi / \text{eV}) \).
- `n_extra_mX` - int
    - Number of 'extra' masses to compute for. Specifing this allows the user to also specify the `mX_2` namelist and add masses that are not uniform logarithmically spaced to the computed mass list.
- `n_med_FF` - int
    - Number of mediator form factor parameters.
- `med_FF_min` - real
    - Minimum mediator form factor parameter, \( \beta \); \( \mathcal{F}_\text{med} = \left( \frac{ \alpha m_e}{q} \right)^\beta \)
- `med_FF_max` - real
    - Maximum mediator form factor parameter, \( \beta \); \( \mathcal{F}_\text{med} = \left( \frac{ \alpha m_e}{q} \right)^\beta \)
- `rhoX_GeV_per_cm3` - real
    - Dark matter density.
    - Units : \( \text{GeV} \; \text{cm}^{-3} \)
- `v0_km_per_sec` - real
    - \( v_0 \) Dark matter velocity distribution parameter.
    - Units : \( \text{km} \; \text{s}^{-1} \)
- `vE_km_per_sec` - real
    - Earth velocity, dark matter velocity distribution parameter.
    - Units : \( \text{km} \; \text{s}^{-1} \)
- `vEsc` - real
    - Escape velocity, dark matter velocity distribution parameter.
    - Units : \( \text{km} \; \text{s}^{-1} \)
- `tff_id` - int
    - Dim : [2]
    - Transition form factor ID.

## bins_scatter

- `n_q` - int
    - Number of bins in momentum deposition.
- `n_E` - int
    - Number of bins in energy deposition.
- `E_width` - real
    - Width of the bins in energy deposition. 
    - Units : \( \text{eV} \)
- `q_width` - real
    - Width of the bins in momentum deposition.
    - Units : \( \text{eV} \)

## numerics_s_vc

Numerics parameters specific to valence to conduction scattering rate calculations.

- `n_val_max` - int
    - Maximum number of valence bands included in the calculation.
- `n_cond_max` - int
    - Maximum number of conduction bands included in the calculation.

## numerics_s_cc

Numerics parameters specific to core to conduction scattering rate calculations.

- `n_FFT_grid` - int
    - Dim : [3]
    - FFT grid size specified on input.
- `n_cond_max` - int
    - Maximum number of conduction bands included in the calculation.
- `n_principal_min` - int
    - Minimum principal quantum number, \( n \), to include in the calculation.
- `n_principal_max` - int
    - Maximum principal quantum number, \( n \), to include in the calculation.

## numerics_s_vf

- `n_val_max` - int
    - Maximum number of valence bands included in the calculation.
- `n_kf_theta` - int
    - Number of \( \theta \) angles in the integration over \( \mathbf{k}_f \)
- `n_kf_phi` - int
    - Number of \( \phi \) angles in the integration over \( \mathbf{k}_f \)
- `n_omega` - int
    - Number of \( \omega \) parameters to compute \( \frac{dR}{d \omega} \) for.
- `Zeff_type` - str
    - Specify how the \( Z_\text{eff} \) parameters are computed. Options : 'one' ( \( Z_\text{eff} = 1 \) ), 'Eb' ( \( Z_\text{eff}\) found from binding energy.  )

## numerics_s_cf

- `Ef_min` - real
    - Minimum final electron energy.
    - Units : \( \text{eV} \)
- `Zeff_type` - str
    - Specify how the \( Z_\text{eff} \) parameters are computed. Options : 'one' ( \( Z_\text{eff} = 1 \) ), 'Eb' ( \( Z_\text{eff}\) found from binding energy.  )
- `ki_min` - real
    - Minimum initial electron momentum to integrate over.
    - Units : \( \text{eV} \)
- `ki_s` - real
    - Scale parameter of maximum initial electron momentum to integrate over, \( k_{i, \text{max}} = k_s Z \alpha m_e \)
- `n_kf_theta` - int
    - Number of \( \theta \) angles in the integration over \( \mathbf{k}_f \)
- `n_kf_phi` - int
    - Number of \( \phi \) angles in the integration over \( \mathbf{k}_f \)
- `n_ki` - int
    - Number of points in \( k_i \) to integrate over.
- `n_ki_theta` - int
    - Number of \( \theta \) angles in the integration over \( \mathbf{k}_i \)
- `n_ki_phi` - int
    - Number of \( \phi \) angles in the integration over \( \mathbf{k}_i \)
- `n_principal_min` - int
    - Minimum principal quantum number, \( n \), to include in the calculation.
- `n_principal_max` - int
    - Maximum principal quantum number, \( n \), to include in the calculation.

## widths

Width parameters, \( \delta = \text{min}( \delta_\text{max}, a + b \omega ) \).

- `n_a` - int
    - Number of \( a \) parameters.
- `n_b` - int
    - Number of \( b \) parameters.
- `n_m` - int
    - Number of \( \delta_\text{max} \) parameters.
- `a_min` - real
    - Minimum value of \( a \).
    - Units : eV
- `a_max` - real
    - Maximum value of \( a \).
    - Units : eV
- `log_b_min` - real
    - Minimum value of \( \text{log10}(b/\text{eV}) \).
- `log_b_max` - real
    - Maximum value of \( \text{log10}(b/\text{eV}) \).
- `m_min` - real
    - Minimum value of \( \delta_\text{max} \).
    - Units : eV
- `m_max` - real
    - Maximum value of \( \delta_\text{max} \).
    - Units : eV
- `sigma` - real
    - Off resonance parameter, only include points with \( |\omega - \Delta \omega| < \sigma \omega \).

## numerics_abs

Numerics parameters specific to DM absorption calculations.

- `n_val_max` - int
    - Maximum number of valence bands included in the calculation.
- `n_cond_max` - int
    - Maximum number of conduction bands included in the calculation.
- `n_v` - int
    - Number of \( v \) points in the integration over \( \mathbf{v} \)
- `n_v_theta` - int
    - Number of \( \theta_v \) points in the integration over \( \mathbf{v} \)
- `n_v_phi` - int
    - Number of \( \phi_v \) points in the integration over \( \mathbf{v} \)
- `save_tran_form` - bool
    - Flag to save the transition form factors.

## bins_dielectric

Binning for the dielectric calculation.

- `n_E` - int
    - Number of \( \omega \) bins for numeric screening.
- `n_q` - int
    - Number of \( q \) bins for numeric screening.
- `n_q_theta` - int
    - Number of \( \theta_q \) bins for numeric screening.
- `n_q_phi` - int
    - Number of \( \phi_q \) bins for numeric screening.
- `E_width` - real
    - Width of bins in \( \omega \).
    - Units : \( \text{eV} \)
- `q_width` - real
    - Width of bins in \( q \).
    - Units : \( \text{eV} \)

## numerics_dielectric

- `n_val_max` - int
    - Maximum number of valence bands included in the calculation.
- `n_cond_max` - int
    - Maximum number of conduction bands included in the calculation.
- `n_k_vec` - int
    - Dim : [3]
    - Number of each dimension of a Monkhorst-Pack sampled 1BZ.

## in_med_scr

- `include_screen` - bool
    - Flag to include screening effects.
- `type` - str
    - Screening type. Options : 'analytic', 'numeric'
- `alpha` - real
    - \( \alpha \) shape parameter for analytic screening. \( \alpha \) in Eq 6 from [https://journals.aps.org/prb/pdf/10.1103/PhysRevB.47.9892](https://journals.aps.org/prb/pdf/10.1103/PhysRevB.47.9892)
- `e0` - real
    - Static dielectric parameter for analytic screening. \( \epsilon(0) \) in Eq 6 from [https://journals.aps.org/prb/pdf/10.1103/PhysRevB.47.9892](https://journals.aps.org/prb/pdf/10.1103/PhysRevB.47.9892)
- `q_tf` - real
    - Thomas Fermi momentum for analytic screening. \( q_\text{TF} \) in Eq 6 from [https://journals.aps.org/prb/pdf/10.1103/PhysRevB.47.9892](https://journals.aps.org/prb/pdf/10.1103/PhysRevB.47.9892)
    - Units : \( \text{eV} \)
- `omega_p` - real
    - Plasma frequency for analytic screening.\( \omega_p \) in Eq 6 from [https://journals.aps.org/prb/pdf/10.1103/PhysRevB.47.9892](https://journals.aps.org/prb/pdf/10.1103/PhysRevB.47.9892)
    - Units : \( \text{eV} \)
