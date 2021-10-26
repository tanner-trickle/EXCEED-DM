title: EXCEED-DM Output File
author: Tanner Trickle
date: v0.2.5

\(\require{mediawiki-texvc}\)

**Filetype**: HDF5

**Extension**: `.h5`, `.hdf5`

All of the output data from an EXCEED-DM run.

# Groups

## control

- `calc_mode` - str
    - Calculation mode.
- `process` - str
    - Physical process to calculate.

## io

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
    - Target material name.
- `band_gap` - real - \( E_g \)
    - Target material band gap.
    - Units: \( \text{eV} \)
- `density` - real - \( \rho_T \)
    - Target material density.
    - Units : \( \text{g} \; \text{cm}^{-3} \)
- `pc_vol` - real - \( \Omega \)
    - Volume of the primitive cell.
    - Units : \( \AA^3 \)

## dm_model

- `particle_type` - str
    - Dark matter particle type, e.g. 'fermion'.
- `mX` - real - \( m_\chi \)
    - Dim : [ `n_mX` ]
    - Dark matter masses.
    - Units : \( \text{eV} \)
- `med_FF` - real - \( \beta \)
    - Dim : [`n_med_FF`]
    - Mediator form factor parameters, \( \beta \); \( \mathcal{F}_\text{med} = \left( \frac{ \alpha m_e}{q} \right)^\beta \)
- `n_mX` - int
    - Number of dark matter masses.
- `n_med_FF` - int
    - Number of mediator form factor parameters.
- `v0` - real - \( v_0 \)
    - Dark matter velocity distribution parameter.
    - Units : \( \text{km} \; \text{s}^{-1} \)
- `vE` - real - \( v_E \)
    - Earth velocity, dark matter velocity distribution parameter.
    - Units : \( \text{km} \; \text{s}^{-1} \)
- `vEsc` - real - \( v_\text{esc} \)
    - Escape velocity, dark matter velocity distribution parameter.
    - Units : \( \text{km} \; \text{s}^{-1} \)
- `tff_id` - int
    - Dim : [2]
    - Transition form factor ID.

## experiment

- `E_threshold` - real - \( E_\text{th} \)
    - Experimental energy deposited threshold.
    - Units : \( \text{eV} \)
- `n_time` - int
    - Number of times of day data is taken.
- `times` - real
    - Dim : [`n_time`]
    - Times of the day data is taken.
    - Units : \( \text{day} \)
- `vE_direction` - real - \( \hat{\mathbf{v}_e}(t) \)
    - Dim: [`n_time`, 3]
    - Direction of the Earth velocity in the galactic frame at all the `times` of the day.
- `m_T` - real
    - Mass of the target.
    - Units : \( \text{kg} \)
- `exposure` - real
    - Exposure time of the target.
    - Units : \( \text{yr} \)

## FFT

- `N` - integer - \( N_\text{FFT} \)
    - Number of points in the FFT.
- `n_grid` - integer
    - Dim : [3]
    - Number of points in the FFT in each dimension.
- `q_max` - real - \( q^\text{FFT}_\text{max} \)
    - Maximum \( q \) the FFT is consistent for.
    - Units : \( \text{eV} \)

## PW_dataset

- `n_k` - int - \( N_\mathbf{k} \)
    - Number of \( \mathbf{k} \) points.
- `n_G` - int \( N_\mathbf{G} \)
    - Number of \( \mathbf{G} \) points.
- `n_val` - int - \( N_\text{val} \)
    - Number of valence bands in the dataset.
- `n_cond` - int - \( N_\text{cond} \)
    - Number of conduction bands in the dataset.
- `n_bands` - int - \( N_\text{bands} \)
    - Total number of bands in the dataset.
- `q_cut` - real - \( q_\text{cut} \)
    - Plane wave expansion parameter, \( q_\text{cut} = \text{max}(|\mathbf{k} + \mathbf{G}|) \).
    - Units : \( \text{eV} \)
- `E_cut` - real - \( E_\text{cut} \)
    - Plane wave expansion parameter, \( E_\text{cut} = \frac{q_\text{cut}^2}{2 m_e} \)
    - Units : \( \text{eV} \)
- `Ef_max` - real - \( E_f^\text{max} \)
    - All conduction bands included in the calculation cross below \( E_f^\text{max} \)
    - Units : \( \text{eV} \)
- `k_red_to_xyz` - real
    - Converts between reciprocal reduced coordinates to physical coordinates. The i'th column is the i'th reciprcal lattice vector, \( \mathbf{b}_i \)
    - Units : \( \text{eV} \)
- `red_to_xyz` - real
    - Converts between reduced coordinates to physical coordinates. The i'th column is the i'th lattice vector, \( \mathbf{a}_i \)
    - Units : \( \text{eV}^{-1} \)
- `spin_degen` - real
    - Spin degeneracy of each states, \( = 2 (1) \) for spin independent (dependent) wave functions.  

## core_electron

- `n_atom` - int
    -  Number of atoms in the primitive cell.
- `n_state` - int
    - Number of core states.
- `eq_pos_red` - real
    - Dim : [`n_atom`, 3]
    - Equilibrium position of each atom in the primitive cell in reduced coordinates.
    - Units : None
- `Z` - int
    - Dim : [`n_atom`] 
    - Proton number of each atom in the primitive cell.
- `energy` - real
    - Dim : [`n_state`]
    - Energy of each of the core electron states.
- `config` - int
    - Dim : [`n_state`, 5]
    - Core electron configuration. A map from a state to the [atom id, n, l, m, n_s] information about each core state, where atom id is the index of the atom in the primitive cell, \( n, l, m \) are the standard quantum numbers, and `n_s` is the number of spin states
- `STO/n_j` - int
    - Dim : [`n_state`]
    - Number of STO functions used to model the core electronic state
- `STO/max_n_j` - int
    - The maximum value of the array at `STO/n_j`
- `STO/n_lj` - int
    - Dim : [`n_state`, `STO/max_n_j`]
    - \( n_{lj} \) parameters used to compute the RHF wave functions from the STO's.
- `STO/Z_lj` - real
    - Dim : [`n_state`, `STO/max_n_j`]
    - \( Z_{lj} \) parameters used to compute the RHF wave functions from the STO's.
- `STO/N0_lj` - real
    - Dim : [`n_state`, `STO/max_n_j`]
    - \( N_{lj} \) parameters used to compute the RHF wave functions from the STO's.
- `STO/C_nlj` - real
    - Dim : [`n_state`, `STO/max_n_j`]
    - \( C_{nlj} \) parameters used to compute the RHF wave functions from the STO's.

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

## numerics_scatter_vc

Numerics parameters specific to valence to conduction scattering rate calculations.

- `n_val_max` - int
    - Maximum number of valence bands included in the calculation.
- `n_cond_max` - int
    - Maximum number of conduction bands included in the calculation.

## numerics_scatter_cc

Numerics parameters specific to core to conduction scattering rate calculations.

- `core_id_list` - int
    - List of the core electron states that the rate was computed using. `core_electron/config(core_id_list(i))` is the configuration of the `i`th core state included in the calculation.
- `n_FFT_grid_goal` - int
    - Dim : [3]
    - FFT grid size specified on input.
- `n_cond_max` - int
    - Maximum number of conduction bands included in the calculation.
- `n_principal_min` - int
    - Minimum principal quantum number, \( n \), to include in the calculation.
- `n_principal_max` - int
    - Maximum principal quantum number, \( n \), to include in the calculation.

## numerics_scatter_vf

- `n_val_max` - int
    - Maximum number of valence bands included in the calculation.
- `n_kf_theta` - int
    - Number of \( \theta \) angles in the integration over \( \mathbf{k}_f \)
- `n_kf_phi` - int
    - Number of \( \phi \) angles in the integration over \( \mathbf{k}_f \)
- `n_omega` - int
    - Number of \( \omega \) parameters to compute \( \frac{dR}{d \omega} \) for.
- `omega_list` - real
    - Dim : [`n_omega']
    - List of \( \omega \) parameters \( \frac{dR}{d \omega} \) is computed for.
    - Units : \( \text{eV} \)
- `Zeff` - real
    - Dim: [`n_state`, `n_k`]
    - \( Z_\text{eff} \) parameters used in the Fermi form factor.
    - Units : None

## numerics_scatter_cf

- `Ef_min` - real
    - Minimum final electron energy.
    - Units : \( \text{eV} \)
- `Zeff` - real
    - Dim : [`n_state`, `n_k`]
    - \( Z_\text{eff} \) parameters used in the Fermi form factor.
    - Units : None
- `core_id_list` - int
    - List of the core electron states that the rate was computed using. `core_electron/config(core_id_list(i))` is the configuration of the `i`th core state included in the calculation.
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
- `omega_list` - real
    - Dim: [`n_omega']
    - List of \( \omega \) parameters \( \frac{dR}{d \omega} \) is computed for.
    - Units: eV
- `n_principal_min` - int
    - Minimum principal quantum number, \( n \), to include in the calculation.
- `n_principal_max` - int
    - Maximum principal quantum number, \( n \), to include in the calculation.

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

## widths

Width parameters, \( \delta = \text{min}( \delta_\text{max}, a + b \omega ) \).

- `n_a` - int
    - Number of \( a \) parameters.
- `n_b` - int
    - Number of \( b \) parameters.
- `n_m` - int
    - Number of \( \delta_\text{max} \) parameters.
- `n` - int
    - Total number of width parameters, \( n = n_a \times n_b \times n_m \).
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
- `info` - real
    - Dim : [`n`, 3]
    - List of the width parameters, each entry lists \( [a, b, \delta_\text{max}] \)

## scatter_rates

Binned and total scattering rates, per cross section, for each time, mediator form factor, and DM mass, indexed by, `t`, `f`, `m` respectively. Rates from a given initial state are inside the `/init_()` subgroup. () should be replaced with integers which fall into the appropriate range, e.g. if `dm_model/n_mX` = 5 then `m_5` is the last subgroup.

- `t_()/f_()/m_()/total` - real
    - Total scattering rate, summed over initial states.
    - Units : \( \text{cm}^{-2} \)
- `t_()/f_()/m_()/total_binned` - real
    - Dim : [`n_q`, `n_E`]
    - Total binned scattering rate, summed over initial states.
    - Units : \( \text{cm}^{-2} \)
- `t_()/f_()/m_()/init_()/binned_i` - real
    - Dim : [`n_q`, `n_E`]
    - Binned scattering rate, from the `i`th initial state.
    - Units : \( \text{cm}^{-2} \)

## abs_rates

Absorption rates, assuming DM-electron coupling, \( g = 1 \) for each time, indexed by `t` and width.

- `t_()/width_()` - real
    - Dim : [`n_mX`]
    - Total absorption rate at each \( \omega = m_\chi \).
    - Units : None

## self_energies

Self energies, \( \Pi_{\mathbf{O}_1, \mathbf{O}_2} \), where the \( \mathcal{O} \)'s are specified as the entry subscripts.

- `width_()/pi_v2_v2_r` - real
    - Dim : [`n_mX`]
    - Real part of \( \Pi_{\bar{v}^2, \bar{v}^2} \).
    - Units : \( \text{eV}^2 \)
- `width_()/pi_v2_v2_c` - real
    - Dim : [`n_mX`]
    - Complex part of \( \Pi_{\bar{v}^2, \bar{v}^2} \).
    - Units : \( \text{eV}^2 \)
- `width_()/pi_1_1_mat_r` - real
    - Dim : [3, 3, `n_mX`]
    - Real part of \( \mathbf{\Pi}_{1,1} \), defined such that \( \frac{\mathbf{q}}{m_e} \cdot \mathbf{\Pi}_{11} \cdot \frac{\mathbf{q}}{m_e} = \Pi_{11} \).
    - Units : \( \text{eV}^2 \)
- `width_()/pi_1_1_mat_c` - real
    - Dim : [3, 3, `n_mX`]
    - Complex part of \( \mathbf{\Pi}_{1,1} \), defined such that \( \frac{\mathbf{q}}{m_e} \cdot \mathbf{\Pi}_{11} \cdot \frac{\mathbf{q}}{m_e} = \Pi_{11} \).
    - Units : \( \text{eV}^2 \)
- `width_()/pi_vi_vj_r` - real
    - Dim : [3, 3, `n_mX`]
    - Real part of \( \mathbf{\Pi}_{v^i,v^j} \).
    - Units : \( \text{eV}^2 \)
- `width_()/pi_vi_vj_c` - real
    - Dim : [3, 3, `n_mX`]
    - Complex part of \( \mathbf{\Pi}_{v^i,v^j} \).
    - Units : \( \text{eV}^2 \)

## abs_tran_form

Transition form factors, \( f_{ii'\mathbf{k}}, \mathbf{f}_{ii'\mathbf{k}}, \widetilde{f}_{ii'\mathbf{k}} \) for absorption calculation, indexed by initial state, `init`, and final state, `fin`.

- `init_()/fin_()/tran_form_v2_c`
    - Dim : [`n_k`] ([`n_k`, 2, 2] if wave functions are spin dependent)
    - Complex part of \( \widetilde{f}_{ii'} \) for all \( \mathbf{k} \).
    - Units : None
- `init_()/fin_()/tran_form_v2_r`
    - Dim : [`n_k`] ([`n_k`, 2, 2] if wave functions are spin dependent)
    - Real part of \( \widetilde{f}_{ii'} \) for all \( \mathbf{k} \).
    - Units : None
- `init_()/fin_()/tran_form_v_c`
    - Dim : [`n_k`] ([`n_k`, 2, 2] if wave functions are spin dependent)
    - Complex part of \( \mathbf{f}_{ii'} \) for all \( \mathbf{k} \).
    - Units : None
- `init_()/fin_()/tran_form_v_r`
    - Dim : [`n_k`] ([`n_k`, 2, 2] if wave functions are spin dependent)
    - Real part of \( \mathbf{f}_{ii'} \) for all \( \mathbf{k} \).
    - Units : None
- `init_()/fin_()/tran_form_1_r`
    - Dim : [`n_k`] ([`n_k`, 2, 2] if wave functions are spin dependent)
    - Real part of \( f_{ii'} \) for all \( \mathbf{k} \). (Note: not written if wave functions are spin independent, since this is trivially 0 in that case.)
    - Units : None
- `init_()/fin_()/tran_form_1_c`
    - Dim : [`n_k`] ([`n_k`, 2, 2] if wave functions are spin dependent)
    - Complex part of \( f_{ii'} \) for all \( \mathbf{k} \). (Note: not written if wave functions are spin independent, since this is trivially 0 in that case.)
    - Units : None

## screening

- `type` - str
    - Screening type.
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

## version

- `number` - str
    - Version number of `EXCEED-DM`.

# Warnings

- `Fortran` is a **column-major order** language, which affects how HDF5 datasets are read in. For example, `python` is a **row-major order** language, and therefore for an \( n \times m \) matrix to be read in to `Fortran`, it must be saved from `python` as a \( m \times n \) matrix.
