title: Dielectric Data File
author: Tanner Trickle
date: v0.2.5

Numerically computed dielectric. Can be used as an input file to screen scattering rate calculations.

# Groups

## dielectric

- `dielectric_r`
    - Dim : [`n_E`, `n_q`, `n_q_theta`, `n_q_phi`]
    - Real part of the dielectric function, \( \varepsilon \) averaged in \( \mathbf{q}, \omega \) bins.
- `dielectric_c`
    - Dim : [`n_E`, `n_q`, `n_q_theta`, `n_q_phi`]
    - Complex part of the dielectric function, \( \varepsilon \) averaged in \( \mathbf{q}, \omega \) bins.

## bins_dielectric

- `n_E` - int
    - Number of \( \omega \) bins.
- `n_q` - int
    - Number of \( q \) bins.
- `n_q_theta` - int
    - Number of \( \theta_q \) bins.
- `n_q_phi` - int
    - Number of \( \phi_q \) bins.
- `E_width` - real
    - Width of bins in \( \omega \).
    - Units : \( \text{eV} \)
- `q_width` - real
    - Width of bins in \( q \)
    - Units : \( \text{eV} \)

## material

- `name` - str
    - Target material name.
- `band_gap` - real - \( E_g \)
    - Target material band gap.
    - Units : eV
- `density` - real - \( \rho_T \)
    - Target material density.
    - Units : \( \text{g} \; \text{cm}^{-2} \)

## numerics_dielectric

Numerics parameters specific to the dielectric calculation.

- `n_val_max` - int
    - Maximum number of valence bands included in the calculation.
- `n_cond_max` - int
    - Maximum number of conduction bands included in the calculation.
- `n_k_vec` - int
    - Dim : [3]
    - Number of each dimension of a Monkhorst-Pack sampled 1BZ.

## version

- `number` - str
    - Version number of `EXCEED-DM`.

# Warnings

- `Fortran` is a **column-major order** language, which affects how HDF5 datasets are read in. For example, `python` is a **row-major order** language, and therefore for an \( n \times m \) matrix to be read in to `Fortran`, it must be saved from `python` as a \( m \times n \) matrix.
