title: Binned Wave Function Coefficients (wfc)
author: Tanner Trickle
date:v0.2.6

{!utilities/binned_wfc/README.md!}

# Output File 

## Groups

### PW_dataset

See `EXCEED-DM Output File` documentation.

## Datasets

- `binned_wfc_FT_sq` - real - \( \langle |\widetilde{u}_i|^2 \rangle(q, \Delta q) \)
    - Dim: [ `n_val` + `n_cond`, `n_q_bins` ]
    - PW wave function coefficients binned in momentum transfer.

# Warnings

- `Fortran` is a **column-major order** language, which affects how HDF5 datasets are read in. For example, `python` is a **row-major order** language, and therefore for an \( n \times m \) matrix to be read in to `Fortran`, it must be saved from `python` as a \( m \times n \) matrix.
