title: Core Electron Configuration File
author: Tanner Trickle
date: v0.2.5

\(\require{mediawiki-texvc}\)

**Filetype**: HDF5

**Extension**: `.h5`, `.hdf5`

**Example**: `./examples/Si/core/Si_core_elec_config.hdf5`

Holds data about the core electron configuration.

# Datasets
- `n_atom` - integer
    - Number of atoms in the primitive cell.
- `n_state` - integer
    - Total number of core electron states. Sum of individual atoms number of core electron states.
- `eq_pos_red` - real 
    - Dim : [`n_atom`, 3]
    - Equilibrium positions, \( \mathbf{x}_0 \), of each atom.
- `energy` - real - Dim : [`n_state`]
    - Binding energy of each state. 
- `Z` - real 
    - Dim : [`n_atom`]
    - Proton number of each atom.
- `config` - integer 
    - [`n_state`, 5]
    - Map between each core state and the atom and standard quantum numbers it corresponds to.
        - `config(i, 1)` : atom ID 
        - `config(i, 2)` : n
        - `config(i, 3)` : l
        - `config(i, 4)` : m
        - `config(i, 5)` : number of spin states

# Warnings

- `Fortran` is a **column-major order** language, which affects how HDF5 datasets are read in. For example, `python` is a **row-major order** language, and therefore for an \( n \times m \) matrix to be read in to `Fortran`, it must be saved from `python` as a \( m \times n \) matrix.

