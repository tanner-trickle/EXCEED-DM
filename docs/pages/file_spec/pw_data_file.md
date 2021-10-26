title: PW Data File
author: Tanner Trickle
date: v0.2.6

\(\require{mediawiki-texvc}\)

**Filetype**: HDF5

**Extension**: `.h5`, `.hdf5`

**Example**: `./examples/Si/dft/Si_2x2x2_AE.hdf5`

Holds data about the Bloch wave function coefficients, in the Fourier basis, \( \widetilde{u}_{i, \mathbf{k}, \mathbf{G}, s} \),
$$\begin{align}
    \psi_{i, \mathbf{k}, s}(\mathbf{x}) = \frac{1}{\sqrt{V}} 
    \sum_\mathbf{G} \widetilde{u}_{i, \mathbf{k}, \mathbf{G}, s} e^{i (\mathbf{k} + \mathbf{G}) \cdot \mathbf{x}}
    \nonumber
\end{align}$$
where \( s \) denotes the component in the \( S_z \) basis. Spin dependence is optional.

# Datasets

- `n_k` - integer - \( N_\mathbf{k} \)
    - Number of \( \mathbf{k} \) points.
- `n_val` - integer - \( N_\mathrm{val} \)
    - Number of valence bands.
- `n_cond` - integer - \( N_\mathrm{cond} \)
    - Number of conduction bands.
- `n_G` - integer
    - Number of \( \mathbf{G} \) points.
- `Ef_max` - real - \( E_f^\mathrm{max} \)
    - Only conduction bands which reach below `Ef_max` are included in the calculation of the band structure. Relative to Fermi surface/valence band maximum.
    - Units : eV
- `G_grid_red` - integer 
    - Dim : [`n_G`, 3]
    - \( \mathbf{G} \) vectors in reduced coordinates.
- `a_vecs_A` - real - \( \mathbf{a}_i \) 
    - Dim : [3, 3]
    - Primitive lattice vectors. Formatted such that `a_vecs_A(i, :)` is the ith primitive lattice vector, \( \mathbf{a}_i \).
    - Units: \( \AA \)
- `b_vecs_A` - real  - \( \mathbf{b}_i \)
    - Dim : [3, 3]
    - Reciprocal lattice vectors. Formatted such that `b_vecs_A(i, :)` is the ith reciprocal lattice vector, \( \mathbf{b}_i \).
    - Units: \( \AA^{-1} \)
- `k_weight` - real - \( w_\mathbf{k} \) 
    - Dim : [`n_k`]
    - Weight associated with each \( \mathbf{k} \) point for integration over the 1BZ.
    - Units: None
- `k_grid_red` - real
    - Dim : [`n_k`, 3]
    - \( \mathbf{k} \) vectors in reduced coordinates.
    - Units: None
- `energy_bands` - real - \( E_{i, \mathbf{k}} \) 
    - Dim : [`n_k`, `n_bands`]
    - `n_bands = n_val + n_cond`
    - Energy of each band, for each \( \mathbf{k} \) point.

# Groups

- `wfc_FT_r`
    - `wfc_FT_r/i_<band number>`
        - `wfc_FT_r/i_<band number>/k_<k number>`
            - Dim : [`n_G`] ([`n_G`, 2] if spin dependent wave functions)
            - Real part of the Bloch wave function coefficients for the `i`th band and `k`th k point, \( \text{Re} \left( \widetilde{u}_{i, \mathbf{k}, \mathbf{G}} \right) \). \( 1 \leq i \leq N_\mathrm{val} + N_\mathrm{cond} \) and \( 1 \leq k \leq N_\mathbf{k} \).

- `wfc_FT_c`
    - `wfc_FT_c/i_<band number>`
        - `wfc_FT_c/i_<band number>/k_<k number>`
            - Dim : [`n_G`] ([`n_G`, 2] if spin dependent wave functions)
            - Imaginary part of the Bloch wave function coefficients for the `i`th band and `k`th k point, \( \text{Im} \left( \widetilde{u}_{i, \mathbf{k}, \mathbf{G}} \right) \). \( 1 \leq i \leq N_\mathrm{val} + N_\mathrm{cond} \) and \( 1 \leq k \leq N_\mathbf{k} \).

# Notes

- Orthogonality of lattice vectors, \( \mathbf{a}_i \cdot \mathbf{b}_j = 2 \pi \, \delta_{ij} \).

# Warnings

- `Fortran` is a **column-major order** language, which affects how HDF5 datasets are read in. For example, `python` is a **row-major order** language, and therefore for an \( n \times m \) matrix to be read in to `Fortran`, it must be saved from `python` as a \( m \times n \) matrix.

