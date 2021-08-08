title: RHF WF Data File
author: Tanner Trickle
date: v0.2.5

\(\require{mediawiki-texvc}\)

**Filetype**: HDF5

**Extension**: `.h5`, `.hdf5`

**Example**: `./examples/RHF_wf_data.hdf5`

RHF wave function coefficients, \( C_{nlj}, Z_{jl}, n_{jl}, n_j \), which semi-analytically model the atomic wave functions,
$$\begin{align*}
    \psi^\text{atom}_{nlm}(\mathbf{x}) & = \sum_{j = 1}^{n_j} C_{nlj} R_\text{STO}(x; Z_{jl}, n_{jl}) Y^m_l( \hat{\mathbf{x}} ) \\
    R_\text{STO}(r; Z, n) & = a_0^{-3/2} \frac{(2 Z)^{ n + \frac{1}{2}} }{\sqrt{(2n)!}} \left( \frac{r}{a_0} \right)^{n - 1} e^{-Zr/a_0}
\end{align*}$$
where \( a_0 \) is the Bohr radius.

The data inside `examples/RHF_wf_data.hdf5` are the tabulated coefficients for all atoms with \( 2 \leq Z \leq 54 \) and can be found [here](https://www.sciencedirect.com/science/article/pii/S0092640X8371003X?via%3Dihub).

# Groups

- `Z_()/n_()/l_()/C_lnj` - real
    - Dim : [`nj`]
    - \( C_{nlj} \) STO coefficients to compute the RHF wave functions for an atom with proton number, \( Z \), and quantum numbers, \( n, l \).
- `Z_()/n_()/l_()/N_lj` - real
    - Dim : [`nj`]
    - \( N_{lj} \) STO coefficients to compute the RHF wave functions for an atom with proton number, \( Z \), and quantum numbers, \( n, l \).
- `Z_()/n_()/l_()/Z_lj` - real
    - Dim : [`nj`]
    - \( Z_{jl} \) STO coefficients to compute the RHF wave functions for an atom with proton number, \( Z \), and quantum numbers, \( n, l \).
- `Z_()/n_()/l_()/n_lj` - real
    - Dim : [`nj`]
    - \( n_{jl} \) STO coefficients to compute the RHF wave functions for an atom with proton number, \( Z \), and quantum numbers, \( n, l \).
- `Z_()/n_()/l_()/nj` - real
    - Dim : [`nj`]
    - \( n_{j} \) number of STO coefficients to compute the RHF wave functions for an atom with proton number, \( Z \), and quantum numbers, \( n, l \).
- `Z_()/n_()/l_()/energy` - real
    - Energy of state with proton number, \( Z \), and quantum numbers, \( n, l \).
    - Units : eV

# Warnings

- `Fortran` is a **column-major order** language, which affects how HDF5 datasets are read in. For example, `python` is a **row-major order** language, and therefore for an \( n \times m \) matrix to be read in to `Fortran`, it must be saved from `python` as a \( m \times n \) matrix.

