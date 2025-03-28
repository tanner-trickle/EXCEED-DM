v1.1.0
---

- Update to compute results shown in: https://arxiv.org/abs/2303.11344 .
    - Added particle_type = 'mdm', 'edm' options for calculation of the absorption rate of DM coupling to electron MDM/EDM.
    - Added a new type of electronic state approximation, 'elec_state_atomic_type', which has two sub-types corresponding to bound states expanded in an STO basis, 'elec_state_atomic_STO_basis_type', and proper continuum, Coloumb wave function states which solve V = -Z/r, 'elec_state_atomic_continuum_type'. These are suitable for describing an 'atomic target', see the above paper for more details.
        - These states are defined by the property that they can be written as \psi(x) = R(r) Y(\theta, \phi), where Y are spherical harmonics. Future subtypes may be incorporated if they provide a 'compute_radial_wf' and 'compute_wf' procedure.
        - These are currently only used in absorption rate calculations.
    - Added a new method for computing transition matrix elements between atomic states, in the q -> 0 limit.
        - Currently, only transitions from atomic states to atomic states are supported, i.e., there is no way to calculation bloch state -> atomic state, or vice versa. Additionally, only initial STO -> final continuum is supported since the STO basis provides a nice analytic formula for the first and second derivatives. Initial atomic STO -> final atomic STO may be added with relative ease, whereas initial continuum -> continuum would be more difficult (although it is unclear what circumstance that would be important).  
    - Added new self-energy, and corresponding absorption rate, calculations for the MDM/EDM models. 
    - Added atomic target and crystal target EDM absorption examples.
    - Updated input file documentation.

v1.0.0
---

- Release of v1.0.0 to accompany the user manual: https://arxiv.org/abs/2210.14917.

v1.0.0-beta
---

- New Electronic Configuration file
    - Contains all the relevant information about the electronic configuration, no more need to separate ‘core’, ‘free’, and ‘valence and conduction’ states
- Extended absorption and dielectric calculations
    - Include all states (e.g. core and free)
- Separated electronic state and calculation dependence
    - Generalized input electronic states, initial/final states can now be specified in different bases, i.e., PW basis, STO basis, or single PW’s (previously ‘free’ states)
        - For example, allows for calculation of STO basis → STO basis transitions
    - Easily add different bases
    - Only need to specify how matrix elements are computed for a specific approximation of the wave functions, the ‘particle physics’ part of the code base doesn’t need to change
- Automatically generated input documentation
- Increased use of CMake and preprocessor directives
    - Version number automatically passed to program and documentation
- Use hdf5_utils for even easier access to data inside an hdf5 file
- Building new documentation with Sphinx
- Improved parallelizability
    - Everything is parallelized over total number of transitions

    - All calculations use MPI_Reduce commands
- Removed dependence on k_weight which only caused confusion (sum to 1 or 2?)
- New input format - Configuration File
    - Major reasons
        - Dynamically sized 1D arrays
            - No need to specify number of elements in an input array
        - Enforced defaults
        - Easily add documentation

        - Split in to as many files as you want
    - Minor reasons
        - No need for ! at end of file
- Added many new examples which all use the new input, and electronic configuration files.

v0.3.0
---

- Updated calculation of some self energies to avoid numerically small numbers. For example, the dielectric computed with either \( \Pi_{1, 1} \) or \( \Pi_{v^i, v^j} \) is identical. 
- Added the calculation of \( \Pi_{\mathbf{v} \cdot \sigma, \mathbf{v} \cdot \sigma} \).
    - Used in generalizing the pseudoscalar absorption rate for spin dependent wave functions.  
- Added a timer type for more granular timing data.
- Complete rework of examples so each example has a specific ID and `readme.md` description file.
    - Removed tests in favor of a script which runs all the examples, `examples/run_all_examples.sh`. 
- Speed improvement of the dielectric calculation in the scattering ( \( q \gg \omega \) ) regime.

v0.2.8
---

- Minor updates to `dm_model` and `PW_dataset` types
- Added routine to find the maximum momentum transfer within the 1BZ, `get_q_max_1BZ`.
- Generalized dielectric inputs, anisotropic dielectrics in each (q, omega) bin are now accepted. See the `Dielectric Data File` documentation for more details. 
- Added `cc_ext` calculation mode, a new way to compute the core -> conduction contribution which allows the q range to extend to the kinematic bound.
- Added a new integration routine when going from differential quantities such as \( \frac{dR}{d \omega} \), \( \frac{dR}{dq} \) to the binned values.
- Added timing information to printed and saved output.
- Updated documentation.

v0.2.7
---

- Fixed bug with `Ef_max` cutoff by shifting valence band maximum to be the zero-point energy.
- Updated build system to `CMake`, removing FoBiS entirely.
- The OpenMPI subroutine, `MPI_Reduce` is used when appropriate verus manual send and recieve calls.
    - Currently implemented in scattering rate and dielectric calculations
- Temporary file is deleted when creating the core electron configuration file with `utilities/core_elec_config/create_elec_config.ipynb`
- Removed `TODO.md` which was not being used (may reintroduce in the future if it becomes useful).
- Updated the `binned_wfc` utility.
- Added documentation for the `binned_wfc` utility.
- Optimized dielectric calculation, improved binned procedure to save memory.
- Update to documentation and refactored main folder.
- Added [home page](https://exceed-dm.caltech.edu) to `README.md`.

v0.2.6
---

- Added the option to build with `cmake`. This will become the default build option in `v1.0.0`. Installation instructions can be found in `install-cmake.md` or on the website.   
- Communication of scattering rate data is now done with an `MPI_Reduce` command versus manual send and recv calls.
- Added a timing subroutine when computing the dielectric.
- Miscellaneous updates to documentation.
- Minor bug fixes: some variables were not saving to the output files. 

v0.2.5
---

- **Major code refactor** and other miscellaneous improvements on the way to v1.0.0
    - Changed valence -> conduction parallelization scheme from {i, i'} -> {i, k}, allowing for significant speedups when ( number of valence bands ) x ( number of conduction bands ) < ( number of processors )
    - Added `n_val_max`, `n_cond_max` options to specify the maximum number of valence and conduction bands to keep in the calculation (valence -> conduction scattering, absorption, and dielectric calculations). `n_val_max` counts down from the Fermi surface and `n_cond_max` counts up from the Fermi surface.
    - Changed core -> conduction parallelization scheme from {i, i'} -> {i, kf}, allowing for significant speedups when ( number of core states ) x ( number of conduction bands ) < ( number of processors )
    - Added `n_prinicipal_min`, `n_principal_max` to set the minimum and maximum principal quantum numbers to include in scattering rate calculations involving core initial states.
    - Changed valence -> free parallelization scheme to {i, k}.
    - Changed absorption parallelization scheme to {i, k}
    - Added option to skip saving transition form factors in absorption calculation.
    - Reworked dielectric calculation, added more namelist input options.
    - Changed dielectric parallelization scheme from {i, i'} -> {i, k}, allowing for significant speedups when ( number of valence bands ) x ( number of conduction bands ) < ( number of processors )
    - Added `dielectric` process which computes just the dielectric.
    - Utilizing `type` structures to make code more modular/reusable. Each type has at least its own `load`, `save`, and `print` procedure.
        - `PW_dataset` - handles the plane wave Bloch wave function coefficients 
        - `dm_model` - dark matter model parameters
        - `material` - collection of target material parameters
        - `expt` - experimental parameters
        - `core_electron` - Core electron configuration and STO wave function coefficients for the core electrons
        - and more!
    - Output scattering rate units are now cm^(-2), for easier conversion to cross section constraints.
    - Output absorption rates are total rates given use specified experimental masses and exposures.
    - `EXCEED-DM` version is now written to output for easy comparison with previous and future versions.
    - Added `k_` subgroup to Bloch coefficient data structure for easy access to the wave function coefficients at a given `i`, `k`.
        - Updated example files in `examples/dft` accordingly.
    - Updated core electron configuration file specification
        - Updated example files in `examples/(Si, Ge)/core` accordingly.
        - Updated `utilities/create_elec_config.ipynb`.
    - Added example input files for Germanium.
    - Improved/standardized output printing with an `info_messages` module.
    - Improved comments inside the code, many variable explanations have associated LaTeX'ed equations which can be read by viewing the documentation in a browser.
    - Major update to documentation. Check it out [here](https://tanner-trickle.github.io) folder.
    - All modules, procedures, and types have some documentation.
    - Specific documentation pages for all input and output files.
    - `/examples/Si/dft/Si_2x2x2_AE_spin.hdf5` are now realistic spin dependent wave function coefficients for Si (just spin-independent ones doubled.).
    - Added preliminary logo, `docs/media/exdm-prelim-logo.png`.

v0.2.4
---

- Reworked the implementation of the vector/pseudoscalar DM absorption calculation.
- Fixed bug in absorption rate calculation introduced in v0.2.3 which overwrote the main processors transition form factors when computing the self energies.

v0.2.3
---

- Calculation of the dielectric for targets with spin-dependent wave functions is now supported.
- Added routine to time the dielectric calculation.
- Added `spin_degen` which accounts for the spin degeneracy factor of the initial (valence) states.
- Added `scalar_LO` absorption calculation mode to just compute the leading order contribution for scalar DM absorption.
- Generalized vector/pseudoscalar DM absorption calculation for anisotropic targets.
  - Self-energy, Pi_vi_vj, is now computed along with the other self energies.
- Implemented routine to rigorously find maximum magnitude of q for which an FFT will give consistent results across meshes. See `find_q_max_FFT` routine in `FFT_util`.
- Updated `install.md`

v0.2.2
---

- Improvements to `absorption` module.
  - Initial support for spin-dependent wave functions.
    - More general transition form factors can be computed.
  - Parallelization of velocity integral
- **LAPACK and BLAS are now required.** Additional installation instructions have been added to the `ubuntu-gnu` build.
- `q_s_FFT` has been removed as a `numerics` input option. This is now computed directly with `LAPACK` routines.

v0.2.1
---

- Added a `dielectric` module which computes the dielectric in the scattering kinematic regime.
  - To screen the rate with a numeric dielectric model, set `screen_type = numeric`.
  - If `load_dielectric_from_file = .FALSE.`, the dielectric will be computed from scratch. Otherwise the screening factor will come from the dielectric matrix in the `dielectric_filename` file.
  - Variables relevant for the calculation are loaded through the `dielectric` namelist. Check out examples in the `examples/` folder.
  - Note: the dielectric here is only used to screen. Currently only v -> c transitions are included in the loop. For now, spin-indpendent wave functions only. 
- Added `ubuntu-gnu-debug` build mode to debug with.

v0.2.0
---

- Calculation of dark matter (scalar, pseudoscalar (axion-like), vector) absorption on electrons!
  - Compute the absorption rate, self-energies, and generalized crystal form factors in the absorption limit.
  - Setting `process = 'absorption'` in the input file switches the calculation to absorption mode. Check out examples in the `examples/` folder.
  - Set a variety of electron lifetime/width parameters, `width = min( a + b omega , width_max )` 
  - See https://arxiv.org/abs/2106.12586 for details of the formulation.
  - (Currently only spin-independent, valence -> conduction transitions are supported.) 

v0.1.3
---

- Partial implementation of general transition form factor in module `transition_form_factor`.
    - Compute scattering rates for interactions that depend on electron spin with spin-dependent electronic wave functions.
    - Support for valence -> conduction transitions.
    - Note : default is to compute spin-independent scattering rates
- Added simple testing routines which will help make sure new additions do not break old functionality.
    - Run test with `FoBiS.py rule -ex tests`.
    - To add new tests just add files to the list inside `tests/run_tests.sh`
- Updated documentation

v0.1.2
---

- Added functionality for spin dependent (two component) wave functions!
  - Automatically detect whether input DFT data has spin dependence.
  - Valence -> conduction transition rates can be computed with spin dependent wave functions.
  - Added example input file for Si which has spin dependent wave functions.

v0.1.1
---

- When timer = .TRUE. a smaller version of the program will be run before the main program, and an estimate of the run time of the full program will be printed.
- Improved output printing.
- Updated Si example input files to use file paths relative to the main folder.

v0.1.0
---

- Initial beta release.
