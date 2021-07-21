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
