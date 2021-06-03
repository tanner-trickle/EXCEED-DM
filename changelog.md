v0.1.3
---

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

---


- Added functionality for spin dependent (two component) wave functions!
  - Automatically detect whether input DFT data has spin dependence.
  - Valence -> conduction transition rates can be computed with spin dependent wave functions.
  - Added example input file for Si which has spin dependent wave functions.

v0.1.1
---

---

- When timer = .TRUE. a smaller version of the program will be run before the main program, and an estimate of the run time of the full program will be printed.
- Improved output printing.
- Updated Si example input files to use file paths relative to the main folder.

v0.1.0
---

---

- Initial beta release.
