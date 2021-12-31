## ID: 5

- Target: Si

- Process: scattering

- Description: Spin independent valence to conduction transition from highest valence band to lowest conduction band. Single dark matter mass and form factor. Screening effects included by numerically computing the dielectric and screening the rate with it.

If `./output/dielectric.hdf5` does not exist the dielectric will be computed from scratch. If the file does exist, e.g. after running the example once, then the file will be loaded and will not be computed from scratch.

Note that the total scattering rate should be smaller than the one computed in example #1.
