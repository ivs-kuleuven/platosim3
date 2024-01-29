* Changelog PlatoSim 3.6.1

** Improvements

*** Redundant time column in `ACS`, `StellarPositions`, and `GhostPositions` have been removed and added as an individual time column (see addition below). This is breaking backward compatible change only for users that do not use the `SimFile.py` class.

*** Small bugfixes to code

*** Small errors in documentation

*** Reduce time it takes to run the validation tests. (Issue #869)

*** Removed declared but unused variables in Sky.cpp





** Bug fixes

*** Bugfix in DatectorWithAnalyticNonGaussian::applyPhotometry() module where not the right mask was used. (Issue #913)

*** Bugfix in generation of distortionmap with mapped PSF (Issue #811)

*** Simple CTI is no longer written to HDF5 file (Issue #863)

*** Fixed sign error in the telescope tilt rotation matrix in both Telescope.cpp and referenceFrames.py (Issue #857)



** New features/functionality

*** New time column is saved to the HDF5 file by default.

*** The PLATOnium toolkit can now simulate the F-CAMs.

*** Added AppleClang (difference from Clang) as a possible compiler.

*** Added non-linear gain.

*** Generally a lot of small bugfixes for the PLATOnium toolkit has been made (Issue #776, #883, #884, #885, #886, #887, #888, #889, #898, #908, #914 #916).

*** Jupyter tutorial notebooks did not work for user without a functional LaTeX installation (configured through `.matplotlibrc`). Now the module `import matplotlibrc` checks if the user has a valid LaTeX installation, and if so, activates the LaTeX rendering, and if not, fall back to normal text rendering.
