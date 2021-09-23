* Release Notes PlatoSim 3.5.0



** Improvements

*** Update `showSim.py` to include biasMapsRight and biasMapsLeft
*** New analytic PSF model and set of parameters for N6000K
*** Include more accurate PSF files for mapped PSF model. The new files can be downloaded from the `Prerequisites` section of the PlatoSim website. 
*** Updated mail.cpp to allow a new log level: 0: only shows errors and no warning. 
*** Update `Simulation::writeInputParametersToHDF5` function.
*** Implemented a new mapped distortion method for the mapped PSF model. The distortion table is included in psf files. 
*** Updated website 




** Bug fixes

*** Fixed issue where timeshift was applied when reading out CCDs for the F-CAMs (GitHub #540)
*** In `hdf5ToFits.py` typecheck before converting to `string` (GitHub #600)
*** Removed fortran dependencies in fftw  install script



    
** New features/functionality

*** Added a more accurate aberration model. Instead of assuming a circular orbit with constant speed around the sun, we can now include the path of the spacecraft in an orbit file to simulate any time-dependent velocity. An accurate orbit file is included in the `inputfiles` directory. 
