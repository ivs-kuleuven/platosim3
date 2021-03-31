* Release Notes PlatoSIM 3.4.0


** Improvements

*** Updated HDF5 library from v10.2 to v12.0

*** Provide package acces to python scipts (GitHub feature request #548)
The python code for PlatoSim3 is now packed in a package called 'platosim'. This makes it easier to combine 
the python code with other projects. To import python modules one can simply ~from platosim import module_name~
or ~import platosim.module_name~.

*** Replaced natural & mechanical vignetting with overall relative transmissivity (GitHub issue #478)
The vignetting (natural + mechanical) effect is replaced by the overall transmissivity. This consists out of the following 
contributions:
- natural vignetting
- mechanical vignetting
- glass absorption + AR coating.

*** Read BFE coefficients from file instead of calculating them (GitHub issue #434)
The Brighter Fatter Effect is now calculated using coefficients that are read from a file. The file's location should be 
specified in the input.yaml file CCD/BFE/CoefficientsFileName.

*** CCD ID new standard (GitHub issue #416)
The CCD ID does now follow a different standard to match the latest convention used in the 
PLATO-OHB-PL-11-0009 i4.0 with respect to the focal plane.

*** Implemented time-dependent CTI (GitHub issue #476)
Implemented time-dependent CTI trap densities for the Short et al. treatment of the CTI.

*** Unique ID needed in ClosedLoopUtility and Log-File (GitHub issue #488)
Implemented a unique identity string for both sockets in ClosedLoopUtility.cpp 
     

** Bug Fixes

*** PSF output contains NaN values (GitHub issue #543)

*** Photon noise applied after CTI (GitHub issue #431)

*** CTI by Short et al. has no effect on the maximum value (GitHub issue #403)

*** CTI Short2013 model has no effect (GitHub issue #446)

*** Incorrect FEE temperature, from file (GitHub issue #413)

*** BFE seems to increase total flux (GitHub issue #397)

*** Simulation with too short external jitter file (GitHub issue #327)

*** Input file of the photometry tutorial needs to be updated (GitHub issue #368)

*** Edge effects introduced by convolution with non-analytical PSF (GitHub issue #357)

*** Full-image smearing strange value for column 2254 (GitHub issue #450)

*** isnan ambiguous in Detector.cpp (GitHub issue #492)

*** Some trouble with the routine createStarCatalogFileFromPixelCoordinates (GitHub issue #506)

*** Analytical Non Gaussian PSF: strange orientation of the PSF on the CCD (GitHub issue #511)

*** Use of a MappedFromFileAsymmetrical PSF with sub-pixel resolution is 1/64 generates an error (GitHub issue #530)

    
** New features/functionality

*** Validation & verification notebooks
Added validation notebooks in /docs/validation/. These notebooks give a graphical validation of functionality of PlatoSim.  

*** Validation & verification scripts
Added validation tests to the /tests/validationTests/scripts/ directory. All of the tests can be ran individually and 
the restult of these tests is then printed out to the terminal. The input and output files of the simulation run are
stored into the /tests/validationTests/ioFiles/ directory. Alternatively, all the test can be run automatically by running 
the ~run.py~ file. The results of these tests is then both printed out to the terminal and saved into the /ioFiles/result.txt
file. 

*** Added a method to create a single point source
Added the file python/platosim/InserPointSource.py. This module contains functionality to create a starcatlog with a single point
located at a specified position. 

*** Record cosmic particle hits in output HDF5 file (GitHub feature request #553)
It is now possible to have PlatoSim record the columns and row pixel values and the flux of the cosmic particles. These values can be 
found in the output HDF5 file. This behavior can be switched on or off using the WriteCosmics paramter in the input file. 

*** Added method to obtain cosmics into SimFile class
Added the method ~getCosmicsCoordinates()~ in the Simfile class. This method returns column, row and flux value of the recorded cosmics. 

*** 6.25s time-shift between CCDs (GitHub issue #401)
Added the  timeshift that exists between the readout of the CCDs. 

*** Using non-rotationally symmetrical PSF from file (GitHub issue #407)
It is possible to obtain the PSF from a file based on the location of the star on the CCD. This is done by using the following configuration
~PSF/Model = "MappedFromFileAsymmetrical"~ and specifying the path location of file with PSF in ~PSF/MappedFromFileAsymmetrical/Filename = "path_to_file.hdf5"~.

*** Access methods for datasets at top level in HDF5 file (GitHub issue #437)
Implemented a function readStringDatasetAttribute(string groupName, string datasetName, string attributeName) for string attribute access in HDF5 files in the file ~HDF5.cpp~. 

*** Conversion between CCD coordinates & field angles (GitHub issue #460)
Added the functions focalPlaneAngles2pixelCoordinates() and pixelCoordinates2FocalPlaneAngles() to referenceFrames.py 

*** Transmission variation as a function of FOV (GitHub issue #415 and GitHub issue #435)
Mechanical vignetting was previously implemented as a sharp blockage of incoming flux beyond the edge of the field of view. However this effect should already
cause a decrease in efficiency at smaller distances from the opptical axis. This gradual decrease in efficiency is now implemented. 

*** Group all contributions to the overall relative transmissivity (GitHub issue #478)
Also mentioned in the Improvements sections.  
Overall transmissivity combines the effect of: 
- natural + mechanical vingnetting
- glass absorption + AR coating

*** Implement ghosts (GitHub issue #515)
Implemented symmetric point-like ghosts and extended ghosts. 

Point-like ghosts that appear on the focal plane at position (xFP, yFP) are caused by stars that 
appear on the focal plane at position (-xFP, -fFP). These ghosts only appear within 8° of the optical axis. 

Extended ghosts that are caused by star at position (xFP, yFP) on the focal plane appear on the position (1.065 * xFP, 1.065 * yFP). The total flux of these ghosts 
is 0.00003% of the total flux of the original star, and they appear all the way to the edge of the field of view. 

*** Add linear irradiance ratio decrease of point-like ghosts (GitHub issue #526)
The flux ratio on-axis is 0.08% and decreases linearly to 0% at 8° after which no point-like ghosts exists.

*** Added getPointLikeGhostCoordinates() and getExtendedGhostCoordinates() methods to SimFile class
These methods are defined in the SimFile class in the /python/platosim/simfile.py file. These methods return the star ID, pixel and focal-plane 
coordinates, flux and in case of the ~getExtendedGhostCoordinates()~ method also the radius of the respective ghosts. 




