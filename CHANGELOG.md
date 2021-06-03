# Changelog for PlatoSim

<!-- 3.4.1 -->
<!-- ***** -->

## 02/06/2021: 3.4.1

## Fixed
* Corrected the implementation of the Jitter (GitHub issue #508)
* Corrected the implementation of the Drift
* Fixed the conda build in Jenkins

## Changed
* Added option to the method to include time in output of `getYawPitchRoll`
* Apply the BFE after full-well saturation (GitHub issue #584)
* Added the option to have conda install for python 3.6, 3.7, 3.8 and 3.9

## Added
* Test that checks the Jitter on different CCD's
* Added method `getYawPitchRollFromDrift` in the simfile.py file
* Option to save diffused PSF in output HDF5 file (GitHub issue #564)
* Added method `setSubfieldAroundPixelRows` to simulation.py
* Added option to (not) include in output file:
    - High resolution PSF (if PSF is not Analytic Gaussian)
    - Star Catalog
    - Platoform Yaw, Pitch, Roll
    - Transmission Efficiency



<!-- 3.4.0 -->
<!-- ***** -->

## 31/03/2021: 3.4.0

## Fixed
* HighResMap matrix is now initialized with zeros (GitHub issue #543)
* Photon noise applied after CTI (GitHub issue #431)
* CTI by Short et al. has no effect on the maximum value (GitHub issue #403)
* CTI Short2013 model has no effect (GitHub issue #446)
* Incorrect FEE temperature, from file (GitHub issue #413)
* BFE seems to increase total flux (GitHub issue #397)
* Simulation with too short external jitter file (GitHub issue #327)
* Input file of the photometry tutorial needs to be updated (GitHub issue #368)
* Edge effects introduced by convolution with non-analytical PSF (GitHub issue #357)
* Full-image smearing strange value for column 2254 (GitHub issue #450)
* isnan ambiguous in Detector.cpp (GitHub issue #492)
* Some trouble with the routine createStarCatalogFileFromPixelCoordinates (GitHub issue #506)
* Analytical Non Gaussian PSF: strange orientation of the PSF on the CCD (GitHub issue #511)
* Use of a MappedFromFileAsymmetrical PSF with sub-pixel resolution is 1/64 generates an error (GitHub issue #530)

## Changed
* Provide package acces to python scipts (GitHub feature request #548)
* Updated HDF5 library from v10.2 to v12.0 	
* Replaced natural & mechanical vignetting with overall relative transmissivity (GitHub issue #478)
* Read BFE coefficients from file instead of calculating them (GitHub issue #434)
* CCD ID new standard (GitHub issue #416)
* Implemented time-dependent CTI (GitHub issue #476)
* Unique ID needed in ClosedLoopUtility and Log-File (GitHub issue #488)
* Implemented time-dependent CTI (GitHub issue #476)
	
## Added
* Validation & verification notebooks
* Validation & verification scripts
* Added a method to create a single point source
* Record cosmic particle hits in output HDF5 file (GitHub feature request #553)
* Added method to obtain cosmics into SimFile class
* 6.25s time-shift between CCDs (GitHub issue #401)
* Using non-rotationally symmetrical PSF from file (GitHub issue #407)
* Access methods for datasets at top level in HDF5 file (GitHub issue #437)
* Conversion between CCD coordinates & field angles (GitHub issue #460)
* Mechanical vignetting < edge of FOV (GitHub issue #435)
* Transmission variation as a function of FOV (GitHub issue #415)
* Group all contributions to the overall relative transmissivity (GitHub issue #478)
* Implement ghosts (GitHub issue #515)
* Add linear irradiance ratio decrease of point-like ghosts (GitHub issue #526)
* Added getPointLikeGhostCoordinates() and getExtendedGhostCoordinates() methods to SimFile class

	


<!-- 3.3.7 -->
<!-- ***** -->

## 25/02/2020: 3.3.7

## Fixed

* Input file of the photometry tutorial needs to be updated (GitHub issue #368)
* Order of the effects (GitHub issue #394)

## Changed

* Input parameter update (GitHub issue #377):
  + wavelength range
  + irradiance over the PLATO wavelength range (PIS only; nothing changes for PlatoSim)
  + added: orientation angle of the solar panels
  + transmission efficiency
  + throughput bandwidth
  + FEE readout noise
  + gain + stability of the CCD gain
  + QE
  + digital saturation for fast cameras (this is now the same as for the normal cameras)
* Dynamic frame transfer times (GitHub issue #369)

## Added

* FEE over/undershoot (GitHub issue #376)
* Dump Analytical PSF within hdf5 output file (GitHub issue #379)
* Jitter from network





<!-- 3.3.6 -->
<!-- ***** -->

## 10/04/2019: 3.3.6

### Fixed

* Bias register map expressed as additional columns (GitHub issue #290)

* Implement partial readout (GitHub issue #285)

* Backward compatibility to the configuration files? (GitHub issue #292)

* Improve documentation on supplementary input files of PlatoSim (GitHub issue #308)

* Migtool error (GitHub issue #307)

* Update HDF5 dependency to 1.10.2 (GitHub issue #322)

* Remove python 3.5 dependency (GitHub issue #318)

* Improved error trapping for CCD/ReadoutMode/ReadoutMode (GitHub issue #302)

* Incorrect sky background level (GitHub issue #325)

* Incorrect flatfield level (GitHub issue #326)

* Segfault on running test harness (GitHub issue #329)

* Open-shutter smearing not accounted for (GitHub issue #339)

* Implement mechanical vignetting (GitHub issue #334)

* Partial-readout parameters not read out correctly (GitHub issue #346)

* Conserving disc space by writing the images as int matrices into the .hdf5 files (GitHub issue #348)

* Bug fix in the createStarCatalogFileFromPixelCoordinates() method.



### Added

* Temperature dependency of the dark current (space environment)



<!-- 3.3.5 -->
<!-- ***** -->

## 26/10/2018: 3.3.5

### Fixed

* Inconsistency in star position output (GitHub issue #294)





<!-- 3.3.4 -->
<!-- ***** -->

## 21/09/2018: 3.3.4

### Fixed

* bug when using demo_fgs.py with field distortion (GitHub issue #280)
* Wrong number of cosmics for small images (GitHub issue #283)



### Changed

* Cosmics can be enabled/disabled per area (image area / bias map / smearing map)

* Random seeds = -1 => use computer time instead (no longer fast-forward random distributions)



### Added

* Configurable log level

* Updated documentation





<!-- 3.3.3 -->
<!-- ***** -->

## 17/05/2018: 3.3.3

### Fixed

* Dark edge seen at bottom of the sub-field due to CTI by Short et al. (GitHub issue #263)

* User-given sky background not multiplied with the transmissivity of the optics (GitHub issue #265)

* Photon flux of stars should be floored instead of rounded (in the jitter steps) (GitHub issue #267)

* Open-shutter smearing outside sub-field should not take numRowsBiasMap into account (GitHub issue #269)





<!-- 3.3.2 -->
<!-- ***** -->

## 27/04/2018: 3.3.2

### Fixed

* Images datasets output dataype and specification (GitHub issue # 167)

### Added

* Documented h5ls and h5get

* Documented output control parameters

* Documented how to install via conda when no pop-up window would appear, asking for the credentials



### Changed

* Parameter values after release of v1.4 of the data package





<!-- 3.3.1 -->
<!-- ***** -->

## 17/04/2018: 3.3.1

### Fixed

* Path to executable set to /build in simulation.py (GitHub issue #257)



<!-- 3.3.0 -->
<!-- ***** -->

## 9/04/2018: 3.3.0

### Added

* Charge diffusion + jitter smoothing

* Documented installation via conda

* Use an external defined star ID in input star catalogues (GitHub issue #229)

* Dark current

* Brighter-fatter effect (BFE)

* Possibility of header lines and custom star IDs in the star catalogue ASCII files
    
* Stellar variability

* Created an <code>Examples</code> folder (in the <code>python</code> directory) where demo scripts show how to use the simulator from Python

* Safety checks in <code>Parameter.h</code> to ensure that the time series from a file has time points in strictly increasing order

* Flag to limit size of HDF5 output files

* Group in the input files (<code>ControlHDF5Content</code>) to control the content of the HDF5 output file

    - <code>WriteSubPixelImages</code> (moved from the <code>CCD</code> group: Boolean flag for writing the sub-pixel images to the HDF5 file [default=no]
    - <code>WriteStarPositions</code> (new): Boolean flag for writing the star positions to the HDF5 file [default=yes]

* Cosmics

* <code>Sky</code> section in the input file with the configuration parameters for the sky background and cosmics

* Added scripts for comparison with PIS

* Time dependency for:

	- PSF sigma (analytic non-Gaussian PSF)
	- focal length 
	- throughput maps



### Fixed

* Bug in jitter from file that caused negative heartbeat intervals in some specific cases

* Kernel dimension restrictions (GitHub issue #211)

* Number of cosmics too high (GitHub issue #206)

* Bug in <code>DetectorWithAnalyticNonGaussianPSF:addFlux()</code> that caused the CCD orientation angle to be ignored when the analytic non-Gaussian PSF was chosen

* Bug in <code>createStarCatalogFileFromPixelCoordinates()</code> so that it now also works when in the input yaml file the telescope group ID and/or CCD position is not "Custom"

* ThroughputMap reset? (GitHub issue # 202)



### Changed

* Project number in the documentation (3.2 -> 3.3)

* Removed relative paths in tests (required for automatic testing in Jenkins)

* Reading in exposure time as double (instead of integer)

* <code>JitterFromFile</code> and <code>ThermoElasticDriftFromFile</code> now only read in the relevant parts of the files

* Improved comments in the default input file

* If <code>UseJitter == no</code>, then the jitter file is no longer read, even when <code>UseJitterFromFile == yes</code> (idem for drift)

* Parameter review -> update of configuration parameter values

* Update of the field distortion polynomial

* Update of the documentation pages:

	- re-structuring
	- added information about new configuration parameters and (optional) input files
	- improved the description of configuration parameters and procedures

* Updated <code>CMakeLists.txt</code> to use C++14 rather than C++11

* Incorporated <code>StarCatalog</code> in <code>Sky</code> so that the former becomes obsolete

* Extended <code>Parameter<T,N></code>  such that it can accommodate arrays of scalars (such as the distortion coefficients)



### Removed

* Unnecessary log messages





<!-- 3.2.1 -->
<!-- ***** -->

## 06/04/2018: 3.2.1

### Added

* Clarification of the <code>SubPixel</code> parameter in the configuration file (GitHub issue #175)

* Added parameters for camera groups and pre-defined CCD positions



### Fixed

* Corrected relative paths in tests (for automatic builds + testing)

* Correction of the application of the FEE and CCD gain + values in the configuration files + documentation

* Throughput made dependent of CCD position (GitHub issue # 193)

* Correction of the angle dependency of QE and polarisation

* Corrected orientation angles for CCD 2 and 4



### Changed

* Updated version number in the documentation

* Using new gain and readout noise in the calculation of the photometry

* CCD code A, B, C, and D were replaced by 3, 2, 4, and 1





<!-- 3.2.0 -->
<!-- ***** -->

## 05/07/2017: 3.2.0

### Fixed

* Flux values calculated in Camera (expressed in photons) are now rounded instead of floored



### Changed

* Flatfield map only generated if <code>IncludeFlatfield == "yes"</code>





<!-- 3.2.0 RC2 -->
<!-- ********* -->

## 16/06/2017: 3.2.0-RC2

### Changed

* Set the <code>CCD/Position</code> to <code>Custom</code> such that the default settings are used for backward compatibility with the previous release of PlatoSim3

* Updated documentation with respect to

	* Reference frames
	* Description of Camera groups and pre-defined CCDs
	
* Updates to tutorials to bring them in-line with the changes in this release





<!-- 3.2.0 RC1 -->
<!-- ********* -->

## 15/05/2017: 3.2.0-RC1

### Added

* More detailed throughput specifications (following PLATO-DLR-PL-RP-001): in addition to vignetting:

	- particulate and molecular contamination,
	- angle-dependent quantum efficiency,
	- and angle-dependent polarisation
	
* Readout noise: contribution of the FEE and CCD added in quadrature

* Gain: contribution from the CCD and the FEE (different for both detector halves and both ADCs)

* Temperature dependency:

	- Implemented for FEE and CCD gain, and for electronic offset
	- Either fixed at the nominal operating temperature of the component or read from a file; similar to jitter and drift

* Quantisation (i.e. combined effect of (1) gain, (2) electronic offset, (3) rounding pixel values, and (4) digital saturation) can be switched on/off

* Performance optimisation:

	- Split time series into chunks
	- Distribution over nodes using Slurm
	- Angle-dependent analytic Gaussian PSF
	- Angle-dependent analytic non-Gaussian PSF
	
* Kinematic aberration:

	- Differential & absolute
	- Baseline: circular Earth orbit

* Ageing: linear degradation implemented for the transmission efficiency

* New dependency: Faddeeva (used by the analytical PSFs) 

* Visualisation of the output with <code>h5ls</code> and <code>h5get</code> (Python functions)

* Easy selection of requested camera group from the input file, use <code>Telescope/GroupID</code> = [1,2,3,4,Fast,Custom]

* Easy selection of requested CCD position, use <code>CCD/Position</code> = [1,2,3,4,Custom]



### Changed

* Updated reference frames, following PLATO-DLR-PL-TN-016 and PLATO-OHB-PL-LI-009 (see PLATO-KUL-PL-TN-001), in particular: sunshield pointing towards the Sun

* Electronic offset: moved from CCD to FEE

* Updated input files, incorporating the changes in configuration parameters

* Updated documentation, describing the new configuration parameters

* Updated tutorials, incorporating the changes in configuration parameters
