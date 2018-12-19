# Changelog for PlatoSim



<!-- 3.3.6 -->
<!-- ***** -->

## ??/??/2018: 3.3.6

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