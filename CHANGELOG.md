# Changelog for PlatoSim

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





## 05/07/2017: 3.2.0

### Fixed

* Flux values calculated in Camera (expressed in photons) are now rounded instead of floored

### Changed

* Flatfield map only generated if <code>IncludeFlatfield == "yes"</code>





## 16/06/2017: 3.2.0-RC2

### Changed

* Set the <code>CCD/Position</code> to <code>Custom</code> such that the default settings are used for backward compatibility with the previous release of PlatoSim3

* Updated documentation with respect to

	* Reference frames
	* Description of Camera groups and pre-defined CCDs
	
* Updates to tutorials to bring them in-line with the changes in this release





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