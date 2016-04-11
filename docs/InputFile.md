# Description of the input file

To configure the Plato Simulator, a large set of input parameters is required.  The input file format use for PlatoSim3 is YAML (see https://learnxinyminutes.com/docs/yaml/). We use only a very limited set of the YAML functionality, enough to allow us to provide input files for different parts of the simulator. 

Any desired simulation can be obtained by modifying the following input:
	* @ref configurationParameters "configuration parameters" (in the YAML file):
		- @ref generalParameters "general parameters"
		- @ref observingParameters "observing parameters"
		- @ref platformParameters "platform parameters"
		- @ref telescopeParameters "telescope parameters"
		- @ref cameraParameters "camera parameters"
		- @ref psfParameters "PSF parameters"
		- @ref ccdParameters "CCD parameters"
		- @ref subFieldParameters "sub-field parameters"
		- @ref seedParameters "seed parameters"
	* file comprising a @ref starCatalogue "star catalogue" of the region of the sky of interest
	* optional file comprising @ref psfFile "pre-computed PSFs"
	* @ref jitterFile "jitter" file (only required when the jitter option has been enabled in the configuration file)
 
In the following sections we describe these parameters for the simulations in detail.
 
@section configurationParameters Configuration Parameters 

The configuration parameters for the simulation are stored in a YAML file, e.g. inputfile.yaml, in the /inputfiles directory. This section describes the parameters in the different blocks of the configuration file. These blocks reflect their function in the simulation.

@subsection generalParameters General Parameters

The general configuration parameters a listed in the <b>General</b> block of the configuration file.  The structure of this block is the following:

\code{.yaml}
General:
   
    ProjectLocation:             ENV['PLATO_PROJECT_HOME']
\endcode


@subsubsection projectLocation ProjectLocation

<i>Allowed values:</i> name of an existing directory on disk or environment variable in the format <dfn>ENV['PLATO_PROJECT_HOME']</dfn>

Full path of the directory in which you have checked out the PlatoSim3 project, or an environment variable, e.g. PLATO_PROJECT_HOME, containing the full path to that directory.  In the latter case, you must make sure you have exported this variable before initiating a simulation:

\code{.unparsed}
 export PLATO_PROJECT_HOME = <full path to the PlatoSim3 directory>
\endcode

@subsection observingParameters Observing Parameters

The <b>ObservingParameters</b> block of the configuration file contains the configuration parameters that are specific to the simulated observation and are not specific for the hardware components of the satellite.  The structure of this block is the following:

\code{.yaml}
ObservingParameters:

	NumExposures:                40              
    ExposureTime:                23              
    RApointing:                  180              
    DecPointing:                 -70             
    Fluxm0:                      1.00238e8       
    SkyBackground:               220.            
    StarCatalogFile:             inputfiles/starcatalog.txt
\endcode

@subsubsection numExposures NumExposures

<i>Allowed values:</i> > 0

Number of exposures to generate for the simulation.

@subsubsection exposureTime ExposureTime

<i>Allowed values:</i> > 0

Integration time of one exposure, expressed in seconds. Note that the total integration time is the sum of the exposure and the readout time:

	\f[ t_{integration} = t_{exposure} + t_{readout}.\f]

@subsubsection  raPointing RApointing

<i>Allowed values:</i>  ∈ [0, 360]

Right ascension of the pointing, expressed in degrees.

@subsubsection  decPointing DecPointing

<i>Allowed values:</i> ∈ [-90, 90]

Declination of the pointing, expressed in degrees.

@subsubsection fluxm0 Fluxm0

<i>Allowed values:</i> > 0

Flux of a star of zero magnitude (\f$ m_{\lambda} = 0 \f$), expressed in photons \f$ \cdot \f$  s<sup>-1</sup> \f$  \cdot \f$  cm<sup>-2</sup> in the passband of the magnitudes listed in the star catalogue.

For an exposure of \f$t_{exp}\f$ seconds, the measured flux \f$F_{phot}\f$ of a star is computed from its catalogue magnitude \f$m_{\lambda}\f$, the effective light-collecting area \f$A\f$ (in cm<sup>2</sup>) of the telescope, the transmission efficiency \f$T_{\lambda}\f$ of the optical system, the quantum efficiency \f$Q\f$ of the detector, and the flux per second \f$F_0\f$ of a star with zero magnitude (\f$m_{\lambda} = 0\f$) from the equation

\f[
F_{phot} = t_{exp} \cdot F_0 \cdot T_{\lambda} \cdot Q \cdot A \cdot 10^{-0.4} m_{\lambda}
\f]

where the \f$\lambda\f$ subscript refers to the wavelength range in which the simulation is performed.

@subsubsection skyBackground SkyBackground

<i>Allowed values:</i> < 0 for automatic calculation, ≥ 0 to use the input value

In case a positive value is given, the sky background (zodiacal + galactic), is set to the given value, expressed in photons \f$ \cdot \f$ s<sup>-1</sup> \f$ \cdot \f$ pixel<sup>-1</sup>.

In case a negative value is given, the sky background is computed automatically from tabular values, interpolated to the central coordinates of the sub-field. A constant sky background is assumed for the whole sub-field. Note that for some regions in the sky the automatic computation of the sky background may fail due to the lack of tabulated values. In that case you can set the sky background manually.

@subsubsection  starCatalogFile StarCatalogFile

Full path to the @ref starCatalogFile "star catalogue file".


@subsection platformParameters Platform Parameters

The <b>Platform</b> block of the configuration file contains all the information that is specific to the platform of the satellite.  The structure of this block is the following:

\code{.yaml}
Platform:

    UseJitter:                   yes             
    UseJitterFromFile:           no              
    JitterYawRms:                2.3             
    JitterPitchRms:              2.3             
    JitterRollRms:               2.3             
    JitterTimeScale:             3600.           
    JitterFileName:              inputfiles/jitter.txt
\endcode

@subsubsection useJitter UseJitter

<i>Allowed values:</i> "yes" and "no"

Indicates whether pointing variations should be taken into account.

@subsubsection useJitterFromFile UseJitterFromFile

<i>Allowed values:</i> "yes" and "no"

Indicates whether the jitter time series must be read from a jitter file ("yes") or the jitter positions must be generated from the jitter parameters ("no").

@subsubsection  jitterYawRms JitterYawRms

<i>Allowed values:</i> ≥ 0, only required if the jitter positions must be generated from jitter parameters (@ref useJitterFromFile = 0)

Standard deviation (expressed in arcsec / s) of the normal distribution (with zero mean) describing the yaw value from one jitter position to the next one.

@subsubsection jitterPitchRms JitterPitchRms

<i>Allowed values:</i> ≥ 0, only required if the jitter positions must be generated from jitter parameters (@ref useJitterFromFile = 0)

Standard deviation (expressed in arcsec / s) of the normal distribution (with zero mean) describing the pitch value from one jitter position to the next one.

@subsubsection jitterRollRms JitterRollRms

<i>Allowed values:</i> ≥ 0, only required if the jitter positions must be generated from jitter parameters (@ref useJitterFromFile = 0)

Standard deviation (expressed in arcsec / s) of the normal distribution (with zero mean) describing the roll value from one jitter position to the next one.

@subsubsection jitterTimeScale JitterTimeScale

<i>Allowed values:</i> > 0

Timescale of the jitter, expressed in seconds.

@subsubsection jitterFileName JiterFileName

Full path of the jitter file. This is only required if the jitter positions must be read from a file (@ref useJitterFromFile = 1).

@subsection telescopeParameters Telescope Parameters

The <b>Telescope</b> block of the configuration file contains all the information that is specific to the telescope.  The structure of this block is the following:


\code{.yaml}
Telescope:
    
    LightCollectingArea:         113.1         
    TransmissionEfficiency:      0.757         
    FOVSquareDegrees:            1072.0        
    DriftYawRms:                 2.3           
    DriftPitchRms:               2.3           
    DriftRollRms:                2.3           
    DriftTimeScale:              3600.         
\endcode

@subsubsection lightCollectingArea LightCollectingArea

<i>Allowed values:</i> > 0

Light-collecting area of one telescope, expressed in cm<sup>2</sup>.

@subsubsection transmissionEfficiency TransmissionEfficiency

<i>Allowed values:</i> ∈ [0,1]

Tranmission efficiency of the optical system, considering the passband and spectral energy distribution of the stars, given the Fluxm0 parameter and the magnitudes in the star catalogue.

@subsubsection  fovSquareDegrees FOVSquareDegrees

<i>Allowed values:</i> > 0

Sky area covered by the FOV of one telescope, expressed in degrees<sup>2</sup>.

@subsubsection driftYawRms DriftYawRms

<i>Allowed values:</i> ≥ 0

Standard deviation (expressed in arcsec /  min) of the normal distribution (with zero mean) describing the yaw value from one thermo-elastic drift position to the next one.

@subsubsection driftPitchRms DriftPitchRms

<i>Allowed values:</i> ≥ 0

Standard deviation (expressed in arcsec / min) of the normal distribution (with zero mean) describing the pitch value from one thermo-elastic drift position to the next one.

@subsubsection driftRollRms DriftRollRms

<i>Allowed values:</i> ≥ 0

Standard deviation (expressed in arcsec / min) of the normal distribution (with zero mean) describing the roll value from one thermo-elastic drift position to the next one.

@subsubsection driftTimeScale DriftTimeScale

<i>Allowed values:</i> > 0

Timescale of the thermo-elastic drift, expressed in seconds.

@subsection cameraParameters Camera Parameters

The <b>Camera</b> block of the configuration file contains all the information that is specific to the camera.  The structure of this block is the following:

\code{.yaml}
Camera:
    
    FocalPlaneOrientation:       0.0             
    PlateScale:                  0.8333          
    FocalLength:                 0.24712595      
    ThroughputBandwidth:         400             
    ThroughputLambdaC:           600             
    IncludeFieldDistortion:      yes             
    FieldDistortion:                             
        Type:                    Polynomial1D
        Degree:                  3
        Coefficients:            [-0.0036696919678, 1.0008542317, -4.12553764817e-05, 5.7201219949e-06]
        InverseCoefficients:     [-0.00458067036444, 1.00110311283, -5.61136295937e-05, -4.311925329e-06]
\endcode

@subsubsection focalPlaneOrientation FocalPlaneOrientation

<i>Allowed values:</i> any

Orientation angle of the focal plane, expressed in degrees. For an angle of 0°, the y-axis of the CCD (with an orientation angle of 0°) points towards the North. A positive angle corresponds to a counterclockwise rotation. Have a look at Fig. 1 for more details.

@subsubsection plateScale PlateScale

<i>Allowed values:</i> > 0

Nominal plate scale in arcsec / micron. This value affects the visible FOV of the CCD.

@subsubsection focalLength FocalLength

<i>Allowed values:</i> > 0

Focal length as recovered from the Zemax model, expresse in m.

@subsubsection throughputBandwidth ThroughputBandwidth

<i>Allowed values:</i> > 0

FWHM of the throughput passband, expressed in nm.

@subsubsection throughputLambdaC ThroughputLambdaC

<i>Allowed values:</i> > 0

Central wavelength of the throughput passband, expressed in nm.

@subsubsection includeFieldDistortion  IncludeFieldDistortion

<i>Allowed values:</i> "yes" and "no"

Indicates whether or not the field distortion must be taken into account.

@subsubsection fieldDistortionType FieldDistortion: Type

<i>Allowed values:</i> "Polynomial1D" or "Polynomial2D"

Indicates that the field distortion is calculated by means of either a 1D or a 2D polynomial.

A 1D polynomial of degrees \f$n\f$ can be written as 

\f[P(x) = c_{0} + c_{1} \cdot x + c_{2} \cdot x^{2} + ... + c_{n} \cdot x^{n}\f]

and a d polynomial as \f[P(x, y) = c_{00} + c_{10} \cdot x + ... + c_{n0} \cdot x^{n} + c_{01} \cdot y      + ... + c_{0n} \cdot y^{n} + c_{11} \cdot x \cdot y + c_{12} \cdot x \cdot y^{2}      + ... + c_{1(n - 1)} \cdot x \cdot y^{n - 1} + ... + c_{(n - 1)1} \cdot x^{n - 1} \cdot y \f]

@subsubsection fieldDistortionDegree FieldDistortion: Degree

<i>Allowed values:</i> > 0; for the 2D polynominal: ≤ 3

Degree \f$n\f$ of the polynomial describing the field distortion.

@subsubsection fieldDistortionCoefficients FieldDistortion: Coefficients

Coefficients of the polynomial describing the field distortion.  For a 1D polynomial the coefficients are specified as \f$ [c_0, c_1,..., c_n] \f$, whilst for a 2D polynomial as \f$ [[c_{00}, c_{01},..., c_{0n}], [c_{10}, c_{11},..., c_{1n}],..., [c_{n0}, c_{n1},..., c_{nn}]] \f$.

@subsubsection fieldDistortionInverseCoefficients FieldDistortion: InverseCoefficients

Coefficients for inverse polynomial of the polynomial describing the field distortion.

@subsection psfParameters PSF Parameters

The <b>PSF</b> block of the configuration file contains all the information that is specific to the PSF.  The structure of this block is the following:


\code{.yaml}
PSF:

    Model:                       Gaussian 
    Gaussian:                             
      Sigma:                     0.25     
      NumberOfPixels:            8        
    FromFile:                             
      Filename:                  inputfiles/psf.hdf5 
      DistanceToOA:              10       
      RotationAngle:             45         
      NumberOfPixels:            8          
\endcode

@subsubsection  psfModel Model
<i>Allowed values:</i> "Gaussian" and "FromFile"

Indicates whether to use a Gaussian PSF or to read the PSF from an HDF5 file.

@subsubsection  gaussSigma Gaussian: Sigma

<i>Allowed values:</i> > 0, only required if a Gaussian PSF must be used (@ref psfModel = Gaussian).

Width (σ) of the two-dimensional Gaussian PSF, expressed in pixels.

@subsubsection  gaussNumPixels Gaussian: NumberOfPixels

<i>Allowed values:</i> > 0, only required if a Gaussian PSF must be used (@ref psfModel = Gaussian).

Number of pixels (in both directions) for which the Gaussian PSF must be generated.

@subsubsection  psfFilename FromFile: Filename

<i>Allowed values:</i> only required if a pre-computed PSF must be used (@ref psfModel = FromFile).

Full path to the file holding the location independent @ref psfFile "pre-computed PSF".

@subsubsection  psfDistance FromFile: DistanceToOA

<i>Allowed values:</i> -1 for automatic calculation, ≥ 0 to use the input value; only required if a pre-computed PSF must be used (@ref psfModel = FromFile).

In case a positive value is given the input value will be used for the angular distance to the optical axis.

In case a negative value is given, the angular distance to the optical axis will be calculated automatically.

@subsubsection psfRotation FromFile: RotationAngle

<i>Allowed values:</i> Any, only required if a pre-computed PSF must be used (@ref psfModel = FromFile).

Arbitrary rotation angle of the PSF, expressed in degrees and measured counterclockwise.

@subsubsection  psfNumPixels FromFile: NumberOfPixels

<i>Allowed values:</i> > 0, only required if a pre-computed PSF must be used (@ref psfModel = FromFile).

Number of pixels (in both directions) for which the PSF was generated.

@subsection ccdParameters CCD Parameters

The <b>CCD</b> block of the configuration file contains all the information that is specific to the CCD.  The structure of this block is the following:

\code{.yaml}
CCD:

    OriginOffsetX:               0         
    OriginOffsetY:               0         
    Orientation:                 0         
    NumColumns:                  4510      
    NumRows:                     4510      
    PixelSize:                   18        
    Gain:                        16             
    QuantumEfficiency:           0.8745         
    FullWellSaturation:          1000000        
    DigitalSaturation:           65535          
    ReadoutNoise:                28             
    ElectronicOffset:            100            
    ReadoutTime:                 2              
    FlatfieldPtPNoise:           0.016          
    CTEMean:                     0.99999        
    IncludeFlatfield:            no             
    IncludePhotonNoise:          yes            
    IncludeReadoutNoise:         yes            
    IncludeCTIeffects:           yes            
    IncludeOpenShutterSmearing:  yes            
    IncludeVignetting:           yes             
    IncludeConvolution:          yes            
    IncludeFullWellSaturation:   yes            
    IncludeDigitalSaturation:    yes            
    WriteSubPixelImagesToHDF5:   no              
\endcode

@subsubsection originOffsetX OriginOffsetX

<i>Allowed values:</i> Any

Offset of the CCD origin from the centre of the optical plane (i.e. the intersection of the optical axis with the focal plane) in the x-direction, expressed in mm. The origin of the CCD is defined as the point where the readout register is located. See Fig. 1 for more details.

@subsubsection originOffsetY OriginOffsetY

<i>Allowed values:</i> Any

Offset of the CCD origin from the centre of the optical plane (i.e. the intersection of the optical axis with the focal plane) in the y-direction, expressed in mm. The origin of the CCD is defined as the point where the readout register is located. See Fig. 1 for more details.

@subsubsection ccdOrientation Orientation

<i>Allowed values:</i> Any

Orientation angle of the CCD w.r.t. the orientation of the focal plane, measured counterclockwise and expressed in degrees. This rotation is performed around the offset origin of the CCD. See Fig. 1 for more details.

@subsubsection ccdNumColumns NumColumns

<i>Allowed values:</i> > 0

Number of pixels of the CCD in the x-direction (i.e. number of columns).

@subsubsection ccdNumRows NumRows

<i>Allowed values:</i> > 0

Number of pixels of the CCD in the y-direction (i.e. number of rows).

@subsubsection pixelSize PixelSize
<i>Allowed values:</i> > 0

Nominal pixel size, expressed in micron.
        
@subsubsection Gain             
<i>Allowed values:</i> > 0

CCD gain, expressed in e<sup>-</sup> / ADU and assumed to be constant throughout a simulation. This parameter relates the number of electrons per pixel to the number of counts (i.e. ADU) per pixel.


@subsubsection QuantumEfficiency   
<i>Allowed values:</i> ∈ [0,1]

Quantum efficiency of the detector, considering the passband and the spectral energy distribution of the stars given the @ref fluxm0 parameter and the magnitude of the stars in the @ref starCatalogue. This is the ratio of the number of collected electrons to the number of incident photons.

@subsubsection fullWellSaturation FullWellSaturation
<i>Allowed values:</i> > 0
     
Full-well saturation limit of a single CCD pixel, expressed in e<sup>-</sup> / pixel. If a pixels receives more electrons than its full-well saturation limit, the additional electrons flow evenly distributed in positive and negative charge-transfer direction, a phenomenon called <i>blooming</i>. The electrons reaching the edge of the CCD will not be detected.
  
@subsubsection digitalSaturation DigitalSaturation
<i>Allowed values:</i> > 0

Digital saturation limit of the CCD to which pixel values are topped off, expressed in ADU / pixel. This value depends on the A/D convertor of the detector. For a 16-bit convertor, the digital saturation limit is 65536 ADU.

The @ref gain of the detector should be such that the full-well saturation results in values below the digital saturation limit.
     
@subsubsection readoutNoise ReadoutNoise   
<i>Allowed values:</i> ≥ 0

Mean readout noise of the detector, expressed in e<sup>-</sup>.

Readout noise occurs due to the imperfect nature of the CCD amplifiers. When the electrons are transferred to the amplifier, the induced voltage is measured. However, this measurement is not perfect, but gives a value which is on average too high by an amount of the readout noise, with the squareroot of the readout noise as standard deviation.
       
@subsubsection electronicOffset ElectronicOffset
<i>Allowed values:</i> ≥ 0

Electronic offset or bias level, expressed in ADU, that is added to the digital signal in order to avoid negative readout values. The electronic offset can be measured in a pre-scan strip, which essentially consists of a few additional rows of the CCD. These rows only contain the electronic offset and the readout noise. This pre-scan strip consisting of @ref piasPreScanRows rows will be stored in the output file.
        
@subsubsection readoutTime ReadoutTime

<i>Allowed values:</i> ≥ 0

Time required to read out an entire CCD working in frame transfer mode, expressed in seconds. Because of the absence of a shutter (which is common in space-based instruments), the CCD still receives light during frame transfer. The flux of each sub-pixel is affected by the flux of the sub-pixels in the same column. Because the CCD is exposed during the whole readout and multiple exposures are created, also the sub-pixels further away from the readout register have their influence.

For non-frame-transfer CCDs the readout time should be set to zero.

@subsubsection flatfieldPtPNoise FlatfieldPtPNoise
<i>Allowed values:</i> ∈ [0,1]

Fractional peak-to-peak amplitude of the pixel non-uniform sensitivity response.
  
@subsubsection CTEMean cteMean
Allowed values: ∈ [0,1]

Mean charge-transfer efficiency (CTE) of the detector.

Because of detector defects, electrons can get trapped in the readout process. The trapped charge ends up getting dissociated from its original pixel and eventually gets released into another pixel. The result is that the original image gets smeared out in the direction away from the readout amplifier (visible in the appearance of "charge trails"). This is known as imperfect CTE (Charge-Transfer Efficiency) or alternatively as CTI (Charge-Transfer Inefficiency). The fraction of the charge that is successfully transferred from one row to the next row is expressed by this parameter.

The charge trails impact photometry, noise, and astrometry of sources. CTI removes flux from the central pixel and thus degrades the expected S/N for an observation. CTI trails bias measurements of source along the trail direction, which can severely impact high-precision astrometry.

@subsubsection inclFlatfield IncludeFlatfield
<i>Allowed values:</i> "yes" and "no"

Indicates whether or not to include the flatfield.

@subsubsection inclPhotonNoise IncludePhotonNoise
<i>Allowed values:</i> "yes" and "no"

Indicates whether or not to include photon noise.

@subsubsection inclReadoutNoise IncludeReadoutNoise
<i>Allowed values:</i> "yes" and "no"

Indicates whether or not to include readout noise.

@subsubsection inclCTI IncludeCTIeffects
<i>Allowed values:</i> "yes" and "no"

Indicates whether or not to include CTI effects.

@subsubsection inclOpenShutterSmearing IncludeOpenShutterSmearing
<i>Allowed values:</i> "yes" and "no"

Indicates whether or not to include open-shutter smearing effects.

@subsubsection inclVignetting IncludeVignetting          
<i>Allowed values:</i> "yes" and "no"

Indicates whether or not to include brightness attenuation towards the edge of the FOV due to vignetting.
        
@subsubsection inclConvolution IncludeConvolution
<i>Allowed values:</i> "yes" and "no"

Indicates whether or not the sub-pixel map must be convolved with the PSF.

@subsubsection inclFullWellSaturation IncludeFullWellSaturation
<i>Allowed values:</i> "yes" and "no"

Indicates whether or not to apply full-well saturation.

@subsubsection inclDigitalSaturation IncludeDigitalSaturation
<i>Allowed values:</i> "yes" and "no"

Indicates whether or not to apply digital saturation.

@subsubsection writeSubPixelImages WriteSubPixelImagesToHDF5
<i>Allowed values:</i> "yes" and "no"

Indicates whether or not the sub-pixel images must be written to the HDF5-file.  Use this for a limited number of exposures, as it takes a lot of space.

@subsection subFieldParameters Sub-Field Parameters

The <b>SubField</b> block of the configuration file contains all the information about the sub-field of the CCD that is modelled by the simulation.  The structure of this block is the following:

\code{.yaml}
SubField:    

    ZeroPointRow:                0
    ZeroPointColumn:             0
    NumColumns:                  10
    NumRows:                     10
    NumBiasPrescanRows:          5
    NumSmearingOverscanRows:     5
    SubPixels:                   4
\endcode

@subsubsection zeroPointRow ZeroPointRow
<i>Allowed values:</i> > 0

Row of the origin of the sub-field in the detector, expressed in pixels.

@subsubsection zeroPointColumn ZeroPointColumn
<i>Allowed values:</i> > 0

Column of the origin of the sub-field in the detector, expressed in pixels.

@subsubsection  numColumns NumColumns
<i>Allowed values:</i> ≥ 8

Number of columns in the sub-field, expressed in pixels.

@subsubsection numRows NumRows
<i>Allowed values:</i> ≥ 8

Number of rows in the sub-field, expressed in pixels.

@subsubsection  numPreScanRows NumBiasPrescanRows
<i>Allowed values:</i> ≥ 0

Number of rows in the pre-scan strip (to determine the bias), expressed in pixels.

@subsubsection numOverScanRows NumSmearingOverscanRows
<i>Allowed values:</i> ≥ 0

Number of rows in the over-scan strip (to determine the smearing), expressed in pixels.

@ubsubsection numSubPixels SubPixels
<i>Allowed values:</i> power of 2 (≤ 128)

Number of sub-pixels per pixel in both directions.

@subsection seedParameters Seed Parameters

The <b>RandomSeeds</b> block of the configuration file contains all the seeds for random-number generation in the simulator.  The structure of this block is the following:

\code{.yaml}
    ReadOutNoiseSeed:            1424949740 
    PhotonNoiseSeed:             1433320336 
    JitterSeed:                  1433320381 
    FlatFieldSeed:               1425284070 
    CTESeed:                     1424949740 
    DriftSeed:                   1433429158 
\endcode

@subsubsection ReadOutNoiseSeed
<i>Allowed values:</i> > 0

Seed for the random-number generator used for the readout noise.

@subsubsection PhotonNoiseSeed
<i>Allowed values:</i> > 0

Seed for the random-number generator used for the photon noise.

@subsubsection JitterSeed
<i>Allowed values:</i> > 0

Seed for the random-number generator used for the jitter.

@subsubsection FlatFieldSeed
<i>Allowed values:</i> > 0

Seed for the random-number generator used for the flatfield.

@subsubsection CTESeed
<i>Allowed values:</i> > 0

Seed for the random-number generator used for the CTE.

@subsubsection DriftSeed
<i>Allowed values:</i> > 0

Seed for the random number generator used for the drift.

@section starCatalogue Star Catalogue

A star catalogue must be provided in a file in ASCII format. This file should contain three columns, separated by a space, holding the following information:
* right ascension of the stars [degrees]
* declination of the stars [degrees]
* stellar magnitude

The path of this file, relative to the @ref projectLocation "project location", must be provided via the @ref starCatalogFile parameter in the configuration file (under observing parameters). 

@section psfFile PSF File (Optional)

Unless you indicate you want to generate a Gaussian PSF, pre-computed normalised PSFs must be provided in the form of an HDF5-file (@ref useGauss = 0). The path of this file, relative to the @ref projectLocation "project location", is specified via the @ref psfFilename parameter.

The simulator will automatically select the PSF for which the angular distance to the optical axis matches best for the simulated sub-field.

@section jitterFile Jitter File (Optional)


