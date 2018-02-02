# Description of the Input File {#InputFileDescription}

To configure the Plato Simulator, a large set of input parameters is required.  The input file format use for PlatoSim3 is <a href="https://learnxinyminutes.com/docs/yaml/">YAML</a>. We use only a very limited set of the YAML functionality, enough to allow us to provide input files for different parts of the simulator. 

Any desired simulation can be obtained by modifying the following input:
	*  [configuration parameters](#configurationParameters) (in the YAML file):
		- [general parameters](#generalParameters)
		- [observing parameters](#observingParameters)
		- [platform parameters](#platformParameters)
		- [telescope parameters](#telescopeParameters)
		- [camera parameters](#cameraParameters)
		- [PSF parameters](#psfParameters)
		- [FEE parameters](#feeParameters)
		- [CCD parameters](#ccdParameters)
		- [sub-field parameters](#subFieldParameters)
		- [seed parameters](#seedParameters)
	* file comprising a [star catalogue](#starCatalogue) of the region of the sky of interest
	* optional file comprising [pre-computed PSFs](#psfFile)
	* [jitter](#jitterFile) file (only required when the jitter option has been enabled in the configuration file)

Additionally, there are two blocks that hold pre-defined settings (which you should NOT alter):
 	* [camera group 1, 2, 3, and 4, and fast cameras](#cameraGroups)
 	* [CCD 1, 2, 3, and 4](#ccdPositions)
 
In the following sections we describe these parameters for the simulations in detail.

For more details on the reference frames we are using (for the spacecraft, telescope, focal plane, and CCD), radial dependency of the PSF, and rotation angles for platform jitter and telescope drift, please, have a look at technical note [PLATO-KUL-PL-TN-0001](../technicalnotes/PLATO-KUL-PL-TN-0001.pdf).



<!-- ************************ -->
<!-- Configuration Parameters -->
<!-- ************************ -->

## <a name="configurationParameters"></a>Configuration Parameters

The configuration parameters for the simulation are stored in a YAML file, e.g. <code>inputfile.yaml</code>, in the <code>/inputfiles</code> directory. This section describes the parameters in the different blocks of the configuration file. These blocks reflect their function in the simulation.



<!-- General Parameters -->
<!-- ****************** -->

### <a name="generalParameters"></a>General Parameters

The general configuration parameters a listed in the <b>General</b> block of the configuration file.  The structure of this block is the following:

\code{.yaml}
General:
   
    ProjectLocation:             ENV['PLATO_PROJECT_HOME']
\endcode



#### <a name="projectLocation"></a>ProjectLocation

<i>Allowed values:</i> name of an existing directory on disk or environment variable, in the format <dfn>ENV['PLATO_PROJECT_HOME']</dfn>.

Full path of the directory in which you have checked out the PlatoSim3 project, or an environment variable, e.g. PLATO_PROJECT_HOME, containing the full path to that directory.  In the latter case, you must make sure you have exported this variable before initiating a simulation:

\code{.unparsed}
 export PLATO_PROJECT_HOME=<full path to the PlatoSim3 directory>
\endcode



### <a name="observingParameters"></a>Observing Parameters

The <b>ObservingParameters</b> block of the configuration file contains the configuration parameters that are specific to the simulated observation and are not specific for the hardware components of the satellite.  The structure of this block is the following:

\code{.yaml}
ObservingParameters:

	MissionDuration:             6.0
	BeginExposureNr:             0
	NumExposures:                40              
    ExposureTime:                23              
    RApointing:                  180              
    DecPointing:                 -70             
    Fluxm0:                      1.00238e8       
    SkyBackground:               220.            
    StarCatalogFile:             inputfiles/starcatalog.txt
\endcode




#### <a name="missionDuration"></a>MissionDuration
<i>Allowed values:</i> > 0

Total duration of the mission (from BOL til EOL), expressed in years.  This will be used to model parameter degradation over time.



#### <a>beginExposureNr</a>BeginExposureNr
<i>Allowed values:</i> ≥ 0

Sequential number of the first exposure. Useful for <a href="https://en.wikipedia.org/wiki/Slurm_Workload_Manager">Slurm</a> parallellisation.  In that case, long simulations (i.e. with a large number of exposures) will be chopped up into smaller simulations, covering [a small number of exposures](#NumExposures) (see Fig. 1).

@image html /images/chopUpSimulation.png "Figure 1: Long simulations will be chopped up into smaller simulations that can be executed in parallel."


#### <a name="numExposures"></a>NumExposures
<i>Allowed values:</i> > 0

Number of exposures to generate in the simulation.



#### <a name="exposureTime"></a>ExposureTime
<i>Allowed values:</i> > 0

Integration time of one exposure, expressed in seconds. Note that the total integration time is the sum of the exposure and the readout time:

	\f[ t_{integration} = t_{exposure} + t_{readout}.\f]



#### <a name="raPointing"></a>RApointing
<i>Allowed values:</i>  ∈ [0, 360]

Right ascension of the pointing, expressed in degrees.



#### <a name="decPointing"></a>DecPointing
<i>Allowed values:</i> ∈ [-90, 90]

Declination of the pointing, expressed in degrees.



#### <a name="fluxm0"></a>Fluxm0
<i>Allowed values:</i> > 0

Flux of a star of zero magnitude (\f$ m_{\lambda} = 0 \f$), expressed in photons \f$ \cdot \f$  s<sup>-1</sup> \f$  \cdot \f$  cm<sup>-2</sup> in the passband of the magnitudes that are listed in the [star catalogue](#starCatalogue).

For an exposure of \f$t_{exp}\f$ seconds, the measured flux \f$F_{phot}\f$ of a star, expressed in photons, is computed from its catalogue magnitude \f$m_{\lambda}\f$, the [effective light-collecting area](#lightCollectingArea) \f$A\f$ (in cm<sup>2</sup>) of the telescope, the  [transmission efficiency](#transmissionEfficiency) \f$T_{\lambda}\f$ of the optical system, the [quantum efficiency](#quantumEfficiency) \f$Q\f$ of the detector, and the flux per second \f$F_0\f$ of a star with zero magnitude (\f$m_{\lambda} = 0\f$) from the equation

\f[F_{phot} = t_{exp} \cdot F_0 \cdot T_{\lambda} \cdot Q \cdot A \cdot 10^{-0.4 \cdot m_{\lambda}}\f]

where the \f$\lambda\f$ subscript refers to the wavelength range in which the simulation is performed.



#### <a name="skyBackground"></a>SkyBackground
<i>Allowed values:</i> < 0 for automatic calculation, ≥ 0 to use the input value

In case a positive value is given, the sky background (zodiacal + galactic), is set to the given value, expressed in photons \f$ \cdot \f$ s<sup>-1</sup> \f$ \cdot \f$ pixel<sup>-1</sup>.

In case a negative value is given, the sky background is computed automatically from tabular values, interpolated to the central coordinates of the sub-field. A constant sky background is assumed for the whole sub-field. Note that for some regions in the sky the automatic computation of the sky background may fail due to the lack of tabulated values. In that case you can set the sky background manually.



####<a name="starCatalogFile"></a> StarCatalogFile

Path to the [star catalogue file](#starCatalogFile), relative to the [project location](#projectLocation).





<!-- Platform Parameters -->
<!-- ******************* -->

## <a name="platformParameters"></a>Platform Parameters


The <b>Platform</b> block of the configuration file contains all the information that is specific to the platform of the satellite.  The structure of this block is the following:

\code{.yaml}
Platform:

    UseJitter:                   yes             
    UseJitterFromFile:           no              
    JitterYawRms:                2.3             
    JitterPitchRms:              2.3             
    JitterRollRms:               2.3             
    JitterTimeScale:             3600.           
    JitterFileName:              /inputfiles/jitter.txt
\endcode



#### <a name="useJitter"></a>UseJitter
<i>Allowed values:</i> "yes" and "no"

Indicates whether pointing variations should be taken into account.

The Plato Simulator can also account for pointing variations of the spacecraft, so-called jitter. A time series of pointing displacement, expressed in Euler angles (yaw, pitch, roll), either has to be provided as a jitter file or will be generated based on the given jitter parameters (see further).

To ensure a realistic modelling of the jitter, the [time step of the jitter time series](#jitterTimeScale) must be smaller than the [exposure time](#exposureTime).

The configuration of the jitter axes is depicted below.  The Euler angles that characterise the jitter are defined w.r.t. to the spacecraft coordinate system (see Fig. 2).  The origin of this coordinate system is the geometric centre of the interface between the bottom of the optical bench and the service module.  The positive roll axis z<sub>SC</sub> points towards the operator-given mean payload line-of-sight, given by the equatorial coordinates ([RApointing](#raPointing), [DecPointing](#decPointing)).

The angles are defined such that they increase with a clockwise rotation, when looking along the positive axes. First a roll rotation is done around the z<sub>SC</sub> axis, then a pitch rotation is done around the rotated y<sub>SC</sub> axis, and finally a yaw rotation is done around the twice-rotated x<sub>SC</sub> axis.

@image html /images/jitterConfiguration.png "Figure 2: Configuration of the jitter axes for the Plato Simulator, defined w.r.t. the spacecraft coordinate system (xSC, ySC, z<sub>SC</sub>).  The origin of this coordinate system is the geometric centre of the interface between the bottom of the optical bench and the service module.  The positive zSC axis points towards the operator-given pointing coordinates. The xSC axis points in the direction of the highest point of the sunshield."



#### <a name="useJitterFromFile"></a>UseJitterFromFile
<i>Allowed values:</i> "yes" and "no"

Indicates whether the jitter time series must be read from a jitter file ("yes") or the jitter positions must be generated from the jitter parameters ("no").



#### <a name="jitterYawRms"></a>JitterYawRms
<i>Allowed values:</i> ≥ 0, only required if the jitter positions must be generated from jitter parameters ([useJitterFromFile](#useJitterFromFile) = no)

Standard deviation (expressed in arcsec) of the normal distribution (with zero mean) describing the yaw value from one jitter position to the next one.



#### <a name="jitterPitchRms"></a>JitterPitchRms
<i>Allowed values:</i> ≥ 0, only required if the jitter positions must be generated from jitter parameters ([useJitterFromFile](#useJitterFromFile) = no)

Standard deviation (expressed in arcsec) of the normal distribution (with zero mean) describing the pitch value from one jitter position to the next one.



#### <a name="jitterRollRms"></a>JitterRollRms
<i>Allowed values:</i> ≥ 0, only required if the jitter positions must be generated from jitter parameters ([useJitterFromFile](#useJitterFromFile) = no)

Standard deviation (expressed in arcsec) of the normal distribution (with zero mean) describing the roll value from one jitter position to the next one.



#### <a name="jitterTimeScale"></a>JitterTimeScale
<i>Allowed values:</i> > 0

Timescale of the jitter (i.e. time between two subsequent jitter positions), expressed in seconds.



#### <a name="jitterFileName"></a>JitterFileName

Path of the jitter file, relative to the [project location](#projectLocation). This is only required if the jitter positions must be read from a file ([UseJitterFromFile](#useJitterFromFile) = yes).





<!-- Telescope Parameters -->
<!-- ******************** -->

### <a name="telescopeParameters"></a>Telescope Parameters

The <b>Telescope</b> block of the configuration file contains all the information that is specific to the telescope.  The structure of this block is the following:

\code{.yaml}
Telescope:
    
    GroupID:                     Custom
    AzimuthAngle:                0.0
    TiltAngle:                   0.0
    LightCollectingArea:         113.1         
    TransmissionEfficiency:      
        BOL:                     0.800
        EOL:                     0.757
    UseDrift:                    yes
    UseDriftFromFile:            no      
    DriftYawRms:                 2.3           
    DriftPitchRms:               2.3           
    DriftRollRms:                2.3           
    DriftTimeScale:              3600.
    DriftFileName:               /inputfiles/drift.txt         
\endcode



#### <a name="groupID"></a>GroupID
<i>Allowed values:</i> ∈ [1, 2, 3, 4, Fast, Custom]

The telescope group identifier can be used to select a telescope group. There are four groups that have a tilt angle of 9.2º from the optical axis of the satellite, and one group for the fast camera's which is aligned with the satellite Z-axis. When you specify GroupID=Custom, the [tilt angle](#tiltAngle) and [azimuth angle](#azimuthAngle) below the GroupID in the inputfile are used, otherwise the angles are taken from pre-defined parameters in the [CameraGroups](#cameraGroups) block of the configuration file.

@image html /images/telescopeGroups.png "Figure: Field of View for the different telescope groups"

#### <a name="tiltAngle"></a>TiltAngle
<i>Allowed values:</i> > 0

Tilt angle of the telescope, expressed in degrees. This angle, together with the [azimuth angle](#azimuthAngle), characterises the orientation of the telescope pointing (i.e. telescope optical axis) w.r.t. the spacecraft/platform pointing. 

The tilt angle is the offset between the telescope optical axis and the platform pointing, i.e. the angle between the telescope line-of-sight (positive z<sub>telescope</sub>)-axis and the positive z<sub>PLM</sub>-axis (see Figs. 3 and 4).

This parameter is only used when the [GroupID](#groupID)=Custom.

#### <a name="azimuthAngle"></a>AzimuthAngle
<i>Allowed values:</i> Any

Azimuth angle of the telescope, expressed in degrees. This angle, together with the [tilt angle](#tiltAngle), characterises the orientation of the telescope pointing (i.e. telescope optical axis) w.r.t. the spacecraft/platform pointing. 

The azimuth angle is the position angle of the rotation of the telescope around the positive z<sub>PLM</sub>-axis (see Figs. 3 and 4).

This parameter is only used when the [GroupID](#groupID)=Custom.

@image html /images/tiltAzimuth.png "Figure 3: Tilt and azimuth of a telescope."



#### <a name="lightCollectingArea"></a>LightCollectingArea
<i>Allowed values:</i> > 0

Light-collecting area of one telescope, expressed in cm<sup>2</sup>.



#### <a name="transmissionEfficiencyBOL"></a>TransmissionEfficiency: BOL
<i>Allowed values:</i> ∈ [0,1]

Tranmission efficiency of the optical system, considering the passband and spectral energy distribution of the stars, given the Fluxm0 parameter and the magnitudes in the [star catalogue](#starCatalogue), at the beginning of the mission (beginning-of-life).  This parameter is used to model the (linear) degradation in transmission efficiency over the [mission lifetime](#missionDuration).



#### <a name="transmissionEfficiencyEOL"></a>TransmissionEfficiency: EOL
<i>Allowed values:</i> ∈ [0,1]

Tranmission efficiency of the optical system, considering the passband and spectral energy distribution of the stars, given the Fluxm0 parameter and the magnitudes in the [star catalogue](#starCatalogue), at the end of the mission (end-of-life).  This parameter is used to model the (linear) degradation in transmission efficiency over the [mission lifetime](#missionDuration).



#### <a name="useDrift"></a>UseDrift
<i>Allowed values:</i> "yes" or "no"

Indicates whether the thermo-elastic drift of the telescope (w.r.t. the platform) should be taken into account.

Similar to the [UseJitter](#useJitter) parameter for the platform jitter.

The Plato Simulator can also account for the thermo-elastic drift, of the telescope (w.r.t. the platform). A time series of displacement, expressed in Euler angles (yaw, pitch, roll), either has to be provided as a drift file or will be generated based on the given drift parameters (see further).

The Euler angles (yaw, pitch, roll) are defined as the rotation angles around the z<sub>SC</sub>, y'<sub>SC</sub> and z<sub>telescope</sub> = z<sub>FP</sub> axes (see Fig. 4), such that the anges increase with a clockwise rotation when looking along the positive axes.



@image html /images/TelescopeCoordinateSystem.png "Figure 4: The optical axis zFP can be obtained from the spacecraft/platform pointing axis zSC by first rotating the (xSC, ySC) plane around the pointing axis zSC over the azimuth angle (left-hand side) nad then rotating the resulting zSC' axis over the tilt angle (right-hand side)."


#### <a name="useDriftFromFile"></a>UseDriftFromFile
<i>Allowed values:</i> "yes" or "no"

Indicates whether the thermo-elastic drift of the telescope (w.r.t. the platform) should be taken into account.

Similar to the [UseJitterFromFile](#useJitterFromFile) parameter for the platform jitter.




#### <a name="driftYawRms"></a>DriftYawRms
<i>Allowed values:</i> ≥ 0

Standard deviation (expressed in arcsec) of the normal distribution (with zero mean) describing the yaw value from one thermo-elastic drift position to the next one.

Similar to the [JitterYawRms](#jitterYawRms) parameter for the platform jitter.



#### <a name="driftPitchRms"></a>DriftPitchRms
<i>Allowed values:</i> ≥ 0

Standard deviation (expressed in arcsec) of the normal distribution (with zero mean) describing the pitch value from one thermo-elastic drift position to the next one.

Similar to the [JitterPitchRms](#jitterPitchRms) parameter for the platform jitter.



#### <a name="driftRollRms"></a>DriftRollRms
<i>Allowed values:</i> ≥ 0

Standard deviation (expressed in arcsec) of the normal distribution (with zero mean) describing the roll value from one thermo-elastic drift position to the next one.

Similar to the [JitterRollRms](#jitterRollRms) parameter for the platform jitter.



#### <a name="driftTimeScale"></a>DriftTimeScale
<i>Allowed values:</i> > 0

Timescale of the thermo-elastic drift (i.e. time between two subsequent drift positions), expressed in seconds.


#### <a name="driftFileName"></a>DriftFileName
Path of the drift file, relative to the [project location](#projectLocation). This is only required if the drift positions must be read from a file ([UseDriftFromFile](#useDriftFromFile) = yes).




<!-- Camera Parameters -->
<!-- ***************** -->

## <a name="cameraParameters"></a>Camera Parameters

The <b>Camera</b> block of the configuration file contains all the information that is specific to the camera.  The structure of this block is the following:

\code{.yaml}
Camera:
    
    FocalPlaneOrientation:       0.0             
    PlateScale:                  0.8333          
    FocalLength:                 0.24712595      
    ThroughputBandwidth:         400             
    ThroughputLambdaC:           600        
    IncludeAberrationCorrection: yes     
    AberrationCorrection:
        Type:                    differential
    IncludeFieldDistortion:      yes             
    FieldDistortion:                             
        Type:                    Polynomial1D
        Degree:                  3
        Coefficients:            [-0.0036696919678, 1.0008542317, -4.12553764817e-05, 5.7201219949e-06]
        InverseCoefficients:     [-0.00458067036444, 1.00110311283, -5.61136295937e-05, -4.311925329e-06]
\endcode




#### <a name="focalPlaneOrientation"></a>FocalPlaneOrientation
<i>Allowed values:</i> Any

Orientation angle of the focal plane, expressed in degrees. For an angle of 0°, the y-axis of the CCD (with an orientation angle of 0°) points towards the North. A positive angle corresponds to a counterclockwise rotation. Have a look at Fig. 2 for more details.

@image html /images/FocalPlaneCoordinateSystem.png "Figure 2: A schematic overview of the focal plane with 4 CCDs. The optical axis zFP is the blue dot in the middle of the 4 CCDs and points in the positive direction towards the reader. The jitter roll axis zSC is the purple dot, and also points in the positive direction towards the reader.  The focal plane is rotated by the angle γFP w.r.t. to the North direction. The origin of the CCD in the focal plane is defined by its offset (ΔxCCD, ΔyCCD) in mm from the centre of the focal plane. It is then rotated by the angle γCCD round its origin."








#### <a name="plateScale"></a>PlateScale
<i>Allowed values:</i> > 0

Nominal plate scale in arcsec / micron. This value affects the visible FOV of the CCD.




#### <a name="focalLength"></a>FocalLength
<i>Allowed values:</i> > 0

Focal length as recovered from the Zemax model, expressed in m.




#### <a name="throughputBandwidth"></a>ThroughputBandwidth
<i>Allowed values:</i> > 0

FWHM of the throughput passband, expressed in nm.




#### <a name="throughputLambdaC"></a>ThroughputLambdaC
<i>Allowed values:</i> > 0

Central wavelength of the throughput passband, expressed in nm.




#### <a name="includeFieldDistortion"></a>IncludeFieldDistortion
<i>Allowed values:</i> "yes" and "no"

Indicates whether or not the field distortion must be taken into account.



#### <a name=includeAberrationCorrection></a>IncludeAberrationCorrection
<i>Allowed values:</i> "yes" and "no"

Indicates whether or not to apply the aberration correction to all star positions in the [star catalogue](#starCatalogue).  This calculation is an approximation based on a circular earth orbit around the sun and does not take the Lissajous orbit of the satellite around L2 into account. We do calculate the aberration, however, which takes into account the aberration correction done for the spacecraft pointing.



#### <a name=aberrationCorrectionType></a>IncludeAberrationCorrection: Type
<i>Allowed values:</i> "differential" and "absolute"

Indicates whether to apply either differential or absolute aberration correction (if ([IncludeAbberationCorrection](#includeAberrationCorrection) = yes).



#### <a name="fieldDistortionType"></a>FieldDistortion: Type
<i>Allowed values:</i> "Polynomial1D" or "Polynomial2D"

Indicates that the field distortion is calculated by means of either a 1D or a 2D polynomial.

A 1D polynomial of degrees \f$n\f$ can be written as 

\f[P(x) = c_{0} + c_{1} \cdot x + c_{2} \cdot x^{2} + ... + c_{n} \cdot x^{n}\f]

and a d polynomial as \f[P(x, y) = c_{00} + c_{10} \cdot x + ... + c_{n0} \cdot x^{n} + c_{01} \cdot y      + ... + c_{0n} \cdot y^{n} + c_{11} \cdot x \cdot y + c_{12} \cdot x \cdot y^{2}      + ... + c_{1(n - 1)} \cdot x \cdot y^{n - 1} + ... + c_{(n - 1)1} \cdot x^{n - 1} \cdot y \f]




#### <a name="fieldDistortionDegree"></a>FieldDistortion: Degree
<i>Allowed values:</i> > 0; for the 2D polynominal: ≤ 3

Degree \f$n\f$ of the polynomial describing the field distortion.




#### <a name="fieldDistortionCoefficients"></a>FieldDistortion: Coefficients

Coefficients of the polynomial describing the field distortion.  For a 1D polynomial the coefficients are specified as \f$ [c_0, c_1,..., c_n] \f$, whilst for a 2D polynomial as \f$ [[c_{00}, c_{01},..., c_{0n}], [c_{10}, c_{11},..., c_{1n}],..., [c_{n0}, c_{n1},..., c_{nn}]] \f$.




#### <a name="fieldDistortionInverseCoefficients"></a>FieldDistortion: InverseCoefficients

Coefficients for inverse polynomial of the polynomial describing the field distortion.





<!-- PSF Parameters -->
<!-- ************** -->

### <a name="psfParameters"></a>PSF Parameters

The <b>PSF</b> block of the configuration file contains all the information that is specific to the PSF.  The structure of this block is the following:


\code{.yaml}
PSF:

    Model:                       MappedGaussian 
    MappedGaussian:                             
      Sigma:                     0.25     
      NumberOfPixels:            8        
    MappedFromFile:                             
      Filename:                  inputfiles/psf.hdf5 
      DistanceToOA:              10       
      RotationAngle:             45         
      NumberOfPixels:            8
    AnalyticGaussian:
      Sigma00:                   1.0
      SigmaX18:                  5.0
      SigmaY18:                  2.0
    AnalyticNonGaussian:
      ParameterFileName:         inputfiles/parameters.txt
\endcode




#### <a name="psfModel"></a>Model
<i>Allowed values:</i> "MappedGaussian", "MappedFromFile", "AnalyticGaussian", "AnalyticNonGaussian

Indicates whether to use a Gaussian PSF, to read the PSF from an HDF5 file, or to use an analytical model (Gaussian or non-Gaussian):

- MappedGaussian: the PSF is a circular Gaussian, the size of which does not change over the FOV;
- MappedFromFile: the PSF is selected from an HDF5 file with pre-computed PSFs, based on the angular distance to the optical axis;
- AnalyticGaussian: the PSF is an elongated Gaussian (the symmetry axes being parallel to the x- and y-axis), for which the width and the height are given at the centre of the FOV and at 18 degrees from the optical axis;
- AnalyticNonGaussian: the PSF is an analytical non-Gaussian model, the parameters of which are stored in a separate file.




#### <a name="gaussSigma"></a>MappedGaussian: Sigma
<i>Allowed values:</i> > 0, only required if a Gaussian PSF must be used ([psfModel](#Model) = MappedGaussian).

Width (σ) of the two-dimensional Gaussian PSF, expressed in pixels.  This Gaussian PSF does not vary in size over the FOV.




#### <a name="gaussNumPixels"></a>MappedGaussian: NumberOfPixels
<i>Allowed values:</i> > 0, only required if a Gaussian PSF must be used ([Model](#psfModel) = MappedGaussian).

Number of pixels (in both directions) for which the Gaussian PSF must be generated.




#### <a name="psfFilename"></a>MappedFromFile: Filename
<i>Allowed values:</i> only required if a pre-computed PSF must be used ([psfModel](#Model) = MappedFromFile).

Path to the file, relative to the [project location](#projectLocation), holding the location independent [pre-computed PSF](#psfFile).



#### <a name="psfDistance"></a>MappedFromFile: DistanceToOA
<i>Allowed values:</i> -1 for automatic calculation, ≥ 0 to use the input value; only required if a pre-computed PSF must be used ([Model](#psfModel) = MappedFromFile).

In case a positive value is given the input value will be used for the angular distance to the optical axis.

In case a negative value is given, the angular distance to the optical axis will be calculated automatically.




#### <a name="psfRotation"></a>MappedFromFile: RotationAngle
<i>Allowed values:</i> Any, only required if a pre-computed PSF must be used ([Model](#psfModel) = MappedFromFile).

Arbitrary rotation angle of the PSF, expressed in degrees and measured counterclockwise.




#### <a name="psfNumPixels"></a>MappedFromFile: NumberOfPixels
<i>Allowed values:</i> > 0, only required if a pre-computed PSF must be used ([Model](#psfModel) = MappedFromFile).

Number of pixels (in both directions) for which the PSF was generated.



#### <a name=sigma00></a>AnalyticGaussian: Sigma00
<i>Allowed values:</i> > 0, only required if a pre-computed PSF must be used ([Model](#psfModel) = AnalyticGaussian).

Standard deviation of the analytical Gaussian PSF in the x- and y-direction at the optical axis, expressed in pixels.



#### <a name=sigmaX18></a>AnalyticGaussian: SigmaX18
<i>Allowed values:</i> > 0, only required if a pre-computed PSF must be used ([Model](#psfModel) = AnalyticGaussian).

Standard deviation of the analytical PSF in the x-direction at 18 degrees from the optical axis, expressed in pixels.



#### <a name=sigmaY18></a>AnalyticGaussian: SigmaY18
<i>Allowed values:</i> > 0, only required if a pre-computed PSF must be used ([Model](#psfModel) = AnalyticGaussian).

Standard deviation of the analytical PSF in the y-direction at 18 degrees from the optical axis, expressed in pixels.



#### <a name=analyticPsfFile></a>AnalyticNonGaussian: ParameterFileName
<i>Allowed values:</i> only required if a  ([psfModel](#Model) = AnalyticNonGaussian).

Path to the file, relative to the [project location](#projectLocation), holding the parameters characterising the analytical model.





<!-- FEE Parameters -->

## <a name="feeParameters"></a>FEE Parameters

The <b>FEE</b> block of the configuration file contains all the information that is specific to the front-end electronics (FEE).  The structure of this block is the following:

\code{.yaml}
FEE:

    NominalOperatingTemperature: 210.15
    Temperature:                 Nominal
    TemperatureFileName:         inputfiles/feeTemperature.txt     
    ReadoutNoise:                40.5         
    Gain:          
    		RefValue:            0.0278
    		ThreeSigma:          0.0    
    		Stability:           -100    		
    ElectronicOffset:           
    		RefValue:            100
    		Stability:           18.8875 
\endcode




#### <a name="nominalTempFEE"></a>NominalOperatingTemperature
<i>Allowed values:</i> > 0

Nominal operating temperature of the FEE, expressed in Kelvin.



#### <a name="tempFEE"></a>Temperature
<i>Allowed values:</i>"Nominal" or "FromFile"

Indicates whether the temperature of the FEE should be fixed at the nominal operating temperature or temperature variations should be read from a file.



#### <a name=tempFileFEE></a>TemperatureFileName
Path to the file, relative to the [project location](#projectLocation), holding the location of the file with the temperature variations of the FEE.



#### <a name=readoutNoiseFEE></a>ReadoutNoise

<i>Allowed values:</i> ≥ 0

Mean readout noise of the FEE, expressed in e<sup>-</sup>/pixel.  This is the same for both ADCs.



#### <a name=gainRefValueFEE></a>Gain: RefValue

<i>Allowed values:</i> > 0

Reference value of the gain of the FEE at its [nominal operating temperature](#nominalTempFEE), expressed in ADU/µV.  The actually gain for the FEE will be different for both ADCs.  We generate a normal distribution, centred at [the reference value](#gainRefValueFEE) and with a width characterised by the [ThreeSigma](#gain3SigmaFEE) parameter, and make two random draws from this distribution.  The outcome will act as gain for ADC1 and ADC2 resp.



#### <a name=gain3SigmaFEE></a>Gain: ThreeSigma

<i>Allowed values:</i> ∈ [0,100]

Percentage of the [reference value for the gain](#gainRefValueFEE) that will act as 3σ for the normal distribution, centred around the (reference value)[#gainRefValueFEE], from which to draw the gain for both ADCs.


#### <a name=gainStabilityFEE></a>Gain: Stability

<i>Allowed values:</i> Any

Change in gain (for both ADCs) with temperature deviations from the nominal operating temperature, expressed in ADU/µV/K.



#### <a name=electronicOffset></a>ElectronicOffset: RefValue

<i>Allowed values:</i> ≥ 0

Electronic offset or bias level at the nominal operating temperature of the FEE, expressed in ADU, that is added to the digital signal in order to avoid negative readout values. The electronic offset can be measured in a pre-scan strip, which essentially consists of a few additional rows of the CCD. These rows only contain the electronic offset and the readout noise. This pre-scan strip consisting of [NumPreScanRows](#numPreScanRows) rows will be stored in the output file.  This is the same for both ADCs.



#### <a name=electronicOffsetStability></a>ElectronicOffset: Stability

<i>Allowed values:</i> Any

Change in electronic offset (for both ADCs) with temperature deviations from the nominal operating temperature, expressed in ADU/pixel/K.




<!-- CCD Parameters -->

## <a name="ccdParameters"></a>CCD Parameters

The <b>CCD</b> block of the configuration file contains all the information that is specific to the CCD.  The structure of this block is the following:

\code{.yaml}
CCD:

    Position:                    Custom

    OriginOffsetX:               0         
    OriginOffsetY:               0         
    Orientation:                 0         
    NumColumns:                  4510      
    NumRows:                     4510      
    FirstRowExposed:             0

    PixelSize:                   18        
    Gain:                        
    		RefValue:            1.80
    		ThreeSigma:          15.0   
    		Stability:           -0.004    
    QuantumEfficiency:           
    		Efficiency:          0.925
    		RefAngle:            45.0
    		ExpectedValue:       0.993    
    Polarization:           
    		Efficiency:          0.978
    		RefAngle:            18.8875
    		ExpectedValue:       0.989      
    Vignetting:
    		ExpectedValue:       0.945 
    Contamination:
    		ParticulateContaminationEfficiency:  0.98
    		MolecularContaminationEfficiency:    0.0566
    FullWellSaturation:          1000000        
    DigitalSaturation:           65535          
    ReadoutNoise:                28             
    ElectronicOffset:            100            
    ReadoutTime:                 2              
    FlatfieldPtPNoise:           0.016      
    CTI:
    		Model:			 Simple
    		Simple:
    		   CTEMean:        0.99999
    	      Short2013:
    	          Beta:		0.37
    	          Temperature:	203.0
    	          NumTrapSpecies:[9.8, 3.31, 1.56, 13.24]
    	          TrapDensity:   [2.46e-20, 1.74e-22, 7.05e-23, 2.45e-23]
    	          ReleaseTime:   [2.37e-4, 2.43e-2, 2.03e-3, 1.40e-1]
    NominalOperatingTemperature: 203.15
    Temperature:                 Nominal
    TemperatureFileName:         inputfiles/ccdTemperature.txt     
    IncludeFlatfield:                 no             
    IncludePhotonNoise:               yes            
    IncludeReadoutNoise:              yes            
    IncludeCTIeffects:                yes            
    IncludeOpenShutterSmearing:       yes            
    IncludeVignetting:                yes   
    IncludePolarization:              yes
    IncludeParticulateContamination:  yes
    IncludeMolecularContamination:    yes
    IncludeQuantumEfficiency:         yes
    IncludeConvolution:               yes            
    IncludeFullWellSaturation:        yes            
    IncludeDigitalSaturation:         yes      
    IncludeQuantisation:              yes      
    WriteSubPixelImagesToHDF5:        no              
\endcode




#### <a name="position"></a>Position
<i>Allowed values:</i> ∈ [1, 2, 3, 4, Custom]

The CCD position can be used to select a specific pre-defined CCD or a custom one.

The pre-defined CCD positions are shown in the figures below.

@image html "/images/CCD Array Configuration - Normal Camera.png" "Figure: Layout of the CCDs for the normal camera's."
@image html "/images/CCD Array Configuration - Fast Camera.png" "Figure: Layout of the CCDs for the fast camera's."

Note that we now use 1, 2, 3, and 4 rather than A, B, C, D, for the normal cameras as well as for the fast ones.

|In the past|Now |
|---|---|
| A  | 3  |
| B  | 2  |
| C  | 4  |
| D  | 1  |

When you specify [Position](#position)=Custom, the origin offset ([OriginOffsetX](#originOffsetX) and [OriginOffsetY](#originOffsetY)), the [orientation](#ccdOrientation), [number of rows](#ccdNumRows) and [columns](#ccdNumColumns), and the [first exposed row](#firstRowExposed) of the CCD are read from the configuration parameters in the CCD block.

In case a pre-defined position is used, these configuration parameters are read from the [CCDPositions](#ccdPositions) block (see below).




#### <a name="orginOffsetX"></a>OriginOffsetX
<i>Allowed values:</i> Any

Offset of the CCD origin from the centre of the optical plane (i.e. the intersection of the optical axis with the focal plane) in the x-direction, expressed in mm. The origin of the CCD is defined as the point where the readout register is located. See Fig. 2 for more details (Δx<sub>CCD</sub>).

This parameter is only used when the [Position](#position)=Custom.




#### <a name="originOffsetY"></a>OriginOffsetY
<i>Allowed values:</i> Any

Offset of the CCD origin from the centre of the optical plane (i.e. the intersection of the optical axis with the focal plane) in the y-direction, expressed in mm. The origin of the CCD is defined as the point where the readout register is located. See Fig. 2 for more details (Δy<sub>CCD</sub>).

This parameter is only used when the [Position](#position)=Custom.


#### <a name="ccdOrientation"></a>Orientation
<i>Allowed values:</i> Any

Orientation angle of the CCD w.r.t. the orientation of the focal plane, measured counterclockwise and expressed in degrees. This rotation is performed around the offset origin of the CCD. See Fig. 2 for more details (γ<sub>CCD</sub>).

This parameter is only used when the [Position](#position)=Custom.




#### <a name="ccdNumColumns"></a>NumColumns
<i>Allowed values:</i> > 0

Number of pixels of the CCD in the x-direction (i.e. number of columns).

This parameter is only used when the [Position](#position)=Custom.




#### <a name="ccdNumRows"></a>NumRows
<i>Allowed values:</i> > 0

Number of pixels of the CCD in the y-direction (i.e. number of rows).

This parameter is only used when the [Position](#position)=Custom.



#### <a name="firstRowExposed"></a>FirstRowExposed
<i>Allowed values:</i> > 0

Row index of the first row in the CCD that is illuminated (the row closest to the readout register is row 0).

This parameter is only used when the [Position](#position)=Custom.




#### <a name="pixelSize"></a>PixelSizeS<i>Allowed values:</i> > 0

Nominal pixel size, expressed in micron.


        
#### <a name=gainRefValueCCD></a>Gain: RefValue

<i>Allowed values:</i> > 0

Reference value of the gain of the CCD at its [nominal operating temperature](#nominalTempCCD), expressed in µV/e<sup>-</sup>.  The actually gain for the CCD will be different for both CCD halves.  We generate a normal distribution, centred at [the reference value](#gainRefValueCCD) and with a width characterised by the [ThreeSigma](#gain3SigmaCCD) parameter, and make two random draws from this distribution.  The outcome will act as gain for the left and right CCD half resp.



#### <a name=gain3SigmaCCD></a>Gain: ThreeSigma

<i>Allowed values:</i> ∈ [0,100]

Percentage of the [reference value for the gain](#gainRefValueCCD) that will act as 3σ for the normal distribution, centred around the (reference value)[#gainRefValueCCD], from which to draw the gain for both CCD halves.


#### <a name=gainStabilityCCD></a>Gain: Stability

<i>Allowed values:</i> Any

Change in gain (for both CCD halves) with temperature deviations from the nominal operating temperature, expressed in µV/e<sup>-</sup>/K.



#### <a name="quantumEfficiency"></a>QuantumEfficiency: Efficiency
<i>Allowed values:</i> ∈ [0,1]

Throughput efficiency due to the quantum efficiency at the given reference angle, considering the passband and the spectral energy distribution of the stars given the [Fluxm0](#fluxm0) parameter and the magnitude of the stars in the [star catalogue](#starCatalogue). This is the ratio of the number of collected electrons to the number of incident photons.



#### <a name="quantumEfficiencyRefAngle"></a>QuantumEfficiency: RefAngle
<i>Allowed values:</i> Any

Reference angle for the throughput efficiency due to the quantum efficiency, expressed in degrees.



#### <a name="quantumEfficiencyExpectedValue"></a>QuantumEfficiency: ExpectedValue
<i>Allowed values:</i> ∈ [0,1]

Expected value of the throughput efficiency due to quantum efficiency (i.e. the mean over all pixels of one detector).



#### <a name="polarizationEfficiency"></a>Polarization: Efficiency
<i>Allowed values:</i> ∈ [0,1]

Throughput efficiency due to the polarisation at the given reference angle.



#### <a name="PolarizationRefAngle"></a>Polarization: RefAngle
<i>Allowed values:</i> Any

Reference angle for the throughput efficiency due to the polarisation, expressed in degrees.



#### <a name="polarizationExpectedValue"></a>Polarization: ExpectedValue
<i>Allowed values:</i> ∈ [0,1]

Expected value of the throughput efficiency due to polarisation (i.e. the mean over all pixels of one detector).




#### <a name="vignettingExpectedValue"></a>Vignetting: ExpectedValue
<i>Allowed values:</i> ∈ [0,1]

Expected value of the throughput efficiency due to vignetting (i.e. the mean over all pixels of one detector).




#### <a name="particulateContamination"></a>Contamination: ParticulateContaminationEfficiency
<i>Allowed values:</i> ∈ [0,1]

Throughput efficiency due to particulate contamination.



#### <a name="molecularContamination"></a>Contamination: MolecularContaminationEfficiency
<i>Allowed values:</i> ∈ [0,1]

Throughput efficiency due to molecular contamination.





#### <a name="fullWellSaturation"></a>FullWellSaturation
<i>Allowed values:</i> > 0
     
Full-well saturation limit of a single CCD pixel, expressed in e<sup>-</sup> / pixel. If a pixels receives more electrons than its full-well saturation limit, the additional electrons flow evenly distributed in positive and negative charge-transfer direction, a phenomenon called <i>blooming</i>. The electrons reaching the edge of the CCD will not be detected..



#### <a name="digitalSaturation"></a>DigitalSaturation
<i>Allowed values:</i> > 0

Digital saturation limit of the CCD to which pixel values are topped off, expressed in ADU / pixel. This value depends on the A/D convertor of the detector. For a 16-bit convertor, the digital saturation limit is 65536 ADU.

The gain of the front-end electronics and detector should be such that the [full-well saturation](#fullWellSaturation) results in values below the digital saturation limit.


     

#### <a name="readoutNoise"></a>ReadoutNoise
<i>Allowed values:</i> ≥ 0

Mean readout noise of the detector, expressed in e<sup>-</sup>.

Readout noise occurs due to the imperfect nature of the CCD amplifiers. When the electrons are transferred to the amplifier, the induced voltage is measured. However, this measurement is not perfect, but gives a value which is on average too high by an amount of the readout noise, with the squareroot of the readout noise as standard deviation (we add the readout noise of the FEE and the CCD in quadrature).
       
       
        
#### <a name="readoutTime"></a>ReadoutTime
<i>Allowed values:</i> ≥ 0

Time required to read out an entire CCD working in frame transfer mode, and the pre-scan and the over-scan strips (to estimate the bias and the smearing resp.), expressed in seconds. Because of the absence of a shutter (which is common in space-based instruments), the CCD still receives light during frame transfer. The flux of each sub-pixel is affected by the flux of the sub-pixels in the same column. Because the CCD is exposed during the whole readout and multiple exposures are created, also the sub-pixels further away from the readout register have their influence.

For non-frame-transfer CCDs the readout time should be set to zero.




#### <a name="flatfieldPtPNoise"></a>FlatfieldPtPNoise
<i>Allowed values:</i> ∈ [0,1]]

Fractional peak-to-peak amplitude of the pixel non-uniform sensitivity response..



  
#### <a name="CTImodel"></a>CTI: Model

<i>Allowed values:</i> "Simple" and "Short2013"

Because of detector defects, electrons can get trapped in the readout process. The trapped charge ends up getting dissociated from its original pixel and eventually gets released into another pixel. The result is that the original image gets smeared out in the direction away from the readout amplifier (visible in the appearance of "charge trails"). This is known as imperfect CTE (Charge-Transfer Efficiency) or alternatively as CTI (Charge-Transfer Inefficiency).

The charge trails impact photometry, noise, and astrometry of sources. CTI removes flux from the central pixel and thus degrades the expected S/N for an observation. CTI trails bias measurements of source along the trail direction, which can severely impact high-precision astrometry.

PlatoSim3 offers two implementation of the CTI:

<ul>
	<li>The simple implementation ("Simple") assumes that for each row transfer (in the direction of the readout register) a fraction of the charge is not transferred to the next row, but stays behind.  It will be released by later row transfers.</li>
	<li>A more sophisticated implementation ("Short2013") is based on <a href="http://mnras.oxfordjournals.org/content/430/4/3078.full.pdf">Short et al., MNRAS 430, 3078-3085 (2013)</a>, in which only parallel readout is taken into account.</li>
</ul>



  
#### <a name="MeanCTE"></a>CTI: Simple: MeanCTE
<i>Allowed values:</i> ∈ [0,1]]

Mean charge-transfer efficiency (CTE) of the detector.  The fraction of the charge that is successfully transferred from one row to the next row is expressed by this parameter.



#### <a name="beta"></a>CTI: Short2013: Beta

Exponent β in Eq. (1) of Short et al. 2013 describing the relationship between the volume of the charge cloud (<i>V<sub>c</sub></i>), the number of electrons in a pixel (<i>N<sub>e</sub></i>), the full-well capacity in electrons (<i>FWC</i>), and the assumed maximum geometrical volume that electrons can occupy within a pixel (<i>V<sub>g</sub></i>):

\f[\frac{V_c}{V_g} = \left( \frac{N_e}{FWC} \right)^{\beta}.\f]




#### <a name="temperature"></a>CTI: Short2013: Temperature

<i>Allowed values:</i> ≥ 0

Temperature <i>T</i> that is used to calculated the thermal velocity <i>v<sub>t</sub></i> of the electrons:

\f[v_t = \frac{3kT}{m_e^{\ast}},\f]

where <i>k</i> is the Boltzmann constant and <i>m<sub>e</sub><sup>*</sup></i> is the effective electron mass in silicon, which we approximate by half the free electron rest mass.




#### <a name="numTrapSpecies"></a>CTI: Short2013: NumTrapSpecies

<i>Allowed values:</i> > 0

Number of trap species that is used in the CTI model by Short et al. 2013.




#### <a name="trapDensity"></a>CTI: Short2013: TrapDensity

<i>Allowed values:</i> Array holding one non-negative entry per trap species. 

Array holding the trap density <i>n<sub>t</sub></i> for each of the considered trap species, expressed in number of traps per pixel.  This is used to calculate the γ-value in Eq. (22) of Short et al. 2013.




#### <a name="trapCaptureCrossSection"></a>CTI: Short2013: TrapCaptureCrossSection

<i>Allowed values:</i> Array holding one non-negative entry per trap species.

Array holding the trap capture cross-section <i>σ</i> for each of the considered trap species, expressed in m<sup>2</sup>.  This is used to calculated the α-value in Eq. (22) of Short et al. 2013.  In this formula, the charge transfer time is used as value for <i>t</i>.




#### <a name="releaseTime"></a>CTI: Short2013: ReleaseTime 

<i>Allowed values:</i> Array holding one non-negative entry per trap species.

Array holding the trap release time constants <i>τ<sub>r</sub></i> for each of the considered trap species, expressed in seconds.



#### <a name="nominalTempCCD"></a>NominalOperatingTemperature

<i>Allowed values:</i> > 0

Nominal operating temperature of the CCD, expressed in Kelvin.



#### <a name="tempCCD"></a>Temperature
<i>Allowed values:</i>"Nominal" or "FromFile"

Indicates whether the temperature of the CCD should be fixed at the nominal operating temperature or temperature variations should be read from a file.



#### <a name=tempFileCCD></a>TemperatureFileName
Path to the file, relative to the [project location](#projectLocation), holding the location of the file with the temperature variations of the CCD.





#### <a name="inclFlatfield"></a>IncludeFlatfield
<i>Allowed values:</i> "yes" and "no"

Indicates whether or not to include the flatfield.




#### <a name="inclPhotonNoise"></a>IncludePhotonNoise
<i>Allowed values:</i> "yes" and "no"

Indicates whether or not to include photon noise.




#### <a name="inclReadoutNoise"></a>IncludeReadoutNoise
<i>Allowed values:</i> "yes" and "no"

Indicates whether or not to include readout noise.




#### <a name="inclCTI"></a>IncludeCTIeffects
<i>Allowed values:</i> "yes" and "no"

Indicates whether or not to include CTI effects.




#### <a name="inclOpenShutterSmearing"></a>IncludeOpenShutterSmearing
<i>Allowed values:</i> "yes" and "no"

Indicates whether or not to include open-shutter smearing effects.




#### <a name="inclVignetting"></a>IncludeVignetting
<i>Allowed values:</i> "yes" and "no"

Indicates whether or not to include brightness attenuation towards the edge of the FOV due to vignetting.


        

#### <a name="inclPolarization"></a>IncludePolarization
<i>Allowed values:</i> "yes" and "no"

Indicates whether or not to include loss of throughput efficiency due to polarisation.


        

#### <a name="inclParticulateContamination"></a>IncludeParticulateContamination
<i>Allowed values:</i> "yes" and "no"

Indicates whether or not to include loss of throughput efficiency due to particulate contamination.


        


#### <a name="inclMolecularContamination"></a>IncludeMolecularContamination
<i>Allowed values:</i> "yes" and "no"

Indicates whether or not to include loss of throughput efficiency due to molecular contamination.


        

#### <a name="inclQuantumEfficiency"></a>IncludeQuantumEfficiency
<i>Allowed values:</i> "yes" and "no"

Indicates whether or not to include loss of throughput efficiency due to quantum efficiency.


        

#### <a name="inclConvolution"></a>IncludeConvolution
<i>Allowed values:</i> "yes" and "no"

Indicates whether or not the sub-pixel map must be convolved with the PSF.




#### <a name="inclFullWellSaturation"></a>IncludeFullWellSaturation
<i>Allowed values:</i> "yes" and "no"

Indicates whether or not to apply full-well saturation.




#### <a name="inclDigitalSaturation"></a>IncludeDigitalSaturation
<i>Allowed values:</i> "yes" and "no"

Indicates whether or not to apply digital saturation.



#### <a name="inclQuantisation"></a>IncludeQuantisation
<i>Allowed values:</i> "yes" and "no"

Indicates whether or not to apply quantisation.  This includes:
	* applying gain (FEE and CCD), hence converting from electrons to ADUs
	* adding electronic offset
	* forcing the ADU values to be integers
	* applying digital saturation (can be [switched off separately](#inclDigitalSaturation))



#### <a name="writeSubPixelImages"></a>WriteSubPixelImages
<i>Allowed values:</i> "yes" and "no"

Indicates whether or not the sub-pixel images must be written to the HDF5-file.  Use this for a limited number of exposures, as it takes a lot of space.





<!-- Sub-Field Parameters -->
<!-- ******************** -->

## <a name="subFieldParameters"></a>Sub-Field Parameters

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




@image html /images/subField.png "Figure 3: Schematic presentation of the modelled sub-field and the parameters required to define it. The pixel coordinates of the origin of the sub-field relative to the CCD are xs and ys."




#### <a name="zeroPointRow"></a>ZeroPointRow
<i>Allowed values:</i> > 0

Row of the origin of the sub-field in the detector, expressed in pixels.  See Fig. 3 for more details (y<sub>s</sub>).




#### <a name="zeroPointColumn"></a>ZeroPointColumn
<i>Allowed values:</i> > 0

Column of the origin of the sub-field in the detector, expressed in pixels.  See Fig. 3 for more details (x<sub>s</sub>).




#### <a name="numColumns"></a>NumColumns
<i>Allowed values:</i> ≥ 8

Number of columns in the sub-field, expressed in pixels.




#### <a name="numRows"></a>NumColumns
<i>Allowed values:</i> ≥ 8

Number of rows in the sub-field, expressed in pixels.




#### <a name="numPreScanRows"></a>NumBiasPrescanRows
<i>Allowed values:</i> ≥ 0

Number of rows in the pre-scan strip (see Fig. 3), expressed in normal pixel units.   This strip is located at the bottom of the sub-field that is modelled in detail and contains the electronic offset and readout noise.




#### <a name="numOverScanRows"></a>NumSmearingOverscanRows
<i>Allowed values:</i> ≥ 0

Number of rows in the over-scan strip (see Fig. 3), expressed in normal pixel units. This strip is located at the top of the sub-field that is modelled in detail and contains the star smearing due to the absence of a shutter. This flux in this strip is also affected by the electronic offset, readout noise, and shot noise. Not included are the PRNU, cosmic hits, and charge-transfer efficiency (CTE).




#### <a name="numSubPixels"></a>SubPixels
<i>Allowed values:</i> power of 2 (≤ 128)

Number of sub-pixels per pixel in both directions.

If you want a pixel of 256 x 256 = 65536 sub-pixels, you should specify in the configuration file 256. The total number of sub-pixels per pixel will then be 65536. If you want 256 = 16 x 16 sub-pixels per pixel, then you should specify 16 in the configuration file.





<!-- Seed Parameters -->
<!-- *************** -->

## <a name="seedParameters"></a>Seed Parameters

The <b>RandomSeeds</b> block of the configuration file contains all the seeds for random-number generation in the simulator.  The structure of this block is the following:

\code{.yaml}
RandomSeeds:

    ReadOutNoiseSeed:            1424949740 
    PhotonNoiseSeed:             1433320336 
    JitterSeed:                  1433320381 
    FlatFieldSeed:               1425284070 
    CTESeed:                     1424949740 
    DriftSeed:                   1433429158 
    FeeGainSeed:                 1429485030
    CcdGainSeed:                 1420450443
\endcode



#### ReadOutNoiseSeed
<i>Allowed values:</i> > 0

Seed for the random-number generator used for the readout noise.



#### PhotonNoiseSeed
<i>Allowed values:</i> > 0

Seed for the random-number generator used for the photon noise.



#### JitterSeed
<i>Allowed values:</i> > 0

Seed for the random-number generator used for the jitter.



#### FlatFieldSeed
<i>Allowed values:</i> > 0

Seed for the random-number generator used for the flatfield.



#### CTESeed
<i>Allowed values:</i> > 0

Seed for the random-number generator used for the CTE.



#### DriftSeed
<i>Allowed values:</i> > 0

Seed for the random number generator used for the drift.


#### FeeGainSeed
<i>Allowed values:</i> > 0

Seed for the random-number generator used for the FEE gain.


#### CcdGainSeed
<i>Allowed values:</i> > 0

Seed for the random-number generator used for the CCD gain.




<!-- Camera groups -->
<!-- ************* -->

### <a name="cameraGroups"></a>Camera groups

The <b>CameraGroups</b> block in the configuration file is used in case a pre-defined camera groups (1, 2, 3, 4, or Fast) was selected via the [GroupID](#groupID) parameter in the [Telescope](#telescopeParameters) block in the configuration file.  The structure of this block is the following:

\code{.yaml}
CameraGroups:

    AzimuthAngle:            [45.0, 135.0, -135.0, -45.0, 0.0] 
    TiltAngle:               [9.2, 9.2, 9.2, 9.2, 0.0] 
\endcode

Mind you, you are NOT supposed to alter this section of the configuration file!



#### AzimuthAngle

Azimuth angle, expressed in degrees, for camera group 1, 2, 3 and 4, and for the fast camera.  Depending on the value of the [GroupID](#groupID) parameter in the [Telescope](#telescopeParameters) block, the appropriate value will be selected from the list.

#### TiltAngle

Tilt angle, expressed in degrees, for camera group 1, 2, 3 and 4, and for the fast camera.  Depending on the value of the [GroupID](#groupID) parameter in the [Telescope](#telescopeParameters) block, the appropriate value will be selected from the list.






<!-- CCD Positions -->
<!-- ************* -->
### <a name="ccdPositions"></a>CCD Positions

The <b>CCDPositions</b> block in the configuration file is used in case a pre-defined CCD position (1, 2, 3, or 4) was selected via the [Position](#position) parameter in the [CCD](#ccdParameters) block in the configuration file.  The structure of this block is the following:

\code{.yaml}
CCDPositions:

    OriginOffsetX:                   [-1, -1, -1, -1]
    OriginOffsetY:                   [82.18, 82.18, 82.18, 82.18]
    Orientation:                     [0, 90, 180, 270]
    NumColumns:                      [4510, 4510, 4510, 4510]
    NumRows:                         [4510, 4510, 4510, 4510]
    FirstRowForNormalCamera:         [0, 0, 0, 0]
    FirstRowForFastCamera:           [2255, 2255, 2255, 2255]       
\endcode

Mind you, you are NOT supposed to alter this section of the configuration file!

#### OriginOffsetX

Offset of the CCD origin from the centre of the optical plane in the x-direction, expressed in mm, for CCD positions 1, 2, 3, and 4.  Depending on the value of the [Position](#position) parameter in the [CCD](#ccdParameters) block, the appropriate value will be selected from the list.

#### OriginOffsetY

Offset of the CCD origin from the centre of the optical plane in the y-direction, expressed in mm, for CCD positions 1, 2, 3, and 4.  Depending on the value of the [Position](#position) parameter in the [CCD](#ccdParameters) block, the appropriate value will be selected from the list.


#### Orientation

Orientation angle of the CCD w.r.t. the orientation of the focal plane, measured counterclockwise and expressed in degrees, for CCD positions 1, 2, 3, and 4.  Depending on the value of the [Position](#position) parameter in the [CCD](#ccdParameters) block, the appropriate value will be selected from the list.

#### NumColumns

Number of pixels of the CCD in the x-direction (i.e. number of columns), for CCD positions 1, 2, 3, and 4.  Depending on the value of the [Position](#position) parameter in the [CCD](#ccdParameters) block, the appropriate value will be selected from the list.

#### NumRows

Number of pixels of the CCD in the y-direction (i.e. number of rows), for CCD positions 1, 2, 3, and 4.  Depending on the value of the [Position](#position) parameter in the [CCD](#ccdParameters) block, the appropriate value will be selected from the list.


#### FirstRowForNormalCamera

Row index of the first row in the CCD that is illuminated (the row closest to the readout register is row 0), for CCD positions 1, 2, 3, and 4, in case of a normal camera ([GroupID](#groupID)=1, 2, 3, or 4).  For normal cameras, the whole CCD is illuminated.  Depending on the value of the [Position](#position) parameter in the [CCD](#ccdParameters) block, the appropriate value will be selected from the list.



#### FirstRowForFastCamera

Row index of the first row in the CCD that is illuminated (the row closest to the readout register is row 0), for CCD positions 1, 2, 3, and 4, in case of a fast camera ([GroupID](#groupID)=Fast).  For fast cameras, only the upper half of the CCD is illuminated.  Depending on the value of the [Position](#position) parameter in the [CCD](#ccdParameters) block, the appropriate value will be selected from the list.



<!-- ************** -->
<!-- Star Catalogue -->
<!-- ************** -->

# <a name="starCatalogue"></a>Star Catalogue

A star catalogue must be provided in a file in ASCII format. This file should contain three columns, separated by a space, holding the following information:
* right ascension of the stars [degrees]
* declination of the stars [degrees]
* stellar magnitude

The path of this file, relative to the [project location](#projectLocation), must be provided via the [StarCatalogFile](#starCatalogFile) parameter in the configuration file (under observing parameters). 





<!-- ******** -->
<!-- PSF File -->
<!-- ******** -->

# <a name="psfFile"></a>PSF File (Optional)

Unless you indicate you want to generate a Gaussian PSF, pre-computed normalised PSFs must be provided in the form of an HDF5-file. The path of this file, relative to the [project location](#projectLocation), is specified via the [Filename](#psfFilename) parameter.

The simulator will automatically select the PSF for which the angular distance to the optical axis matches best for the simulated sub-field.





<!-- *********** -->
<!-- Jitter File -->
<!-- *********** -->

# <a name="jitterFile"></a>Jitter File (Optional)

If required ([UseJitterFromFile](#useJitterFromFile) = "yes"), a jitter time series must be provided in a file in ASCII format. This file should contain four columns, separated by a space, holding the following information:

* time [s]
* yaw [arsec]
* pitch [arcsec]
* roll [arcsec]

The path of this file, relative to the [project location](#projectLocation), must be provided via the  [JitterFileName](#jitterFileName) parameter in the configuration file.

To ensure a realistic modelling of the jitter, the time step in the jitter time series must be smaller than the [exposure time](#exposureTime).


