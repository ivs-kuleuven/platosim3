# Configuration Parameters {#ConfigurationParameters}

To configure the Plato Simulator, a large set of configuration parameters is required.  The input file format used for PlatoSim3 is <a href="https://learnxinyminutes.com/docs/yaml/">YAML</a>, e.g. <code>inputfile.yaml</code>, in the <code>/inputfiles</code> directory.  The different blocks in the configuration files reflect their function in the simulation.

We use only a very limited set of the YAML functionality, enough to allow us to provide input files for different parts of the simulator. 

Any desired simulation can be obtained by modifying the following input:
		* [general parameters](#generalParameters)
		* [observing parameters](#observingParameters)
		* [sky parameters](#skyParameters)
		* [platform parameters](#platformParameters)
		* [telescope parameters](#telescopeParameters)
		* [camera parameters](#cameraParameters)
		* [PSF parameters](#psfParameters)
		* [FEE parameters](#feeParameters)
		* [CCD parameters](#ccdParameters)
		* [sub-field parameters](#subFieldParameters)
        * [photometry parameters](#photometryParameters)
		* [seed parameters](#seedParameters)
        * [control TCP connection parameters](#controlTcpConnection)
		* additionally, there are two blocks that hold pre-defined settings (which you should NOT alter):
			- [camera group 1, 2, 3, and 4, and fast cameras](#cameraGroups)
			- [CCD 1, 2, 3, and 4](#ccdPositions)
 
In the following sections we describe these parameters for the simulations in detail.

For more details on the reference frames we are using (for the spacecraft, telescope, focal plane, and CCD), radial dependency of the PSF, and rotation angles for platform jitter and telescope drift, please, have a look at technical note [PLATO-KUL-PL-TN-0001](../technicalnotes/PLATO-KUL-PL-TN-0001.pdf).

---





<!-- ****************** -->
<!-- General Parameters -->
<!-- ****************** -->

## <a name="generalParameters"></a>General Parameters

The general configuration parameters are listed in the <b>General</b> block of the configuration file.  The structure of this block is the following:

\code{.yaml}
General:
   
    ProjectLocation:             ENV['PLATO_PROJECT_HOME']
\endcode



### <a name="projectLocation"></a>ProjectLocation

<i>Allowed values:</i> name of an existing directory on disk or environment variable, in the format <code>ENV['PLATO_PROJECT_HOME']</code>.

Full path of the directory in which you have checked out the PlatoSim3 project, or an environment variable, e.g. <code>PLATO_PROJECT_HOME</code>, containing the full path to that directory.  In the latter case, you must make sure you have exported this variable before initiating a simulation:

\code{.unparsed}
 export PLATO_PROJECT_HOME=<full path to the PlatoSim3 directory>
\endcode

---





<!-- ******************** -->
<!-- Observing Parameters -->
<!-- ******************** -->

## <a name="observingParameters"></a>Observing Parameters

The <b>ObservingParameters</b> block of the configuration file contains the configuration parameters that are specific to the simulated observation and are not specific for the hardware components of the satellite.  The structure of this block is the following:

\code{.yaml}
ObservingParameters:

	MissionDuration:             6.0
	BeginExposureNr:             0
	NumExposures:                40              
    CycleTime:                   25              
    RApointing:                  180              
    DecPointing:                 -70             
    Fluxm0:                      1.00179e8 
    StarCatalogFile:             inputfiles/starcatalog.txt
\endcode




### <a name="missionDuration"></a>MissionDuration
<i>Allowed values:</i> > 0

Total duration of the mission (from BOL till EOL), expressed in years.  This will be used to model parameter degradation over time.



### <a name="beginExposureNr"></a>BeginExposureNr
<i>Allowed values:</i> \f$\ge \f$ 0

Sequential number of the first exposure. Useful for <a href="https://en.wikipedia.org/wiki/Slurm_Workload_Manager">Slurm</a> parallelisation.  In that case, long simulations (i.e. with a large number of exposures) will be chopped up into smaller simulations, covering [a small number of exposures](#NumExposures) (see Fig. 1).

@image html /images/chopUpSimulation.png "Figure 1: Long simulations will be chopped up into smaller simulations that can be executed in parallel."



### <a name="numExposures"></a>NumExposures
<i>Allowed values:</i> > 0

Number of exposures to generate in the simulation.



### <a name="cycleTime"></a>CycleTime
<i>Allowed values:</i> > 0

Image cycle time, expressed in seconds.  This is the sum of the integration time of one exposure and the duration of the readout of one exposure before the next exposure start:

	\f[ t_{cycle} = t_{exposure} + t_{readout, before}.\f]

For the normal cameras, the latter is the total readout time; for the fast cameras, it is the time for the frame transfer (i.e. to transfer the content of the upper CCD half to the lower CCD half).



### <a name="raPointing"></a>RApointing
<i>Allowed values:</i>  \f$\in \f$ [0, 360]

Right ascension of the pointing, expressed in degrees.



### <a name="decPointing"></a>DecPointing
<i>Allowed values:</i> \f$\in \f$ [-90, 90]

Declination of the pointing, expressed in degrees.



### <a name="fluxm0"></a>Fluxm0
<i>Allowed values:</i> > 0

Flux of a star of zero magnitude (\f$ m_{\lambda} = 0 \f$), expressed in photons \f$ \cdot \f$  s<sup>-1</sup> \f$  \cdot \f$  cm<sup>-2</sup> in the passband of the magnitudes that are listed in the [star catalogue](#starCatalogue).

For an exposure of \f$t_{exp}\f$ seconds, the measured flux \f$F_{phot}\f$ of a star, expressed in photons, is computed from its catalogue magnitude \f$m_{\lambda}\f$, the [effective light-collecting area](#lightCollectingArea) \f$A\f$ (in cm<sup>2</sup>) of the telescope, the [transmission efficiency](#transmissionEfficiency) \f$T_{\lambda}\f$ of the optical system, the [quantum efficiency](#quantumEfficiency) \f$Q\f$ of the detector, and the flux per second \f$F_0\f$ of a star with zero magnitude (\f$m_{\lambda} = 0\f$) from the equation

\f[F_{phot} = t_{exp} \cdot F_0 \cdot T_{\lambda} \cdot Q \cdot A \cdot 10^{-0.4 \cdot m_{\lambda}}\f]

where the \f$\lambda\f$ subscript refers to the wavelength range in which the simulation is performed.



### <a name="starCatalogFile"></a> StarCatalogFile

Path to the star catalogue file, relative to the [project location](#projectLocation).

---





<!-- ************** -->
<!-- Sky Parameters -->
<!-- ************** -->

## <a name="skyParameters"></a>Sky Parameters

The <b>Sky</b> block of the configuration file contains all the information that is specific to the sky, i.e. sky background and cosmics.  The structure of this block is the following:

\code{.yaml}
Sky:

	SkyBackground:               342.
	IncludeVariableSources:      no
    	VariableSourceList:          inputfiles/varsource.txt
	IncludeCosmicsInSubField:        yes
    IncludeCosmicsInSmearingMap:     yes
    IncludeCosmicsInBiasMap:         yes    
	Cosmics:
		CosmicHitRate:                      10
		TrailLength:                   [0, 15]
		Intensity:               [2000, 40000] 
\endcode





### <a name="skyBackground"></a>SkyBackground
<i>Allowed values:</i> < 0 for automatic calculation, \f$\ge \f$ 0 to use the input value

In case a positive value is given, the sky background (zodiacal + galactic), is set to the given value, expressed in photons \f$ \cdot \f$ s<sup>-1</sup> \f$ \cdot \f$ pixel<sup>-1</sup>.  Note that this value has not been multiplied with the tranmission efficiency yet.

In case a negative value is given, the sky background is computed automatically from tabular values, interpolated to the central coordinates of the sub-field. A constant sky background is assumed for the whole sub-field. Note that for some regions in the sky the automatic computation of the sky background may fail due to the lack of tabulated values. In that case you can set the sky background manually.



### <a name="inclVariability"></a>IncludeVariableSources
<i>Allowed values:</i> "yes" and "no"

Indicates whether or not stellar variability must be included.


### <a name="variableSourceList"></a>VariableSourceList
<i>Allowed values:</i> only required if stellar variability is to be included ([IncludeVariableSources](#inclVariability) = "yes")

Path to the file, relative to the [project location](#projectLocation), indicating how the magnitude of the sources varies over time.



### <a name="inclCosmicsInSubField"></a>IncludeCosmicsInSubField
<i>Allowed values:</i> "yes" and "no"

Indicates whether or not cosmics must be added to the pixel map.




### <a name="inclCosmicsInBiasMap"></a>IncludeCosmicsInBiasMap
<i>Allowed values:</i> "yes" and "no"

Indicates whether or not cosmics must be added to the bias register map.




### <a name="inclCosmicsInSmearingMap"></a>IncludeCosmicsInSmearingMap
<i>Allowed values:</i> "yes" and "no"

Indicates whether or not cosmics must be added to the smearing map.




### <a name="cosmics"></a>Cosmics

The configuration parameters in the <b>Cosmics</b> section are the parameters characterising the cosmic hits.  The excess electrons of the cosmic hits are distributed over a trail, characterised by a decay function of the form:

\f[ f(t) = \exp{\left(\frac{-t^2}{2 * \sigma^2}\right)}. \f]




#### <a name="cosmicHitRate"></a>Cosmics: CosmicHitRate
<i>Allowed values:</i> >= 0

Mean cosmic hit rate, expressed in events / cm<sup>2</sup> / s.  The actual cosmic hit rate for any exposure is sampled randomly from a Poisson distribution with the mean cosmic hit rate as mean.  The number of cosmic hits in the simulated sub-field is calculated by multiplying the actual cosmic hit rate with the size of the sub-field (expressed in cm<sup>2</sup>) and the [cycle time](#cycleTime) (expressed in s).



#### <a name="cosmicTrailLength"></a>Cosmics: TrailLength

Interval for the allowed length of the cosmic trails, expressed in pixels.


#### <a name="cosmicIntensity"></a>Cosmics: Intensity

Interval for the allowed number of electrons comprised in a cosmic hit (that are to be spread over the trail).

---





<!-- ******************* -->
<!-- Platform Parameters -->
<!-- ******************* -->

## <a name="platformParameters"></a>Platform Parameters


The <b>Platform</b> block of the configuration file contains all the information that is specific to the platform of the satellite.  The structure of this block is the following:

\code{.yaml}
Platform:

    SolarPanelOrientation:       0
    UseJitter:                   yes             
    JitterSource:                FromRedNoise
    JitterYawRms:                1.0             
    JitterPitchRms:              1.0             
    JitterRollRms:               1.0             
    JitterTimeScale:             3600.           
    JitterFileName:              /inputfiles/jitter.txt
\endcode



### <a name="solarPanelOrientation"></a>SolarPanelOrientation
<i>Allowed values:</i> 0, 90, 180, and 270

Orientation angle of the solar panel, expressed in degrees.  This is the roll angle of the platform, which enables orienting the solar panels towards the Sun each quarter, i.e. at the beginning of each quarter the roll angle must be increased by 90 degrees.  By convention the roll angle must be set to zero degrees at the beginning of the first quarter.

Note that - to properly account for the re-orientation of the solar panels - simulations must be chopped up in chunks of maximum three months.



### <a name="useJitter"></a>UseJitter
<i>Allowed values:</i> "yes" and "no"

Indicates whether pointing variations should be taken into account.

The Plato Simulator can also account for pointing variations of the spacecraft, so-called jitter. A time series of pointing displacement, expressed in Euler angles (yaw, pitch, roll), either has to be provided as a jitter file or will be generated based on the given jitter parameters (see further).

To ensure a realistic modelling of the jitter, the [time step of the jitter time series](#jitterTimeScale) must be smaller than the [exposure time](#cycleTime).

The configuration of the jitter axes is depicted below.  The Euler angles that characterise the jitter are defined w.r.t. to the spacecraft coordinate system (see Fig. 2).  The origin of this coordinate system is the geometric centre of the interface between the bottom of the optical bench and the service module.  The positive roll axis \f$z_{\rm SC} \f$ points towards the operator-given mean payload line-of-sight, given by the equatorial coordinates ([RApointing](#raPointing), [DecPointing](#decPointing)).

The angles are defined such that they increase with a clockwise rotation, when looking along the positive axes. First a roll rotation is done around the \f$z_{\rm SC} \f$ axis, then a pitch rotation is done around the rotated \f$y_{\rm SC} \f$ axis, and finally a yaw rotation is done around the twice-rotated \f$x_{\rm SC} \f$ axis.

@image html /images/jitterConfiguration.png "Figure 2: Configuration of the jitter axes for the Plato Simulator, defined w.r.t. the spacecraft coordinate system (\f$xSC \f$, \f$y_{\rm SC} \f$, \f$z_{\rm SC} \f$).  The origin of this coordinate system is the geometric centre of the interface between the bottom of the optical bench and the service module.  The positive \f$z_{\rm SC} \f$ axis points towards the operator-given pointing coordinates. The xSC axis points in the direction of the highest point of the sunshield."

@image html /images/jitterConfiguration.png "Figure 2: Configuration of the jitter axes for the Plato Simulator, defined w.r.t. the spacecraft coordinate system (x<sub>SC</sub>, y<sub>SC</sub>, z<sub>SC</sub>).  The origin of this coordinate system is the geometric centre of the interface between the bottom of the optical bench and the service module.  The positive z<sub>SC</sub> axis points towards the operator-given pointing coordinates. The x<sub>SC</sub> axis points in the direction of the highest point of the sunshield."



### <a name="jitterSource"></a>JitterSource
<i>Allowed values:</i> "yes" and "no"

Indicates whether the jitter time series must be read from a jitter file ("yes") or the jitter positions must be generated from the jitter parameters ("no").
FromFile, FromRedNoise, or FromNetwork

<i>Allowed values:</i> "FromRedNoise", "FromFile", and "FromNetwork"

Indicates from where to read the jitter:

- FromRedNoise: the jitter positions must be generator from the jitter parameters;
- FromFile: the jitter time series must be read from a jitter file;
- FromNetwork: the jitter positions must be read from a network (which is configured in the [ControlTcpConnection](#controlTcpConnection) block).


### <a name="jitterYawRms"></a>JitterYawRms
<i>Allowed values:</i> \f$\ge \f$ 0, only required if the jitter positions must be generated from jitter parameters ([JitterSource](#jitterSource) = FromRedNoise)

Standard deviation (expressed in arcsec) of the normal distribution (with zero mean) describing the yaw value from one jitter position to the next one.



### <a name="jitterPitchRms"></a>JitterPitchRms
<i>Allowed values:</i> \f$\ge \f$ 0, only required if the jitter positions must be generated from jitter parameters ([JitterSource](#jitterSource) = FromRedNoise)

Standard deviation (expressed in arcsec) of the normal distribution (with zero mean) describing the pitch value from one jitter position to the next one.



### <a name="jitterRollRms"></a>JitterRollRms
<i>Allowed values:</i> \f$\ge \f$ 0, only required if the jitter positions must be generated from jitter parameters ([JitterSource](#jitterSource) = FromRedNoise)

Standard deviation (expressed in arcsec) of the normal distribution (with zero mean) describing the roll value from one jitter position to the next one.



### <a name="jitterTimeScale"></a>JitterTimeScale
<i>Allowed values:</i> > 0

Timescale of the jitter (i.e. time between two subsequent jitter positions), expressed in seconds.



### <a name="jitterFileName"></a>JitterFileName

Path of the jitter file, relative to the [project location](#projectLocation). This is only required if the jitter positions must be read from a file ([JitterSpirce](#jitterSource) = FromFile).

---





<!-- ******************** -->
<!-- Telescope Parameters -->
<!-- ******************** -->

## <a name="telescopeParameters"></a>Telescope Parameters

The <b>Telescope</b> block of the configuration file contains all the information that is specific to the telescope.  The structure of this block is the following:

\code{.yaml}
Telescope:
    
    GroupID:                     Custom
    AzimuthAngle:                0.0
    TiltAngle:                   0.0
    LightCollectingArea:         113.1         
    TransmissionEfficiency:      
        BOL:                     0.7191
        EOL:                     0.7191
    UseDrift:                    yes
    UseDriftFromFile:            no      
    DriftYawRms:                 2.3           
    DriftPitchRms:               2.3           
    DriftRollRms:                2.3           
    DriftTimeScale:              3600.
    DriftFileName:               /inputfiles/drift.txt         
\endcode





### <a name="groupID"></a>GroupID
<i>Allowed values:</i> \f$\in \f$ [1, 2, 3, 4, Fast, Custom]

The telescope group identifier can be used to select a telescope group. There are four groups that have a tilt angle of \f$9.2^{\circ} \f$- from the optical axis of the satellite, and one group for the fast camera's which is aligned with the satellite Z-axis. When you specify GroupID=Custom, the [tilt angle](#tiltAngle) and [azimuth angle](#azimuthAngle) below the GroupID in the inputfile are used, otherwise the angles are taken from pre-defined parameters in the [CameraGroups](#cameraGroups) block of the configuration file.

@image html /images/telescopeGroups.png "Figure 3: Field of View for the different telescope groups"



### <a name="tiltAngle"></a>TiltAngle
<i>Allowed values:</i> > 0

Tilt angle of the telescope, expressed in degrees. This angle, together with the [azimuth angle](#azimuthAngle), characterises the orientation of the telescope pointing (i.e. telescope optical axis) w.r.t. the spacecraft/platform pointing. 

The tilt angle is the offset between the telescope optical axis and the platform pointing, i.e. the angle between the telescope line-of-sight (positive \f$z_{\rm telescope} \f$-axis and the positive \f$z_{\rm PLM} \f$-axis (see Figs. 4 and 5).

This parameter is only used when the [GroupID](#groupID)=Custom.



### <a name="azimuthAngle"></a>AzimuthAngle
<i>Allowed values:</i> Any

Azimuth angle of the telescope, expressed in degrees. This angle, together with the [tilt angle](#tiltAngle), characterises the orientation of the telescope pointing (i.e. telescope optical axis) w.r.t. the spacecraft/platform pointing. 

The azimuth angle is the position angle of the rotation of the telescope around the positive \f$z_{\rm PLM} \f$-axis (see Figs. 3 and 4).

This parameter is only used when the [GroupID](#groupID)=Custom.

@image html /images/tiltAzimuth.png "Figure 4: Tilt and azimuth of a telescope."



### <a name="lightCollectingArea"></a>LightCollectingArea
<i>Allowed values:</i> > 0

Light-collecting area of one telescope, expressed in cm<sup>2</sup>.



### <a name="transmissionEfficiency"></a>TransmissionEfficiency
The transmission efficiency of the optical system, considering the passband and the spectral energy distribution of the stars, given the [Fluxm0](#fluxm0) parameter and the magnitudes in the [star catalogue](#starCatalogue) degrades linearly over the [mission lifetime](#missionDuration).



#### <a name="transmissionEfficiencyBOL"></a>TransmissionEfficiency: BOL

<i>Allowed values:</i> \f$\in \f$ [0,1]

Tranmission efficiency of the optical system, considering the passband and spectral energy distribution of the stars, given the [Fluxm0](#fluxm0) parameter and the magnitudes in the [star catalogue](#starCatalogue), at the beginning of the mission (beginning-of-life).  This parameter is used to model the (linear) degradation in transmission efficiency over the [mission lifetime](#missionDuration).



#### <a name="transmissionEfficiencyEOL"></a>TransmissionEfficiency: EOL
<i>Allowed values:</i> \f$\in \f$ [0,1]

Tranmission efficiency of the optical system, considering the passband and spectral energy distribution of the stars, given the Fluxm0 parameter and the magnitudes in the [star catalogue](#starCatalogue), at the end of the mission (end-of-life).  This parameter is used to model the (linear) degradation in transmission efficiency over the [mission lifetime](#missionDuration).



### <a name="useDrift"></a>UseDrift
<i>Allowed values:</i> "yes" or "no"

Indicates whether the thermo-elastic drift of the telescope (w.r.t. the platform) should be taken into account.

Similar to the [UseJitter](#useJitter) parameter for the platform jitter.

The Plato Simulator can also account for the thermo-elastic drift, of the telescope (w.r.t. the platform). A time series of displacement, expressed in Euler angles (yaw, pitch, roll), either has to be provided as a drift file or will be generated based on the given drift parameters (see further).

The Euler angles (yaw, pitch, roll) are defined as the rotation angles around the \f$z_{\rm SC} \f$, \f$y'_{\rm SC} \f$ and 
\f$z_{\rm telescope} = z_{\rm FP} \f$ axes (see Fig. 5), such that the anges increase with a clockwise rotation when looking along the positive axes.



@image html /images/TelescopeCoordinateSystem.png "Figure 5: The optical axis zFP can be obtained from the spacecraft/platform pointing axis z<sub>SC</sub> by first rotating the (x<sub>SC</sub>, y<sub>SC</sub>) plane around the pointing axis z<sub>SC</sub><sub>SC</sub> over the azimuth angle (left-hand side) nad then rotating the resulting zSC' axis over the tilt angle (right-hand side)."


### <a name="useDriftFromFile"></a>UseDriftFromFile
<i>Allowed values:</i> "yes" or "no"

Indicates whether the thermo-elastic drift of the telescope (w.r.t. the platform) should be taken into account.

Similar to the [UseJitterFromFile](#useJitterFromFile) parameter for the platform jitter.




### <a name="driftYawRms"></a>DriftYawRms
<i>Allowed values:</i> \f$\ge \f$ 0

Standard deviation (expressed in arcsec) of the normal distribution (with zero mean) describing the yaw value from one thermo-elastic drift position to the next one.

Similar to the [JitterYawRms](#jitterYawRms) parameter for the platform jitter.



### <a name="driftPitchRms"></a>DriftPitchRms
<i>Allowed values:</i> \f$\ge \f$ 0

Standard deviation (expressed in arcsec) of the normal distribution (with zero mean) describing the pitch value from one thermo-elastic drift position to the next one.

Similar to the [JitterPitchRms](#jitterPitchRms) parameter for the platform jitter.



### <a name="driftRollRms"></a>DriftRollRms
<i>Allowed values:</i> \f$\ge \f$ 0

Standard deviation (expressed in arcsec) of the normal distribution (with zero mean) describing the roll value from one thermo-elastic drift position to the next one.

Similar to the [JitterRollRms](#jitterRollRms) parameter for the platform jitter.



### <a name="driftTimeScale"></a>DriftTimeScale
<i>Allowed values:</i> > 0

Timescale of the thermo-elastic drift (i.e. time between two subsequent drift positions), expressed in seconds.


### <a name="driftFileName"></a>DriftFileName
Path of the drift file, relative to the [project location](#projectLocation). This is only required if the drift positions must be read from a file ([UseDriftFromFile](#useDriftFromFile) = yes).

---





<!-- ***************** -->
<!-- Camera Parameters -->
<!-- ***************** -->

## <a name="cameraParameters"></a>Camera Parameters

The <b>Camera</b> block of the configuration file contains all the information that is specific to the camera.  The structure of this block is the following:

\code{.yaml}
Camera:
    
    FocalPlaneOrientation:
        Source:                      ConstantValue
        ConstantValue:               0.0 
        FromFile:                    inputfiles/fporientation.txt           
    PlateScale:                  0.8333          
    FocalLength: 
        Source:                      FromFile 
        ConstantValue:               0.24752
        FromFile:                    inputfiles/focallength.txt   
    ThroughputBandwidth:         550             
    ThroughputLambdaC:           638        
    IncludeAberrationCorrection: yes     
    AberrationCorrection:
        Type:                    differential
    IncludeFieldDistortion:      yes             
    FieldDistortion:
        Type:                        Polynomial1D
        Source:                      FromFile 
        ConstantCoefficients:        [0.316257210577,  0.066373219688,  0.372589221219]
        ConstantInverseCoefficients: [-0.317143032936, 0.242638513347, -0.459260203502]
        CoefficientsFromFile:        inputfiles/distortioncoefficients.txt
        InverseCoefficientsFromFile: inputfiles/distortioninversecoefficients.txt
\endcode




### <a name="focalPlaneOrientation"></a> FocalPlaneOrientation



The orientation of the focal plane can either be kept constant over the simulations or vary, according to the values in a file provided by the user.

For an angle of 0°, the y-axis of the CCD (with an orientation angle of 0°) points towards the North. A positive angle corresponds to a counterclockwise rotation. Have a look at Fig. 6 for more details.

@image html /images/FocalPlaneCoordinateSystem.png "Figure 6: A schematic overview of the focal plane with 4 CCDs. The optical axis zFP is the blue dot in the middle of the 4 CCDs and points in the positive direction towards the reader. The jitter roll axis zSC is the purple dot, and also points in the positive direction towards the reader.  The focal plane is rotated by the angle \f$\gamma \f$FP w.r.t. to the North direction. The origin of the CCD in the focal plane is defined by its offset (\f$\Delta \f$xCCD, \f$\Delta \f$yCCD) in mm from the centre of the focal plane. It is then rotated by the angle \f$\gamma \f$CCD round its origin."

#### <a name="focalPlaneOrientationSource"></a>FocalPlaneOrientation: Source
<i>Allowed values:</i> "ConstantValue" and "FromFile".

Indicates whether the value of the focal-plane orientation angle must be constant or is allowed to vary over time, according to the values in a user-provided file.



#### <a name="focalPlaneOrientationConstant"></a>FocalPlaneOrientation: ConstantValue
<i>Allowed values:</i> Any

Orientation angle of the focal plane, expressed in degrees, in case the focal-plane orientation must remain the same over the duration of the whole simulation ([FocalPlaneOrientation: Source](#focalPlaneOrientationSource) = ConstantValue).



#### <a name="focalPlaneOrientationFromFile"></a>FocalPlaneOrientation: FromFile
Path of the file with the focal-plane orientation angle as it varies over time, in case the focal-plane orientation angle must be read from a file ([FocalPlaneOrientation: Source](#focalPlaneOrientationSource) = FromFile).





### <a name="plateScale"></a>PlateScale
<i>Allowed values:</i> > 0

Nominal plate scale in arcsec / micron. This value affects the visible FOV of the CCD.




### <a name="focalLength"></a>FocalLength
The focal length can either be kept constant over the simulations or vary, according to the values in a file provided by the user.



#### <a name="focalLengthSource"></a>FocalLength: Source

<i>Allowed values:</i> "ConstantValue" and "FromFile".

Indicates whether the value of the focal length must be constant or is allowed to vary over time, according to the values in a user-provided file.



#### <a name="focalLengthConstant"></a>FocalLength: ConstantValue
<i>Allowed values:</i> > 0

Focal length as recovered from the Zemax model, expressed in m.

Orientation angle of the focal plane, expressed in degrees, in case the focal length must remain the same over the duration of the whole simulation ([FocalLength: Source](#focalLengthSource) = ConstantValue).



#### <a name="focalLengthFromFile"></a>FocalLength: FromFile
Path of the file with the focal plane as it varies over time, in case the focal length must be read from a file ([FocalLength: Source](#focalLengthSource) = FromFile).



### <a name="throughputBandwidth"></a>ThroughputBandwidth
<i>Allowed values:</i> > 0

FWHM of the throughput passband, expressed in nm.



### <a name="throughputLambdaC"></a>ThroughputLambdaC
<i>Allowed values:</i> > 0

Central wavelength of the throughput passband, expressed in nm.




### <a name="includeFieldDistortion"></a>IncludeFieldDistortion
<i>Allowed values:</i> "yes" and "no"

Indicates whether or not the field distortion must be taken into account.



### <a name=includeAberrationCorrection></a>IncludeAberrationCorrection
<i>Allowed values:</i> "yes" and "no"

Indicates whether or not to apply the aberration correction to all star positions in the [star catalogue](#starCatalogue). 



### <a name=aberrationCorrection></a>AberrationCorrection

The calculation of the aberration correction is an approximation based on a circular earth orbit around the sun and does not take the Lissajous orbit of the satellite around L2 into account. We do calculate the aberration, however, which takes into account the aberration correction done for the spacecraft pointing.



#### <a name=aberrationCorrectionType></a>IncludeAberrationCorrection: Type
<i>Allowed values:</i> "differential" and "absolute"

Indicates whether to apply either differential or absolute aberration correction (if ([IncludeAbberationCorrection](#includeAberrationCorrection) = yes).



### <a name="fieldDistortion"></a>FieldDistortion

The field distortion is represented by either a 1D or a 2D polynomial, the coefficients of which must be kept constant over the simulations or vary, according to the values in a file provided by the user.


#### <a name="fieldDistortionType"></a>FieldDistortion: Type
<i>Allowed values:</i> "Polynomial1D"

Indicates that the field distortion is calculated by means of a 1D polynomial.

Assumed is that the distortion is radial and only depends on the radial distance _r_ in the focal plane.  This type of distortion can be modelled with a pincushion or barrel distortion:

\f[P(r) = k_1 \cdot r^3 + k_2 \cdot r^5 + k_3 \cdot x^7.\f]


#### <a name="fieldDistortionSource"></a>FieldDistortion: Source

Indicates whether the coefficients for the polynomial describing the field distortion must be constant or are allowed to vary over time, according to the values in a user-provided file.



#### <a name="fieldDistortionCoefficientsConstant"></a>FieldDistortion: ConstantCoefficients

Coefficients for the 1D polynomial that converts the normalised undistorted pixel coordinates (i.e. pixel coordinates divided by the focal length in pixels) to the distortion, expressed in normalised pixel coordinates, in case the coefficients must remain the same over the duration of the whole simulation ([FieldDistortion: Source](#fieldDistortionSource) = ConstantValue).



#### <a name="fieldDistortionInverseCoefficientsConstant"></a>FieldDistortion: ConstantInverseCoefficients

Inverse coefficients for the 1D polynomial that converts the normalised distorted pixel coordinates (i.e. pixel coordinates divided by the focal length in pixels) to the (negative) distortion, expressed in normalised pixel coordinates, in case the inverse coefficients must remain the same over the duration of the whole simulation ([FieldDistortion: Source](#fieldDistortionSource) = ConstantValue).



#### <a name="fieldDistortionCoefficientsFromFile"></a>FieldDistortion: CoefficientsFromFile

Coefficients for inverse polynomial of the polynomial describing the field distortion, in case the coefficients must be read from a file ([FieldDistortion: Source](#fieldDistortionSource) = FromFile).



#### <a name="fieldDistortionInverseCoefficientsConstant"></a>FieldDistortion: InverseCoefficientsFromFile

Inverse coefficients for inverse polynomial of the polynomial describing the field distortion, in case the coefficients must be read from a file ([FieldDistortion: Source](#fieldDistortionSource) = ConstantValue).

---




<!-- ************** -->
<!-- PSF Parameters -->
<!-- ************** -->

## <a name="psfParameters"></a>PSF Parameters

The <b>PSF</b> block of the configuration file contains all the information that is specific to the PSF.  The structure of this block is the following:


\code{.yaml}
PSF:

    Model:                       MappedGaussian 
    MappedGaussian:                             
      Sigma:                     0.638     
      NumberOfPixels:            8   
      ChargeDiffusionStrength:     0.2
      IncludeChargeDiffusion:      no
      IncludeJitterSmoothing:      no
    MappedFromFile:                             
      Filename:                  inputfiles/psf.hdf5 
      DistanceToOA:              -1       
      RotationAngle:             -1         
      NumberOfPixels:            8
      ChargeDiffusionStrength:     0.2
      IncludeChargeDiffusion:      no
      IncludeJitterSmoothing:      no
    AnalyticGaussian:
      Sigma00:                   1.0
      SigmaX18:                  5.0
      SigmaY18:                  2.0
    AnalyticNonGaussian:
      ParameterFileName:         inputfiles/parameters.txt
      Sigma:  
            Source:                  ConstantValue 
            ConstantValue:           0.5
            FromFile:                inputfiles/sigmaPSF.txt
\endcode




### <a name="psfModel"></a>Model
<i>Allowed values:</i> "MappedGaussian", "MappedFromFile", "AnalyticGaussian", and "AnalyticNonGaussian

Indicates whether to use a Gaussian PSF, to read the PSF from an HDF5 file, or to use an analytical model (Gaussian or non-Gaussian):

- MappedGaussian: the PSF is a circular Gaussian, the size of which does not change over the FOV;
- MappedFromFile: the PSF is selected from an HDF5 file with pre-computed PSFs, based on the angular distance to the optical axis;
- AnalyticGaussian: the PSF is an elongated Gaussian (the symmetry axes being parallel to the x- and y-axis), for which the width and the height are given at the centre of the FOV and at 18 degrees from the optical axis;
- AnalyticNonGaussian: the PSF is an analytical non-Gaussian model, the parameters of which are stored in a separate file.



### <a name="mappedGaussian"></a>MappedGaussian

The PSF is a circular Gaussian, the size of which does not change over the FOV.

#### <a name="gaussSigma"></a>MappedGaussian: Sigma
<i>Allowed values:</i> > 0, only required if a Gaussian PSF must be used ([psfModel](#Model) = MappedGaussian)

Width (\f$\sigma \f$) of the two-dimensional Gaussian PSF, expressed in pixels.  This Gaussian PSF does not vary in size over the FOV.




#### <a name="gaussNumPixels"></a>MappedGaussian: NumberOfPixels
<i>Allowed values:</i> > 0, only required if a Gaussian PSF with a fixed size over the FOV must be used ([Model](#psfModel) = MappedGaussian)

Number of pixels (in both directions) for which the Gaussian PSF must be generated.



#### <a name="chargeDiffusionStrength"></a>MappedGaussian: ChargeDiffusionStrength

<i>Allowed values:</i> \f$\ge \f$ 1 / [SubPixels](#numSubPixelsl)

Charge diffusion has been modelled by a convolution with a Gaussian diffusion kernel, of which this is the standard deviation.



#### <a name="inclChargeDiffusion"></a>MappedGaussian: IncludeChargeDiffusion

<i>Allowed values:</i> "yes" and "no"

Indicates whether or not to include charge diffusion.


#### <a name="inclJitterSmoothing"></a>MappedGaussian: IncludeJitterSmoothing

<i>Allowed values:</i> "yes" and "no"

Indicates whether or not to include jitter smoothing.  This is implemented as charge diffusion with a kernel width of 0.5 sub-pixels and alleviates the problem of jitter discontinuities.  Only applicable in case [IncludeChargeDiffusion](#inclChargeDiffusion) = False.



### <a name="mappedFromFile"></a>MappedFromFile

The PSF is selected from an HDF5 file with pre-computed PSFs, based on the angular distance to the optical axis.



#### <a name="psfFilename"></a>MappedFromFile: Filename
<i>Allowed values:</i> only required if a pre-computed PSF must be used ([psfModel](#Model) = MappedFromFile)

Path to the file, relative to the [project location](#projectLocation), holding the location independent [pre-computed PSF](#psfFile).



#### <a name="psfDistance"></a>MappedFromFile: DistanceToOA
<i>Allowed values:</i> -1 for automatic calculation, \f$\ge \f$ 0 to use the input value; only required if a pre-computed PSF must be used ([Model](#psfModel) = MappedFromFile)

In case a positive value is given the input value will be used for the angular distance to the optical axis.

In case a negative value is given, the angular distance to the optical axis will be calculated automatically.




#### <a name="psfRotation"></a>MappedFromFile: RotationAngle
<i>Allowed values:</i> Any, only required if a pre-computed PSF must be used ([Model](#psfModel) = MappedFromFile)

Arbitrary rotation angle of the PSF, expressed in degrees and measured counterclockwise.




#### <a name="psfNumPixels"></a>MappedFromFile: NumberOfPixels
<i>Allowed values:</i> > 0, only required if a pre-computed PSF must be used ([Model](#psfModel) = MappedFromFile)

Number of pixels (in both directions) for which the PSF was generated.



#### <a name="chargeDiffusionStrengthFile"></a>MappedFromFile: ChargeDiffusionStrength

<i>Allowed values:</i> \f$\ge \f$ 1 / [SubPixels](#numSubPixels)

Charge diffusion has been modelled by a convolution with a Gaussian diffusion kernel, of which this is the standard deviation.



#### <a name="inclChargeDiffusionFile"></a>MappedFromFile: IncludeChargeDiffusion

<i>Allowed values:</i> "yes" and "no"

Indicates whether or not to include charge diffusion.




#### <a name="inclJitterSmoothingFile"></a>MappedGaussian: IncludeJitterSmoothing

<i>Allowed values:</i> "yes" and "no"

Indicates whether or not to include jitter smoothing.  This is implemented as charge diffusion with a kernel width of 0.5 sub-pixels and alleviates the problem of jitter discontinuities.  Only applicable in case [IncludeChargeDiffusion](#inclChargeDiffusionFile) = False.




### <a name="analyticGaussian"></a>AnalyticGaussian

The PSF is an elongated Gaussian (the symmetry axes being parallel to the x- and y-axis), for which the width and the height are given at the centre of the FOV and at 18 degrees from the optical axis.



#### <a name=sigma00></a>AnalyticGaussian: Sigma00
<i>Allowed values:</i> > 0, only required if a analytic Gaussian PSF must be used ([Model](#psfModel) = AnalyticGaussian)

Standard deviation of the analytical Gaussian PSF in the x- and y-direction at the optical axis, expressed in pixels.



#### <a name=sigmaX18></a>AnalyticGaussian: SigmaX18
<i>Allowed values:</i> > 0, only required if a analytic Gaussian PSF must be used ([Model](#psfModel) = AnalyticGaussian)

Standard deviation of the analytical PSF in the x-direction at 18 degrees from the optical axis, expressed in pixels.



#### <a name=sigmaY18></a>AnalyticGaussian: SigmaY18
<i>Allowed values:</i> > 0, only required if a analytic Gaussian PSF must be used ([Model](#psfModel) = AnalyticGaussian)

Standard deviation of the analytical PSF in the y-direction at 18 degrees from the optical axis, expressed in pixels.


### <a name="analyticNonGaussian"></a>AnalyticNonGaussian

The PSF is an analytical non-Gaussian model, the parameters of which are stored in a separate file.



#### <a name=analyticPsfFile></a>AnalyticNonGaussian: ParameterFileName
<i>Allowed values:</i> only required if a analytic non-Gaussian PSF must be used ([Model](#psfModel) = AnalyticNonGaussian)

Path to the file, relative to the [project location](#projectLocation), holding the parameters characterising the analytical model.


#### <a name=analyticPsfSigma></a>AnalyticNonGaussian: Sigma

The width of the analytic non-Gaussian PSF, equal to sigma for a Gaussian PSF (expressed in pixels), can either be kept constant over the simulations or vary, according to the values in a file provided by the user. 


##### <a name="analyticPsfSigmaSource"></a>AnalyticNonGaussian: Sigma: Source

<i>Allowed values:</i> "ConstantValue" or "FromFile", only required if a analytic non-Gaussian PSF must be used ([Model](#psfModel) = AnalyticNonGaussian)

Indicates whether the value of the width of the analytic non-Gaussian PSF must be constant or is allowed to vary over time, according to the values in a user-provided file.
 

##### <a name="analyticPsfSigmaConstant"></a>AnalyticNonGaussian: Sigma: ConstantValuee

<i>Allowed values:</i> > 0, only required if a analytic non-Gaussian PSF must be used ([Model](#psfModel) = AnalyticNonGaussian) and [Sigma: Source] = ConstantValue

Width of the analytic non-Gaussian PSF, equal to sigma for a Gaussian PSF (expressed in pixels), in case it must remain the same over the duration of the whole simulation ([Sigma: Source] = ConstantValue).



##### <a name="analyticPsfSigmaConstant"></a>AnalyticNonGaussian: Sigma: FromFile

<i>Allowed values:</i> only required if a analytic non-Gaussian PSF must be used ([Model](#psfModel) = AnalyticNonGaussian) and [Sigma: Source] = ConstantValue

Path to the file, relative to the [project location](#projectLocation), with the width of the analytic non-Gaussian PSF, in case it must be read from a file ([Sigma: Source] = FromFile).

Charge diffusion can be switched on/off by reading the appropriate parameter file.  One will be provided with the parameters with charge diffusion and one without. 

---





<!-- ************** -->
<!-- FEE Parameters -->
<!-- ************** -->

## <a name="feeParameters"></a>FEE Parameters

The <b>FEE</b> block of the configuration file contains all the information that is specific to the front-end electronics (FEE).  The structure of this block is the following:

\code{.yaml}
FEE:

    NominalOperatingTemperature: 210.15
    Temperature:                 Nominal
    TemperatureFileName:         inputfiles/feeTemperature.txt     
    ReadoutNoise:                40.5         
    Gain:          
    		RefValueLeft:        11.1
    		RefValueRight:       11.1
    		ThreeSigma:          0.0    
    		AllowedDifference:   -100    		
    ElectronicOffset:           
    		RefValue:            100
    		Stability:           1
    OverAndUnderShoot:
        Strength:                    0.003867
        DecaySpeed:                  0.755
        DecayRate:                   1.277
        Range:                       5
    IncludeOverAndUnderShoot:        no

\endcode




### <a name="nominalTempFEE"></a>NominalOperatingTemperature
<i>Allowed values:</i> > 0

Nominal operating temperature of the FEE, expressed in Kelvin.



### <a name="tempFEE"></a>Temperature
<i>Allowed values:</i>"Nominal" or "FromFile"

Indicates whether the temperature of the FEE should be fixed at the nominal operating temperature or temperature variations should be read from a file.



### <a name=tempFileFEE></a>TemperatureFileName
Path to the file, relative to the [project location](#projectLocation), holding the location of the file with the temperature variations of the FEE.



### <a name=readoutNoiseFEE></a>ReadoutNoise

<i>Allowed values:</i> \f$\ge \f$ 0

Mean readout noise of the FEE, expressed in e<sup>-</sup>/pixel.  This is the same for both ADCs.



### <a name=gainFEE></a>Gain

The actual gain for the FEE will be different for both ADCs.  A reference value is given for ADC1 and ADC2, and the difference should not exceed the [specified allowed difference](#gainAllowedDiffFEE).



#### <a name=gainRefValueLeftFEE></a>Gain: RefValueLeft

<i>Allowed values:</i> > 0

Reference value of the gain of ACD1 of the FEE at its [nominal operating temperature](#nominalTempFEE), expressed in ADU/µV. 


#### <a name=gainRefValueRightFEE></a>Gain: RefValueRight

<i>Allowed values:</i> > 0

Reference value of the gain of ACD2 of the FEE at its [nominal operating temperature](#nominalTempFEE), expressed in ADU/µV.  



#### <a name=gainAllowedDiffFEE></a>Gain: AllowedDifference

<i>Allowed values:</i> \f$\in \f$ [0,100]

Percentage of the reference values for the gain of ADC1 and ADC2 that indicates the maximum allowed difference between these gain values.



#### <a name=gainStabilityFEE></a>Gain: Stability

<i>Allowed values:</i> Any

Change in gain (for both ADCs) with temperature deviations from the nominal operating temperature, expressed in ADU/µV/K.



### <a name=electronicOffset></a>ElectronicOffset

The electronic offset or bias level is added  to the digital signal in order to avoid negative readout values. The electronic offset can be measured in a pre-scan strip, which essentially consists of a few additional rows of the CCD. These rows only contain the electronic offset and the readout noise. This pre-scan strip consisting of [NumPreScanRows](#numPreScanRows) rows will be stored in the output file.  This is the same for both ADCs.


#### <a name=electronicOffsetRefValue></a>ElectronicOffset: RefValue

<i>Allowed values:</i> \f$\ge \f$ 0

Electronic offset or bias level at the nominal operating temperature of the FEE, expressed in ADU.



#### <a name=electronicOffsetStability></a>ElectronicOffset: Stability

<i>Allowed values:</i> Any

Change in electronic offset (for both ADCs) with temperature deviations from the nominal operating temperature, expressed in ADU/pixel/K.



### <a name="overAndUnderShoot"></a>OverAndUnderShoot

Over-/undershoot has been noticed in F-FEE measurements.  Looking at the content of a readout register at any given time, the charges in any pixels will affect the next pixels in the readout register, further away from the readout register, e.g. at a distance \f$\Delta x \f$.  For a difference in signal between two such pixels, \f$\Delta S \f$, the
induced over-/undershoot will be

\f[a \cdot \Delta S \cdot \exp{(- \lambda \cdot \Delta x^b)}.\f]

Both detector halves are treated independently.



#### <a name=""></a>OverAndUnderShoot: Strength

<i>Allowed values:</i> > 0

Parameter \f$a \f$ in the formula above.



#### <a name=""></a>OverAndUnderShoot: DecayRate

<i>Allowed values:</i> > 0

Parameter \f$\lambda \f$ in the formula above.



#### <a name=""></a>OverAndUnderShoot: DecaySpeed

<i>Allowed values:</i> > 0

Parameter \f$b \f$ in the formula above.



#### <a name=""></a>OverAndUnderShoot: Range

<i>Allowed values:</i> > 0

Maximum distance \f$\Delta x \f$ over which a pixel in the readout register exerts over-/undershoot on other pixels in the readout register, further away from the readout electronics.



### <a name="inclOverAndUnderShoot"></a>MappedGaussian: IncludeOverAndUnderShoot

<i>Allowed values:</i> "yes" and "no"

Indicates whether or not to include F-FEE over-/undershoot.  Only applicable in case [GroupID](#groupID) = Fast.

---





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
    TimeShift:                   0.0
    PixelSize:                   18      
    BFE:
        CoefficientsFileName:
    Gain:                        
        RefValueLeft:        1.80
        RefValueRight:       1.80
        AllowedDifference:   15.0   
        Stability:           -0.004    
    QuantumEfficiency:
        MeanQuantumEfficiency:       0.5985
        MeanAngleDependency:         1.01
    Polarization:
    		ExpectedValue:       0.989      
    RelativeTransmissivity:
        Coefficients:           [4.18e-2, -5.65e-5, 2.37e-7]
        RadiusFOV:              18.8908
        ExpectedValue:          0.920
    Contamination:
    		ParticulateContaminationEfficiency:  0.98
    		MolecularContaminationEfficiency:    0.0566
    DarkSignal:
      		DarkCurrent:                  1.2
      		DSNU:                         10.0
      		Stability:                    5.0
    FullWellSaturation:          1000000        
    DigitalSaturation:           65535          
    ReadoutNoise:                28
    SerialTransferTime:              340
    ParallelTransferTime:            110
    ParallelTransferTimeFast:        90
    ReadoutMode:
       ReadoutMode:                  Nominal
       Partial:
          FirstRowReadout:           0
          NumRowsReadout:            4510
    ElectronicOffset:            100 
    FlatfieldNoiseRMS:           0.010      
    CTI:
    		Model:			 Simple
    		Simple:
    		   CTEMean:        0.99999
    	      Short2013:
    	          Beta:		0.37
    	          Temperature:	203.0
    	          NumTrapSpecies:[9.8, 3.31, 1.56, 13.24]
    	          TrapDensity:   
                        BOL:    [0.0, 0.0, 0.0, 0.0]
                        EOL:    [2.46e-20, 1.74e-22, 7.05e-23, 2.45e-23]
    	          ReleaseTime:   [2.37e-4, 2.43e-2, 2.03e-3, 1.40e-1]
    NominalOperatingTemperature: 203.15
    Temperature:                 Nominal
    TemperatureFileName:         inputfiles/ccdTemperature.txt     
    IncludeFlatfield:                 no            
    IncludeDarkSignal:                yes 
    IncludePhotonNoise:               yes            
    IncludeReadoutNoise:              yes            
    IncludeCTIeffects:                yes            
    IncludeOpenShutterSmearing:       yes            
    IncludeRelativeTransmissivity:    yes   
    IncludePolarization:              yes
    IncludeParticulateContamination:  yes
    IncludeMolecularContamination:    yes
    IncludeQuantumEfficiency:         yes
    IncludeConvolution:               yes            
    IncludeFullWellSaturation:        yes            
    IncludeDigitalSaturation:         yes      
    IncludeQuantisation:              yes      
\endcode




### <a name="position"></a>Position
<i>Allowed values:</i> \f$\in \f$ [1, 2, 3, 4, Custom]

The CCD position can be used to select a specific pre-defined CCD or a custom one.

The pre-defined CCD positions are shown in the figures below.

@image html "/images/CCD Array Configuration - Normal Camera.png" "Figure 7: Layout of the CCDs for the normal camera's."
@image html "/images/CCD Array Configuration - Fast Camera.png" "Figure 8: Layout of the CCDs for the fast camera's."

<!-- Note that we now use 1, 2, 3, and 4 rather than A, B, C, D, for the normal cameras as well as for the fast ones.

|In the past|Now |
|---|---|
| A  | 3  |
| B  | 2  |
| C  | 4  |
| D  | 1  | -->

When you specify [Position](#position)=Custom, the origin offset ([OriginOffsetX](#originOffsetX) and [OriginOffsetY](#originOffsetY)), the [orientation](#ccdOrientation), [number of rows](#ccdNumRows) and [columns](#ccdNumColumns), and the [first exposed row](#firstRowExposed) of the CCD are read from the configuration parameters in the CCD block.

In case a pre-defined position is used, these configuration parameters are read from the [CCDPositions](#ccdPositions) block (see below).




### <a name="orginOffsetX"></a>OriginOffsetX
<i>Allowed values:</i> Any

Offset of the CCD origin from the centre of the optical plane (i.e. the intersection of the optical axis with the focal plane) in the x-direction, expressed in mm. The origin of the CCD is defined as the point where the readout register is located. See Fig. 6 for more details (\f$\Delta x_{\rm CCD} \f$).

This parameter is only used when the [Position](#position)=Custom.




### <a name="originOffsetY"></a>OriginOffsetY
<i>Allowed values:</i> Any

Offset of the CCD origin from the centre of the optical plane (i.e. the intersection of the optical axis with the focal plane) in the y-direction, expressed in mm. The origin of the CCD is defined as the point where the readout register is located. See Fig. 6 for more details (\f$\Delta y_{\rm CCD} \f$).

This parameter is only used when the [Position](#position)=Custom.


### <a name="ccdOrientation"></a>Orientation
<i>Allowed values:</i> Any

Orientation angle of the CCD w.r.t. the orientation of the focal plane, measured counterclockwise and expressed in degrees. This rotation is performed around the offset origin of the CCD. See Fig. 6 for more details (\f$\gamma_{\rm CCD} \f$).

This parameter is only used when the [Position](#position)=Custom.




### <a name="ccdNumColumns"></a>NumColumns
<i>Allowed values:</i> > 0

Number of pixels of the CCD in the x-direction (i.e. number of columns).

This parameter is only used when the [Position](#position)=Custom.




### <a name="ccdNumRows"></a>NumRows
<i>Allowed values:</i> > 0

Number of pixels of the CCD in the y-direction (i.e. number of rows).

This parameter is only used when the [Position](#position)=Custom.



### <a name="firstRowExposed"></a>FirstRowExposed
<i>Allowed values:</i> > 0

Row index of the first row in the CCD that is illuminated (the row closest to the readout register is row 0).

This parameter is only used when the [Position](#position)=Custom.



### <a name="timeShift"></a>TimeShift

<i>Allowed values:</i> > 0

Time shift between the readout of the CCDs [s].  Will only be used if [Position](#position)=Custom.




### <a name="pixelSize"></a>PixelSize

<i>Allowed values:</i> > 0

Nominal pixel size, expressed in micron.

### <a name=BFE></a>BFE

The brighter-fatter effect (BFE) is modelled following the method proposed in [Guyonnet et al. 2015](https://arxiv.org/abs/1501.01577).

#### <a name=coefficientsBFE></a>BFE: CoefficientsFileName

Path to the HDF5 file comprising the coefficients _a_ for the BFE.

### <a name=gainFEE></a> Gain

The actual gain for the CCD will be different for both detector halves.  A reference value is given for both halves, and the difference should not exceed the [specified allowed difference](#gainAllowedDiffCCD).



#### <a name=gainRefValueLeftCCD></a>Gain: RefValueLeft

<i>Allowed values:</i> > 0

Reference value of the gain of the left- and right-hand side of the CCD at its [nominal operating temperature](#nominalTempCCD), expressed in µV/e<sup>-</sup>.


#### <a name=gainRefValueRightFEE></a>Gain: RefValueRight

<i>Allowed values:</i> > 0

Reference value of the gain of ACD2 of the FEE at its [nominal operating temperature](#nominalTempFEE), expressed in µV/e<sup>-</sup>.



#### <a name=gainAllowedDiffCCD></a>Gain: AllowedDifference

<i>Allowed values:</i> \f$\in \f$ [0,100]

Percentage of the reference values for the gain of the left- and right-hand side of the CCD that indicates the maximum allowed difference between these gain values.


#### <a name=gainStabilityCCD></a>Gain: Stability

<i>Allowed values:</i> Any

Change in gain (for both CCD halves) with temperature deviations from the nominal operating temperature, expressed in µV/e<sup>-</sup>/K.



### <a name=quantumEfficiency></a>QuantumEfficiency

Quantum efficiency is the ratio of the number of collected electrons to the number of incident photons, considering the passband and the spectral energy distribution of the stars given the [Fluxm0](#fluxm0) parameter and the magnitude of the stars in the [star catalogue](#starCatalogue).



<!-- #### <a name="quantumEfficiencyRelRefEfficiency"></a>QuantumEfficiency: RelativeRefEfficiency
<i>Allowed values:</i> \f$\in \f$ [0,1]

Relative efficiency due to angle dependency of the quantum efficiency at the reference angle.



#### <a name="quantumEfficiencyRefAngle"></a>QuantumEfficiency: RefAngle
<i>Allowed values:</i> Any

Reference angle for the throughput efficiency due to the quantum efficiency, expressed in degrees. -->



#### <a name="quantumEfficiencyMean"></a>QuantumEfficiency: MeanQuantumEfficiency
<i>Allowed values:</i> \f$\in \f$ [0,1]

Mean throughput efficiency due to quantum efficiency (i.e. the mean over all pixels of one detector).



#### <a name="angleDependencyQE"></a>QuantumEfficiency: MeanAngleDependency
<i>Allowed values:</i> > 0

Mean efficiency caused by the angle dependency of the quantum efficiency.



### <a name=polarization></a>Polarization

Optical elements induce a preferred direction for the propagation of light.  This effect is called polarisation.


<!-- #### <a name="polarizationEfficiency"></a>Polarization: Efficiency
<i>Allowed values:</i> \f$\in \f$ [0,1]

Throughput efficiency due to the polarisation at the given reference angle.



#### <a name="PolarizationRefAngle"></a>Polarization: RefAngle
<i>Allowed values:</i> Any

Reference angle for the throughput efficiency due to the polarisation, expressed in degrees. -->



#### <a name="polarizationExpectedValue"></a>Polarization: ExpectedValue
<i>Allowed values:</i> \f$\in \f$ [0,1]

Expected value of the throughput efficiency due to polarisation (i.e. the mean over all pixels of one detector).  Currently no information on the angle dependency of polarisation is available and hence this value will be used for the whole FOV, until further notice.



### <a name="relativeTransmissivity"></a>RelativeTransmissivity

On top of the (time-dependent) [transmission efficiency](#transmissionEfficiency), the overall relative transmissivity should be taken into account.  This decrease in efficiency with distance to the optical axis, comprises the following contributions:

* natural vignetting (brightness attenuation towards the edges of the FOV, introduced by the view factor of the entrance pupil);
* mechanical vignetting (due to the undersized mask at the entrance pupil), incl. total blockage of all incoming radiation beyond the edge of the FOV;
* glass absorption + anti-reflective coating;



#### <a name="relativeTransmissivityCoefficients"></a>RelativeTransmissivity: Coefficients
<i>Allowed values:</i> > 0

Coefficients \f$k_1, k_2, k_3 \f$ for the polynomial that converts the distance from the optical axis, \f$\theta \f$ (expressed in degrees), to the variation in the overall relative transmissivity (expressed in percentage):

\f[P(\theta) = k_1 \cdot \theta^2 + k_2 \cdot \theta^4 + k_3 \cdot \theta^6.\f]



#### <a name="mechanicalVignettingRadiusFOV"></a>RelativeTransmissivity: RadiusFOV
<i>Allowed values:</i> > 0

Radius of the FOV, expressed in degrees.  Beyond this radius all incoming flux (apart from the cosmic hits) is shielded off.



#### <a name="relativeTransmissivityExpectedValue"></a>RelativeTransmissivity: ExpectedValue
<i>Allowed values:</i> \f$\in \f$ [0,1]

Expected value of the throughput efficiency due to the overall relative transmissivity (i.e. the mean over all pixels of one detector, within the FOV).


### <a name=contamination></a>Contamination

The contribution to contamination is two-fold:

* particulate contamination is the unintended presence of particles on (optical) surfaces, which leads to straylight and influences the effciency;
* molecular contamination is caused chiefly by outgassing of materials in the first phase of the mission and affects all surfaces (the down side of L2 and the CCD will be affected the most).

#### <a name="particulateContamination"></a>Contamination: ParticulateContaminationEfficiency
<i>Allowed values:</i> \f$\in \f$ [0,1]

Throughput efficiency due to particulate contamination.



#### <a name="molecularContamination"></a>Contamination: MolecularContaminationEfficiency
<i>Allowed values:</i> \f$\in \f$ [0,1]

Throughput efficiency due to molecular contamination.



### <a name="darkSignal"></a>DarkSignal

Dark signal is the relatively small electric current that is generated in the CCD when no outside radiation is entering the device.

#### <a name="darkCurrent"></a>DarkSignal: DarkCurrent
<i>Allowed values:</i> > 0

Dark current, expressed in e<sup>-</sup> / s.  This is the nominal value of the dark signal.



#### <a name="dsnu"></a>DarkSignal: DSNU
<i>Allowed values:</i> \f$\in \f$ [0,100]

Dark signal non-uniformity, expressed as a percentage of the [dark current](#darkCurrent).  This is the systematic (fixed-pattern) deviation of a pixel's dark current from its nominal value.



#### <a name="darkCurrentStability"></a>DarkSignal: Stability
<i>Allowed values:</i> \f$\ge \f$ 0

Temperature stability of the dark current, expressed in in e<sup>-</sup> / K / s.





### <a name="fullWellSaturation"></a>FullWellSaturation
<i>Allowed values:</i> > 0
     
Full-well saturation limit of a single CCD pixel, expressed in e<sup>-</sup> / pixel. If a pixels receives more electrons than its full-well saturation limit, the additional electrons flow evenly distributed in positive and negative charge-transfer direction, a phenomenon called <i>blooming</i>. The electrons reaching the edge of the CCD will not be detected..



### <a name="digitalSaturation"></a>DigitalSaturation
<i>Allowed values:</i> > 0

Digital saturation limit of the CCD to which pixel values are topped off, expressed in ADU / pixel. This value depends on the A/D convertor of the detector. For a 16-bit convertor, the digital saturation limit is 65536 ADU.

The gain of the front-end electronics and detector should be such that the [full-well saturation](#fullWellSaturation) results in values below the digital saturation limit.


     

### <a name="readoutNoise"></a>ReadoutNoise
<i>Allowed values:</i> \f$\ge \f$ 0

Mean readout noise of the detector, expressed in e<sup>-</sup>.

Readout noise occurs due to the imperfect nature of the CCD amplifiers. When the electrons are transferred to the amplifier, the induced voltage is measured. However, this measurement is not perfect, but gives a value which is on average too high by an amount of the readout noise, with the squareroot of the readout noise as standard deviation (we add the readout noise of the FEE and the CCD in quadrature).



### <a name="serialTransferTime"></a>SerialTransferTime
<i>Allowed values:</i> ≥ 0

Time required to shift the content of the readout register over one pixel, towards the output node.  This is not only relevant for the image area but also for the serial pre-scan (i.e. bias register map) and the serial over-scan (which has not been implemented).



### <a name="parallelTransferTime"></a>ParallelTransferTime
<i>Allowed values:</i> \f$\ge \f$ 0

Time required to shift the charges one row down (towards the readout register) in case the readout register will be read out by the FEE.

The difference with [ParallelTransferTimeFast](#parallelTransferTimeFast) is due to two delay parameters recommended by Teledyne e2v for clock settling.  Settling times are to ensure one clock has reached its low level before another clock starts to rise.  It ensures good charge transfer, as without it you could cause charge to be lost.  This is particularly the case when dealing with the (parallel) transfer into the readout register, as the rise and fall times of the image clocks are a lot slower than those of the register clocks.



### <a name="parallelTransferTimeFast"></a>ParallelTransferTimeFast
<i>Allowed values:</i> \f$\ge \f$ 0

Time required to shift the charges one row down (towards the readout register) in case the readout register will not be read out by the FEE.  In that case clock settling is not needed, hence the difference with [ParallelTransferTime](#parallelTransferTime).




### <a name="flatfieldPtPNoise"></a>FlatfieldNoiseRMS
<i>Allowed values:</i> \f$\ge \f$ 0

Local PRNU (Pixel Response Non-Uniformity), defined as the standard deviation in the signal level (\f$S_i \f$), divided by the mean signal (\f$ \overline{S}\f$):

\f[PRNU = \frac{\sqrt{\frac{1}{N} \sum_{i=1}^{N} (S_i - \overline{S})^2}}{\overline{S}}\f]



### <a name=cti></a>CTI

Because of detector defects, electrons can get trapped in the readout process. The trapped charge ends up getting dissociated from its original pixel and eventually gets released into another pixel. The result is that the original image gets smeared out in the direction away from the readout amplifier (visible in the appearance of "charge trails"). This is known as imperfect CTE (Charge-Transfer Efficiency) or alternatively as CTI (Charge-Transfer Inefficiency).

The charge trails impact photometry, noise, and astrometry of sources. CTI removes flux from the central pixel and thus degrades the expected S/N for an observation. CTI trails bias measurements of source along the trail direction, which can severely impact high-precision astrometry.


  
#### <a name="CTImodel"></a>CTI: Model

<i>Allowed values:</i> "Simple" and "Short2013"

PlatoSim3 offers two implementation of the CTI:

<ul>
	<li>The simple implementation ("Simple") assumes that for each row transfer (in the direction of the readout register) a fraction of the charge is not transferred to the next row, but stays behind.  It will be released by later row transfers.</li>
	<li>A more sophisticated implementation ("Short2013") is based on <a href="http://mnras.oxfordjournals.org/content/430/4/3078.full.pdf">Short et al., MNRAS 430, 3078-3085 (2013)</a>, in which only parallel readout is taken into account.</li>
</ul>



#### <a name="simpleCTI"></a>CTI: Simple
The simple implementation ("Simple") assumes that for each row transfer (in the direction of the readout register) a fraction of the charge is not transferred to the next row, but stays behind.  It will be released by later row transfers.

  
##### <a name="MeanCTE"></a>CTI: Simple: MeanCTE
<i>Allowed values:</i> \f$\in \f$ [0,1]

Mean charge-transfer efficiency (CTE) of the detector.  The fraction of the charge that is successfully transferred from one row to the next row is expressed by this parameter.


#### <a name=ShortCTI></a>CTI: Short
A more sophisticated implementation ("Short2013") is based on <a href="http://mnras.oxfordjournals.org/content/430/4/3078.full.pdf">Short et al., MNRAS 430, 3078-3085 (2013)</a>, in which only parallel readout is taken into account.


##### <a name="beta"></a>CTI: Short2013: Beta

Exponent \f$\beta \f$ in Eq. (1) of Short et al. 2013 describing the relationship between the volume of the charge cloud ( \f$V_c \f$), the number of electrons in a pixel (\f$N_e \f$), the full-well capacity in electrons (\f$FWC \f$), and the assumed maximum geometrical volume that electrons can occupy within a pixel (\f$V_g \f$):

\f[\frac{V_c}{V_g} = \left( \frac{N_e}{FWC} \right)^{\beta}.\f]




##### <a name="temperature"></a>CTI: Short2013: Temperature

<i>Allowed values:</i> \f$\ge \f$ 0

Temperature \f$T \f$ that is used to calculated the thermal velocity \f$v_t \f$ of the electrons:

\f[v_t = \frac{3kT}{m_e^{\ast}},\f]

where \f$k \f$ is the Boltzmann constant and \f$m_e^* \f$ is the effective electron mass in silicon, which we approximate by half the free electron rest mass.




##### <a name="numTrapSpecies"></a>CTI: Short2013: NumTrapSpecies

<i>Allowed values:</i> > 0

Number of trap species that is used in the CTI model by Short et al. 2013.



##### <a name="trapDensity"></a>CTI: Short2013: TrapDensity

We assume the trap density for each of the considered trap densities to increase linearly over time (in absence of charge injection).



###### <a name="trapDensityBOL"></a>CTI: Short2013: TrapDensity: BOL

<i>Allowed values:</i> Array holding one non-negative entry per trap species. 

Array holding the trap density \f$n_t \f$ at BOL for each of the considered trap species, expressed in number of traps per pixel.  This is used to calculate the \f$\gamma \f$-value in Eq. (22) of Short et al. 2013.



###### <a name="trapDensityEOL"></a>CTI: Short2013: TrapDensity: EOL

<i>Allowed values:</i> Array holding one non-negative entry per trap species. 

Array holding the trap density \f$n_t \f$ at EOL for each of the considered trap species, expressed in number of traps per pixel.  This is used to calculate the \f$\gamma \f$-value in Eq. (22) of Short et al. 2013.




##### <a name="trapCaptureCrossSection"></a>CTI: Short2013: TrapCaptureCrossSection

<i>Allowed values:</i> Array holding one non-negative entry per trap species.

Array holding the trap capture cross-section \f$\sigma\f$ for each of the considered trap species, expressed in m<sup>2</sup>.  This is used to calculated the \f$\alpha\f$-value in Eq. (22) of Short et al. 2013.  In this formula, the charge transfer time is used as value for \f$t \f$.




##### <a name="releaseTime"></a>CTI: Short2013: ReleaseTime 

<i>Allowed values:</i> Array holding one non-negative entry per trap species.

Array holding the trap release time constants \f$\tau_r\f$ for each of the considered trap species, expressed in seconds.



### <a name="nominalTempCCD"></a>NominalOperatingTemperature

<i>Allowed values:</i> > 0

Nominal operating temperature of the CCD, expressed in Kelvin.



### <a name="tempCCD"></a>Temperature
<i>Allowed values:</i>"Nominal" or "FromFile"

Indicates whether the temperature of the CCD should be fixed at the nominal operating temperature or temperature variations should be read from a file.



### <a name=tempFileCCD></a>TemperatureFileName
Path to the file, relative to the [project location](#projectLocation), holding the location of the file with the temperature variations of the CCD.





### <a name="inclFlatfield"></a>IncludeFlatfield
<i>Allowed values:</i> "yes" and "no"

Indicates whether or not to include the flatfield.



### <a name="inclDarkSignal"></a>IncludeDarkSignal
<i>Allowed values:</i> "yes" and "no"

Indicates whether or not to include dark signal.



### <a name="inclPhotonNoise"></a>IncludePhotonNoise
<i>Allowed values:</i> "yes" and "no"

Indicates whether or not to include photon noise.




### <a name="inclReadoutNoise"></a>IncludeReadoutNoise
<i>Allowed values:</i> "yes" and "no"

Indicates whether or not to include readout noise.




### <a name="inclCTI"></a>IncludeCTIeffects
<i>Allowed values:</i> "yes" and "no"

Indicates whether or not to include CTI effects.




### <a name="inclOpenShutterSmearing"></a>IncludeOpenShutterSmearing
<i>Allowed values:</i> "yes" and "no"

Indicates whether or not to include open-shutter smearing effects.




### <a name="inclRelativeTransmissivity"></a>IncludeRelativeTransmissivity
<i>Allowed values:</i> "yes" and "no"

Indicates whether or not to include the overall relative transmissivity.

        

### <a name="inclPolarization"></a>IncludePolarization
<i>Allowed values:</i> "yes" and "no"

Indicates whether or not to include loss of throughput efficiency due to polarisation.


        

### <a name="inclParticulateContamination"></a>IncludeParticulateContamination
<i>Allowed values:</i> "yes" and "no"

Indicates whether or not to include loss of throughput efficiency due to particulate contamination.


        


### <a name="inclMolecularContamination"></a>IncludeMolecularContamination
<i>Allowed values:</i> "yes" and "no"

Indicates whether or not to include loss of throughput efficiency due to molecular contamination.


        

### <a name="inclQuantumEfficiency"></a>IncludeQuantumEfficiency
<i>Allowed values:</i> "yes" and "no"

Indicates whether or not to include loss of throughput efficiency due to quantum efficiency.


        

### <a name="inclConvolution"></a>IncludeConvolution
<i>Allowed values:</i> "yes" and "no"

Indicates whether or not the sub-pixel map must be convolved with the PSF.  This applies only to the Gaussian and the pre-computed PSF.  When using the analytic PSF, the PSF is always applied (irrespective of the value of this configuration parameter)!




### <a name="inclFullWellSaturation"></a>IncludeFullWellSaturation
<i>Allowed values:</i> "yes" and "no"

Indicates whether or not to apply full-well saturation.




### <a name="inclDigitalSaturation"></a>IncludeDigitalSaturation
<i>Allowed values:</i> "yes" and "no"

Indicates whether or not to apply digital saturation.



### <a name="inclQuantisation"></a>IncludeQuantisation
<i>Allowed values:</i> "yes" and "no"

Indicates whether or not to apply quantisation.  This includes:
	* applying gain (FEE and CCD), hence converting from electrons to ADUs,
	* adding electronic offset,
	* forcing the ADU values to be integers,
	* and applying digital saturation (can be [switched off separately](#inclDigitalSaturation)).



### <a name="writeSubPixelImages"></a>WriteSubPixelImages
<i>Allowed values:</i> "yes" and "no"

Indicates whether or not the sub-pixel images must be written to the HDF5-file.  Use this for a limited number of exposures, as it takes a lot of space.

---





<!-- ******************** -->
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
    NumBiasPrescanRows:          10
    NumBiasPrescanColumns:       5
    NumSmearingOverscanRows:     5
    SubPixels:                   4
\endcode




@image html /images/subField.png "Figure 9: Schematic presentation of the modelled sub-field and the parameters required to define it. The pixel coordinates of the origin of the sub-field relative to the CCD are xs and ys."





### <a name="zeroPointRow"></a>ZeroPointRow
<i>Allowed values:</i> > 0

Row of the origin of the sub-field in the detector, expressed in pixels.  See Fig. 9 for more details ( \f$y_s \f$).




### <a name="zeroPointColumn"></a>ZeroPointColumn
<i>Allowed values:</i> > 0

Column of the origin of the sub-field in the detector, expressed in pixels.  See Fig. 9 for more details (\f$x_s \f$).




### <a name="numColumns"></a>NumColumns
<i>Allowed values:</i> \f$\ge \f$ 8

Number of columns in the sub-field, expressed in pixels.




### <a name="numRows"></a>NumColumns
<i>Allowed values:</i> \f$\ge \f$ 8

Number of rows in the sub-field, expressed in pixels.



### <a name="numPreScanRows"></a>NumBiasPrescanRows
<i>Allowed values:</i> \f$\ge \f$ 0

Number of rows in the pre-scan strip (see Fig. 9), expressed in normal pixel units.  There are two such strips, on either side of the detector image, and they contain the electronic offset and readout noise of the adjacent detector half.

This parameter is configurable (and not fixed to the number of rows in the detector or the sub-field), because we want (1) to avoid the bias maps to take up too much space in the output file and (2) be able to do accurate bias correction for the photometric reduction (thus want to avoid small noisy bias maps).




### <a name="numPreScanColumns"></a>NumBiasPrescanColumns
<i>Allowed values:</i> \f$\ge \f$ 0

Number of columns in the pre-scan strip (see Fig. 9), expressed in normal pixel units.  There are two such strips, on either side of the detector image, and they contain the electronic offset and readout noise of the adjacent detector half.




### <a name="numOverScanRows"></a>NumSmearingOverscanRows
<i>Allowed values:</i> \f$\ge \f$ 0

Number of rows in the over-scan strip (see Fig. 9), expressed in normal pixel units. This strip is located at the top of the sub-field that is modelled in detail and contains the star smearing due to the absence of a shutter. This flux in this strip is also affected by the electronic offset, readout noise, and shot noise. Not included are the PRNU, cosmic hits, and charge-transfer efficiency (CTE).




### <a name="numSubPixels"></a>SubPixels
<i>Allowed values:</i> power of 2 (\f$\le \f$ 128)

Number of sub-pixels per pixel in both directions.

If you want a pixel of 256 x 256 = 65536 sub-pixels you should specify in the configuration file 256. The total number of subpixels per pixel will then be 65536. If you want 256 = 16 x 16 sub-pixels per pixel, then you should specify 16 in the configuration file.

---




<!-- ********************* -->
<!-- Photometry Parameters -->
<!-- ********************** -->

## <a name="photometryParameters"></a>Photometry Parameters

The <b>Photometry</b> block of the configuration file contains all information to run photometry.  The structure of this block is the following:

\code{.yaml}
Photometry:

    IncludePhotometry:               no
    ContaminationRadius:             4
    MaskUpdateInterval:              14.0
    TargetFileName:                  inputfiles/photometryTargets.txt
\endcode

Note that photometry can only be performed if [Model](#psfModel)=AnalyticNonGaussian.


### IncludePhotometry
<i>Allowed values:</i> "yes" and "no"

Indicates whether or not photometry should be run.  Photometry can only be performed if [Model](#psfModel)=AnalyticNonGaussian.



### ContaminationRadius
<i>Allowed values:</i> > 0

Radius [pixels] around a target within which sources are considered contaminants when calculating the photometry for that target.



### MaskUpdateInterval
<i>Allowed values:</i> > 0

Update interval [days] to update the photometry mask.



### TargetFileName

Path of the file comprising the list of targets identifiers (as listed in the [star catalogue](#starCatalogue)) for which to calculate the photometry, relative to the [project location](#projectLocation).

---





<!-- *************** -->
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
    CosmicSeed:                  1494750830
    DarkSignalSeed:              1468838669
\endcode





### ReadOutNoiseSeed
<i>Allowed values:</i> > 0 or -1

Seed for the random-number generator used for the readout noise.

In case a value of -1 is given as input, the computer time at the start of the simulation will be used instead.  That way, the fast-forward of the random generator when using Slurm is no longer needed (which is better for performance reasons).  The actual value that is used, will be written to the output HDF5 file.



### PhotonNoiseSeed
<i>Allowed values:</i> > 0 or -1

Seed for the random-number generator used for the photon noise.

In case a value of -1 is given as input, the computer time at the start of the simulation will be used instead.  That way, the fast-forward of the random generator when when chopping up the simulation in chunks (Slurm) is no longer needed (which is better for performance reasons).



### JitterSeed
<i>Allowed values:</i> > 0 or -1

Seed for the random-number generator used for the jitter.

In case a value of -1 is given as input, the computer time at the start of the simulation will be used instead.  That way, the fast-forward of the random generator when when chopping up the simulation in chunks (Slurm) is no longer needed (which is better for performance reasons).  The actual value that is used, will be written to the output HDF5 file.  To avoid jumps in the power spectrum when using auto-generated jitter values, it is advised to generate the jitter values for the whole simulation beforehand, write these values to a file, and reading in that file when simulating the different chunks.



### FlatFieldSeed
<i>Allowed values:</i> > 0 or -1

Seed for the random-number generator used for the flatfield.

In case a value of -1 is given as input, the computer time at the start of the simulation will be used instead.  That way, the fast-forward of the random generator when when chopping up the simulation in chunks (Slurm) is no longer needed (which is better for performance reasons).



### CTESeed
<i>Allowed values:</i> > 0 or -1

Seed for the random-number generator used for the CTE.

In case a value of -1 is given as input, the computer time at the start of the simulation will be used instead.  That way, the fast-forward of the random generator when when chopping up the simulation in chunks (Slurm) is no longer needed (which is better for performance reasons).



### DriftSeed
<i>Allowed values:</i> > 0 or -1

Seed for the random number generator used for the drift.

In case a value of -1 is given as input, the computer time at the start of the simulation will be used instead.  That way, the fast-forward of the random generator when when chopping up the simulation in chunks (Slurm) is no longer needed (which is better for performance reasons).  To avoid jumps in the power spectrum when using auto-generated drift values, it is advised to generate the drift values for the whole simulation beforehand, write these values to a file, and reading in that file when simulating the different chunks.


### CosmicSeed
<i>Allowed values:</i> > 0 or -1

Seed for the random-number generators for the cosmics.

In case a value of -1 is given as input, the computer time at the start of the simulation will be used instead.  That way, the fast-forward of the random generator when when chopping up the simulation in chunks (Slurm) is no longer needed (which is better for performance reasons).


### DarkSignalSeed
<i>Allowed values:</i> > 0 or -1

Seed for the random-number generators for the dark signal.

In case a value of -1 is given as input, the computer time at the start of the simulation will be used instead.  That way, the fast-forward of the random generator when when chopping up the simulation in chunks (Slurm) is no longer needed (which is better for performance reasons).

---





<!-- ************************** -->
<!-- Content Control Parameters -->
<!-- *************************** -->

## <a name="contentControlParameters"></a>Content Control Parameters

The <b>ControlHDF5Content</b> block of the configuration file contains all the seeds for random-number generation in the simulator.  The structure of this block is the following:

\code{.yaml}
ControlHDF5Content:

    WritePixelMaps:                  yes
    WriteBiasMaps:                   yes
    WriteSmearingMaps:               yes
    WriteThroughputMaps:             yes
    WriteFlatfieldMap:               yes
    WriteSubPixelImages:             no
    WriteStarPositions:              yes
\endcode



### WritePixelImages
<i>Allowed values:</i> "yes" and "no"

Indicates whether or not the pixel maps must be stored in the output file.



### WriteBiasMaps
<i>Allowed values:</i> "yes" and "no"

Indicates whether or not the bias register maps must be stored in the output file.



### WriteSmearingMaps
<i>Allowed values:</i> "yes" and "no"

Indicates whether or not the smearing maps must be stored in the output file.



### WriteThroughputMaps
<i>Allowed values:</i> "yes" and "no"

Indicates whether or not the throughput maps must be stored in the output file.



### WriteFlatfieldMap
<i>Allowed values:</i> "yes" and "no"

Indicates whether or not the flatfield maps must be stored in the output file.



### WriteSubPixelImages
<i>Allowed values:</i> "yes" and "no"

Indicates whether or not the sub-pixel maps must be stored in the output file.  Only do this for a small number of exposures, as this takes up a lot of space.  Note that this only takes into account the effects that are applied until right before the re-binning (see @ref SimulationSteps "here").


### WriteStarPositions
<i>Allowed values:</i> "yes" and "no"

Indicates whether or not the star positions should be stored in pixel and focal plane coordinates in the output file.  This scales with the number of exposures and the number of stars in the sub-field.

---





<!-- ********************************* -->
<!-- Control TCP connection parameters -->
<!-- ********************************* -->

## <a name="controlTcpConnection"></a>ControlTcpConnection

The <b>ControlTcpConnection</b> block in the configuration file is used in case the jitter is read from a network ([JitterSource](#jitterSource) = FromNetwork).  The structure of this block is the following:

\code{.yaml}
ControlTcpConnection:

    SendImagettesToClients:        no
    GetWindowPositionsFromServer:  no

    WindowPositionServerAddress: tcp://localhost:5558
    JitterServerAddress:         tcp://localhost:5559
    ImagetteClientAddress:       tcp://localhost:5560

    WindowPositionSocketTimeout:    100
    JitterSocketTimeout:            100
\endcode



### <a name="sendImagettesToClients"></a>SendImagettesToClients

<i>Allowed values:</i> "yes" and "no"

Indicates whether or not the simulated imagettes should be sent to the client.



### <a name="getWindowPositionsFromServer"></a>GetWindowPositionsFromServer

<i>Allowed values:</i> "yes" and "no"

Indicates whether or not the window positions should be taken from a server and be updated in upcoming simulations.



### <a name="windowPositionServerAddress"></a>WindowPositionServerAddress

Address from which to read the window positions in case this is requested ([GetWindowPositionsFromServer](#getWindowPositionsFromServer) = yes).



### <a name="jitterServerAddress"></a>JitterServerAddress

Address from which to read the jitter positions.



### ImagetteClientAddress

Client address to which to send the simulated imagettes in case this is requested ([SendImagettesToClients](#sendImagettesToClients) = yes).



### WindowPositionSocketTimeout

<i>Allowed values:</i> > 0

Number of seconds of not receiving window positions (from the [window position server](#windowPositionServerAddress)), after which the connection to that server is regarded as stalled / broken.



### JitterSocketTimeout

<i>Allowed values:</i> > 0

Number of seconds of not receiving jitter positions (from the [jitter server](#jitterServerAddress)), after which the connection to that server is regarded as stalled / broken.

---





<!-- ************* -->
<!-- Camera groups -->
<!-- ************* -->

## <a name="cameraGroups"></a>Camera groups

The <b>CameraGroups</b> block in the configuration file is used in case a pre-defined camera groups (1, 2, 3, 4, or Fast) was selected via the [GroupID](#groupID) parameter in the [Telescope](#telescopeParameters) block in the configuration file.  The structure of this block is the following:

\code{.yaml}
CameraGroups:

    AzimuthAngle:            [45.0, 135.0, 225.0, 315.0, 0.0] 
    TiltAngle:               [9.2, 9.2, 9.2, 9.2, 0.0] 
\endcode

Mind you, you are NOT supposed to alter this section of the configuration file!



### AzimuthAngle

Azimuth angle, expressed in degrees, for camera group 1, 2, 3 and 4, and for the fast camera.  Depending on the value of the [GroupID](#groupID) parameter in the [Telescope](#telescopeParameters) block, the appropriate value will be selected from the list.



### TiltAngle

Tilt angle, expressed in degrees, for camera group 1, 2, 3 and 4, and for the fast camera.  Depending on the value of the [GroupID](#groupID) parameter in the [Telescope](#telescopeParameters) block, the appropriate value will be selected from the list.

---





<!-- ************* -->
<!-- CCD Positions -->
<!-- ************* -->

## <a name="ccdPositions"></a>CCD Positions

The <b>CCDPositions</b> block in the configuration file is used in case a pre-defined CCD position (1, 2, 3, or 4) was selected via the [Position](#position) parameter in the [CCD](#ccdParameters) block in the configuration file.  The structure of this block is the following:

\code{.yaml}
CCDPositions:

    OriginOffsetX:                   [-1, -1, -1, -1]
    OriginOffsetY:                   [82.18, 82.18, 82.18, 82.18]
    Orientation:                     [180, 270, 0, 90]
    NumColumns:                      [4510, 4510, 4510, 4510]
    NumRows:                         [4510, 4510, 4510, 4510]
    FirstRowForNormalCamera:         [0, 0, 0, 0]
    FirstRowForFastCamera:           [2255, 2255, 2255, 2255]
    TimeShift:                       [0.0, 6.25, 12.5, 18.75]
\endcode

Mind you, you are NOT supposed to alter this section of the configuration file!



### OriginOffsetX

Offset of the CCD origin from the centre of the optical plane in the x-direction, expressed in mm, for CCD positions 1, 2, 3, and 4.  Depending on the value of the [Position](#position) parameter in the [CCD](#ccdParameters) block, the appropriate value will be selected from the list.



### OriginOffsetY

Offset of the CCD origin from the centre of the optical plane in the y-direction, expressed in mm, for CCD positions 1, 2, 3, and 4.  Depending on the value of the [Position](#position) parameter in the [CCD](#ccdParameters) block, the appropriate value will be selected from the list.



### Orientation

Orientation angle of the CCD w.r.t. the orientation of the focal plane, measured counterclockwise and expressed in degrees, for CCD positions 1, 2, 3, and 4.  Depending on the value of the [Position](#position) parameter in the [CCD](#ccdParameters) block, the appropriate value will be selected from the list.



### NumColumns

Number of pixels of the CCD in the x-direction (i.e. number of columns), for CCD positions 1, 2, 3, and 4.  Depending on the value of the [Position](#position) parameter in the [CCD](#ccdParameters) block, the appropriate value will be selected from the list.



### NumRows

Number of pixels of the CCD in the y-direction (i.e. number of rows), for CCD positions 1, 2, 3, and 4.  Depending on the value of the [Position](#position) parameter in the [CCD](#ccdParameters) block, the appropriate value will be selected from the list.



### FirstRowForNormalCamera

Row index of the first row in the CCD that is illuminated (the row closest to the readout register is row 0), for CCD positions 1, 2, 3, and 4, in case of a normal camera ([GroupID](#groupID)=1, 2, 3, or 4).  For normal cameras, the whole CCD is illuminated.  Depending on the value of the [Position](#position) parameter in the [CCD](#ccdParameters) block, the appropriate value will be selected from the list.



### FirstRowForFastCamera

Row index of the first row in the CCD that is illuminated (the row closest to the readout register is row 0), for CCD positions 1, 2, 3, and 4, in case of a fast camera ([GroupID](#groupID)=Fast).  For fast cameras, only the upper half of the CCD is illuminated.  Depending on the value of the [Position](#position) parameter in the [CCD](#ccdParameters) block, the appropriate value will be selected from the list.


### TimeShift

Time shift [s] of the readout of the individual CCDs, w.r.t. to readout of CCD1.  Will only be used if [Position](#position)=1, 2, 3, or 4.