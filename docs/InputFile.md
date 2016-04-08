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

<i>Allowed values:</i> \f$ \in [0, 360]\f$

Right ascension of the pointing, expressed in degrees.

@subsubsection  decPointing DecPointing

<i>Allowed values:</i> \f$ \in [-90, 90]\f$

Declination of the pointing, expressed in degrees.

@subsubsection fluxm0 Fluxm0

<i>Allowed values:</i> > 0

Flux of a star of zero magnitude (\f$ m_{\lambda} = 0 \f$), expressed in \f$ photons \cdot s^{-1} \cdot cm^{-2}\f$ in the passband of the magnitudes listed in the star catalogue.

For an exposure of \f$t_{exp}\f$ seconds, the measured flux \f$F_{phot}\f$ of a star is computed from its catalogue magnitude \f$m_{\lambda}\f$, the effective light-collecting area \f$A\f$ (in \f$cm^2\f$) of the telescope, the transmission efficiency \f$T_{\lambda}\f$ of the optical system, the quantum efficiency \f$Q\f$ of the detector, and the flux per second \f$F_0\f$ of a star with zero magnitude (\f$m_{\lambda} = 0\f$) from the equation

\f[
F_{phot} = t_{exp} \cdot F_0 \cdot T_{\lambda} \cdot Q \cdot A \cdot 10^{-0.4} m_{\lambda}
\f]

where the \f$\lambda\f$ subscript refers to the wavelength range in which the simulation is performed.

@subsubsection skyBackground SkyBackground

<i>Allowed values:</i> ≥ 0

The sky background (zodiacal + galactic), expressed in \f$photons \cdot s^{-1} \cdot pixel^{-1}\f$.

@subsubsection  starCatalogFile StarCatalogFile

Full path to the @ref starCatalogFile "star catalogue file".

@subsection platformParameters Platform Parameters

@subsection telescopeParameters Telescope Parameters

@subsection cameraParameters Camera Parameters

@subsection psfParameters PSF Parameters

@subsection ccdParameters CCD Parameters

@subsection subFieldParameters Sub-Field Parameters

@subsection seedParameters Seed Parameters

@section starCatalogue Star Catalogue

@section psfFile PSF File (Optional)

@section jitterFile Jitter File (Optional)


