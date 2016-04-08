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
	ProjectLocation:	ENV['PLATO_PROJECT_HOME']
\endcode


@subsubsection projectLocation ProjectLocation

<i>Allowed values:</i> name of an existing directory on disk or environment variable in the format <dfn>ENV['PLATO_PROJECT_HOME']</dfn>

Full path of the directory in which you have checked out the PlatoSim3 project, or an environment variable, e.g. PLATO_PROJECT_HOME, containing the full path to that directory.  In the latter case, you must make sure you have exported this variable before initiating a simulation:

\code{.unparsed}
 export PLATO_PROJECT_HOME = <full path to the PlatoSim3 directory>
\endcode

@subsection observingParameters Observing Parameters

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
 
The YAML input file consists of the following blocks:

* General

   Specifies general parameters that are used throughout the simulator code. The most important parameter is the *ProjectLocation* which defines the path where input- and output for the simulator run will be stored. The ProjectLocation also serves as the prefix to be used for relative paths used for other input paramaters.

* ObservingParameters
* Platform
* Telescope
* Camera
* CCD
* SubField
* RandomSeeds

