# Description of the input file

The input file format use for PlatoSim3 is YAML (see https://learnxinyminutes.com/docs/yaml/). We use just a very limited set of the YAML functionality, enough to allow us to provide input files for different parts of the simulator. 

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

