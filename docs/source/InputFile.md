# Description of the Input Files {#InputFileDescription}

## Configuration Parameters

To configure the Plato Simulator, a large set of input parameters is required.  We distinguish between the following types of configuration parameters:

* general parameters,
* observing parameters,
* sky parameters,
* platform parameters,
* telescope parameters,
* camera parameters,
* PSF parameters,
* FEE parameters,
* CCD parameters,
* sub-field parameters,
* seed parameters,
* content control parameters,
* and pre-defined settings for the camera groups and CCD positions (specific for PLATO).

These are described in more detail @ref ConfigurationParameters "here".

---





## Supplementary Files

Additionally, depending on the configuration, some additional files may be required:

* a file comprising a star catalogue of the region of the sky of interest,
* a file comprising time series for the jitter angles (yaw, pitch, roll),
* a file comprising time series for the drift angles (yaw, pitch, roll),
* a file comprising a time series for the focal-plane orientation,
* a file comprising a time series for the focal length,
* two files  comprising the time series for the field distortion coefficients and their inverse,
* a file comprising the [pre-computed PSFs](#precomputedPsfFile),
* a file comprising the parameters characterising the [variation in shape of the analytic non-Gaussian PSF](#analyticNonGaussianPsfFile),
* a file comprising a time series for the [width of the analytic non-Gaussian PSF](#analyticNonGaussianPsfFile),
* a file comprising a time series for the operating temperature of the FEE,
* and a file comprising a time series for the operating temperature of the CCD.

These are described @ref SupplementaryFiles "here".




## Migration of Input Files

When new versions of PlatoSim are released, there might be some changes in the configuration or input files. New parameters might have been added or parameters or groups might have been restructured. In the PlatoSim distribution we provide a helper tool (``migtool.py``) to migrate your YAML inputfiles to the new format. This tool will compare your (old) YAML inputfile with the default YAML inputfile that is on your local system, i.e. in the folder ``$PLATO_PROJECT_HOME/inputfiles``. The command can be used as follows (asuming your are located in the PlatoSim project home folder):

```
usage: python python/migtool.py [-hv] [-o outputFilename] inputFilename
```

The ``-h`` option prints an instructive help message. The ``inputFilename`` is your (old) YAML inputfile. The ``-o`` option let's you specify the name of the output file in which the migrated configuration will be saved. When no output file is given, the result will be printed on the screen (``stdout``).

The ``-v`` option prints the changes that will be applied on the screen. That might be useful, because it will signal changes in parameter values. An example is shown below. ``CHECK`` means you will have to check manually. The value is your input file was ``0.016`` and the value in the new input file is ``0.01``. Your value will be retained after the migration. The reason is that you probably have good reasons to have this value different from the default value, and you don't want to loose that during the migration.

```
CHECK - Value changed for CCD.FlatfieldPtPNoise from 0.016 to 0.01
```

Other possible output when using the ``-v`` option is signaled with ``DONE``, i.e. when a new key or subkey has been added or when an obsolete key has been removed.

Unfortunately, the python module that we use to parse the YAML files does not retain the comments that are in the files. Those comments will be lost after migration.

