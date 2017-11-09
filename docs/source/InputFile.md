# Description of the Input Files {#inputFileDescription}

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
* and pre-defined settings for the camera groups and CCD positions (specific for PLATO).

These are described in more detail @ref configurationParameters "here".

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

These are described @ref supplementaryFiles "here".
 