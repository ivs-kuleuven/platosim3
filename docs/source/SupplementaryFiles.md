# Supplementary Files {#SupplementaryFiles}

Depending on you @ref ConfigurationParameters "configuration", some additional files may be required:

* a file comprising a [star catalogue](#starCatalogue) of the region of the sky of interest,
* a file comprising time series for the [jitter angles (yaw, pitch, roll)](#jitterFile),
* a file comprising time series for the [drift angles (yaw, pitch, roll)](#driftFile),
* a file comprising a time series for the [focal-plane orientation](#focalPlaneOrientationFile),
* a file comprising a time series for the [focal length](#focalLengthFile),
* two files  comprising the time series for the [field distortion coefficients and their inverse](#fieldDistortionFiles),
* a file comprising the [pre-computed PSFs](#precomputedPsfFile),
* a file comprising the parameters characterising the [variation in shape of the analytic non-Gaussian PSF](#analyticNonGaussianPsfFile),
* a file comprising a time series for the [width of the analytic non-Gaussian PSF](#analyticNonGaussianPsfFile),
* a file comprising a time series for the [operating temperature of the FEE](#temperatureFeeFile),
* and a file comprising a time series for the [operating temperature of the CCD](#temperatureCcdFile).

In the following sections we describe these files in more detail.

---





<!-- ************** -->
<!-- Star Catalogue -->
<!-- ************** -->

## <a name="starCatalogue"></a>Star Catalogue

A star catalogue must be provided in a file in ASCII format. This file should contain three columns, separated by a space, holding the following information:
* right ascension of the stars [degrees]
* declination of the stars [degrees]
* stellar magnitude

The path of this file, relative to the project location, must be provided via the <code>StarCatalogFile</code> parameter in the <code>ObservingParameters</code> block in the @ref ConfigurationParameters "configuration file". 

---





<!-- *********** -->
<!-- Jitter File -->
<!-- *********** -->

## <a name="jitterFile"></a>Jitter File (Optional)

If required (<code>UseJitterFromFile</code> = <code>yes</code>), a jitter time series must be provided in a file in ASCII format. This file should contain four columns, separated by a space, holding the following information:

* time [s]
* yaw [arsec]
* pitch [arcsec]
* roll [arcsec]

The path of this file, relative to the project location, must be provided via the <code>JitterFileName</code> parameter in the <code>PlatForm</code> block in the @ref ConfigurationParameters "configuration file".
To ensure a realistic modelling of the jitter, the time step in the jitter time series must be smaller than the exposure time.

---





<!-- ********** -->
<!-- Drift file -->
<!-- ********** -->

## <a name="driftFile"></a>Drift File (Optional)

If required (<code>UseDriftFromFile</code> = <code>yes</code>), a drift time series must be provided in a file in ASCII format. This file should contain four columns, separated by a space, holding the following information:

* time [s],
* yaw [arsec],
* pitch [arcsec],
* and roll [arcsec].

The path of this file, relative to the project location, must be provided via the  <code>DriftFileName</code> parameter in the <code>Telescope</code> block in the @ref ConfigurationParameters "configuration file".

---




<!-- **************************** -->
<!-- Focal-Plane Orientation File -->
<!-- **************************** -->

## <a name="focalPlaneOrientationFile"></a>Focal-Plane Orientation File (Optional)

If required (<code>FocalPlaneOrientation: Source</code> = <code>FromFile</code>), a focal-plane orientation time series must be provided in a file in ASCII format.  This file should contain two columns, separated by a space, holding the following information:

* time [s]
* and focal-plane orientatation [degrees].

For an angle of 0°, the y-axis of the CCD (with an orientation angle of 0°) points towards the North. A positive angle corresponds to a counterclockwise rotation. Have a look at Fig. 6 @ref ConfigurationParameters "here".

The path of this file, relative to the project location, must be provided via the  <code>FocalPlaneOrientation: FromFile</code> parameter in the <code>Camera</code> in the @ref ConfigurationParameters "configuration file".

---

<!-- ***************** -->
<!-- Focal-Length File -->
<!-- ***************** -->

## <a name="focalLengthFile"></a>Focal-Length File (Optional)

If required (<code>FocalLength: Source</code> = <code>FromFile</code>), a focal-length time series must be provided in a file in ASCII format.  This file should contain two columns, separated by a space, holding the following information:

* time [s]
* and focal length [m].

The path of this file, relative to the project location, must be provided via the  <code>FocalLength: FromFile</code> parameter in the <code>Camera</code> block in the @ref ConfigurationParameters "configuration file".

---





<!-- ********************** -->
<!-- Field Distortion Files -->
<!-- ********************** -->

## <a name="fieldDistortionFiles"></a>Field Distortion Files (Optional)

If required (<code>FieldDistortion: Source</code> = <code>FromFile</code>), time series must be provided for the field distortion coefficients and their inverse in two files in ASCII format.  This file should contain columns, separated by a space, holding the following information:

* time [s]
* and the field distortion coefficients (or their inverse), all separated by spaces.

The path of these files, relative to the project location, must be provided via the <code>FieldDistortion: CoefficientsFromFile</code> and <code>FieldDistortion: InverseCoefficientsFromFile</code> parameters in the @ref ConfigurationParameters "configuration file".

---





<!-- ******** -->
<!-- PSF File -->
<!-- ******** -->

## <a name="psfFiles"></a>PSF Files (Optional)

### <a name="precomputedPsfFile"></a>Pre-Computed PSF File

If case you want to use pre-computed normalised PSFs (<code>Model</code> = <code>MappedFromFile</code>in the <code>PSF</code> block), you must provide these in the form of an HDF5-file.  The path of this file, relative to the project location, is specified via the <code>MappedFromFile: Filename</code> parameter in the <code>PSF</code> block in the configuration file.  The most recent version of the pre-computed PSFs can be downloaded from our FTP server, as described @ref ReqsRun "here". 

The simulator will automatically select the PSF for which the angular distance to the optical axis matches best for the simulated sub-field.



### <a name="analyticNonGaussianPsfFile"></a>Analytic Non-Gaussian PSF Files (Optional)

In case you want to use an analytic non-Gaussian PSF (<code>Model</code> = <code>AnalyticGaussian</code> in the <code>PSF</code> block), you must provide the parameters characterising it in a file in ASCII format.  The most recent values for these parameters can be found in the file <code>psfallv3.txt</code> file in the <code>/inputfiles</code> directory.

Additionally, if required (<code>AnalyticNonGaussian: Sigma: Source</code> = <code>FromFile</code> in the <code>PSF</code> block), a time series for the width of the analytic non-Gaussian PSF must be provided in a file in ASCII format.  This file should contain columns, separated by a space, holding the following information:

* time [s]
* and the width of the PSF [pixels].

The path of this file, relative to the project location, must be provided via the  <code>AnalyticNonGaussian: Sigma: FromFile</code> parameter in the <code>Camera</code> block in the @ref ConfigurationParameters "configuration file".

---





<!-- **************************** -->
<!-- Temperature File for the FEE -->
<!-- **************************** -->

## <a name="temperatureFeeFile"></a>Temperature File for the FEE (Optional)

If required (<code>Temperature</code> = <code>FromFile</code> in the <code>FEE</code> block), a temperature time series for the FEE must be provided in a file in ASCII format.  This file should contain two columns, separated by a space, holding the following information:

* time [s]
* and operating temperature of the FEE [K].

The path of this file, relative to the project location, must be provided via the  <code>TemperatureFileName</code> parameter in the </code>FEE</code> block in the @ref ConfigurationParameters "configuration file".

---





<!-- **************************** -->
<!-- Temperature File for the CCD -->
<!-- **************************** -->

## <a name="temperatureCcdFile"></a>Temperature File for the CCD (Optional)

If required (<code>Temperature</code> = <code>FromFile</code> in the <code>CCD</code> block), a temperature time series for the CCD must be provided in a file in ASCII format.  This file should contain two columns, separated by a space, holding the following information:

* time [s]
* and operating temperature of the CCD [K].

The path of this file, relative to the project location, must be provided via the  <code>TemperatureFileName</code> parameter in the </code>CCD</code> block in the @ref ConfigurationParameters "configuration file".