* Release Notes PlatoSim 3.5.2



** Improvements

*** Changed field distortion for mapped PSF to deal with more general mapped psf files.

*** Changed field distoriton for analytic PSF from Radial model to Wang model. (GitHub #652)

*** Changed 'distortedToUndistortedFocalPlaneCoordinates' and 'undistortedToDistortedFocalPlaneCoordinates' in 
python/platosim/referenceFrames.py

*** Orientation angle received via network is propagated correctly to the detector (Github #660)

*** Changed deprecated 'append' method for pandas dataframe in MappedGaussianPSF validation test into 'concat' method

*** Changed pixelToSkyCoordinates function in referenceFrames.py to work better with F-cameras





** Bug fixes

*** Corrected bug in validationtest for Cosmics.

*** Corrected bug in python mapped distortion functions in python/platosim/referenceFrames.py. (GitHub #659)

*** Fixed bug in Camera.cpp. Previously distortion for mapped PSF would only be taken into account when 
includeFieldDistoritions was set to True. Now, mapped distoriton always happens independent of that value. 






    
** New features/functionality

*** Added inhomogenous trap density (GitHub #639)

*** Added validationtest for Short2013 CTI

*** Added validationtest for Short2013fromfile CTI

*** Added 'distortioncoefficients.txt' and 'distortioninversecoefficients.txt'

*** Added a first version of the PlatoSim license

*** Added metallic shield around CCD for F-Cameras

*** Added validation test for metallic shield around CCDs for F-Cameras


