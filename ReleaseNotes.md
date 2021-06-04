* Release Notes PlatoSim 3.4.1



** Improvements


*** Added an option to the method `getYawPitchRoll` to also obtain the time (GitHub issue #508)

*** Apply the BFE after full-well saturation (GitHub issue #584)
The current implementation of the BFE fails for very large pixel values, which results in negative pixel values in
the pixelmap. By first applying full-well saturation, these large values will not occur.

*** The conda install should now work for python versions 3.6, 3.7, 3.8 and 3.9.


** Bug fixes

*** Corrected the implementation of the jitter/drift from red noise
The previous implementation of the jitter and drift was not consistent between different CCDs. It is now implemented
to be consistent between different CCDs if the jitter/drift seed is the same.

*** Fixed conda build of PlatoSim in Jenkins


    
** New features/functionality

*** Validation test for Jitter on different CCDs

*** Add method `getYawPitchRollFromDrift` in simfile.py
Added a method to extract the Yaw, Pitch, Roll and time from the output file.

*** Option to add the diffused PSF in the output file. (GitHub Issue #564)

*** Add option `setSubfieldAroundPixelRows` (GitHub Issue #587)

*** Added the option to (not) include in the output file:
    - High resolution PSF (if PSF is not Analytic Gaussian)
    - Star Catalog
    - Platform Yaw, Pitch, Roll
    - Transmission Efficiency

