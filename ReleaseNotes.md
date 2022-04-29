* Release Notes PlatoSim 3.5.3



** Improvements

*** Improved stability of validations tests for drift/jitter from file

*** Renamed previous starcatalog

*** F-Camera can now be simulated with custom CCD





** Bug fixes

*** When cosmics are added for F-Cams, we make sure they can not fall into covered part

*** Corrected bug 'calculateSubfieldAroundCoordinates' in 'referenceFrames.py' where middel pixel of subfield was rounded up/down by 1.

*** orbit.txt file is now only read if we include aberration in the simulation




    
** New features/functionality

*** Added custom inputfile for F-Camera

*** Added new starcatalog

*** Added CTI in SmearingMaps for "Short2013" CTI model


