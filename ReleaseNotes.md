* Release Notes PlatoSim 3.5.4



** Improvements

*** Changed the test for the gain, to correctly test that the left/right CCD have different gains.

*** Changed the start time for the orbit.txt file.

*** changed brighterFatterEffect.py validation test to be able to deal with the different gain values (on different CCDs).





** Bug fixes

*** Fixed bug where incorrect gain was applied to the right CCD.

*** Fixed a bug when applying open shutter smearing when we include relative transmissivity.

*** Fixed issue with CTI (Short2013 model) where we did not use the correct dwell time.





** New features/functionality

*** Added a sensible estimate of the number of occupied traps when applying "Short2013" CTI model.

*** Added comparison between new and old Cosmic ray model to the validation tests.

*** New plot module to plot a star in the CCD focal plane.

*** Added the platonium package into the PlatoSim python files.


