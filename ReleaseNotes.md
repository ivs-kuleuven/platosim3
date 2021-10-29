* Release Notes PlatoSim 3.5.1



** Improvements

*** The dependencies python install files now check that the `Installs` directory exists and creates this directory if it doesn't. 

*** Made the log files for `Camera::makeStarCatalogSelection` clearer. 

*** Changed the python functions in `referenceFrames.py`, `plot.py` and `simulation.py` to deal with mapped distortion.

*** Mapped distortion now uses a continuous approximation insead of the previous (crude) method of one-to-one fitting of closed point.

*** Changed cosmics intensity from uniform to skew-normal. (GitHub #638)

*** Renamed `getCosmicsCoordinates()` to `getCosmicsAffectedPixels()` in `simfile.py`

*** Updated website 




** Bug fixes

*** The diffused PSF that was saved to the output HDF5 is now rotated with respect to the CCD it falls on. (GitHub #627)

*** Fixed bug where the star coordinates where written to the output HDF5 file without taking field distortion into account. (GitHub #631)



    
** New features/functionality

*** Added an option to individually switch on/off extended or pointlike ghosts. 

*** Added `getCosmicsInfo()` method to extract the entry position, entry angle and the trail length of all cosmics


