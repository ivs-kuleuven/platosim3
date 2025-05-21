* Changelog PlatoSim 3.7.0

** Improvements

*** Added more tweaks for l1 pipeline

*** Refactored aberration in Sky.cpp

*** Changed names in technical notebook

*** Cosmetic changes in Detector.cpp

*** Modified applyPhotometry() in DetectorWithAnalyticNonGaussianPSF

*** simfile.getApertureMask() now also return SPR

*** Update the plot.py method used in plot.py

*** Update sigma-clipping lower bound threshold for lightcurve.py

*** Speedup of lowess stitching algorithm of lightcurve.py

*** Upgrade picsim --vizier feature

*** Lowess-Theil-Sen rebust detrending method added to lightcurve.py

*** Move auto sigma selection for clip into platomium

*** Update F-Camera YAML parameters in simulation.py

*** Bugfix for platonium animation calling simfile.getApertureMask()

*** Implemented a passband correction in varsim.py

*** Added parameters from the input yaml file to the HDF5 file

*** Removed the use of "normal" in skyToPixelCoordinates in referenceFrames.py

*** Update of pulsation amplitude distribution in varsource.py

*** Changed how the number of traps for CTI are calculated.

*** Interpolate aberration in between timesteps.

*** Update the plato.py and patonium.py

*** Import cumtrapz or cumulative_trapezoid based on scipy version

*** Changed the comment for the gain nonlinearity in inputfile

*** Update color plotting for Aitoff galactic sky projections

*** getInfoOfStarWithID() in Sky.cpp now returns a tuple of doubles instead of tuple of ints

*** In CMakeLists.txt changed doc string from C++11 to C++14

*** Added option to disable O-C plot in detrend and clip plots

*** Added function to query a target star based on its Gaia DR3 ID in starquery.py

*** Update version of ligo.skymap to 2.0.0

*** Cosmetic changes in referenceFrames.py

*** Update starquery.py to include Gaia DR3 proper motion.




** Bug fixes

*** Bugfix in varsource.py for using wrong time column

*** Bugfix in lightcurve.py of NSR functions

*** Bugfix a conflict in pyproject.toml

*** Fixed missing return value in ClosedLoopDetectorWithAnalyticNonGaussianPSF

*** Bugfix in simulation.py when dealing with F-Camera

*** Bugfix in lightcurve.py using 'single'

*** Bugfix to varsim and varsource

*** Bugfix to platonium using verbosity level 1

*** Bugfix for fetching cadence in lightcurve.py

*** Bugfix for platonium.py and payload.py

*** Bugfix for model selection in detrending algorithm

*** Bugfix of analysis script for PLATO-CS in platonium

*** Bugfix of varsim's argparse string making it fail when consulting help

*** Bugfix for parsing argument ted_ampl in payload.py

*** Bugfix in Camera.cpp so that now star position has a value for every exposure

*** Bugfix of using seed correctly for stellar flares

*** Bugfix of varsim for mocka usage

*** Bugfix of LPV and bCep model in varsim.py

*** Bugfix in RR Lyrae and Cepheid model in varsource.py

*** Bugfix in simulation.py createStarCatalogFileFromPixelCoordinates() when dealing with quaternions.

*** Bugfix in getFlux in simfile.py, this output of this function was too low.

*** Bugfix in platonium.py setting --cadence and using sim.useNormalCamera()


** New features/functionality

*** Added straylight implementation for moon

*** Added validation test for straylight

*** Added magnitude dependent sigma-clipping value to platonium post-processing

*** Added notebook for MOCKA

*** Added validation test for SPR

*** Solar flare model implemented in varsim

*** Added script for running platonium and varsim in parallel

*** Added new analysis script to merge PLATOnium light curves

*** Added lightcurve.gaps method

*** Added new Lesia L1 script to platonium

*** Added script to generate a time series from pulsation modes of varsim.py
