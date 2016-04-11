@mainpage Documentation for the PLATO Simulator

<!-- ************************************ -->
<!-- Welcome to the Plato Simulator Pages -->
<!-- ************************************ -->

@section intro Welcome to the Plato Simulator Pages

The Plato Simulator is an end-to-end software tool, designed to perform realistic simulations of the expected observations of the Plato mission. It can, however, easily be adapted to similar types of missions.

Our simulator models and simulates time series of CCD images by including models of the CCD and its electronics, the telescope optics, the stellar field, the jitter movements of the spacecraft, and all important natural sources of noise.

Many aspects concerning the design trade-off of a space-based instrument and its performance can best be tackled through realistic simulations of the expected observations. The complex interplay of various noise sources in the course of the observations make such simulations an indispensable part of the assessment and design study of any space-based mission.





<!-- ************ -->
<!-- Installation -->
<!-- ************ -->

@section installation Installation

The PlatoSim3 software is being distributed via <a href="https://github.com/IvS-KULeuven/PlatoSim3">GitHub</a> (see screenshot below). 

If you are interested in contributing to the software, you must <code>fork</code> PlatoSim3.  If you are only interested in using it, it suffices to download the ZIP-file.  At a later stage, releases will be distributed.

@image html /images/gitHub.png ""

After you have downloaded the Plato Simulator code, you must first take care of the dependencies before you can actually build and run the Plato Simulator. The next sub-sections describe the requirements and the procedures to install the dependencies and to build the Plato Simulator code.



<!-- Requirements -->
<!-- ************ -->

@subsection requirements Requirements

To be able to install the dependencies and build the code, the following software must be installed on your computer:

* Python: for the installation of the dependencies
*  <a href="https://cmake.org/">CMake</a>: cross-platform open-source build system to control the software compilation process (using simple platform and compiler independent configuration files)


<!-- Dependencies -->
<!-- ************ -->

@subsection dependencies Dependencies

Everything concerning the dependencies for PlatoSim3 can be found in the <code>/dependencies</code> directory.  The <code>/dependencyDownloads</code> contains the tarball or zipball file of all required packages.  In the <code>/installscripts</code> sub-directory you can find Python scripts that help you to unzip or untar these files into the <code>/dependencyInstalls</code> directory.  You can do this for one dependency at a time, like this:

 \code{.py}
python ./installscripts/install_hdf5.py
python ./dependencies/installscripts/install_googletest.py
python ./dependencies/installscripts/install_yaml.py
python ./dependencies/installscripts/install_armadillo.py
python ./dependencies/installscripts/install_fftw.py
\endcode

Alternatively, you can also install the required dependencies and build the code in one go by typing: 

\code{.py}
./install.sh
\endcode



<!-- Build & Install -->
<!-- *************** -->

@subsection build Build & Install

The first time you want to run the Plato Simulator or each time you changed something in the code, the software must be built again. Just <code>cd</code> to the <code>/build</code> directory (the first time you will have to create this directory yourself) and type the following commands:

\code
cmake..
(make clean)
make -j 4
\endcode

This will create two executables :
* <code>platosim</code> to run simulations (see below)
* and <code>testplatosim</code> to run the test harnesses (without arguments).



<!-- Running PlatoSim3 -->
<!-- ***************** -->

@subsection run Running PlatoSim3

After the installation of the software, the Plato Simulator can be run from the <code>/build</code> directory you have created. For the simulation itself, only one input file with configuration parameters (e.g. <code>/inputfiles/inputfile.xml</code>) is required as input. 

To initiate a simulation, <code>cd</code> to the <code>/build</code> directory and type:

\code
./platosim < input file > <output file> [<log file>]
\endcode

The structure of the input file and the meaning of the individual configuration parameters are described @ref docs/InputFile.md "here".





<!-- ********* -->
<!-- Reference -->
<!-- ********* -->

@section reference Reference

We kindly ask you to refer to <a href="http://arxiv.org/abs/1404.1886">this work</a> in any publication the Plato Simulator software contributes to.