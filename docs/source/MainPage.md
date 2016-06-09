@mainpage Documentation for the Plato Simulator

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

PlatoSim has been used on recent OS X and recent Linux (Fedora, Ubuntu, Debian) systems. We cannot support  Windows systems at the moment, so we advise Windows users to install a Virtual Machine, and run PlatoSim inside the VM.

The PlatoSim3 software is being distributed via <a href="https://github.com/IvS-KULeuven/PlatoSim3">GitHub</a> (see screenshot below), so you'll be needing a (free) <a href="https://github.com/join">GitHub account</a>.  As the PlatoSim3 repository is kept private for now, access must be granted by the PlatoSim3 team explicitly.

If you are interested in contributing to the software, you must <code>fork</code> PlatoSim3, using the GitHub web interface.  Just press the <code>"Fork"</code> button on the GitHub web interface (see Fig. 1).

If you are only interested in using PlatoSim3, it suffices to download the ZIP file, by pressing the <code>"Download ZIP"</code> button in the GitHub web interface (see Fig. 1). The disadvantage is that - every time you want to update to a more recent version of the software - you must download a new ZIP file and re-install the dependencies (see below).

If you want to be able to update the software (without having to re-install the dependencies each time), it is better to <code>clone</code> PlatoSim3 by executing the following command in a designated directory (you have to do this only once):

\code git clone https://github.com/IvS-KULeuven/PlatoSim3.git .\endcode

Mind the dot at the end of the command!

You can then update the software by executing the command

\code git pull origin master \endcode 

in the directory in which you installed PlatoSim3.

However, this will only work smoothly if you did not change any of the PlatoSim3 files or added files to the PlatoSim3 folders. The only exception is the <code>/inputfiles</code> folder, where you can add files but should not change the original files. We recommend that you copy the <code>inputfile.yaml</code> file to your own version.

At a later stage, releases will be distributed.  We will then extend the documentation accordingly.

@image html /images/gitHub.png "Figure 1: Screenshot of the GitHub web interface."

After you have downloaded the PlatoSim3 code, you must first take care of the dependencies before you can actually build and run the Plato Simulator. The next sub-sections describe the requirements and the procedures to install the dependencies and to build the PlatoSim3 code.



<!-- Requirements -->
<!-- ************ -->

@subsection requirements Requirements

To be able to install the dependencies and build the code, the following software must be installed on your computer:

* Python: for the installation of the dependencies
* gcc v5.1 or more recent, or clang v3.3 or more recent
*  <a href="https://cmake.org/">CMake</a>: cross-platform open-source build system to control the software compilation process (using simple platform and compiler independent configuration files)
* <a href="http://www.openblas.net">BLAS</a> and <a href = "http://www.netlib.org/lapack/">LAPACK</a>. Without these, the simulator will likely be slower. These libraries    come pre-installed on Mac OS X (so Mac users do not have to do anything). Many Linux distributions also standardly have these libraries installed, or offer a package manager to easily install them.
        
It is also possible to run simulation from within Python.  We recommend downloading a Python distribution such as <a href="https://docs.continuum.io/anaconda/install">Anaconda</a>.

<!-- Dependencies -->
<!-- ************ -->

@subsection dependencies Dependencies

PlatoSim3 relies on a number of other dependencies, which are all included in the packages for your convenience.  Everything concerning the dependencies can be found in the <code>/dependencies</code> directory.  The <code>/dependencyDownloads</code> sub-directory contains the tarball or zipball file of all required packages.  In the <code>/installscripts</code> sub-directory you can find Python scripts that help you to unzip or untar these files into the <code>/dependencyInstalls</code> directory.  You can do this for one dependency at a time, like this:

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

@section run Running PlatoSim3

After the installation of the software, the Plato Simulator can be run from the <code>/build</code> directory you have created. For the simulation itself, only one input file with configuration parameters (e.g. <code>/inputfiles/inputfile.xml</code>) is required as input. 

To initiate a simulation, <code>cd</code> to the <code>/build</code> directory and type:

\code
./platosim <input file> <non-existing output file> [<log file>]
\endcode

The structure of the input file and the meaning of the individual configuration parameters are described in a separate document: @ref InputFileDescription

Note that - before running PlatoSim3 - you must have an environment variable <code>PLATO_PROJECT_HOME</code>, set to the base folder of PlatoSim3.  For example, you can put in your <code>.bash_profile</code> the following:

\code
PLATO_PROJECT_HOME="<path to PlatoSim3>"
export PLATO_PROJECT_HOME
\endcode

Note that you have to use an absolute path!

If you want to use realistic PSF models instead of a Gaussian, you can download these from <a href="ftp://ftp.ster.kuleuven.be/dist/Plato/PlatoSim3/psf.hdf5">our FTP server</a>.  A convenient place to store this file is in the <code>/inputfiles</code> directory.

<!-- Accessing the Output -->
<!-- ********************* -->

@subsection output Accessing the Output

PlatoSim3 writes its output to an HDF5 file. HDF stands for Hierarchical Data Format, and is a next generation file format that was specifically designed to store and organise large amounts of data.

HDF5 behaves much likes a Unix-like folder structure, but where folders are called groups.  Each group can contain other groups, array datasets, and scalar attributes. For example, the first subfield image is located in <code>/Images/image000000</code> in the HDF5 file.

To quickly list the contents of the group structure of an HDF5 file on the command line, make sure that your PATH environment variable includes <code>dependencies/Installs/hdf5-1.8.16/bin/</code>, so that you can execute

\code
$ h5ls myOutputfile.hdf5
\endcode

or e.g.

\code
$ h5ls myOutputfile.hdf5/StarCatalog
\endcode

@subsubsection python Python

For Python users, we provided a <code>simfile.py</code> module in the <code>/python</code> folder, with convenient tools to extract and plot the Simulator output. For example, one can plot a subfield image using

\code{.py}
from simfile import *
myFile = SimFile("myOutputfile.hdf5")
myFile.showImage(0)
\endcode

The top of <code>simfile.py</code> contains documentation with several examples.


You can also access the HDF5 file using the <code>pytables</code> module. For example (using the latest version of <code>PyTables</code>):

\code{.py}
import tables as tbl
myFile = tbl.open_file("myOutputfile.hdf5", "r")
image = myFile.root.Images.image000000
imshow(image, interpolation="nearest", origin="lower")
myFile.root.InputParameters.CCD._v_attrs
...
print(myFile.root.InputParameters.CCD._v_attrs.Gain)
\endcode


@subsubsection idl IDL
IDL user can access the HDF5 file using, for example, 

\code{.idl}
path = FILEPATH(“Simul01.hdf5")
file = H5F_OPEN(path)
contents = H5_PARSE(path)
help, contents, /STRUCTURE
...
help, contents.Images, /STRUCTURE
...
dataset = H5D_OPEN(file,'/Images/image000000') 
image = H5D_READ(dataset)
print, size(image)
\endcode




<!-- ********* -->
<!-- Reference -->
<!-- ********* -->

@section reference Reference

We kindly ask you to refer to <a href="http://arxiv.org/abs/1404.1886">this work</a> in any publication the Plato Simulator software contributes to.


