@mainpage <!-- Documentation for the PLATO Simulator -->





<!-- ************************************ -->
<!-- Welcome to the PLATO Simulator Pages -->
<!-- ************************************ -->

## <a name="welcome"></a>Welcome to the PLATO Simulator Pages

The PLATO Simulator is an end-to-end software tool, designed to perform realistic simulations of the expected observations of the PLATO mission. It can, however, easily be adapted to similar types of missions.

Our simulator models and simulates time series of CCD images by including models of the CCD and its electronics, the telescope optics, the stellar field, the jitter movements of the spacecraft, and all important natural sources of noise.

Many aspects concerning the design trade-off of a space-based instrument and its performance can best be tackled through realistic simulations of the expected observations. The complex interplay of various noise sources in the course of the observations make such simulations an indispensable part of the assessment and design study of any space-based mission.





<!-- ************ -->
<!-- Installation -->
<!-- ************ -->

## <a name="installation"></a>Installation

The PlatoSim3 software is being distributed via <a href="https://github.com/">GitHub</a>, so you'll be needing a (free) <a href="https://github.com/join">GitHub account</a>.  As the 
<a href="https://github.com/IvS-KULeuven/PlatoSim3">PlatoSim3 repository</a> is kept private for now, access must be granted by the PlatoSim3 team explicitly.

PlatoSim3 has been tested on recent Mac OS X and recent Linux (Fedora, Ubuntu, Debian) systems. We cannot support Windows systems at the moment, so we advise Windows users to install a Virtual Machine (VM), and run the PLATO Simulator inside the VM.

### <a name="retrieve&update"></a>Retrieving & Updating the Software

To download the software, go to the <a href="https://github.com/IvS-KULeuven/PlatoSim3">PlatoSim3 repository on GitHub</a>. The GitHub web interface to download the PlatoSim3 software is shown in Fig. 1.

If you are interested in contributing to the software, you must <code>fork</code> PlatoSim3, using the GitHub web interface.  Just press the <code>"Fork"</code> button and follow the instructions.

Although GitHub allows downloading the code as a ZIP file (by pressing the <code>"Download ZIP"</code> button in the GitHub web interface), we strongly discourage this, as this makes the process of updating to a more recent version of the software more complex and tedious.

If you want to be able to update the software (without having to re-install the [dependencies](#dependencies) each time), it is better to <code>clone</code> PlatoSim3 by executing the following command in a designated directory (you have to do this only once!):

\code git clone https://github.com/IvS-KULeuven/PlatoSim3.git .\endcode

Mind the dot at the end of the command!

You can then update the software by executing the command

\code git pull origin master \endcode 

in the directory in which you installed PlatoSim3.

However, this will only work smoothly if you did not change any of the PlatoSim3 files or added files to the PlatoSim3 folders. The only exceptions are the <code>/inputfiles</code> and the <code>/build</code> folder, where you can add files.  Please, do not modify the original files in the <code>/inputfiles</code> folder, as this might cause problems when updating the software.  We recommend that you copy the <code>inputfile.yaml</code> file and modify the copy rather than the original file.

At a later stage, releases will be distributed.  We will then extend the documentation accordingly.

Please note that you have to re-build the code each time you fetch software changes. How to do this is explained [here](#build&install).

@image html /images/gitHub.png "Figure 1: Screenshot of the GitHub web interface."

After you have downloaded the PlatoSim3 code, you first have to install a few packages (so-called [dependencies](#dependencies)) before you can actually build and run the PLATO Simulator.  The next sub-sections describe the [requirements](#requirements), and the procedures to install the [dependencies](#dependencies) and to [build](#build&install) the PlatoSim3 code.



<!-- Requirements -->
<!-- ************ -->

### <a name="requirements"></a>Requirements

To be able to install the dependencies and build the code, the following software must be installed on your computer:

* Python: for the installation of the dependencies
* <a href="https://gcc.gnu.org/">gcc</a> v5.1 or more recent, or <a href="http://clang.llvm.org/">clang</a> v3.3 or more recent
*  <a href="https://cmake.org/">CMake</a>: cross-platform open-source build system to control the software compilation process (using simple platform and compiler independent configuration files)
* <a href="http://www.openblas.net">BLAS</a> and <a href = "http://www.netlib.org/lapack/">LAPACK</a>. Without these, the simulator will likely be slower. These libraries come pre-installed on Mac OS X (so Mac users do not have to do anything). Many Linux distributions also standardly have these libraries installed or offer a package manager to easily install them.
        
It is also possible to run the PLATO Simulator from within Python.  We recommend downloading a Python distribution such as <a href="https://docs.continuum.io/anaconda/install">Anaconda</a>.  More information on setting up your Python distribution and scripting PlatoSim3 using Python can be found @ref ScriptingInPython "here".

<!-- Dependencies -->
<!-- ************ -->

### <a name="dependencies"></a>Dependencies

PlatoSim3 relies on a number of other packages, which are all included in the PlatoSim3 distribution for your convenience.  Everything concerning the dependencies can be found in the <code>/dependencies</code> directory.  The <code>/dependencyDownloads</code> sub-directory contains the tarball or zipball file of all required packages.  In the <code>/installscripts</code> sub-directory you can find Python scripts that help you to unzip or untar these files into the <code>/dependencyInstalls</code> directory.  You can do this for one dependency at a time, like this (in the directory in which PlatoSim3 was downloaded):

 \code{.py}
python ./dependencies/installscripts/install_hdf5.py
python ./dependencies/installscripts/install_googletest.py
python ./dependencies/installscripts/install_yaml.py
python ./dependencies/installscripts/install_armadillo.py
python ./dependencies/installscripts/install_fftw.py
\endcode

Alternatively, you can also install the required dependencies and build the code in one go by typing (also in the directory in which PlatoSim3 was downloaded): 

\code{.py}
./install.sh
\endcode

If problems would arise when executing this command, it may be useful to tried to install the dependencies one-by-one to pinpoint the problem.


<!-- Build & Install -->
<!-- *************** -->

### <a name="build&install"></a>Build & Install

The first time you want to run the PLATO Simulator or each time you changed something in the code, the software must be built again. Just <code>cd</code> to the <code>/build</code> directory (the first time you will have to create this directory yourself) and type the following commands (after emptying the <code>/build</code> directory):

\code
cmake ..
(make clean)
make -j 4
\endcode

This will create two executables :
* <code>platosim</code> to run simulations (see below)
* and <code>testplatosim</code> to run the test harnesses (without arguments).

<!-- Before Running PlatoSim3 -->
<!-- ************************ -->

## <a name="beforeRunningPlatoSim3">Before Running PlatoSim3

After the installation of the software, the PLATO Simulator can be run from the <code>/build</code> directory you have created, but first you have to set a couple of environment variables


### <a name="environmentVariables">Setting up your Environment Variables

To avoid having to hardcode any path in configuration files, tutorials, etc., you must export three environment variables:

- <code>PLATO_PROJECT_HOME</code>: to refer to the directory in which PlatoSim3 was installed,
- <code>PLATO_WORKDIR</code>: to refer to the directory you want to write output or in which you want to store your own configuration files (preferably not within <code>/PlatoSim3</code>),
- and <code>PYTHONPATH</code>: to refer to the directory in which our Python scripts can be found.

This can be done as follows:

\code
PLATO_PROJECT_HOME=<full path to /PlatoSim3>
export PLATO_PROJECT_HOME

PLATO_WORKDIR=<full path to a preferred working directory>
export PLATO_WORKDIR

PYTHONPATH=$PYTHONPATH:<full path to /PlatoSim3/python>
export PYTHONPATH
\endcode

If you want, you can copy this code to make your own little script to set up your environment (e.g. <code>setPlatoEnvironment</code>), or add it to your <code>.bash_profile</code>.

To check the content of these variables (to see whether they are set to the proper location), type

\code
echo $PLATO_PROJECT_HOME

echo $PLATO_WORKDIR

echo $PYTHONPATH
\endcode

### <a name="yourOwnFiles">Where to Store your own Files?

To avoid problems when updating the PlatoSim3 software, it is best to store your own input and output files in a designated working directory, (preferably) outside the installation directory of PlatoSim3.  You can (but should not) add your input files to the <code>/inputfiles</code> directory, but under no circumstances change the original files in that directory!

The path to the designated working directory must be exported as the <code>PLATO_WORKDIR</code> environment variable, as described above.

<!-- Running PlatoSim3 -->
<!-- ***************** -->

## <a name="runningPlatoSim3"></a>Run PlatoSim3

 For the simulation itself, only one input file with configuration parameters (e.g. <code>/inputfiles/inputfile.yaml</code>) is required as input. 

To initiate a simulation, <code>cd</code> to the <code>/build</code> directory and type:

\code
./platosim <input file> <non-existing output file> [<log file>]
\endcode

The structure of the input file and the meaning of the individual configuration parameters are described @ref InputFileDescription "here".

Note that - before running PlatoSim3 - you must have an environment variable <code>PLATO_PROJECT_HOME</code>, set to the base folder of PlatoSim3.

If you want to use realistic PSF models instead of a Gaussian, you can download these from <a href="ftp://plato:miSotalP@ftp.ster.kuleuven.be/psf.hdf5
">our FTP server</a>.  A convenient place to store this file is in the <code>/inputfiles</code> directory.

<!-- Accessing the Output -->
<!-- ********************* -->

### <a name="output"></a>Output

The output of the simulations with PlatoSim3 is stored as an HDF5 file.  The structure of such files and how to access (and visualise) the content is described @ref OutputFileDescription "here".

Note that the x-axis is defined along the serial readout register and corresponds to the columns of the detector (and therefore of the pixel and sub-pixel map).  The y-axis corresponds to the rows of the detector (and therefore of the pixel and sub-pixel map).  In Python, for example, the <code>imshow()</code> method of <code>matplotlib</code> transposes the image as it tries to mimic Matlab.

The <code>Armadillo</code> arrays (that are used internally to store the maps) are column-major rather than row-major.





<!-- ******************* -->
<!-- In Case of Problems -->
<!-- ******************* -->

## <a name=inCaseOfProblems>In Case of Problems

In case you would come across problems you cannot solve yourself, please, let us know!  We would like you to use the issue tracking on GitHub rather than sending an email to the developers, as this helps to better keep track of the issues and their status.  How raise issues (and which information you must provide us with) is desribed @ref IssueTracking "here".

<!-- ********* -->
<!-- Reference -->
<!-- ********* -->

## <a name="reference"></a>Reference

We kindly ask you to refer to <a href="http://arxiv.org/abs/1404.1886">this work</a> in any publication the PLATO Simulator software contributes to.


