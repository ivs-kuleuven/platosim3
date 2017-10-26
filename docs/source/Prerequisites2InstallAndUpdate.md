# Prerequisites {#reqsInstallUpdate}

## GitHub

The PlatoSim3 software is being distributed via <a href="https://github.com/">GitHub</a>, so you'll be needing a (free) <a href="https://github.com/join">GitHub account</a>.  As the 
<a href="https://github.com/IvS-KULeuven/PlatoSim3">PlatoSim3 repository</a> is kept private for now, access must be granted by the PlatoSim3 team explicitly.

## Operating System

PlatoSim3 has been tested on recent Mac OS X and recent Linux (Fedora, Ubuntu, Debian) systems. We cannot support Windows systems at the moment, so we advise Windows users to install a Virtual Machine (VM), and run the PLATO Simulator inside the VM.

## Installed Software Packages

To be able to install the dependencies and build the code, the following software must be installed on your computer:

* Python: for the installation of the dependencies
* <a href="https://gcc.gnu.org/">gcc</a> v5.1 or more recent, or <a href="http://clang.llvm.org/">clang</a> v3.3 or more recent
*  <a href="https://cmake.org/">CMake</a>: cross-platform open-source build system to control the software compilation process (using simple platform and compiler independent configuration files)
* <a href="http://www.openblas.net">BLAS</a> and <a href = "http://www.netlib.org/lapack/">LAPACK</a>. Without these, the simulator will likely be slower. These libraries come pre-installed on Mac OS X (so Mac users do not have to do anything). Many Linux distributions also standardly have these libraries installed or offer a package manager to easily install them.
        
It is also possible to run the PLATO Simulator from within Python (please, use Python 3).  We recommend downloading a Python distribution such as <a href="https://docs.continuum.io/anaconda/install">Anaconda</a>.  More information on setting up your Python distribution and scripting PlatoSim3 using Python can be found @ref ScriptingInPython "here".  We have also provided tutorials in the form of Python notebooks.  More information about those can be found @ref Tutorials "here".