# Prerequisites for Installing and Updating PlatoSim3 {#dev-prerequisites}

## Git

The PlatoSim code is under version control in GitHub. To be able to get the latest version of the code on your local machine and to share possible contributions with other people in the project, you need to install Git on your computer. Installation instructions can be found in the <a href="https://git-scm.com/book/en/v2/Getting-Started-Installing-Git">Git reference documentation</a>.

![](https://git-scm.com/images/logo@2x.png)

---

## GitHub

The PlatoSim code is located in a [private repository in GitHub](https://github.com/IvS-KULeuven/PlatoSim3).  To be able to access it, we have to grant you access explicitly. Please, send your GitHub username to the development team and you will be granted read access to the repository.  You will then get an invitation by email to join.

If you have not done so already, you can make an [account on GitHub](https://github.com/join) for free.

![](images/octocat.png)

---

## Operating System

PlatoSim3 has been tested on recent Mac OS X and recent Linux (Fedora, Ubuntu, Debian) systems. We cannot support Windows systems at the moment, so we advise Windows users to install a Virtual Machine (VM), and run the PLATO Simulator inside the VM.  Alternatively, if you're using Windows 10, you can use [Windows Subsystems for Linux (WSL)](https://docs.microsoft.com/en-us/windows/wsl/install-win10) instead.

---

## Installed Software Packages

To be able to install the dependencies and build the code, the following software must be installed on your computer:

* Python: for the installation of the dependencies
* <a href="https://gcc.gnu.org/">gcc</a> v5.1 or more recent, or <a href="http://clang.llvm.org/">clang</a> v3.3 or more recent
*  <a href="https://cmake.org/">CMake</a>: cross-platform open-source build system to control the software compilation process (using simple platform and compiler independent configuration files)
* <a href="http://www.openblas.net">BLAS</a> and <a href = "http://www.netlib.org/lapack/">LAPACK</a>. Without these, the simulator will likely be slower. These libraries come pre-installed on Mac OS X (so Mac users do not have to do anything). Many Linux distributions also standardly have these libraries installed or offer a package manager to easily install them.
        
It is also possible to run the PLATO Simulator from within Python (please, use Python 3).  We recommend downloading a Python distribution such as <a href="https://docs.continuum.io/anaconda/install">Anaconda</a>.  More information on setting up your Python distribution and scripting PlatoSim3 using Python can be found @ref ScriptingInPython "here".  We have also provided tutorials in the form of Python notebooks.  More information about those can be found @ref Tutorials "here".


