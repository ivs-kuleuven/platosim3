# Building PlatoSim3 {#Building}

## First-Time Installation

PlatoSim3 relies on a number of other packages, which are all included in the PlatoSim3 distribution for your convenience.  Everything concerning the dependencies can be found in the <code>/dependencies</code> directory.  The <code>/dependencyDownloads</code> sub-directory contains the tarball or zipball file of all required packages.  In the <code>/installscripts</code> sub-directory you can find Python scripts that help you to unzip or untar these files into the <code>/dependencyInstalls</code> directory.  You can do this for one dependency at a time, like this (in the directory in which PlatoSim3 was downloaded):

 \code{.py}
python ./dependencies/installscripts/install_hdf5.py
python ./dependencies/installscripts/install_googletest.py
python ./dependencies/installscripts/install_yaml.py
python ./dependencies/installscripts/install_armadillo.py
python ./dependencies/installscripts/install_fftw.py
python ./dependencies/installscripts/install_faddeeva.py
\endcode

Alternatively, you can also install the required dependencies and build the code in one go by typing (also in the directory in which PlatoSim3 was downloaded): 

\code
./install.sh
\endcode

If problems would arise when executing this command, it may be useful to tried to install the dependencies one-by-one to pinpoint the problem.

Building the software will create two executables :

* <code>platosim</code> to run simulations (see below)
* and <code>testplatosim</code> to run the test harnesses (without arguments).


---


## Software Changes

In case of code changes (after retrieving the latest version from GitHub or after introducing changes to the code yourself)

In case you have updated the PlatoSim3 code but the dependencies remain unchanged, you only have to re-build the software but not resolve the dependencies again.  This saves you a tremendous amount of time.  The way to do this, is:

\code
cmake ..
(make clean)
make -j 4
\endcode


---


## Updated Dependencies

At some stage, we will want to update (some of) the dependencies.  You will be notified by the developer team in case this happens.  You will then have to:

* clear the <code>/dependencies/Installs</code> directory
* and run the install script again (as described above for the first-time installation).


---

## Running the Test Harnesses

In order to be able to run the test harnesses, you must first build the code (see above) and export the required environment variables, as explained @ref ReqsRun "here".

The actual command to run the tests must be executed in the <code>/build</code> directory:

\code
./testplatosim
\endcode


---

## Troubleshooting

### Not Using the System Default Compiler

If you want to use a different compiler than the system default to execute the steps described above, you have to export two additional environment variables, <code>CCX</code> and <code>CC</code>, as follows:

\code
export CXX=g++-5
export CC=gcc-5
\endcode

where you may want to adapt the right-hand side of these two lines to the compiler (version) of your choice.

### Still Experiencing Problems?

If you are still experiencing problems following the instructions above, please, @ref IssueTracking "tell us about it"!  It's convenient if you can send us the full error log, which you can get hold of as follows:

\code
./install.sh > output.txt 2> errors.txt
\endcode

---

## Switching between Branches

As you can read @ref Branching "here", we no longer only use the master branch.  If you switch to another branch and want to run simulations with the current branch, you will have to re-build the software.