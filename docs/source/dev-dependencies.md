# Dependencies{#dev-dependencies}

PlatoSim3 relies on a number of other packages, which are all included in the PlatoSim3 distribution for your convenience.  Everything concerning the dependencies can be found in the <code>/dependencies</code> directory.  The <code>/dependencyDownloads</code> sub-directory contains the tarball or zipball file of all required packages.  In the <code>/installscripts</code> sub-directory you can find Python scripts that help you to unzip or untar these files into the <code>/dependencyInstalls</code> directory.  You can do this for one dependency at a time, like this (in the directory in which PlatoSim3 was downloaded):

    $ python ./dependencies/installscripts/install_hdf5.py
    $ python ./dependencies/installscripts/install_googletest.py
    $ python ./dependencies/installscripts/install_yaml.py
    $ python ./dependencies/installscripts/install_armadillo.py
    $ python ./dependencies/installscripts/install_fftw.py
    $ python ./dependencies/installscripts/install_faddeeva.py

Alternatively, you can also install the required dependencies and build the code in one go by typing (also in the directory in which PlatoSim3 was downloaded):

    $ ./install.sh

If problems would arise when executing this command, it may be useful to tried to install the dependencies one-by-one to pinpoint the problem.

Running this script will create two executables :

* <code>platosim</code> to run simulations (see below)
* and <code>testplatosim</code> to run the test harnesses (without arguments).
