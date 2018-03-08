# Installing PlatoSim3 Via Conda {#InstallViaConda}

Before you install another version of PlatoSim3, you must activate the desired conda environment, as described @ref ReqsInstallViaConda "here".  It is not necessary to create a new conda environment every time you install a different version of the software, unless you want to use multiple version in parallel.

The installation procedure will automatically detect which operating system your are running and will install the appropriate packages for you.

To install the latest successfully built version of the <code>master</code>, type:

\code
conda config --add channels conda-forge
conda install -c http://www.miricle.org/platosim/ platosim
\endcode

For the develop branch, the latter command must be replaced by

\code
conda config --add channels conda-forge
conda install -c http://www.miricle.org/platosim.devel/ platosim
\endcode

Please, contact the developer team for the username and password.

We will, at a later stage, at information on how to install a specific version (rather than the latest one).