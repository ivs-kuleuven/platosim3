# Installing (and Updating) PlatoSim3 Via Conda {#InstallViaConda}

Before you install another version of PlatoSim3, you must activate the desired conda environment, as described @ref ReqsInstallViaConda "here".  It is not necessary to create a new conda environment every time you install a different version of the software, unless you want to use multiple version in parallel.

The installation procedure will automatically detect which operating system your are running and will install the appropriate packages for you.

Before you install PlatoSim via conda for the first time in this environment, type:

\code
conda config --add channels conda-forge
\endcode

To install the latest successfully built version of the <code>master</code>, type:

\code
conda install -c http://www.miricle.org/platosim/ platosim
\endcode

For the <code>develop</code> branch, the latter command must be replaced by

\code
conda install -c http://www.miricle.org/platosim.devel/ platosim
\endcode

To install a specific version (only for the <code>master</code> branch), just append <code>=<version></code> to this command.

To update an installed version to the latest one, replace <code>conda install</code> by <code>conda update</code>.

Please, contact the developer team for the username and password.

If no pop-up window, asking for the credentials, would appear, you can adapt the <code>conda install</code> commands from above, by placing <code><username>:<password></code>@ between the <code>http://</code> and the <code>www</code>. 

