# Prerequisites for Installing PlatoSim3: Conda Environments {#ReqsInstallViaConda}

## Creating Conda Environments

To be able to install PlatoSim via conda, you have to have the Python distribution <a href="https://docs.continuum.io/anaconda/install">Anaconda</a> installed.  You have to create an <a href ="https://conda.io/docs/user-guide/tasks/manage-environments.html">conda environments</a>, here called <code>platosim</code>, as follows:

\code
conda create -n platosim python=3.5 anaconda
\endcode

but it is advisable to use multiple conda environments if you want to be able to switch between version and/or branches in a smooth way (e.g. <code>platosimMaster</code> and <code>platosimDevel</code>).

To get an overview of all your conda environments, type

\code
conda env list
\endcode

The active environments will be marked with *.

## Activating and De-activating Conda Environments

Each time you install a different version of PlatoSim or if you want to use a specific version of PlatoSim you have installed (via conda) on your system, you have to activate the appropriate conda environment, which is done as follows:

\code
source activate <environment name>
\endcode

To de-activate the environment you are currently on, type

\code
source deactivate
\endcode