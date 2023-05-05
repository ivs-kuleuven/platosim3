# Prerequisites for Installing PlatoSim3: Conda Environments {#user-prerequisites}

## Creating Conda Environments

To be able to install PlatoSim via conda, you have to have the Python distribution <a href="https://docs.continuum.io/anaconda/install">Anaconda</a> installed.  You have to create an <a href ="https://conda.io/docs/user-guide/tasks/manage-environments.html">conda environments</a>, here called <code>platosim</code>, as follows:

\code
conda create -n platosim python=<Python version>
\endcode

but it is advisable to use multiple conda environments if you want to be able to switch between version and/or branches in a smooth way (e.g. <code>platosimMaster</code> and <code>platosimDevel</code>).  Supported Python versions are 3.5, 3.6, and 3.7.  Please, note that - when you switch to a different version of Python - it is advised to create a new conda environment rather than trying to update your existing one.  It will save you a lot of trouble if you do it like this.


To get an overview of all your conda environments, type

    $ conda env list

The active environments will be marked with *.

---

## Activating and De-activating Conda Environments

Each time you install a different version of PlatoSim or if you want to use a specific version of PlatoSim you have installed (via conda) on your system, you have to activate the appropriate conda environment, which is done as follows:

    $ conda activate <environment name>

To de-activate the environment you are currently on, type

    $ conda deactivate
