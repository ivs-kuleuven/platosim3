Prerequisites
=============

We advise you to install `Anaconda <https://docs.continuum.io/anaconda/install/>`_ as Python distribution, as this already has most of the required packages installed. It comes with a GUI, called Spyder, from which you can run your Python scripts.




.. raw:: html

   <hr>
   
.. _python_users:

*For users*
-----------


Users only have to source the conda environment in which the desired version of the software has been installed in. Also they are advised to create a new conda environment when they want to switch to a different version of Python. After a conda update or a conda install of the PlatoSim software, everything will be configured appropriately, without further action.




.. raw:: html

   <hr>
   
.. _python_developers:

*For developers*
----------------

For developers, to be able to run the PlatoSim and inspect the output, a few Python libraries are needed. Follow the steps below to successfully install all Python libraries needed to run PLATOnium. The blue colored admonition boxes elaborate further steps to be taken if installing on a computing cluster:

1. **Install PlatoSim3 from source**

   Install PlatoSim following :ref:`PlatoSim's developer branche <install_source>`.

   .. admonition:: Installing PlatoSim on a cluster

      Computing clusters typically have a high standard of version control for software and hence the user cannot install software directly, but needs to load pre-installed software **modules**.
	 
2. **Create a Conda environment**
   
   Again make sure that the relevant modules are loaded beforehand. We strongly recommend to install Poetry through a `Conda <https://docs.conda.io/en/latest/>`_ environment, since computing clusters typically only have limited amount of Python versions installed by default. While using Conda all Python versions can be installed and thus an exact freeze of the poetry installation can be made. Thus first install Conda (or `Miniconda <https://docs.conda.io/en/latest/miniconda.html>`_) and then create a new environment e.g. called ``platosim``: 

   .. code-block:: shell
		
      conda create -n platosim python=<Python version>

   Supported Python versions are 3.8 and 3.9. Please, note that when you switch to a different version of Python, it is advised to create a new conda environment rather than trying to update your existing one. This means you (as a developer) will have to re-install all Python packages (see below), but is will save you a lot of trouble if you do it like this.
Now activate your new conda environment:

   .. code-block:: shell

      conda activate platosim
      
3. **Install Python packages with Poetry**

   Since the force of PLATOnium comes from it usage to running simulations on a computing cluster we use `Poetry <https://python-poetry.org/>`_ to handle the package managing (which is far better than conda) as it manage your project installation easily over multiple platforms in a deterministic way.

   If not done so, first install Poetry:

   .. code-block:: shell
  
      curl -sSL https://install.python-poetry.org | python -
   
   .. admonition:: Poetry prerequiste on cluster
	 
      Verify that the python location is correct with ``which python``. Note that clusters typically don't support much disk space for ``$HOME``, and hence it is better to install Poetry on your ``$DATA`` disk storage instead. You can parse the installation location of Poetry (``</path/to/poetry>``) with:

      .. code-block:: shell

	 curl -sSL https://install.python-poetry.org | POETRY_HOME=</path/to/poetry> python -

      Verify that poetry was installed successfully by typing ``poetry --version`` and verify the installation location with ``which poetry``. Next change the installation location of the virtual Poetry environment to:
      
      .. code-block:: shell

	 poetry config virtualenvs.path </path/to/poetry>/virtualenvs

      In order for Poetry to be available from any compute node, you need to include the following path to your ``~/.bashrc`` file:
   
      .. code-block:: shell

	 POETRY=<path/to/poetry>/bin
	 export POETRY

   From the base of the PlatoSim repository, simply install the Python packages with:

   .. code-block:: shell

      poetry install



		
