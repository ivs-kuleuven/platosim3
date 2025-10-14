Quickstart
==========

.. warning::

   PlatoSim is currently only available to memebers of the PLATO Mission Consortium (with a signed Non-Disclosure Agreement: NDA). Thus, to install (or download) the PlatoSim software, please contact one of the PlatoSim developers (see our `GitHub page <https://github.com/IvS-KULeuven/PlatoSim3>`_) for the credentials: **<username>** and **<password>**.

The easiest way to install PlatoSim is through our :doc:`Conda installation <install_conda>` suing the following set of terminal commands:

.. code-block:: shell
		
   conda create -n platosim python=3.10
   conda activate platosim
   conda install -c https://<username>:<password>@jenkins.miricle.org/platosim/ platosim

Create a working directory for which you want to store all your future PlatoSim projects:

.. code-block:: shell

   mkdir -p ~/platosim_workdir/quickstart
   
Export the following paths (and add them to your bash file for future sessions):

.. code-block:: shell
		
   export PLATO_PROJECT_HOME=$CONDA_PREFIX
   export PLATO_WORKDIR=~/plato_workdir
   export PYTHONPATH=$PYTHONPATH:$PLATO_PROJECT_HOME/python
   
PlatoSim also requires a PSF file (which is not included in the default distribution because of its size), which you can download from our FTP site:  

.. code-block:: shell
		
   wget ftp://<username>:<password>@ftp.ster.kuleuven.be/PSF_Focus_0mu.hdf5 $PLATO_PROJECT_HOME/inputfiles

A simple test simulation from the command line looks like:

.. code-block:: shell

   platosim $PLATO_PROJECT_HOME/inputfiles/inputfile.yaml $PLATO_WORKDIR/quickstart/test.hdf5

A simple test simulation from Python looks like:

.. code-block:: python

   import os
   import platosim.simulation as Simulation
   odir = os.getenv('PLATO_WORKDIR') + '/quickstart'
   sim = Simulation('test', outputDir=odir)
   sim.run()
   
We also recommend using :doc:`Platonium <platonium_overview>`, a user-friendly PlatoSim toolkit with lots of extra functionalities that help you create and run simulation projects very efficiently. Simply download this `pyproject.toml <https://github.com/IvS-KULeuven/PlatoSim3/blob/master/pyproject.toml>`_ file and follow :doc:`Platonium's prerequisites <platonium_prerequisites>`.
