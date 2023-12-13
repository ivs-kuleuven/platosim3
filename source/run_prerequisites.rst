Prerequisites
=============

.. raw:: html

   <hr>
   
.. _run_prerequisites_export:
   
*Export environment variables*
------------------------------

To avoid having to hardcode any paths in configuration files, tutorials, etc., you must export three environment variables:

* ``PLATO_PROJECT_HOME``: refer to the directory in which PlatoSim3 was installed
* ``PLATO_WORKDIR``: refer to your working directory ``<path/to/plato_workdir>`` were you should store your own configuration files and save your simulations (preferably not within ``PlatoSim3/``)
* ``PYTHONPATH``: refer to the directory in which our Python scripts can be found

Before exporting any environment variables, first make sure to activate your Conda environment. If you have installed PlatoSim as a :doc:`user (via Conda) <install_conda>`, the home location of the PlatoSim project needs to be exported as:

.. code-block:: shell

   export PLATO_PROJECT_HOME=$CONDA_PREFIX

On the other hand, if you installed PlatoSim as a :doc:`developer (from Source) <install_source>`, this environment variable needs to be exported as:
   
.. code-block:: shell

   export PLATO_PROJECT_HOME=<full/path/to/PlatoSim3>

Next, choose a PlatoSim working directory and export the two following environment variables (which are independent on if you are a user or a developer):
   
.. code-block:: shell
		
   export PLATO_WORKDIR=<full/path/to/working/directory>
   export PYTHONPATH=$PYTHONPATH:$PLATO_PROJECT_HOME/python
   
Lastly, check that the exported paths are set to the proper location with:

.. code-block:: shell

   echo $PLATO_PROJECT_HOME
   echo $PLATO_WORKDIR
   echo $PYTHONPATH
 
We recommend that you place these three environment variables in a ``.bash_profile`` script that can be read by your system at all times.

.. note::

   You should now be able to run PlatoSim from the command line:

   .. code-block:: shell

      platosim -h

   and from a Python shell:

   .. code-block:: python

      import platosim

.. raw:: html

   <hr>






.. _run_prerequisites_storage:
   
*File storage*
--------------

To avoid problems when updating the PlatoSim software, it is best to store your own input and output files in a designated working directory -- i.e. ``PLATO_WORKDIR`` as described above. For developers we recommend to place this directory outside the ``PlatoSim3/`` directory. You can (but should not) add your input files to the ``inputfiles/`` directory, but under no circumstances change the original files in that directory!

In case - as a user - you install a new version of the software in an existing Conda environment, all changes will be overwritten and newly added files will be removed from the intallation. So better store them outside of the installation. For all users that want to run PlatoSim from the terminal, we recommend that you as a developer copy the ``inputfile.yaml`` file and modify the copy rather than the original file. For users you can simply download the ``inputfile.yaml`` directly from GitHub.
