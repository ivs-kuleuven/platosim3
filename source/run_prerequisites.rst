Prerequisites
=============

.. raw:: html

   <hr>
   
.. _run_prerequisites_export:
   
*Export Environment Variables*
------------------------------

To avoid having to hardcode any paths in configuration files, tutorials, etc., you must export three environment variables:

* ``PLATO_PROJECT_HOME``: refer to the directory in which PlatoSim3 was installed
* ``PLATO_WORKDIR``: refer to the working directory ``<path/to/plato_workdir>`` were you should store your own configuration files and save your simulations (preferably not within ``PlatoSim3/``)
* ``PYTHONPATH``: refer to the directory in which our Python scripts can be found

We recommend that you place the following path-exports in a ``.bash_profile`` script that can be read by your system at all times:
  
.. code-block:: shell

   export PLATO_PROJECT_HOME=<full/path/to/PlatoSim3>
   export PLATO_WORKDIR=<full/path/to/working/directory>
   export PYTHONPATH=$PYTHONPATH:$PLATO_PROJECT_HOME/python

In case you've installed PlatoSim via Conda, first activate your Conda environment, and export the former environment variable as:

.. code-block:: shell

   export PLATO_PROJECT_HOME=$CONDA_PREFIX

To check the content of these variables (to see whether they are set to the proper location), type:

.. code-block:: shell

   echo $PLATO_PROJECT_HOME
   echo $PLATO_WORKDIR
   echo $PYTHONPATH
 
Note that developers can take advantage of the :doc:`setup script for PLATOnium <platonium_overview>`.

   
.. raw:: html

   <hr>
   
.. _run_prerequisites_storage:
   
*File storage*
--------------

To avoid problems when updating the PlatoSim software, it is best to store your own input and output files in a designated working directory -- i.e. ``PLATO_WORKDIR`` as described above. For developers we recommend to either place your outside the installation directory of PlatoSim3. You can (but should not) add your input files to the ``inputfiles/`` directory, but under no circumstances change the original files in that directory!

In case - as a user - you install a new version of the software in an existing conda environment, all changes will be overwritten and newly added files will be removed from the intallation. So better store them outside of the installation. For all users that want to run PlatoSim from the terminal, we recommend that you as a developer copy the ``inputfile.yaml`` file and modify the copy rather than the original file. For users you can simply download the ``inputfile.yaml`` directly from GitHub.
