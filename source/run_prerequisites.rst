Prerequisites
=============

.. raw:: html

   <hr>



   
.. _run_prerequisites_export:
   
*Export Environment Variables*
------------------------------

To avoid having to hardcode any path in configuration files, tutorials, etc., you must export three environment variables:

* ``PLATO_PROJECT_HOME``: refer to the directory in which PlatoSim3 was installed
* ``PLATO_WORKDIR``: refer to the working directory ``<path/to/plato_workdir>`` were you should store your own configuration files and save your simulations (preferably not within ``PlatoSim3/``)
* ``PYTHONPATH``: refer to the directory in which our Python scripts can be found

For simplicity you can easily export the above environment variables using the ``platosim`` executeable from the terminal. From the PlatoSim base directory simply run:

.. code-block:: shell

   platosim --setup <path/to/plato_workdir>

The script will finalize your setup doing the following:

* Export the environment variable ``$PLATO_PROJECT_HOME``
* Export the environment variable ``$PLATO_WORKDIR``
* Export the environment variable ``$PYTHONPATH``
* Making all PLATOnium scripts globally executable
* Create a ``.bash_profile`` file within your ``$PLATO_WORKDIR`` directory
  
If you want to change your working directory at a later stage, simply remove the ``.bash_profile`` file and run above command again. More information on how to use the scripts of PLATOnium can be found the :doc:`Tutorials <platonium_tutorials>`.

.. admonition:: Troubleshooting

    If you experience any problems with the environment variables, you can try to check if the paths are set correctly. A correct export looks like:

    .. code-block:: shell

       export PLATO_PROJECT_HOME=<full/path/to/PlatoSim3>
       export PLATO_WORKDIR=<full/path/to/working/directory>
       export PYTHONPATH=$PYTHONPATH:$PLATO_PROJECT_HOME/python

    In case you've installed PlatoSim via Conda, the former environment variable should be exported as:

    .. code-block:: shell

       export PLATO_PROJECT_HOME=$CONDA_PREFIX

    (the ``CONDA_PREFIX`` environment variable is automatically known when you activate the appropriate conda environment).

    To check the content of these variables (to see whether they are set to the proper location), type:

    .. code-block:: shell

       echo $PLATO_PROJECT_HOME
       echo $PLATO_WORKDIR
       echo $PYTHONPATH

.. raw:: html

   <hr>




   
.. _run_prerequisites_storage:
   
*File storage*
--------------

To avoid problems when updating the PlatoSim software, it is best to store your own input and output files in a designated working directory -- i.e. ``PLATO_WORKDIR`` as described above. For developers we recommend to either place your outside the installation directory of PlatoSim3. You can (but should not) add your input files to the ``inputfiles/`` directory, but under no circumstances change the original files in that directory!

In case - as a user - you install a new version of the software in an existing conda environment, all changes will be overwritten and newly added files will be removed from the intallation. So better store them outside of the installation. For all users that want to run PlatoSim from the terminal, we recommend that you as a developer copy the ``inputfile.yaml`` file and modify the copy rather than the original file. For users you can simply download the ``inputfile.yaml`` directly from GitHub.
