Run simulations
===============

To run simulations with PlatoSim, you will have to feed a :ref:`Configuration file <run_input_parameters>` as input to the simulator, possibly together with a few additional :ref:`Supplementary files <run_input_files>`. Launching simulations can be done either from the :ref:`command line <run_command_line>` or in :ref:`Python <run_python>` as explained on this page.

.. raw:: html

   <hr>




.. _run_simulations_command:

*Command line*
--------------
   
To run a simulation from the command line, simply type:

.. code-block:: shell
		
   platosim <inputFile> <outputFile> [<logFile> <logLevel>]

The log level is an integer value that allows to configure what kind of log messages is written to the log file:

1. messages: ``error``, ``warning`` (default);
2. messages: ``error``, ``warning``, ``info``;
3. messages: ``error``, ``warning``, ``info``, ``debug``.
	 
A simple test simulation looks like:

.. code-block:: shell
		
   platosim inputfile.yaml out.hdf5

This will generate two files: the HDF5 output file containing the simulation and a ascii log file with the log messages from PlatoSim.
   
.. raw:: html

   <hr>



   
.. _run_simulations_python:

*Python*
--------

`A full suite of Jupyter notebook tutorials are made available <https://github.com/IvS-KULeuven/PlatoSim3/tree/master/docs/tutorials>`_. Below we provide a minimal example of how to get started running PlatoSim using our Python interface. If you're planning on creating your own scripts in Python, then besides reading our Jupyter notebooks, we recommend spending a bit of time investigating the different built in `Python classes (and their utilities) <https://github.com/IvS-KULeuven/PlatoSim3/tree/master/python/platosim>`_.

We note that the :doc:`PLATOnium toolkit <platonium_overview>` can be used to setup a large scale simulation of imagttes/photometric time series within minutes. This toolkit can be used to generate stellar catalogue from the PIC, generate variable signals, generate instrumental effects, and takes care of the configuration of payload inline with the future observations of PLATO. This toolkit is meant for plug-and-play and in particalar to be used on a computing cluster.

.. raw:: html

   <hr>
   
**Create a simulation object**

To run PlatoSim from within Python, you have to import the simulation module:

.. code-block:: python
		
   import platosim.simulation as Simulation

and create a new Simulation object (the name ``sim`` is arbitrary):

.. code-block:: python
		
   sim = Simulation(<runName>,
                    configurationFile = <full path to the YAML configuration file>,
                    outputDir = <full path to the output directory>)

Note that you can use the working directory (i.e. the directory specified in the ``PLATO_WORKDIR``) environment variable to specify the path to the input configuration file or the output directory. To read the content of this environment variable in Python, use:

.. code-block:: python
		
   import os
   workDir = os.getenv("PLATO_WORKDIR")


**Input configuration file**

The ``configurationFile`` parameter is optional and contains the full path to the configuration file that will be used as a reference. A copy of this file will be stored in the output directory (with the name ``<outputFile>.yaml``) and only the copy (and not the original file) will be modified. We will show later how to do this.

It is possible to also set the path to the reference input file as follows:

.. code-block:: python
		
   sim.readConfigurationFile(<full path to the configuration file>)

In that case, all changes to the configuration parameters you may have applied, will be overwritten.


**Output directory**

Also the ``outputDir`` parameter is optional. This contains the full path to the output directory. This must be specified before you can start the simulation. The output directory can be set can also be set as follows:

.. code-block:: python
		
   sim.outputDir = <full path to the output directory>

In this directory, the following information will be stored when running the simulation:

  * ``<outputFile>.yaml``: copy of the input file with the chosen configuration parameters (copied from the reference configuration file and possibly modified, as shown below);
  * ``<outputFile>.hdf5``: resulting exposures (images, PSF, etc.);
  * ``<outputFile>.log``: log file (to report any problems).

   
**Modify configuration parameters**

Individual configuration parameters can be modified in the copy of the configuration file (which will be stored with the name runName.yaml in the output directory) as follows:

.. code-block:: python
		
   sim["<block>/<parameter>"] = <new value>

for individual parameters in the block of the YAML file (e.g. ``sim["Camera/FocalPlaneOrientation"] = 0``), or:

.. code-block:: python
		
   sim["<block>/<sub-block>/<parameter>"] = <new value>

for parameters in sub-block in the YAML file (e.g. ``sim["Platform/Orientation/Angles/RAPointing"] = 90.0``). Note that the original YAML file is not modified. Only the copy of this fill gets modified. In case you change the configuration file of the Simulation object, the copy ``runName.yaml`` file will be overwritten and all changes to the configuration parameters you have applied, will be lost.

**Run the simulation**

The simulation can be launched with the following command:

.. code-block:: python

   sim.run()

By setting the optional input parameter ``removeOutputFile`` to ``True``, any previous simulatedata with the same name within a given folder will be overwritten:

.. code-block:: python
		
   sim.run(removeOutputFile = True)

.. raw:: html

   <hr>



   


.. _run_simulations_error:

*Error codes*
-------------
		
If PlatoSim raises an error message upon execution, please have a look the log messages before raising a new GitHub issue. The two most frequent *error codes* are:

- **Error code -1**: The HDF5 file already exists, hence remove it and try again
- **Error code 6**: Wrong formatting or typesetting of the YAML file or a supplementary file
