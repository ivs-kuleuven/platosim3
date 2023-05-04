Tutorials
=========

`A full suite of Jupyter notebook tutorials are ready for the user to explore <https://github.com/IvS-KULeuven/PlatoSim3/tree/master/docs/tutorials>`_. Below we provide a minimal example of how to get started running PlatoSim using our Python interface.

If you're planning on creating your own scripts in Python, then besides reading our Jupyter notebooks, we recommend spending a bit of time investigating the different built in `Python classes (and their utilities) <https://github.com/IvS-KULeuven/PlatoSim3/tree/master/python/platosim>`_.

We note that the :ref:`PLATOnium toolkit <python_platonium>` in the next (and last) section can be used to setup a large scale simulation of imagttes/photometric time series within minutes. This toolkit takes care of a realistic PIC catalogue, generating variable signals, generating instrumental effects, and takes care of the configuration of payload inline with the future observations of PLATO. This toolkit is meant for plug-and-play and in particalar running simulations on a computing cluster. Furthermore, PLATOnium can be used in combination with the PLATO pipeline reduction chain (L1 pipeline).


.. raw:: html

   <hr>
   
*Running PlatoSim3 in Python*
-----------------------------

**Creating a Simulation Object**

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
   workDir = os.environ["PLATO_WORKDIR"]


**Input Configuration File**

The configurationFile parameter is optional and contains the full path to the configuration file that will be used as a reference. A copy of this file will be stored in the output directory (with the name ``runName.yaml``) and only the copy (and not the original file) will be modified. We will show later how to do this.

If not specified, the inputfile.yaml configuration file from the ``/inputfiles`` directory will be used (provided you have exported the ``PLATO_PROJECT_HOME`` environment variable). It is possible to also set the path to the reference input file as follows:

.. code-block:: python
		
   sim.readConfigurationFile(<full path to the configuration file>)

In that case, all changes to the configuration parameters you may have applied, will be overwritten.


**Output Directory**

Also the outputDir parameter is optional. This contains the full path to the output directory. This must be specified before you can start the simulation. The output directory can be set can also be set as follows:

.. code-block:: python
		
   sim.outputDir = <full path to the output directory>

In this directory, the following information will be stored when running the simulation:

  * ``runName.yaml``: copy of the input file with the chosen configuration parameters (copied from the reference configuration file and possibly modified, as shown below);
  * ``runName.hdf5``: resulting exposures (images, PSF, etc.);
  * ``runName.log``: log file (to report any problems).

   
**Modifying Configuration Parameters**

Individual configuration parameters can be modified in the copy of the configuration file (which will be stored with the name runName.yaml in the output directory) as follows:

.. code-block:: python
		
   sim["<block>/<parameter>"] = <new value>

for individual parameters in the block of the YAML file (e.g. ``sim["Camera/FocalPlaneOrientation"] = 0``), or:

.. code-block:: python
		
   sim["<block>/<sub-block>/<parameter>"] = <new value>

for parameters in sub-block in the YAML file (e.g. ``sim["PSF/Gaussian/Sigma"] = 0.50``). Note that the original YAML file is not modified. Only the copy of this fill gets modified. In case you change the configuration file of the Simulation object, the copy ``runName.yaml`` file will be overwritten and all changes to the configuration parameters you have applied, will be lost.

**Running the Simulation**

The simulation can be executed with the following command:

.. code-block:: python

   sim.run()

By setting the optional input parameter ``removeOutputFile`` to ``True``, any files with the same names in the specified output directory will be deleted first:

.. code-block:: python
		
   sim.run(removeOutputFile = True)

   
