Migration tool
==============

At the moment the comments behind the configuration parameters are unfortunately not conserved. When new versions of PlatoSim are released, there might be some changes in the configuration or input files. New parameters might have been added, or parameters or groups might have been restructured. While using PlatoSim's Python interface such changes are completely transparent if you simply use the default ``inputfile.yaml`` file and simply modify this through your simulation object in your Python scripts.

If you run PlatoSim from the command line and have already existing YAML configuration files, you will potentially have issues (``errorcode 6`` will show up). In the PlatoSim distribution we provide a helper tool to migrate your YAML inputfiles to the new format. The tool `migtool.py <https://github.com/IvS-KULeuven/PlatoSim3/blob/develop/python/platosim/script/migtool.py>`_ in the ``/python/platosim/script`` folder of the PlatoSim installation. This tool will compare your (old) YAML input file with the default YAML input file that is on your local system, i.e. ``/inputfiles/inputfile.yaml``. The command can be used as follows (assuming your are located in the PlatoSim project home folder):

.. code-block::
   
   migtool.py [-hv] [-o outputFilename] inputFilename

The ``-h`` option prints an instructive help message. The ``inputFilename`` is your (old) YAML inputfile. The ``-o`` option lets you specify the name of the output file in which the migrated configuration will be saved. When no output file is given, the result will be printed on the screen (stdout). The ``-v`` option prints the changes that will be applied on the screen. That might be useful, because it will signal changes in parameter values. An example is shown below.

``CHECK`` means you will have to check manually. The value is your input file was 0.016 and the value in the new input file is 0.01. Your value will be retained after the migration. The reason is that you probably have good reasons to have this value different from the default value, and you don't want to loose that during the migration.

.. code-block::
   
   CHECK - Value changed for CCD.FlatfieldPtPNoise from 0.016 to 0.01

Other possible output when using the ``-v`` option is signaled with ``DONE``, i.e. when a new key or sub-key has been added or when an obsolete key has been removed. Unfortunately, the Python module that we use to parse the YAML files does not retain the comments that are in the files. Those comments will be lost after migration.
