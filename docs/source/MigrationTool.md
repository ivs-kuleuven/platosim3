# Migration Tool {#MigrationTool}

At the moment the comments behind the configuration parameters are unfortunately not conserved.When new versions of PlatoSim are released, there might be some changes in the configuration or input files.  New parameters might have been added, or parameters or groups might have been re-structured.  In the PlatoSim distribution we provide a helper tool (<code>migtool.py</code> in the <code>/python</code> folder of the PlatoSim installation) to migrate your YAML inputfiles to the new format.  This tool will compare your (old) YAML input file with the default YAML input file that is on your local system, i.e. <code>/inputfiles/inputfile.yaml</code>.  The command can be used as follows (assuming your are located in the PlatoSim project home folder):

\code
python python/migtool.py [-hv] [-o outputFilename] inputFilename
\endcode

The <code>-h</code> option prints an instructive help message. The <code>inputFilename</code> is your (old) YAML inputfile. The <code>-o</code> option lets you specify the name of the output file in which the migrated configuration will be saved. When no output file is given, the result will be printed on the screen (<code>stdout</code>).
 The <code>-v</code> option prints the changes that will be applied on the screen. That might be useful, because it will signal changes in parameter values. An example is shown below. 
 
 <code>CHECK</code> means you will have to check manually. The value is your input file was <code>0.016</code> and the value in the new input file is <code>0.01</code>. Your value will be retained after the migration. The reason is that you probably have good reasons to have this value different from the default value, and you don't want to loose that during the migration.
 
\code 
CHECK - Value changed for CCD.FlatfieldPtPNoise from 0.016 to 0.01
\endcode

 Other possible output when using the <code>-v</code> option is signaled with <code>DONE</code>, i.e. when a new key or sub-key has been added or when an obsolete key has been removed.
 Unfortunately, the Python module that we use to parse the YAML files does not retain the comments that are in the files. Those comments will be lost after migration.