# Scripting in Python {#ScriptingInPython}

The basic functionality of the C++ code is to model a part (i.e. sub-field) of one of the four CCDs (A, B, C, or D) that are attached to the camera of one telescope on the platform.  You can initiate Plato simulations from within Python, e.g. to automatically select the correct CCD or to combine multiple telescopes on the platform.

Under <code>/docs/tutorials/ScriptingInPython</code> you can find a basic tutorial on scripting PlatoSim in Python.





## Setting your Python Path


The Python functionality that is offered in PlatoSim can be found in the <code>python</code>.  Make sure to add this directory to your Python path, as follows:

\code{python}
PYTHONPATH=$PYTHONPATH:<full path to the /Platosim3/python directory>
\endcode



## Creating a <code>simulation</code> Object

To run PlatoSim from within Python, you have to import the <code>simulation</code> module:

\code{python}
import simulation
\endcode

and create a new <code>Simulation object</code> (the name <code>sim</code> is arbitrary):

\code{python}
sim = Simulation(<runName>, configurationFile = <full path to the configuration file>, outputDir = <full path to the output directory>)
\endcode





### Input Configuration File

The <code>configurationFile</code> parameter is optional and contains the full path to the configuration file that will be used as a reference.  A copy of this file will be stored in the output directory (with the name <code>runName.yaml</code>) and only the copy (and not the original file) will be modified.  We will show later how to do this.

If not specified, the <code>inputfile.yaml</code> configuration file from the <code>/inputfiles</code> directory will be used (provided you have exported the <code>PLATO_PROJECT_HOME</code> evironment variable).  It is possible to also set the path to the reference input file as follows:

\code{python}
sim.readConfigurationFile(<full path to the configuration file>)
\endcode

### Output

Also the <code>outputDir</code> parameter is optional.  This contains the full path to the output directory.  If this is not specified, the current directory will be used to store the output files.  The output directory can be set can also be set as follows:

\code{python}
sim.outputDir = <full path to the output directory>
\endcode

In this directory, the following information will be stored when running the simulation:

- <code>runName.yaml</code>: copy of the input file with the chosen configuration parameters (copied from the reference configuration file and possibly modified, as shown below);
- <code>runName.hdf5</code>: resulting exposures (images, PSF, etc.);
- <code>runName.log</code>: log file (to report any problems).

## Modifying Configuration Parameters

Individual configuration parameters can be modified in the copy of the configuration file (which will be stored with the name <code>runName.yaml</code> in the output directory) as follows:

\code{python}
sim["<block>/<parameter>"] = <new value>
\endcode

for individual parameters in the block of the YAML file (e.g. <code>sim["Camera/FocalPlaneOrientation"] = 0</code>), or

\code{python}
sim["<block>/<sub-block>/<parameter>"] = <new value>
\endcode

for parameters in sub-block in the YAML file (e.g. <code>sim["PSF/Gaussian/Sigma"] = 0.50</code>).

A detailed description of the configuration parameters can be found @ref InputFileDescription "here".

Note that the original YAML file is not modified.  Only the copy of this fill gets modified.

## Running the Simulation

The simulation can be executed with the following command:

\code{python}
sim.run()
\endcode

By setting the optional input parameter <code>removeOutputFile</code> to <code>True</code>, any files with the same names in the specified output directory will be deleted first:

\code{python}
sim.run(removeOutputFile = True)
\endcode
