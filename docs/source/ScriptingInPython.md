# Scripting in Python {#ScriptingInPython}

The basic functionality of the C++ code is to model a part (i.e. sub-field) of one of the four CCD that are attached to the normal camera (A, B, C, or D) or to the fast camera (AF, BF, CF, or DF) of one of the telescopes on the platform.  You can initiate PLATO simulations from within Python, e.g. to automatically select the correct CCD or to combine multiple telescopes on the platform.

Under <code>/docs/tutorials/ScriptingInPython</code> you can find a basic tutorial on scripting PlatoSim3 in Python.

We advise you to install <a href="https://docs.continuum.io/anaconda/install">Anaconda</a> as Python distribution, as this already has most of the required packages installed.  It comes with a GUI, called Spyder, from which you can run your Python scripts.


## Configuring Python

### Creating a Python Environment

With <a href="http://conda.pydata.org/docs/index.html">conda</a> (which is included in <a href="https://docs.continuum.io/anaconda/install">Anaconda</a>), you can create, export, list, remove, and update environments that have different versions of Python and/or packages installed in them. Switching or moving between environments is called activating the environment.

To create and activate such an environment, type

\code
source /anaconda/bin/activate root
conda create -n platoSim python=3.5 anaconda
\endcode

More information on managing environment with conda can be found <a href="http://conda.pydata.org/docs/using/envs.html#">here</a>.

### Installing Extra Packages

To be able to run the PLATO Simulator and inspect the output, you must install the following packages:

- <a href="http://www.numpy.org/">numpy</a>, a fundamental package for scientific computing in Python,
- <a href="http://matplotlib.org/">matplotlib</a>, a 2D plotting package, 
- <a href="https://pypi.python.org/pypi/pyaml">pyaml</a>, to deal with the configuration files in <a href="https://learnxinyminutes.com/docs/yaml/">YAML</a> format,
- and <a href="http://www.pytables.org/">pytables</a>, to inspect the output files in <a href="https://www.hdfgroup.org/HDF5/">HDF5</a> format.

These packages can be installed with the following commands:

\code
conda install numpy
conda install matplotlib
conda install pytables
pip install pyaml
\endcode

As it comes to the choice between <code>conda install</code> and <code>pip install</code>: first try <code>conda install</code>, and - if that doesn't work - then try <code>pip install</code> (or <code>pip3 install</code> if you're using Python3).

### Managing your Python Path

The Python functionality that is offered in PlatoSim3 can be found in the <code>python</code>.  Make sure to add this directory to your Python path, as follows:

\code{.py}
PYTHONPATH=$PYTHONPATH:<full path to the /Platosim3/python directory>
export PYTHONPATH
\endcode

In Spyder, click on the "Python Path Manager" button to add the path to this directory to your Python path.


## Running PlatoSim3 in Python

### Creating a <code>Simulation</code> Object

To run PlatoSim from within Python, you have to import the <code>simulation</code> module:

\code{.py}
import simulation
\endcode

and create a new <code>Simulation object</code> (the name <code>sim</code> is arbitrary):

\code{.py}
sim = Simulation(<runName>, configurationFile = <full path to the configuration file>, outputDir = <full path to the output directory>)
\endcode

Note that you can use the working directory (i.e. the directory specified in the <code>PLATO_WORKDIR</code>) environment variable to specify the path to the input configuration file or the output directory.  To read the content of this environment variable in Python, use

\code{.py}
import os
workDir = os.environ["PLATO_WORKDIR"]
\endcode


#### Input Configuration File

The <code>configurationFile</code> parameter is optional and contains the full path to the configuration file that will be used as a reference.  A copy of this file will be stored in the output directory (with the name <code>runName.yaml</code>) and only the copy (and not the original file) will be modified.  We will show later how to do this.

If not specified, the <code>inputfile.yaml</code> configuration file from the <code>/inputfiles</code> directory will be used (provided you have exported the <code>PLATO_PROJECT_HOME</code> environment variable).  It is possible to also set the path to the reference input file as follows:

\code{.py}
sim.readConfigurationFile(<full path to the configuration file>)
\endcode

In that case, all changes to the configuration parameters you may have applied, will be overwritten.

#### Output Directory

Also the <code>outputDir</code> parameter is optional.  This contains the full path to the output directory.  This must be specified before you can start the simulation.  The output directory can be set can also be set as follows:

\code{.py}
sim.outputDir = <full path to the output directory>
\endcode

In this directory, the following information will be stored when running the simulation:

- <code>runName.yaml</code>: copy of the input file with the chosen configuration parameters (copied from the reference configuration file and possibly modified, as shown below);
- <code>runName.hdf5</code>: resulting exposures (images, PSF, etc.);
- <code>runName.log</code>: log file (to report any problems).

 
### Modifying Configuration Parameters

Individual configuration parameters can be modified in the copy of the configuration file (which will be stored with the name <code>runName.yaml</code> in the output directory) as follows:

\code{.py}
sim["<block>/<parameter>"] = <new value>
\endcode

for individual parameters in the block of the YAML file (e.g. <code>sim["Camera/FocalPlaneOrientation"] = 0</code>), or

\code{.py}
sim["<block>/<sub-block>/<parameter>"] = <new value>
\endcode

for parameters in sub-block in the YAML file (e.g. <code>sim["PSF/Gaussian/Sigma"] = 0.50</code>).  Note that the original YAML file is not modified.  Only the copy of this fill gets modified.  In case you change the configuration file of the <code>Simulation</code> object, the copy <code>runName.yaml</code> file will be overwritten and all changes to the configuration parameters you have applied, will be lost.

A detailed description of the configuration parameters can be found @ref InputFileDescription "here".



### Running the Simulation

The simulation can be executed with the following command:

\code{.py}
sim.run()
\endcode

By setting the optional input parameter <code>removeOutputFile</code> to <code>True</code>, any files with the same names in the specified output directory will be deleted first:

\code{.py}
sim.run(removeOutputFile = True)
\endcode
