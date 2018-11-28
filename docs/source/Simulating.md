# Running Simulations with PlatoSim3 {#Simulating}

To run simulations with PlatoSim3, you will have to feed one @ref ConfigurationParameters "configuration file" as input to the simulator, possibly together with a few @ref SupplementaryFiles "additional files".  Launching simulations can be done either on the [command line](#commandLine) or in [Python](#pythonWrapper).

Note that - before running PlatoSim3 - you must have an environment variable <code>PLATO_PROJECT_HOME</code>, set to the base folder of PlatoSim3, as described @ref ReqsRun "here".

---





## <a name="commandLine"></a>Command Line

To initiate a simulation on the command line, <code>cd</code> to the <code>/build</code> directory and type:

\code
./platosim <input file> <non-existing output file> [<log file> <log level>]
\endcode

The log level is an integer value that allows to configure what kind of log messages is written to the log file:

* 1: only error and warning messages;
* 2: error, warning, and info messages;
* 3: all messages (error, warning, info, and debug).



---

## <a name="pythonWrapper"></a>Python

Alternatively, you can make use of our Python wrapper to configure and run PlatoSim3.  How to do this and how to set up you Python environment is described @ref ScriptingInPython "here".

Under <code>/docs/tutorials</code> you can find a number of tutorials (in the form of <a href="http://jupyter.org/">Jupyter Notebooks</a>) to get you going.  Have a look @ref Tutorials "here" to find out how to configure your Python environment to be able to get started.
