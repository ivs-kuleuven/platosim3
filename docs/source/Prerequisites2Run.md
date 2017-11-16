# Prerequisites for Running PlatoSim3 {#ReqsRun}

## Data Package

If you want to use realistic PSF models instead of a Gaussian, you can download these from <a href="ftp://plato:miSotalP@ftp.ster.kuleuven.be/psf.hdf5
">our FTP server</a>.  A convenient place to store this file is in the <code>/inputfiles</code> directory.

---

## Environment Variables

After the installation of the software, the PLATO Simulator can be run from the <code>/build</code> directory (you will have to create).  To avoid having to hardcode any path in configuration files, tutorials, etc., you must export three environment variables:

- <code>PLATO_PROJECT_HOME</code>: to refer to the directory in which PlatoSim3 was installed,
- <code>PLATO_WORKDIR</code>: to refer to the directory you want to write output or in which you want to store your own configuration files (preferably not within <code>/PlatoSim3</code>),
- and <code>PYTHONPATH</code>: to refer to the directory in which our Python scripts can be found.

This can be done as follows:

\code
PLATO_PROJECT_HOME=<full path to /PlatoSim3>
export PLATO_PROJECT_HOME

PLATO_WORKDIR=<full path to a preferred working directory>
export PLATO_WORKDIR

PYTHONPATH=$PYTHONPATH:<full path to /PlatoSim3/python>
export PYTHONPATH
\endcode

If you want, you can copy this code to make your own little script to set up your environment (e.g. <code>setPlatoEnvironment</code>), or add it to your <code>.bash_profile</code>.

To check the content of these variables (to see whether they are set to the proper location), type

\code
echo $PLATO_PROJECT_HOME

echo $PLATO_WORKDIR

echo $PYTHONPATH
\endcode

---

## <a name="yourOwnFiles">Where to Store your own Files?

To avoid problems when updating the PlatoSim3 software, it is best to store your own input and output files in a designated working directory, (preferably) outside the installation directory of PlatoSim3.  You can (but should not) add your input files to the <code>/inputfiles</code> directory, but under no circumstances change the original files in that directory!

The path to the designated working directory must be exported as the <code>PLATO_WORKDIR</code> environment variable, as described above.