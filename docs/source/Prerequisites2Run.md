# Prerequisites for Running PlatoSim3 {#ReqsRun}


## Data Package

If you want to use realistic PSF models instead of a Gaussian, you can download these from out FTP server. The default file for most users (in focuss PSF at 6000K) can be found <a href="ftp://plato:miSotalP@ftp.ster.kuleuven.be/PSF_Focus_0mu.hdf5"> here</a>.\n  <!-- A convenient place to store this file is in the <code>/inputfiles</code> directory. -->


### Out off focuss PSF (at 6000K):
| Distortion | download link |
|------------|---------------|
| 10 mu      | <a href="ftp://plato:miSotalP@ftp.ster.kuleuven.be/PSF_Focus_10mu.hdf5"> PSF 10 mu </a>  |
| -10 mu     | <a href="ftp://plato:miSotalP@ftp.ster.kuleuven.be/PSF_Focus_-10mu.hdf5"> PSF -10 mu</a> |
| 20 mu      | <a href="ftp://plato:miSotalP@ftp.ster.kuleuven.be/PSF_Focus_20mu.hdf5"> PSF 20 mu </a>  |
| -20 mu     | <a href="ftp://plato:miSotalP@ftp.ster.kuleuven.be/PSF_Focus_-20mu.hdf5"> PSF -20 mu</a> |
| 40 mu      | <a href="ftp://plato:miSotalP@ftp.ster.kuleuven.be/PSF_Focus_40mu.hdf5"> PSF 40 mu </a>  |
| -40 mu     | <a href="ftp://plato:miSotalP@ftp.ster.kuleuven.be/PSF_Focus_-40mu.hdf5"> PSF -40 mu</a> |
| 60 mu      | <a href="ftp://plato:miSotalP@ftp.ster.kuleuven.be/PSF_Focus_60mu.hdf5"> PSF 60 mu </a>  |
| -60 mu     | <a href="ftp://plato:miSotalP@ftp.ster.kuleuven.be/PSF_Focus_-60mu.hdf5"> PSF -60 mu</a> |
| 80 mu      | <a href="ftp://plato:miSotalP@ftp.ster.kuleuven.be/PSF_Focus_80mu.hdf5"> PSF 80 mu </a>  |
| -80 mu     | <a href="ftp://plato:miSotalP@ftp.ster.kuleuven.be/PSF_Focus_-80mu.hdf5"> PSF -80 mu</a> |
| 100 mu      | <a href="ftp://plato:miSotalP@ftp.ster.kuleuven.be/PSF_Focus_100mu.hdf5"> PSF 100 mu </a>  |
| -100 mu     | <a href="ftp://plato:miSotalP@ftp.ster.kuleuven.be/PSF_Focus_-100mu.hdf5"> PSF -100 mu</a> |



### PSF at different temperatures:
| Temperature | download link |
|-------------|---------------|
| 4000K       | <a href="ftp://plato:miSotalP@ftp.ster.kuleuven.be/FF_0mu_N4000K_4224stars.hdf5"> PSF 4000K </a>    |
| 5000K       | <a href="ftp://plato:miSotalP@ftp.ster.kuleuven.be/FF_0mu_N5000K_4224stars.hdf5"> PSF 5000K </a>    |
| 6000K       | <a href="ftp://plato:miSotalP@ftp.ster.kuleuven.be/FF_0mu_N6000K_4224stars.hdf5"> PSF 6000K </a>    |
| 6500K       | <a href="ftp://plato:miSotalP@ftp.ster.kuleuven.be/FF_0mu_N6500K_4224stars.hdf5"> PSF 6500K </a>    |



<!-- | 40 mu       | <a href="ftp://plato:miSotalP@ftp.ster.kuleuven.be/FF_0mu_B6000K_1084stars.hdf5"> PSF 40 mu </a>    | -->

---

## Source the Conda Environment (Users Only)

After the installation of the software, the PLATO Simulator can be run.  Developers will have built the code in the <code>/build</code> directory and run it from there.  For users it is sufficient to activate the appropriate @ref user-prerequisites "conda environment" (i.e. the conda environment in which they have installed the version of the software they want to use).

---

## Environment Variables {#EvironVar}

To avoid having to hardcode any path in configuration files, tutorials, etc., you must export three environment variables:

- <code>PLATO_PROJECT_HOME</code>: to refer to the directory in which PlatoSim3 was installed,
- <code>PLATO_WORKDIR</code>: to refer to the directory you want to write output or in which you want to store your own configuration files (preferably not within <code>/PlatoSim3</code>),
- and <code>PYTHONPATH</code>: to refer to the directory in which our Python scripts can be found.

This can be done as follows:

    $ PLATO_PROJECT_HOME=<full path to /PlatoSim3>
    $ export PLATO_PROJECT_HOME

    $ PLATO_WORKDIR=<full path to a preferred working directory>
    $ export PLATO_WORKDIR

    $ PYTHONPATH=$PYTHONPATH:$PLATO_PROJECT_HOME/python
    $ export PYTHONPATH

In case you've installed PlatoSim3 @ref user-install "via conda", the former environment variable should be exported as

    $ export PLATO_PROJECT_HOME=$CONDA_PREFIX

(the <code>CONDA_PREFIX</code> environment variable is automatically known when you activate the appropriate conda environment)

If you want, you can copy this code to make your own little script to set up your environment (e.g. <code>setPlatoEnvironment</code>), or add it to your <code>.bash_profile</code>.

To check the content of these variables (to see whether they are set to the proper location), type

    $ echo $PLATO_PROJECT_HOME

    $ echo $PLATO_WORKDIR
    
    $ echo $PYTHONPATH

---

## <a name="yourOwnFiles">Where to Store your own Files?

To avoid problems when updating the PlatoSim3 software, it is best to store your own input and output files in a designated working directory, (preferably) outside the installation directory of PlatoSim3.  You can (but should not) add your input files to the <code>/inputfiles</code> directory, but under no circumstances change the original files in that directory!  In case - as a user - you install a new version of the software in an existing conda environment, all changes will be overwritten and newly added files will be removed from the intallation.  So better store them outside of the installation! 

The path to the designated working directory must be exported as the <code>PLATO_WORKDIR</code> environment variable, as described above.
