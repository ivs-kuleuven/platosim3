Tutorials
=========

We here provide a quicklook of the important scripts and their use cases. PLATOnium is PlatoSim's user friendly command-line interface to quickly setup a large scale multi-camera and multi-quarter simulation using the PLATO Input Catalogue (PIC) and realistic stellar variability input. Thus, almost all simulations will involve the three following steps provided by the software:

#. ``picsim``    : Tool to generate a PIC star catalogue
#. ``varsim``    : Tool to generate noise-less light curves
#. ``payload``   : Tool to generate instrumental systematics
#. ``platonium`` : Tool to run a PlatoSim multi-camera and multi-quarter simulation

Both ``picsim``, ``varsim``, ``payload`` creates immediate userable input files to run ``platonium`` in an easy and generic way.

Notice all script share the same architecture for default input- and output argument parsing:

#. ``-i`` : Input file/folder 
#. ``-o`` : Output file/folder
#. ``-v`` : Verbosity (a.k.a log level)
#. ``-p`` : Flag for plotting

Both ``varsim`` and ``platonium`` have adjustable verbosity levels (i.e. log level), which becomes important when running simulations on a cluster. The verbosity is also identical to PlatoSim usage:

#. ``-v 0`` : Cluster mode: Disabling print and warnings, and no log files are saved
#. ``-v 1`` : Default mode: Print details to bash but do not save log files (Default)
#. ``-v 3`` : Debug mode  : Print details to bash and saves all log files


   

   
.. raw:: html

   <hr>
	   
.. _picsim:

picsim
------

**Create a customized PIC**

To speed up the process of running realistic simulations a script called ``picsim`` (a.k.a. PIC of Destiny) is made available for the user. This script draw target and contaminant stars from the `asPIC 1.1.0 <https://archive.stsci.edu/hlsp/aspic>`_  created by `Montalto et al. (2021) <https://arxiv.org/abs/2108.13712>`_. 

Since ``picsim`` relies on the old version of asPIC currently stars can only be drawn from the old North PLATO Field (NPF; :math:`l=65^{\circ}`, :math:`b=+30^{\circ}`) and South PLATO Field (SPF; :math:`l=253^{\circ}`, :math:`b=-30^{\circ}`). A total number of 302,743 PLATO targets compose the PIC samples (P1, P2, P4, and P5), of which 8,587,898 photometric contaminants are catalogued within a 60 arcsec radial distance from each target. The figure below shows all P1 sample stars from both the NPF and SPF colour coded by their magnitude.

.. image:: ../figures/platonium_P1SampleAllsky.png
   :align: center
   :width: 800

**Usage function:** To get an overview of its usage simply type:

.. code-block::

   picsim -h
      
**Prelook:** As a simple example, say we want to create a catalogue of 100 P1 stars from the SPF. As a sanity check you can parse the flag ``-t`` which produce a few plots of your *target* star selection. If you also want to see the number statistic for the contaminants simply parse the plotting flag `-p` instead. Hence try to launch the following commands yourself:
   
.. code-block::

   picsim 100 P1 SPF -t
   picsim 100 P1 SPF -p

**General examples:** To elaborate on the usage, we here show a few useful examples:
  
.. code-block::

   picsim 100 P1 NPF --camera 24 --project <project_name>
   picsim 100 P5 SPF --spec G --lum V --mag 9.5-12.2 -o </path/to/outdir>
   picsim all P1 NPF --group 1 --mag_limit 18 --dis_limit 30 -o </path/to/outdir>
   
In the first example we draw 100 P1 stars from the NPF but only visibible by all N-CAMs in camera-group 1 and stars being observable with all 24 N-CAMs. In the second example we draw 100 P5 stars from the SPF but limit our catalogue to only contain G dwarf main sequence stars within the PLATO passband magnitude range of 9.5-11.2. In the last example we select all P1 stars from the NPF but limit here number of contaminants stars by only saving stars brighter than 18 magnitude within a maximum radial distance of 30 arcsec (i.e. within 2 pixels) from their target star. Notice that all examples shown we parse the argument ``-o`` which is needed in order to save the stellar catalogue to the . To lower the I/O writing between software packges the stellar catalogues are saved to the binary format called *feather* (with file extension ``.ftr``).

**Combine or replot catalogue(s):** Notice that ``picsim`` also can be used to combine catalogues or replot an old catalogue produced by itself. Both can simply be done by giving the already exsisting feather binary catalogue files as input to ``picsim`` again using the argument ``-i`` as follows:

.. code-block::

   picsim all all SPF -i </path/to/indir>/starcat**.ftr -o </path/to/outdir> -p

Note we here use the asteriks `**` to specify that all files should start with ``starcat`` (which is the default output prefix of ``picsim``) and all files needs to be in feather format ``.ftr``. To make sure that all targets are selected we have here slected ``all`` for the two first mandatory arguments. Currently it is not possible to combine catalogues from different PLATO pointing fields, hence, we imagine that the previous catalogues all drawed stars from the SPF.





.. raw:: html

   <hr>

.. _varsim:

varsim
------

**Generate variable source files**

In order to help the process of generating variable input light curves for PlatoSim, as part of PLATonium the script ``varsim`` is made available. Given a star and a exoplent this script creates a synthetic stellar and exoplanet variability model directly inline with the cadence of your simulated observations. The script is developed around modelling stellar activity modulations, granulation, and stellar convection driven oscillations (p-modes). Exoplanet transits and phase curve variations (due to Doppler beaming and ellisoidal distorsion) can likewise be included. There is a future plan on including a eclipsing binary model as well. The stellar p-modes and transit eclipses are derived from synthetic `PHOENIX <https://phoenix.astro.physik.uni-goettingen.de/>`_ spectra being convolved with the instrumental bandpass. The figure below shows an example of a generated noise-less light curve of a variable star with a Neptune-sized transiting planet. We refer the reader to the technical note `PLATO-PL-KUL-0020 <https://issues.cosmos.esa.int/platowiki/display/PPWS/Simulated+datasets?preview=%2F36548784%2F56427070%2FPLATO-KUL-PL-TN-0020_MultiQuarterCameraSimulation.pdf>`_ for a more detailed description.

.. image:: ../figures/platonium_varsimExample.png
   :align: center
   :width: 800


**Usage function:** To get an overview of its usage simply type:

.. code-block::

   varsim -h


**General usage:** Seen from the usage function there are several ways to create a noise-less light curve using ``varsim``. The following two examples show the most standard way to run ``varsim``:

.. code-block::

   varsim --star_params 5777 4.5 0.0 --planet_params 100 50 0 90 0 1 1 --time 365 -p
   varsim --star_params 5777 4.5 0.0 --planet_params 100 50 0 90 0 1 1 --quarter 1-8 -o </path/to/varsource.txt>

In both examples we parse the stellar- and planet parameters, respectively. Note that the ``--quarter`` argument represent mission quarters (i.e. one 1 quarter is 30 days) and, if invoked, this argument overwrites the ``--time`` argument (likewise given in units of days). The quarter arguments is especially handy as it allows you to e.g. produce noise-less light curves with a time column suited for simulations starting beyond mission BOL (e.g. use ``--quarter 5-8`` to get a ligth curve starting after one year of mission BOL).

Intuitively the the user always needs to parse a star argument: either ``--star`` or ``--star_params``. By default the granulation and pulsation signals are activated, but all other stellar vaiability signals (i.e. spots for now) needs to be activated by parsing the corresponding flag. You can exclude the granulation and pulsations independetly by parsing ``none`` as argument for their scaling relations:

.. code-block::

   varsim --star_params 6000 4.5 0.0 --gran none --puls none --spot --time 100 -p

This example shows how to include only stellar spot modulation in the final light curve.

   
**Case studies:** To make it easier for the user to quickly fetch a favorit star-planet system, it is possible to respectively save your favorit star and your favorit planet to the file ``source_star.py`` and ``source_plaent.py`` placed in the folder ``$PLATONIUM/platonium/var``. Simply copy one of the existing code block starting with ``source == "<name>"`` and add the star/planet parameters. Note you need to obey the unit convention by `Astropy <https://www.astropy.org/>`_ (i.e. multiplying with ``u.<unit>``). A few systems exists by default, and you can likewise cross-match different systems such as: 

.. code-block::

   varsim --star Sun --planet CoRoT-1b --time 30 -p


**Photometric standards:** Regarding photometric standard stars the following examples are two potential solutions:

.. code-block::

   varsim --quarter 2-4 -p
   varsim --star std --quarter 2-4 -p
   
By parsing nothing else than the time/quarter duration, the first example is simply a constant star. The second example calls a simple script that models a short-period rotationally modulated
A type (Ap) star as simple sinusoidal vairations with occasional contribution of the second harmonic. The rotational period is randomly drawn from a uniform distribution between 1-3 days and we assume that a typically distribution of 10-30 mmag in the Kepler passband (due to current limited knowledge of the typical amplitudes of these stars in the PLATO passband).




.. raw:: html

   <hr>

.. _payload:

payload
-------

**Generate instrumental systematics**

To help introduce more realistic (i.e. more noisy) simulations, we have developed the script called ``payload``. This script will generate the following files:

  - ``instrumentPRE.txt``  : Pointing errors for each mission quarter pointing
  - ``instrumentAPE.txt``  : Alignment errors between cameras and the optical bench
  - ``instrumentTED.txt``  : Long-term camera drift due to Thermo-Elastic Distortion
  - ``instrumentGaps.txt`` : Realistic distribution of spacecraft down times
  - ``run.pbs`` & ``data.pbs`` : Job and parameterisation file for array-like parallelisation
    
Beyond the input files that are possible to parse to PlatoSim (see `here <http://ivs-kuleuven.github.io/PlatoSim3/_input_file_description.html>`_) in PLATOnium the first three files generated by this script are so-called Pointing Error Sources (PES) files. I.e. systematic noise sources that directly impact the source PSF either in position, shape, or both. The output format of these two files are directly in an input PlatoSim can interpret:

  - Pointing Repeatability Error (PRE) format: [Quarter (int), Yaw (:math:`"`), Pitch (:math:`"`) Roll (:math:`"`)]
  - Absolute Pointing Error (APE) format: [Tilt (:math:`^{\circ}`), Azimuth (:math:`^{\circ}`)]

The PRE and APE files are generated by drawing each error component from a Gaussian distribution with an translational and rotational error tolerence of :math:`3\sigma` as illustrated in a pixel displacement below:

.. image:: ../figures/platonium_payloadPES.png
   :align: center
   :width: 800
   	   
Next the script generated a polynomial model of the long-term Thermo-Elastic Drift (TED) of each camera for each mission quarter. Such a model looks like the following:

.. image:: ../figures/platonium_payloadTED.png
   :align: center
   :width: 800

Since data gaps is very important for any type time series analysis, we here provide a script that random generates several gaps from: quarterly rolls, data dowlinks, loss of fine guidance, and safe mode events. We draw these from a known distribution of the *Kepler* mission and one example of gapped time series is seen below:
	   
.. image:: ../figures/platonium_payloadGaps.png
   :align: center
   :width: 800
	   
Note also the lenght of each data gap event is randomly drawn to match closer to the future observations.

**Usage function:** To get an overview of its usage simply type:

.. code-block::

   payload -h
      
**General usage:** As a simple example type:
   
.. code-block::

   payload SPF 100 --project <project_name> -p

While using the ``--project`` output path, the files will be generated directly into your working directory and will immediately be known and used when running ``platonium`` later.
   

   


.. raw:: html

   <hr>

.. _platonium:
   
platonium
---------

**Run multi-camera PlatoSim simulations**

This script uses the PIC targets and their contaminants (created with ``picsim``) to simulate realistic imagettes or light curves using PlatoSim. It can simulate any of the 24 normal cameras/telescopes, for an arbitary number of mission quarters. Since a nescessary rotation (along the roll axis) of the spacecraft platform is required in order to repoint the solar panels every 90 days, simulations cannot realistically extend beyond a quarter of a year. Given an PIC input sample, the script will automatically simulation all targets being visible by any camera falling on one of the 4 CCDs. The fact that the PLATO spacecraft will be equipped with 4 camera groups consisting of 6 cameras each, constrains the efficient use of node-cores a computing cluster. In order to make the simulation as realistic as possible, random seeds of various intrinsic- and instrumental effects are included, meaning that each camera within each camera group needs to be simulated independently (which indeed are the raw imaging output of PLATO).

.. image:: ../figures/platonium_starInCCDFocalPlane.png
   :align: center
   :width: 800

**General setup:** Before creating your first project and start running simulations we here present some general information that will help you get started:

#. PLATOnium naturally shares the same path environment ``PLATO_WORKDIR`` of PlatoSim. We strongly recommend to set and use this path environment as the working directory for your projects, as it allows faster I/O formatting.

#. To secure that ``platonium`` can parse all the necessary input files to PlatoSim, within any project directory (``</path/to/plato_workdir/project_name>``) it is mandatory to have a folder called ``input`` where you place all inputfiles. In this folder you need to place your ``inputfile.yaml``, ``starcat**.ftr``, etc. Hence a directory tree example will look like:

   .. code-block::

      </path/to/plato_workdir>
          ├── </project_name>
              ├── /input
	          ├── inputfile.yaml
	          ├── starcat**targets.ftr
	          ├── starcat**contaminants.ftr
	          └── varSourceFile.txt

   When set properly, one can simply parse the argument ``--project <project_name>`` to ``platonium`` and all files within ``</path/to/plato_workdir/project_name>`` will be read automatically.  

#. The general rule while using ``platonium`` is that any argument parsed from the terminal will overwrite the configured parameter within your ``inputfile.yaml`` if it exists. E.g. the optional input argument ``--cadence <cycle_time>`` will overwrite the parameter called ``CycleTime`` given in the input YAML file.
 
**Usage function:** To get an overview of its usage simply type:

.. code-block::

   platonium -h

Seen from the usage function, ``platonium`` takes 4 mandatory input parameters being the star ID in your target catalogue (``starcat**targets.ftr``), the camera-group ID, the camera ID, and lastly the quarter number (where 1 is the frist quarter from mission BOL).

**General usage:** It is possible to parse your input project path in the two following ways:
   
.. code-block::

   platonium 1 1 1 1 --project <project_name> -p
   platonium 1 1 1 1 -i </path/to/plato_wordir/project_name> -p

Some frequent changed observational parameters:
   
.. code-block::

   platonium 2 4 6 24 --project <project_name> --cadence 50
   platonium 2 4 6 24 --project <project_name> --seed 1234567
   platonium 2 4 6 24 --project <project_name> --mask 10

In the first example we change default cadence from 25s to 50s. Next example the ``--seed`` argument should be parsed when you want to reproduce your results. We here use the updated version of `Numpy's Random Generator <https://numpy.org/doc/stable/reference/random/generator.html>`_ to bootstrap or configure all random seeds used in PlatoSim. Lastly, if you have turned on the PlatoSim's built-in photometry extraction module within your ``inputfile.yaml`` then you can change the mask-update using the argument ``--mask`` expressed in days.


**Variable sources:** Easy to include variable sources. ``varsim`` provide the correct format for inclusion here:
   
.. code-block::

   platonium 3 2 6 24 --project <project_name> --varfile </path/to/varSourceFile.txt>
   platonium 3 2 6 24 --project <project_name> --varlist </path/to/varSourceList.txt>

Notice the difference between ``varSourceFile`` and ``varSourceList``. The first is the ascii file containing the noise-less light curve, and the lastter is a ascii file containing the star indices as a first column and the full path and filename of each ``varSourceFile`` you want to include. Thus, first example only adds a variable signal for your target stars, whereas the second example can be used to include variable signals for the stellar contaminants as well (particularly important to get a realistic distribution of *false-positive* detections of exoplanet transits).
