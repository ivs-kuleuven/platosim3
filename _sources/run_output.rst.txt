Output Files
============

The output files of PlatoSim is so-called `HDF5 format <https://www.hdfgroup.org/solutions/hdf5/>`_, and the following figure shows the output structure: 

.. figure:: ../figures/platosim_structureOfHDF5.png
   :align: center
   :width: 1000

The blue boxes show the top-level directories, which can contain information about the input parameters (red box and orange boxes), pixel maps (green boxes), or information about the detected stars / information about the payload (purple boxes):

* The first level shows the group ``InputParameters`` contains a copy of the configuration parameters from the YAML file. Every HDF5 output file has the Git version of the software saved should you need it.
* The second level shows all the groups that contain image data for every exposure. These are the ``Images``, ``SmearingMaps``, ``BiasMaps``, and ``ThroughputMaps``.
* The third level shows all the groups that contain single frame image data. These are the ``PSF``, ``BackgroundMap``, ``FlatField``, and the ``CTI``. Notice that the latter contains a map for each CTI trap specy that is included in the used model (hence 4 species for the Short model).
* The fourth level shows all the groups that contains a mix of stellar parameters, payload parameters, and coordinates. The group ``StarCatalog`` contains the sky coordinates, the pixel coordinates, and the focal plane coordinates of all the stars that were detected during any exposure. The ``starIDs`` map the ID from the input starCatalog that is supplied with the configuration. The two groups ``ASC`` and ``Telescope`` contains time series of the displacement of the platform and camera, respectively. The groups ``StarPositions``, ``PointLikeGhostPositions``, and ``ExtendedGhostPositions`` contain information about the different source images produced, and ``Cosmics`` contains information about the cosmic rays detected within the subfield. The output structure of the four latter groups depend on the :ref:`GroupByExposure <run_input_parameters_hdf5>` setting.  

We recommend you to use the command-line tool `argos <https://pypi.org/project/argos/>`_ to inspect the output files for quick-look purposes. When installed (which is part of PlatoSim's default Python suite) this can simply be done with:

.. code-block:: shell

   argos <outputFileName>.hdf5

To extract information from the HDF5 file we strongly recommend you to use PlatoSim's build-in Python functionalities. We show these works in our `Python Tuturials <https://github.com/IvS-KULeuven/PlatoSim3/tree/master/docs/tutorials>`_.

   
