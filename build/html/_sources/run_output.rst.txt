Output Files
============

The output files of PlatoSim is so-called `HDF5 format <https://www.hdfgroup.org/solutions/hdf5/>`_, and the following figure shows the output structure: 

.. figure:: ../figures/platosim_structureOfHDF5.png
   :align: center
   :width: 1000

The first level shows the group ``InputParameters`` contains a copy of the configuration parameters from the YAML file. Every HDF5 output file has the Git version of the software saved should you need it.
   
The second level shows all the groups that contain image data for every exposure. These are the ``Images``, ``SmearingMaps``, ``BiasMaps``, and ``ThroughputMaps``.

The third level shows all the groups that contain single frame image data. These are the ``PSF``, ``BackgroundMap``, ``FlatField``, and the ``CTI``. Notice that the latter contains a map for each CTI trap specy that is included in the used model (hence 4 species for the Short model).

The fourth level shows all the groups that contains a mix of stellar parameters, payload parameters, and coordinates. The group ``StarCatalog`` contains the sky coordinates, the pixel coordinates, and the focal plane coordinates of all the stars that were detected during any exposure. The ``starIDs`` map the ID from the input starCatalog that is supplied with the configuration.

The exact output structure of the groups ``StarPositions``, ``PointLikeGhostPositions``, ``ExtendedGhostPositions``, and ``Cosmics`` depends on the YAML input entry called ``GoupByExposure`` when controlling what needs to be saved. By default every information is saved in folders *per exposure*, however, if this parameter is ``no`` the information is saved *per star ID*. For long baseline simulations (more than :math:`100\,000` exposures), we recommend to use the latter option since this can significantly lower the computational resources (time and memory). Note that the information about the cosmic rays are never saved per star, however, it is saved into subfolders of 1000 exposures each.

We recommend you to use the Python package `argos <https://pypi.org/project/argos/>`_ to inspect the output files for quick-look purposes. When installed (which is part of PlatoSim's default Python suite) this can simply be done with:

.. code-block:: shell

   argos <outputFileName>.hdf5

To extract information from the HDF5 file we strongly recommend you to use PlatoSim's build-in Python functionalities. We show these works in our `Jupyter notebook tuturials <https://github.com/IvS-KULeuven/PlatoSim3/tree/master/docs/tutorials>`_.

   
