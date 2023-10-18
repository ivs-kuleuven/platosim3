Pipeline
========

The following guide explains how to install, setup, and run the LESIA pipeline as an integrated part of PLATOnium.


.. raw:: html

   <hr>

.. _pipeline_run:

   
Installation
------------   

.. attention::

   Before installing the LESIA pipeline, first :ref:`Setup PLATOnium <platonium_overview_setup>`.
   
   In order to have a frozen functional setup between PlatoSim/PLATOnium and the LESIA pipeline, the lastest tested versions are:

   **PlatoSim develop branch:**

   .. code-block:: shell

      git checkout 3.6.0-225-g552813f9
   
   **LESIA SVN pipeline**

   .. code-block:: shell

      svn checkout 5420
   
   
1. The LESIA pipeline is a SVN (subversion) repository, hence, if not installed, first install SVN

   .. code-block:: shell

      sudo apt update		   
      sudo apt install subversion -y

   Confirm the installation by viewing the version: ``svn --version``.

   Next create a directory called **plato** located in the same (software) folder where **PlatoSim3** is located, and enter the directory:

   .. code-block:: shell

      mkdir </path/to/plato>
      cd </path/to/plato>

   Now download the SVN repository by providing your ``username`` and ``password`` (contact the LESIA team for the credentials): 
      
   .. code-block:: shell

      svn checkout https://version-lesia.obspm.fr/repos/PLATO/algorithms/ --username <username> --password <password>

                 
2. Install the LESIA pipeline with:

   .. code-block:: shell

      cd </path/to/plato>/algorithms
      make install


3. Finally, run the ``setup.sh`` script (again) located within the **PlatoSim3** folder. To setup the pipeline, two arguments need to be parsed; the first being the path to your PlatoSim working directory (``</path/to/plato_workdir>``) and the second being the path to your newly created folder **plato** folder (``</path/to/plato>``). Thus, simply run:
   
   .. code-block:: shell

      ./setup </path/to/plato_workdir> </path/to/plato>



   
.. admonition:: Update the LESIA pipeline:

   To update the SVN repository, simply go to ``cd $PLATO/algortihms`` and do:

   .. code-block:: shell

      svn update
      sudo make

   To checkout a specific revision use:

   .. code-block:: shell

      svn checkout <revision>



      
.. raw:: html

   <hr>

.. _pipeline_run:

      
Run the pipeline
----------------

Compared to the previous tutorial on how to run ``platonium``, we now only have to parse the ``--sample <plato_sample>`` to activate the LESIA pipeline. If ``P1`` is parsed, the on-ground pipeline is activated, and if ``P5`` is parsed, the on-board pipeline is activated. E.g. say that we already made two seperate stellar catalogues of respectively P1 and P5 sample stars using ``picsim``. A simple test of 30 exposures for each chain of the pipeline looks like the following:

.. code-block::

   platonium 1 1 1 1 --project <project_name> --sample P1 --nexp 30
   platonium 1 1 1 1 --project <project_name> --sample P5 --nexp 30

Note that a lot of information is printed to bash. When running on a computing cluster, this behavior is typically not desired and hence the flag ``-v 0`` avoid printing to bash.


.. raw:: html

   <hr>

.. _pipeline_extract:

      
Output files
------------

The output files produced by PLATOnium when the LESIA pipeline is activated is in the form of so-called **feather** binary files.
   


.. raw:: html

   <hr>

.. _platonium_pipeline_troubleshooting:


Troubleshooting
---------------

   
.. _pipeline_troubleshooting_errors:
   
*Makefile errors*
.................

A few obstacles can happen writing a code assuming full sudo rights (and assuming that apt get is available). In order to successfully install the L1 pipeline on a cluster, the following ajustments may need to be inforced:

**Disabled environment:** Make sure that the ``$PLATO`` path is defined. If added to your ``.bashrc`` or ``.profile`` file this path should automatically be loaded every time you enter a compute node on any machine/cluster. However, it occured that it was not and the installation will link to a wrong location and fail.

**aswsim3:** This library ``aswsim3`` sometimes cause problems and since we do not need to PLATOnium we can simply skip this library. However, ``aswsim3`` depends on the library ``platosimlib``, hence, simply replace ``aswsim3`` with ``platosimlib`` in the section ``install: all`` listed within in ``Makefile``. Note the problem is typically that GLS packages cannot be found here, like:

.. code-block:: shell
		
   ../src/test_gsl.c:2:24: fatal error: gsl_matrix.h: No such file or directory
   #include <gsl_matrix.h>


**OpenBLAS vs. BLAS/LAPACK:** If using OpenBLAS as a BLAS (and LAPACK) implementation the ``-lblas`` and ``-llapack`` link arguments would be appropriate when linking to the reference BLAS/LAPACK libraries from NetLib, which are typically called ``libblas.so`` and ``liblapack.so``. However, with OpenBLAS both are in the same shared library called ``libopenblas.so``.

1. The first installation error you may encounter is thus from the library ``WP/326100/spline2dbase``. This library uses the ``setup.py`` script within for the installation which link to the ``-lblas`` and ``-llapack``. Thus, if using OpenBLAS you will need to replace ``-lblas`` and ``-llapack`` with ``openblas`` in the following line:

.. code-block:: shell

   libraries=["m","blas","lapack"],

   to

   libraries=["m","openblas"],
   
For the rest of the libraries you need to look into the ``Makefile`` of each package.

2. Within the main installation Makefile (``$PLATO/algorithms/Makefile``) the inversion modules named ``WP/32100/voxel/..`` link to the openblas libray. If you get a message that ``-llapack not found`` then a potential solution is to remove this flag completely since it is already included within openblas. Hence remove ``-llapack`` from the following line in each file you encounter this error:

.. code-block:: shell

   LDLIBS = -lcvxs -lm -llapack -lopenblas -lpthread

   
**Man pages:** Since the user do not generally have sudo right on a computing cluster the inversion module tries to install the man pages for the packges ``WP/32100/voxel/..``. This can result in the following error:

.. code-block:: shell
  
   warning: failed to load external entity "/usr/share/xml/docbook/stylesheet/nwalsh/html/docbook.xsl"
   cannot parse /usr/share/xml/docbook/stylesheet/nwalsh/html/docbook.xsl

These are really not needed for our computations and hence

1. Remove the ``.1`` and ``.txt`` files from ``all :`` that will be installed
2. Comment out the line ``install -m 644 combine-mat.1 $(PLATOMAN1)/..``
3. Run ``make install`` again


*Python errors*
...............

**pymt64:** The package ``pymt64`` can sometimes cause porblems. A well known issue for this pacakage is the following error:

.. code-block:: shell

   import pymt64
   File "pymt64.pyx", line 1, in init pymt64
   ValueError: numpy.ndarray size changed, may indicate binary incompatibility. Expected 96 from C header, got 88 from PyObject

This can be solved by re-installing the library:

.. code-block:: shell

   pip uninstall pymt64
   pip install --no-cache-dir pymt64

   
*Cannot find script*
....................

Sometimes the following scripts: ``combine-mat``, ``invert-rls``, ``spline2real``, and ``find-amp`` cannot be found globally (typically true on a computing cluster by the local node you are running simulations on). By default these scripts are installed into your ``$PLATO/bin`` folder, however, the script ``setup.sh`` should take care of copying these files to your Conda environment **bin** folder. If the above script are not directly executeable from bash (e.g. ``combine-mat -h``) then the PLATOnium setup script most likely failed. Thus, simply copy all of these files to your Conda environment **bin** folder. You can find this location by first activating your Conda environment and thereafter type ``which python``. Copy the files into the **bin** folder and not the **bin/python** folder!
