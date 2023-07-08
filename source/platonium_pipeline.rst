Pipeline
========

.. important::

   The following guide explains how to install PLATOnium in combination with PlatoSim, however, if you're interested in using the L1 pipeline to extract *reduced* on-ground and on-board PLATO light curves you'll need additionally to follow the installation below. We note that PlatoSim/PLATOnium **do allow** fast and generic on-board light curve extraction for the generation of *raw* light curves.


.. attention::

   Note that if you are using the OpenBLAS library then BLAS and occationally LAPACK becomes redundant. Thus, several prerequites may be needed to succesfully install the PLATO LESIA repository containing the script of the provisional L1-pipeline. Please have a look at the :doc:`Troubleshooting <l1_troubleshooting>` page.

   
1. The LESIA repository is a SVN (subversion) repository. Hence first install svn (using ``sudo apt install <subversion>``). Create a directory called ``plato`` (preferably on your ``$DATA`` storage if installing on a cluster) and enter the directory:

   .. code-block:: shell

      mkdir </path/to/plato>
      cd </path/to/plato>

   Now download the SVN repository by providing your ``username`` and ``password``: 
      
   .. code-block:: shell

      svn checkout https://version-lesia.obspm.fr/repos/PLATO/algorithms/ --username <username> --password <password>

      
2. Prior to the installation you need to run setup script called ``/PLATOnium/setup.sh``, also explained in the general :ref:`step-by-step installation instructions <repo_install_step>`. This scirpt takes in fact two arguments; the frist being ``</path/to/plato_workdir>`` and the second being the path to your newly created folder ``</path/to/plato>``. Thus, simply run:
   
   .. code-block:: shell

      platosim --setup </path/to/plato_workdir> </path/to/plato>
   
      
3. First make sure that your PLATOnium (Poetry) Python environment is activated, and then install PLATO by:

   .. code-block:: shell

      cd $PLATO/algorithms
      make install

      
.. admonition:: Update PLATO repository

   To update the repository (whenever a new change is made) simply go to ``cd $PLATO/algortihms`` and do:

   .. code-block:: shell

      svn update
      sudo make





      

Troubleshooting
===============


.. raw:: html

   <hr>

.. _l1_trouble_python:

*Python errors*
---------------

pymt64
......

The package ``pymt64`` can sometimes cause porblems. A well known issue for this pacakage is the following error:

.. code-block:: shell

   import pymt64
   File "pymt64.pyx", line 1, in init pymt64
   ValueError: numpy.ndarray size changed, may indicate binary incompatibility. Expected 96 from C header, got 88 from PyObject


This can be solved by re-installing the library:

.. code-block:: shell

   pip uninstall pymt64
   pip install --no-cache-dir pymt64





   
.. raw:: html

   <hr>
   
.. _l1_trouble_makefile:
   
*Makefile errors*
-----------------

A few obstacles can happen writing a code assuming full sudo rights (and assuming that apt get is available). In order to successfully install the L1 pipeline on a cluster, the following ajustments may need to be inforced:

Disabled environment
....................

Make sure that the ``$PLATO`` environment is loaded. If added to your ``.bashrc`` or ``.profile`` file this path should automatically be loaded every time you enter a compute node on any machine/cluster. However, it occured that it was not and the installation will link to a wrong location and fail.

aswsim3
.......

This library ``aswsim3`` sometimes cause errors and since we do not need to PLATOnium we can simply skip this library. However, ``aswsim3`` depends on the library ``platosimlib`` hence simply replace ``aswsim3`` with ``platosimlib`` in the section ``install: all`` listed within in ``Makefile``. Note the problem is typically that GLS packages cannot be found here, like:

../src/test_gsl.c:2:24: fatal error: gsl_matrix.h: No such file or directory
#include <gsl_matrix.h>

OpenBLAS vs. BLAS/LAPACK
........................

If using OpenBLAS as a BLAS (and LAPACK) implementation the ``-lblas`` and ``-llapack`` link arguments would be appropriate when linking to the reference BLAS/LAPACK libraries from NetLib, which are typically called ``libblas.so`` and ``liblapack.so``. However, with OpenBLAS both are in the same shared library called ``libopenblas.so``.

1. The first installation error you may encounter is thus from the library ``WP/326100/spline2dbase``. This library uses the ``setup.py`` script within for the installation which link to the ``-lblas`` and ``-llapack``. Thus, if using OpenBLAS you will need to replace ``-lblas`` and ``-llapack`` with ``openblas`` in the following line:

.. code-block:: shell

   libraries=["m","blas","lapack"],  # Unix-like specific

   to

   libraries=["m","openblas"],  # Unix-like specific
   
For the rest of the libraries you need to look into the ``Makefile`` of each package.

   
2. Within the main installation Makefile (``$PLATO/algorithms/Makefile``) the inversion modules named ``WP/32100/voxel/..`` link to the openblas libray. If you get a message that ``-llapack not found`` then a potential solution is to remove this flag completely since it is already included within openblas. Hence remove ``-llapack`` from the following line in each file you encounter this error:

.. code-block:: shell

   LDLIBS = -lcvxs -lm -llapack -lopenblas -lpthread

Man pages
.........
		
Since the user do not generally have sudo right on a computing cluster the inversion module tries to install the man pages for the packges ``WP/32100/voxel/..``. This can result in the following error:

.. code-block:: shell
  
   warning: failed to load external entity "/usr/share/xml/docbook/stylesheet/nwalsh/html/docbook.xsl"
   cannot parse /usr/share/xml/docbook/stylesheet/nwalsh/html/docbook.xsl

These are really not needed for our computations and hence

1. Remove the ``.1`` and ``.txt`` files from ``all :`` that will be installed
2. Comment out the line ``install -m 644 combine-mat.1 $(PLATOMAN1)/..``
3. Run ``make install`` again





.. raw:: html

   <hr>
   
*Cannot find script*
--------------------

Sometimes the following scripts: ``combine-mat``, ``invert-rls``, ``spline2real``, and ``find-amp`` cannot be found globally by the local node you are running simulations on. By default these scripts are installed into your ``$PLATO/bin`` folder, thus, simply copy all of these files to your Poetry environment **bin** folder. You can find the Poetry environment folder by first activating your Poetry environment and thereafter type ``which python``. Copy the files into the **bin** folder and not the **bin/python** folder!





.. raw:: html

   <hr>

*Remove warnings*
-----------------

If you attend to run a lot of simulations, it is quite important ro ignore all messages and warning printed to your *standard error output file* (assuming that you have double checked that your scripts run successfully of course). The reason is that you might overflow your memory storage space which can result in that all data is useless from the unset of low data storage (while the job schedular might tell you that the simulations are still running). Thus, always use the argument ``-v 0`` for ``platonium`` and add the following python code to the top of the files ``pproc.py`` and ``jittercorrection.py``:

.. code-block:: python

   import warnings
   warnings.filterwarnings("ignore")
