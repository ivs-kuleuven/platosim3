Overview
========

Welcome to **PLATOsim's Numerical Imaging testbed Utilizing Multi-disciplinary** (PLATOnium) simulations! In short, PLATOnium is a Python wrapper around PlatoSim and thus takes advantages of all the utilities and scripts that continuously are being develop for PlatoSim. This toolkit can speed up the (often) lengthly procedure of generating multi-camera simulations in a highly realistic manner. PLATOnium was initially developed in order to bridge the payload development activities with the (core and complementary) science activities of the PLATO mission.

.. attention:: 

   PLATOnium is only available when cloning the PlatoSim GitHub repository. We strongly recommend to use **Python 3.9**, Conda as your environment manager, and Poetry as your package manager (as we will do in the following) Also note that the setup of PLATOnium and the LESIA pipeline is only tested for Linux users, since these tools should be used on a computing cluster where Linux is the standard. Mac users may have to adapt especially the ``setup.py`` script.

   

.. raw:: html

   <hr>

.. _platonium_overview_setup:

      
Setup PLATOnium
----------------
   
1. First create a new Conda environment and activate it:

   .. code-block:: shell

      conda create -n platonium python=3.9
      conda activate platonium


2. Install all Python dependencies using Poetry from within the **PlatoSim3** folder:

   .. code-block:: shell

      poetry install --with platonium

   If all packages are installed successful, then Poetry will report that ``platosim 3.6.0`` is installed. You can also verify this by listing the installed packages either using Poetry (``poetry show``) or Conda (``conda list``).

3. Lastly, from within the **PlatoSim3** folder run:

   .. code-block:: shell

      ./setup <path/to/plato_workdir>

   The script will finalize your setup doing the following:

   * Export the environment variable ``$PLATO_PROJECT_HOME``
   * Export the environment variable ``$PLATO_WORKDIR``
   * Export the environment variable ``$PYTHONPATH``
   * Make all the PLATOnium scripts globally executable (see :doc:`Tutorials <platonium_tutorials>`)
   * Create a ``.bash_profile`` file within your **PlatoSim3** folder

   If you want to change your working directory at a later stage, simply run ``setup.sh`` again.


.. admonition:: Acronym fun fact
   
   *Platonium, named after the Ancient Greek philosopher Plato, is a hypothetical element with atomic number of 1030. Like all elements after lead (excluding stable isotopium), Platonium (Pto) has no stable isotopes, but its most stable isotope is supposedly 2576Pto.*   
