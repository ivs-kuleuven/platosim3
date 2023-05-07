Overview
========

Welcome to **PLATOsim's Numerical Imaging testbed Utilizing Multi-disciplinary** (PLATOnium) simulations! In short, PLATOnium is a Python wrapper around PlatoSim and thus takes advantages of all the utilities and scripts that continuously are being develop for PlatoSim. This toolkit can speed up the (often) lengthly procedure of generating multi-camera simulations in a highly realistic manner. PLATOnium was initially developed in order to bridge the payload development activities with the (core and complementary) science activities of the PLATO mission.

.. admonition:: Setup script for PLATOnium
  
   **For now, PLATOnium is only available for developers.** Before using this toolkit, first make sure that you have installed the :ref:`required Python packages <install_source_python>` (with Poetry this is simply: ``poetry install --with platonium``). Next you need to go to the PlatoSim base directory and run:

   .. code-block:: shell

      ./setup <path/to/plato_workdir>

   The script will finalize your setup doing the following:

   * Export the environment variable ``$PLATO_PROJECT_HOME``
   * Export the environment variable ``$PLATO_WORKDIR``
   * Export the environment variable ``$PYTHONPATH``
   * Make all the PLATOnium scripts globally executable (see :doc:`Tutorials <platonium_tutorials>`)
   * Create a ``.bash_profile`` file within your ``$PLATO_PROJECT_HOME`` directory

   If you want to change your working directory at a later stage, simply remove the ``.bash_profile`` file and run above command again.

.. admonition:: Acronym fun fact
   
   *Platonium, named after the Ancient Greek philosopher Plato, is a hypothetical element with atomic number of 1030. Like all elements after lead (excluding stable isotopium), Platonium (Pto) has no stable isotopes, but its most stable isotope is supposedly 2576Pto.*

   
