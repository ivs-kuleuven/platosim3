.. PLATOnium documentation master file, created by
   sphinx-quickstart on Mon Jan 10 22:56:41 2022.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Welcome to PlatoSim's documentation!
====================================

The PLATO Simulator (PlatoSim) is an end-to-end software tool, designed to perform realistic simulations of the expected observations of the PLATO mission. Our simulator simulates time series of CCD images by including models of the stellar field, the spacecraft systematics, the jitter movements of the spacecraft, the camera optics, the CCD and its electronics, and all important natural sources of noise.

.. image:: ../figures/logo_platosim.png
   :align: center
   :width: 600

PlatoSim is a highly versatile tool ideal to model the complex interplay of various noise sources needed for the performance assessment and design study of any space-based mission. PlatoSim is widely used by the PLATO Mission Consortium (PMC) to also conduct scientific studies related to the core and complimentary program of the mission. Being an indispensable tool in the preparation of the PLATO mission, we hope that you find PlatoSim a useful tool in your research. Please consult:

- All information needed to install and update PlatoSim is described in the sections: :doc:`Installation <install_overview>` and :doc:`Running PlatoSim <run_overview>`. Please, follow these instructions carefully before running PlatoSim.
- The PlatoSim software package includes a Python interface (with a suite of tools) which we explain in our `Python tutorials <https://github.com/IvS-KULeuven/PlatoSim3/tree/master/docs/tutorials>`_. We highly recommend to use this interface to avoid sudden software errors when updating your code (since our Python interface is kept backward compatible with previous software versions).
- Particularly relevant for simulating multi-camera imagette time series or light curves, the PlatoSim toolkit called :doc:`PLATOnium <platonium_overview>` may be of interest.
- In case of questions or problems, see our :doc:`Troubleshooting <basic_troubleshooting>` section.
- For more information about the building blocks of PlatoSim, consult the :doc:`Architecture <basic_architecture>` section and the :doc:`PlatoSim3 paper <basic_acknowledgements>`.  

 
.. raw:: html

   <hr>


Table of Content
================

.. toctree::
   :maxdepth: 1
   :caption: General:

   basic_quickstart
   basic_architecture
   basic_troubleshooting
   basic_acknowledgements
   
.. toctree::
   :maxdepth: 1
   :caption: Installing & Updating:

   install_overview
   install_conda
   install_source
   
.. toctree::
   :maxdepth: 1
   :caption: Running PlatoSim

   run_overview
   run_prerequisites
   run_simulations
   run_input
   run_output
   
.. toctree::
   :maxdepth: 1
   :caption: PLATOnium toolkit

   platonium_overview
   platonium_prerequisites
   platonium_tutorials
   platonium_pipeline_v2
   .. platonium_cluster

   
	     
