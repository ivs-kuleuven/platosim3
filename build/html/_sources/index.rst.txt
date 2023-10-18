.. PLATOnium documentation master file, created by
   sphinx-quickstart on Mon Jan 10 22:56:41 2022.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Welcome to PlatoSim's documentation!
====================================

The PLATO Simulator (PlatoSim) is an end-to-end software tool, designed to perform realistic simulations of the expected observations of the PLATO mission. Our simulator models and simulates time series of CCD images by including models of the CCD and its electronics, the telescope optics, the stellar field, the jitter movements of the spacecraft, and all important natural sources of noise.

.. image:: ../figures/logo_platosim.png
   :align: center
   :width: 600

Many aspects concerning the design trade-off of a space-based instrument and its performance can best be tackled through realistic simulations of the expected observations. The complex interplay of various noise sources in the course of the observations make such simulations an indispensable part of the assessment and design study of any space-based mission.


.. raw:: html

   <hr>


Table of Content
================

.. toctree::
   :maxdepth: 1
   :caption: General:

   basic_intro
   basic_troubleshooting
   basic_acknowledgements
   basic_architecture
   
.. toctree::
   :maxdepth: 1
   :caption: Installing & Updating:

   install_overview
   install_conda
   install_source
   install_migtool
   
.. toctree::
   :maxdepth: 1
   :caption: Running PlatoSim

   run_overview
   run_prerequisites
   run_input
   run_simulations
   run_output
   
.. toctree::
   :maxdepth: 1
   :caption: PLATOnium toolkit

   platonium_overview
   platonium_tutorials
   platonium_pipeline
   ..
      platonium_cluster

   
	     
