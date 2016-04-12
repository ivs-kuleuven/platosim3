@mainpage Documentation for the PLATO Simulator

@section intro_sec Introduction

The Plato Simulator is an end-to-end software tool, designed to perform realistic simulations of the expected observations of the Plato mission. It can, however, easily be adapted to similar types of missions.

Our simulator models and simulates time series of CCD images by including models of the CCD and its electronics, the telescope optics, the stellar field, the jitter movements of the spacecraft, and all important natural sources of noise.

Many aspects concerning the design trade-off of a space-based instrument and its performance can best be tackled through realistic simulations of the expected
observations. The complex interplay of various noise sources in the course of the observations make such simulations an indispensable part of the assessment and 
design study of any space-based mission.

@section install_sec Installation

@subsection step1 Step 1: Download PlatoSim3 from GitHub

After you have downloaded the Plato Simulator code (either by using git clone or by downloading the zip-ball or tarball file), you must first take care of the
dependencies before you can actually build and run the Plato Simulator. This section describes the requirements and the procedures to install the dependencies 
and to build the Plato Simulator code.

Have a look here for a detailed description of the different components of the Plato Simulator package.

@section input_sec Inputfiles

The complete description of the input files can be found in [a seperate document](@ref InputFileDescription)