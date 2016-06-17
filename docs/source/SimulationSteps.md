# Simulation Steps {#SimulationSteps}

The goal of the PlatoSim C++ code is to model a part (referred to as the "sub-field") of one CCD of 
a single telescope on the platform.  On this page we describe the steps that are executed within a 
simulation in more detail.

The PlatoSim3 has been rewritten from scratch in order to get rid of historical baggage and 
inefficient code constructions that had been built up for years. The PlatoSim v2 became impossible 
to maintain and it was very hard to add new features with confidence.

We therefore started from scratch with an object-oriented design and proper testing. We minimised 
the dependencies to external libraries and replaced inactive libraries with code that is still 
actively maintained.



## Conceptual Design

The PLATO Simulator design started from the concept pictured below (see Fig. 1). The main hardware components 
that make up the PLATO spacecraft are the platform, telescope, camera, and detectors, all of which are represented by a separate class in PlatoSim3. 
Since we are simulating noise features on an image taken from the sky, we added the <code>Sky</code> class to 
the concept. All these components are controlled by a <code>Simulation</code> object.


@image html "/images/PlatoSim3 Conceptual Design.png" "Figure 1: Conceptual design of the PLATO Simulator PlatoSim3."

The more detailed design of the PlatoSim3 is depicted in Fig. 2. The <code>Simulation</code> class in the middle controls the process and knows about all the main components in the system. Each of these main components has its own set of responsibilities:

- <code>Platform</code> provides information about the jitter,
- <code>Telescope</code> controls the thermo-elastic drift,
- <code>Camera</code> provides information on the Point Spread Function (represented by the <code>PSF</code> clas),
- and <code>Detector</code> controls all the different (sub-)pixel maps. 

The blue file icons attached to several boxes indicate that information about its responsibility is provided in the @ref InputFileDescription "YAML configuration file", while the green file icons indicate that the component is writing information into the @ref OutputFileDescription "HDF5 output file". 

@image html "/images/PlatoSim3 Detailed Design.png" "Figure 2: Detailed design of the PLATO Simulator PlatoSim3."


## Control Flow

Figure 3 shows the initialisation steps of the PlatoSim3. 

@image html "/images/PlatoSim3 Initialisation.png" "Figure 3: Initialisation steps for the PLATO Simulator PlatoSim3."

Figure 4 shows the control flow of the PLATO Simulator PlatoSim3.

The simulation is a loop over the number of exposures that is requested in the @ref InputFileDescription "YAML configuration file". The two main steps in the simulation process are the integration of the light on the detector and the readout process. 

The flux of each star in the sub-field is added to the sub-pixel map with the time interval of the pointing jitter and the duration of the exposure time. The background flux is then added and the sub-pixel map is convolved with a rebinned point spread function. After applying the flatfield, the sub-pixel map is rebinned to a pixel map, and the (geometric) vignetting is applied.

During the readout process all the noise features are applied on the pixel map. Many of the features can be switched on or off in the @ref InputFileDescription "YAML configuration file", some features need to get a value of 1.0 or 0.0 in order to be ignored.

After processing each exposure, the pixel map is written to the @ref OutputFileDescription "output HDF5 file" and the option to write also the sub-pixel maps can be enabled in the @ref InputFileDescription "YAML configuration file". Details about detected stars, their positions, and the pointing information is finally written to the @ref OutputFileDescription "output HDF5 file".

@image html "/images/PlatoSim3 Control Flow.png" "Figure 4:  Control Flow for the PLATO Simulator PlatoSim3."

