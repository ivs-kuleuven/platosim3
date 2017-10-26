# Downloading PlatoSim {#downloading}

To download the software, go to the <a href="https://github.com/IvS-KULeuven/PlatoSim3">PlatoSim3 repository on GitHub</a>. The GitHub web interface to download the PlatoSim3 software is shown in Fig. 1.

@image html /images/gitHub.png "Figure 1: Screenshot of the GitHub web interface."

If you are interested in contributing to the software, you must <code>fork</code> PlatoSim3, using the GitHub web interface.  Just press the <code>"Fork"</code> button and follow the instructions.

Although GitHub allows downloading the code as a ZIP file (by pressing the <code>"Download ZIP"</code> button in the GitHub web interface), we strongly discourage this, as this makes the process of updating to a more recent version of the software more complex and tedious.

If you want to be able to update the software (without having to re-install the dependencies each time), it is better to <code>clone</code> PlatoSim3 by executing the following command in a designated directory (you have to do this only once!):

\code git clone https://github.com/IvS-KULeuven/PlatoSim3.git .\endcode

Mind the dot at the end of the command!

After you have downloaded the PlatoSim3 code, you first have to install a few packages (so-called dependencies) before you can actually build and run the PLATO Simulator.  How to do this, is described @ref building "here".