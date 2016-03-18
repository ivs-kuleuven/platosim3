
from simfile import SimFile
from simulation import Simulation

# Specify the absolute paths of some of the input files and the output folder

outputDir   = "/Users/rik/Work/PLATO/Simulation/"

# Set up a Simulation object

sim = Simulation("Simul01", configurationFile="/Users/rik/Work/PLATO/myInputs/inputfile.yaml")
sim.setOutputDir(outputDir)

# Run the simulation without the flux extraction, with an HDF5 file as output

simFile = sim.run(doPhotometry=False)

# Look at an exposure image

simFile.showImage(0)

# Read an imagette from the output HDF5 file

imagette = simFile.getImagette(15209, 0, radius=2)
