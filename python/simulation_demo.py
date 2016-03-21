
from simfile import SimFile
from simulation import Simulation

# Specify the absolute paths of some of the input files and the output folder

outputDir   = "/Users/rik/Work/PLATO/Simulation/"

# Set up a Simulation object

sim = Simulation("Simul01", configurationFile="/Users/rik/Work/PLATO/myInputs/inputfile.yaml")
sim.outputDir = "/Users/rik/Work/PLATO/Simulations"

# Run the simulation

simFile = sim.run()

# Look at the first exposure image

simFile.showImage(0)

# Read an imagette from the output HDF5 file

imagette = simFile.getImagette(15209, 0, radius=2)

