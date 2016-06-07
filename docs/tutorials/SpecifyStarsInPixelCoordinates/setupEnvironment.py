import os, sys

projectDir = '/Users/rik/Git/PlatoSim3'
workDir = '/Users/rik/Work/PLATO'

os.environ['PLATO_PROJECT_HOME'] = projectDir
os.environ['PLATO_WORK_DIR'] = workDir

sys.path.append(projectDir + "/python")
