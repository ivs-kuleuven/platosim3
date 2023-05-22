#!/usr/bin/env bash

# Stop script if an error is encountered
set -e

# Set parsed arguments
PLATO_WORKDIR=$1
PLATO_PIPELINE=$2

# Set global parameters
POETRY=$HOME/.local/bin/poetry   # TODO check install location
PLATO_SETUP=$PWD/.bash_profile
PLATO_PROJECT_HOME=$PWD
PLATO_MANPATH=$PWD/pipeline/man

# If no arguments are given write usage message
if [ -z "$PLATO_WORKDIR" ] || [ "$1" = "-h" ]; then    
    echo ""
    echo "Usage: setup.sh"
    echo ""
    echo "ARGUMENTS:"
    echo "</path/to/plato_workdir>  : Location to PlatoSim working directory"
    echo "</path/to/plato_pipeline> : Location to PLATO LESIA SVN repo (optional)"
    echo ""
    echo "DESCRIPTION:"
    echo "PLATOnium is a PlatoSim toolkit that allows the user to run multi-camera"
    echo "simulations and while also using the L1 pipeline in extention to extract"
    echo "on-ground or on-board data products in accordance to the mission strategy."
    echo "Given the PLATO_WORDIR and PLATO_PIPELINE (optional) as input arguments"
    echo "this script sets up all your path environment so they are globally defined"
    echo "on your system. If the latter argument for the pipeline is provided the"
    echo "script further sets up and installs the L1 pipeline."
    echo ""
    exit 1
    
else

    # Footprint of project setup
    if [ -f "$PLATO_SETUP" ]; then
	echo "------------------------------------------"
	echo ":   PlatoSim has already been set up!    :"
	echo "------------------------------------------"
	echo "To redo the setup, simply remove the file:"
	echo "$PLATO_SETUP"
	echo "and run this setup.sh script again."
	exit 1
    fi

    # Write a bash_profile that will be loaded
    echo "# >>> -- Export all PlatoSim paths -- <<<" >> $HOME/.bashrc
    echo "source $PLATO_SETUP"                       >> $HOME/.bashrc
    echo "# >>> ------------------------------- <<<" >> $HOME/.bashrc

    # Create setup file
    touch $PLATO_SETUP

    # Export paths
    echo "#!/usr/bin/env bash" >> $PLATO_SETUP
    echo "" >> $HOME/.bashrc
    echo "export POETRY=$HOME/.local/bin/poetry" >> $PLATO_SETUP
    echo "export PLATO_PROJECT_HOME=$PWD"        >> $PLATO_SETUP 
    echo "export PLATO_WORKDIR=$PLATO_WORKDIR"   >> $PLATO_SETUP
    
    # Add and export paths to $HOME/.bashrc    
    if [ $# -eq 2 ]; then

	# Define path of pipeline
	# NOTE: PLATO cannot be changed
	echo "export PLATO=$PLATO_PIPELINE" >> $PLATO_SETUP

	# Create folder structure
	mkdir -p $PLATO_PIPELINE/bin
	mkdir -p $PLATO_PIPELINE/lib
	mkdir -p $PLATO_PIPELINE/man
	mkdir -p $PLATO_PIPELINE/man/man1
	mkdir -p $PLATO_PIPELINE/include
	
	# Create env.bash file and write to it
	echo "unset MANPATH"                                  >> $PLATO_SETUP
	echo "export PLATO_MANPATH=$(manpath):$PLATO_MANPATH" >> $PLATO_SETUP
	echo "export PLATO_PYLIB=$PLATO_PIPELINE/lib/python"  >> $PLATO_SETUP
    fi

    # Lastly, export the python path to bashrc
    echo "export PYTHONPATH=${PYTHONPATH}:$CONDA_PREFIX:$PLATO_PROJECT_HOME/python:$PLATO_WORKDIR:$PLATO_PIPELINE:$PLATO_MANPATH:$PLATO_PYLIB" >> $PLATO_SETUP
    echo "export PATH=${PATH}:$CONDA_PREFIX/bin:$PLATO_PROJECT_HOME/build:$PLATO_PIPELINE/bin" >> $PLATO_SETUP
       
    # Add code to global executeables (-i overwrite old files) TODO
    cp -rf $PLATO_PROJECT_HOME/python/platosim/platonium/picsim    $CONDA_PREFIX/bin
    cp -rf $PLATO_PROJECT_HOME/python/platosim/platonium/varsim    $CONDA_PREFIX/bin
    cp -rf $PLATO_PROJECT_HOME/python/platosim/platonium/payload   $CONDA_PREFIX/bin
    cp -rf $PLATO_PROJECT_HOME/python/platosim/platonium/platonium $CONDA_PREFIX/bin
    
    # Reload .bashrc TODO mac
    . $HOME/.bashrc

    ## Lastly, try to install L1 pipeline
    # if [ -f "$PLATO_PIPELINE/algorithms/Makefile" ]; then
    # 	echo "----------------------------"
    # 	echo " Installing the L1 pipeline "
    # 	echo "----------------------------"
    # 	cd $PLATO_PIPELINE/algorithms 
    # 	make install
    ## The following is needed (bug?) to locate these two files
    # cp $PALTO/algorithms/WP/321000/invert/invert_parabolic1_multi $CONDA_PREFIX/bin
    # cp $PALTO/algorithms/WP/321000/microscan/discretize.py $CONDA_PREFIX/bin/discretize
    # chmod 755 $CONDA_PREFIX/bin/discretize
    # fi
    
    # Fix Jupyter-notebook problem "module not found"
    # pip install ipykernel
    # python -m ipykernel install --user

    # Finish with prolog message
    echo "---------------------------"
    echo " PLATOnium has been set up!"
    echo "---------------------------"
    echo "From bash checkout:"
    echo ">> picsim -h"
    echo ">> varsim -h"
    echo ">> payload -h"
    echo ">> platonium -h"
    if [ -f "$PLATO_PIPELINE/algorithms/Makefile" ]; then
	echo "---------------------------"
	echo " Pipeline has been set up !"
	echo "--------------------------"
	echo ">> pproc.py -h"
	echo ">> psffit.py -h"
	echo ">> photometry.py -h"
	echo ">> jittercorrection.py -h"
    fi
fi
