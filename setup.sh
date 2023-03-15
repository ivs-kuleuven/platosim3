#!/usr/bin/env bash

# Stop this script if we encounter an error with one of the packages
set -e

# Set parsed arguments
PLATO_WORKDIR=$1
PLATO_PIPELINE=$2

# Set global parameters
POETRY=$HOME/.local/bin/poetry
PLATO_SETUP=$PWD/.bash_profile
PLATO_PROJECT_HOME=$PWD
PLATO_MANPATH=$PWD/pipeline/man
PLATO_PYLIB=$PLATO_PIPELINE/lib/python

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
    echo "Script to setup the PLATOnium software on your local machine or cluster."
    echo "Note that clusters typically don't support much disk space for HOME, and"
    echo "hence it is better to install Poetry on your DATA disk storage."
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
	echo "export PLATO_PIPELINE=$PLATO_PIPELINE" >> $PLATO_SETUP

	# Create folder structure
	mkdir -p $PLATO_PIPELINE/bin
	mkdir -p $PLATO_PIPELINE/lib
	mkdir -p $PLATO_PIPELINE/man
	mkdir -p $PLATO_PIPELINE/man/man1
	mkdir -p $PLATO_PIPELINE/include
	
	# Create env.bash file and write to it
	echo "unset MANPATH"                                  >> $PLATO_SETUP
	echo "export PLATO_MANPATH=$(manpath):$PLATO_MANPATH" >> $PLATO_SETUP
	echo "export PLATO_PYLIB=$PLATO_PYLIB"                >> $PLATO_SETUP
    fi

    # Lastly, export the python path to bashrc
    echo "export PYTHONPATH=${PYTHONPATH}:$POETRY:$PLATO_PROJECT_HOME/python:$PLATO_WORKDIR:$PLATO_PIPELINE:$PLATO_MANPATH:$PLATO_PYLIB" >> $PLATO_SETUP
    echo "export PATH=${PATH}:$HOME/.local/bin/:$PLATO_PROJECT_HOME/build:$PLATO_PIPELINE/bin" >> $PLATO_SETUP
       
    # Add code to global executeables (-i overwrite old files)
    cp -rf $PLATO_PROJECT_HOME/python/platosim/platonium/platonium $HOME/.local/bin/
    cp -rf $PLATO_PROJECT_HOME/python/platosim/picsim/picsim       $HOME/.local/bin/
    cp -rf $PLATO_PROJECT_HOME/python/platosim/varsim/varsim       $HOME/.local/bin/

    # Reload .bashrc
    . $HOME/.bashrc

    # Lastly, try to install L1 pipeline
    if [ -f "$PLATO_PIPELINE/algorithms/Makefile" ]; then
	echo "----------------------------"
	echo " Installing the L1 pipeline "
	echo "----------------------------"
	cd $PLATO_PIPELINE/algorithms 
	#make install	
    fi
    
    # Fix Jupyter-notebook problem "module not found"
    # pip install ipykernel
    # python -m ipykernel install --user

    # Finish with prolog message
    echo "----------------------------"
    echo " Platosim has been set up!  "
    echo "----------------------------"
    echo "From bash checkout:"
    echo ">> platonium -h"
    echo ">> picsim -h"
    echo ">> varsim -h"
    if [ -f "$PLATO_PIPELINE/algorithms/Makefile" ]; then
	echo ">> pproc.py -h"
	echo ">> jittercorrection.py -h"
	echo ">> photometry.py -h"
	echo ">> psffit.py -h"
    fi
fi
