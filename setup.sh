#!/usr/bin/env bash

# Stop script if an error is encountered
set -e

# Set parsed arguments
PLATO_WORKDIR=$1
PLATO_PIPELINE=$2
REINSTALL_PIPELINE=$3

# Set global parameters
POETRY=$HOME/.local/bin
PLATO_SETUP=$PWD/.bash_profile
PLATO_PROJECT_HOME=$PWD
PLATONIUM=$PLATO_PROJECT_HOME/python/platosim/platonium

# If no arguments are given write usage message
if [ -z "$PLATO_WORKDIR" ] || [ "$1" = "-h" ] || [ "$1" = "--help" ]; then    
    echo ""
    echo "Usage: setup.sh"
    echo ""
    echo "ARGUMENTS:"
    echo "</path/to/plato_workdir>  : Location to PlatoSim working directory"
    echo "</path/to/plato_pipeline> : Location to PLATO LESIA SVN repo (optional)"
    echo "<reinstall_pipeline>      : Flag to force a reinstall of LESIA pipeline"
    echo ""
    echo "DESCRIPTION:"
    echo "PLATOnium is a PlatoSim toolkit that allows the user to run multi-camera"
    echo "simulations and while also using the LESIA pipeline in extention to extract"
    echo "on-ground or on-board data products in accordance to the mission strategy."
    echo "Given the PLATO_WORDIR and PLATO_PIPELINE (optional) as input arguments"
    echo "this script sets up all your path environment so they are globally defined"
    echo "on your system. If the latter argument for the pipeline is provided, the"
    echo "script further sets up and installs the LESIA pipeline."
    echo ""
    exit 1
    
else

    # SETUP PLATONIUM
    
    # Footprint of project setup (overwrite)
    if [ -f "$PLATO_SETUP" ]; then
	rm $PLATO_SETUP
    fi
    
    # Create setup file
    touch $PLATO_SETUP

    # Export paths
    echo "#!/usr/bin/env bash"                   >> $PLATO_SETUP
    echo "export POETRY=$POETRY"                 >> $PLATO_SETUP
    echo "export PLATO_PROJECT_HOME=$PWD"        >> $PLATO_SETUP 
    echo "export PLATO_WORKDIR=$PLATO_WORKDIR"   >> $PLATO_SETUP

    # Add Write a bash_profile that will be loaded
    #if grep -q "source $PLATO_SETUP" "$HOME/.bashrc"; then
    echo ""                                          >> $HOME/.bashrc
    echo "# >>> -- Export all PlatoSim paths -- <<<" >> $HOME/.bashrc
    echo "source $PLATO_SETUP"                       >> $HOME/.bashrc
    echo "# >>> ------------------------------- <<<" >> $HOME/.bashrc
    #fi

    # Add scripts to global executables
    cp -rf $PLATONIUM/picsim.py             $CONDA_PREFIX/bin/picsim
    cp -rf $PLATONIUM/varsim.py             $CONDA_PREFIX/bin/varsim
    cp -rf $PLATONIUM/payload.py            $CONDA_PREFIX/bin/payload
    cp -rf $PLATONIUM/platonium.py          $CONDA_PREFIX/bin/platonium
    cp -rf $PLATONIUM/../scripts/migtool.py $CONDA_PREFIX/bin/migtool
    cp -rf $PLATONIUM/../scripts/pdshow.py  $CONDA_PREFIX/bin/pdshow

    
    # LESIA PIPELINE
    
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

	# Install the LESIA pipeline if no files are in bin or if reinstall
	if [ ! "$(ls -A $PLATO_PIPELINE/bin)" ] || [ $# -eq 3 ]; then
	    cd $PLATO_PIPELINE/algorithms 
	    sudo make install
	fi
	
	# Copy missing files -> Seems like a bug
	PLATO_WP321=$PLATO_PIPELINE/algorithms/WP/321000
	cp $PLATO_WP321/invert/invert_parabolic1_multi $CONDA_PREFIX/bin
	cp $PLATO_WP321/microscan/discretize.py $CONDA_PREFIX/bin/discretize
	chmod 755 $CONDA_PREFIX/bin/discretize

	# Lastly, export the python path to bashrc
	echo "export PYTHONPATH=${PYTHONPATH}:$CONDA_PREFIX:$PLATO_PROJECT_HOME/python:$PLATO_WORKDIR:$PLATO_PIPELINE:$PLATO_MANPATH:$PLATO_PYLIB" >> $PLATO_SETUP
	echo "export PATH=${PATH}:$POETRY:$CONDA_PREFIX/bin:$PLATO_PROJECT_HOME/build:$PLATO_PIPELINE:$PLATO_PIPELINE/bin" >> $PLATO_SETUP
    else
	# Exclude LESIA paths if not requested
	echo "export PYTHONPATH=${PYTHONPATH}:$CONDA_PREFIX:$PLATO_PROJECT_HOME/python:$PLATO_WORKDIR" >> $PLATO_SETUP
	echo "export PATH=${PATH}:$POETRY:$CONDA_PREFIX/bin:$PLATO_PROJECT_HOME/build" >> $PLATO_SETUP
    fi


    # PROLOGUE
        
    # Finish with prolog message
    echo "----------------------------"
    echo " PLATOnium has been set up !"
    echo "----------------------------"
    echo ""
    echo " From bash checkout:"
    echo ""
    echo ">> picsim -h"
    echo ">> varsim -h"
    echo ">> payload -h"
    echo ">> platonium -h"
    echo ""
    echo " See also help functions:"
    echo ""
    echo ">> pdshow -h"    
    echo ">> migtool -h"
    echo ""
    if [ -f "$PLATO_PIPELINE/algorithms/Makefile" ]; then
	echo " Pipeline extention is added:"
	echo ""
	echo ">> pproc.py -h"
	echo ">> psffit.py -h"
	echo ">> photometry.py -h"
	echo ">> jittercorrection.py -h"
	echo ""
    fi

    # Reload .bashrc
    bash
fi
