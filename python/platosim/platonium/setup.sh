#!/usr/bin/env bash

# Set parsed arguments
PLATO_WORKDIR=$1
PLATO_LESIA=$2

# If no arguments are given write usage message
if [ -z "$PLATO_WORKDIR" ] || [ "$1" = "-h" ]; then    
    echo ""
    echo "Usage: setup.sh"
    echo ""
    echo "ARGUMENTS:"
    echo "</path/to/plato_workdir>  : Location to PlatoSim working directory"
    echo "</path/to/plato_lesia>    : Location to PLATO LESIA SVN repo (optional)"
    echo ""
    echo "DESCRIPTION:"
    echo "Script to setup the PLATOnium software on your local machine or cluster."
    echo "Note that clusters typically don't support much disk space for HOME, and"
    echo "hence it is better to install Poetry on your DATA disk storage."
    echo ""
    exit 1
    
else

    # Footprint of project setup
    if [ -f "$PWD/.platonium_init" ]; then
	echo "--------------------------------"
	echo "Project has already been set up!"
	echo "--------------------------------"
	echo "To redo the setup, simply remove the file:"
	echo "$PWD/.platonium_init"
	echo "and run this setup.py script again."
	exit 1
    else
	touch $PWD/.platonium_init
    fi

    # Add and export paths to $HOME/.bashrc
    if [[ ! -z "echo $PLATONIUM" ]]; then

	echo "" >> $HOME/.bashrc
	echo "# >>> Add and export PLATOnium paths <<<" >> $HOME/.bashrc
	echo "" >> $HOME/.bashrc

	# Export PLATO_WORKDIR
	echo "PLATO_WORKDIR=${PLATO_WORKDIR}" >> $HOME/.bashrc
	echo "export PLATO_WORKDIR" >> $HOME/.bashrc
	echo "" >> $HOME/.bashrc

	# Poetry environment
	if [[ -z "echo $POETRY" ]]; then
    	    echo "POETRY=$HOME/.local/bin/poetry" >> $HOME/.bashrc
    	    echo "export POETRY" >> $HOME/.bashrc
	    echo "" >> $HOME/.bashrc
	fi

	# For PLATO-LESIA
	if [[ -z "$PLATO" ]]; then

	    # Export PLATO path
	    echo "export PLATO=</path/to/plato>" >> $HOME/.bashrc
	    echo "source $PLATO/etc/env.bash" >> $HOME/.bashrc
	    echo "" >> $HOME/.bashrc

	    # Create PLATO_LESIA/etc/env.bash file and write to it
	    echo "#!/usr/bin/env bash" >> $PLATO/etc/env.bash
	    echo "" >> $PLATO/etc/env.bash
	    echo "unset MANPATH" >> $PLATO/etc/env.bash
	    echo "export MANPATH=$PLATO/man:$(manpath)" >> $PLATO/etc/env.bash
	    echo "export PATH=$PLATO/bin:${PATH}" >> $PLATO/etc/env.bash
	    echo "export PLATOPY=$PLATO/lib/python" >> $PLATO/etc/env.bash
	    echo "if [ "$PYTHONPATH" != "" ]; then" >> $PLATO/etc/env.bash
	    echo "    export PYTHONPATH=${PLATOPY}:${PYTHONPATH}:$HOME/lib/python/" >> $PLATO/etc/env.bash
	    echo "else" >> $PLATO/etc/env.bash
	    echo "    export PYTHONPATH=${PLATOPY}:$HOME/lib/python/" >> $PLATO/etc/env.bash
	    echo "fi" >> $PLATO/etc/env.bash
	fi

	# Export Python and Path
	echo "export PYTHONPATH=$PYTHONPATH:$POETRY:$PLATONIUM" >> $HOME/.bashrc
	echo "" >> $HOME/.bashrc
	echo "export PATH=$PATH:$POETRY:$PLATONIUM" >> $HOME/.bashrc
	echo "" >> $HOME/.bashrc
	echo "# >>> ----------------------------- <<<" >> $HOME/.bashrc
	echo "" >> $HOME/.bashrc
    fi

    # Reload .bashrc
    . $HOME/.bashrc

    # Add code to global executeables (-i overwrite old files)
    cp -rf $PLATO_PROJECT_HOME/python/platonium/platonium $HOME/.local/bin/
    cp -rf $PLATO_PROJECT_HOME/python/picsim/picsim $HOME/.local/bin/
    cp -rf $PLATO_PROJECT_HOME/python/varsim/varsim $HOME/.local/bin/
    #cp -u $PLATONIUM/quicktools.py $HOME/.local/bin/

    # Finish with prolog message
    echo "----------------------------"
    echo " PLATOnium has been set up! "
    echo "----------------------------"
    if [ -z "$PLATO" ]; then
	echo "L1 pipeline has been set up!"
	echo "----------------------------"
    fi
    echo "From bash checkout:"
    echo ">> platonium -h"
    echo ">> picsim -h"
    echo ">> varsim -h"
    
fi
