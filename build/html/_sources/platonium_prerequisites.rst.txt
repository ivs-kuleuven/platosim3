Prerequisites
=============

This section gives a detailed explaination on how to finalise your setup to use PLATOnium. At this stage you should be able to run PlatoSim both from the command line and from a Python shell.

**Extra Python packages**

To use PLATOnium we only have to add a few more packages to our Conda environment using Poetry. Move to the directory where the ``pyproject.toml`` file is located (if you followed our suggestin during the installation this would be ``$CONDA_PREFIX/python``). Depending on your needs, run the Poetry install command again but now with the additional argument(s):

- ``--with platonium``
- ``--with platonium --with pipeline``

The first add-on command is sufficient for running the :doc:`PLATOnium tutorials <platonium_overview>`. Only add-on the second command if you want to run the (experimental) setup between the PLATOnium toolkit and the :doc:`LESIA prototype pipeline <platonium_pipeline>`.
  
**Script executeables**

Before you continue, verify that the following path variables are defined (preferably in a ``.bash_profile`` file which are read by the configuration file of your shell environment):

.. code-block:: shell

   echo $PLATO_PROJECT_HOME
   echo $PLATO_WORKDIR
   echo $PYTHONPATH

PLATOnium consist of four command line scripts that needs to be copied to the ``bin`` directory of your Conda environment. If you are a user, make sure that ``PLATO_PROJECT_HOME`` points to your ``CONDA_PREFIX`` path. Now define the following path variable:

.. code-block:: shell
		
   PLATONIUM=$PLATO_PROJECT_HOME/python/platosim/platonium

Copy the following scripts into your Conda environment:
   
.. code-block:: shell

    cp $PLATONIUM/picsim.py    $CONDA_PREFIX/bin/picsim
    cp $PLATONIUM/varsim.py    $CONDA_PREFIX/bin/varsim
    cp $PLATONIUM/payload.py   $CONDA_PREFIX/bin/payload
    cp $PLATONIUM/platonium.py $CONDA_PREFIX/bin/platonium

Lastly, make them executable with:

.. code-block:: shell

    chmod +x $CONDA_PREFIX/bin/picsim
    chmod +x $CONDA_PREFIX/bin/varsim
    chmod +x $CONDA_PREFIX/bin/payload
    chmod +x $CONDA_PREFIX/bin/platonium
