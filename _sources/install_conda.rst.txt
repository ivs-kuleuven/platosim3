For users (via Conda)
=====================

To make life easier on the people who want to use PlatoSim without ever wanting to touch the code, we use `Jenkins <https://www.jenkins.io/>`_ to automatically build PlatoSim. Jenkins enables you to download the latest successfully built version of the :ref:`master or develop branch <install_source_brancing>`, or any specific version(s) of these. See :ref:`A word about Jenkins <install_conda_jenkins>` for more information. The flowchart below summarises the steps you have to take.

.. image:: ../figures/flowchart_installPlatoSimViaConda.png
   :align: center
   :width: 1000
   :alt: Alternative text 
	 



.. raw:: html

   <hr>
   
.. _install_conda_prerequisites:

*Prerequisites*
---------------

.. important::
   
   To be able to install PlatoSim via Conda you need to have the Python distribution `Anaconda <https://docs.continuum.io/anaconda/install/>`_ installed. Well tested Jenkins builds are with **Python version 3.8 and 3.9**.

First, create a Conda environment (e.g. called ``platosim``):

.. code-block:: shell
		
   conda create -n <environment name> python=<Python version>

Activate your new conda environment:

.. code-block:: shell

   conda activate <environment name>		
   
It is advised to use multiple Conda environments if you want to be able to switch between versions and/or branches in a smooth way (e.g. between ``platosim_master`` and ``platosim_develop``). We also recommend to create a new Conda environment, rather than trying to update your existing one, when you switch to a different version of Python. Find more information on how to use and manage `Conda environments <https://conda.io/projects/conda/en/latest/user-guide/tasks/manage-environments.html>`_.






.. raw:: html

   <hr>

.. _install_conda_installing:

*Installing* 
------------

.. warning::

   To install the PlatoSim software you need a set of credentials (**username** and **password**). If you have not received these, please contact one of the `PlatoSim developers <https://github.com/IvS-KULeuven/PlatoSim3>`_ and request them. Having the credentials at hand, use them in the installation procedure below by replacing the entries ``<username>`` and ``<password>`` with the given username and password, respectively.

Start by activating your desired Conda environment. The installation procedures below will automatically detect which operating system your are running and will install the appropriate packages for you.

Before you install PlatoSim via conda, for the first time in this environment, type: 

.. code-block:: shell
		
   conda config --add channels conda-forge
   
We recommend to install the ``master`` branch, unless you are interested in a feature that only exists on the ``develop`` branch. To install the latest version of the ``master`` branch use:

.. code-block:: shell

   conda install -c https://<username>:<password>@jenkins.miricle.org/platosim/ platosim

To install the latest version of the ``develop`` branch, use:

.. code-block:: shell

   conda install -c https://<username>:<password>@jenkins.miricle.org/platosim.devel/ platosim

To install a specific version of either the ``master`` or ``develop`` branch, simply append ``<version>=`` to the above commands. 







.. raw:: html

   <hr>

.. _install_conda_updating:

*Updating* 
----------

You may wish to update your installation when a new PlatoSim release becomes available. To update the ``master`` branch, simply use:

.. code-block:: shell

   conda update --force-reinstall -c  https://<username>:<password>@jenkins.miricle.org/platosim/ platosim
      
Accordingly, if you want to update the ``develop`` branch, type:

.. code-block:: shell

   conda update --force-reinstall -c  https://<username>:<password>@jenkins.miricle.org/platosim.devel/ platosim








.. raw:: html

   <hr>
   
.. _install_conda_python:

*Python packages*
-----------------

Conda, for which Jenkins depends on, is unfortunately not capable of installing the full suite of Python packages needed to use PlatoSim's built-in Python interface. Thus, minimally you need to install the packages listed under the section ``tool.poetry.dependencies`` in the file `pyproject.toml <https://github.com/IvS-KULeuven/PlatoSim3/blob/master/pyproject.toml>`_. As an alternative to installing each of the Python libraries one-by-one yourself using Conda/PyPi, in the following we show how to use `Poetry <https://python-poetry.org/>`_ to handle the package managing and install all packages with one command. Poetry manages your project installation easily over multiple platforms in a deterministic way, which has a clear advantage when you want to install an exact copy of your PlatoSim setup on a different machine (e.g. on a computing cluster).

**Install Poetry**

First install Poetry by following `Poetry's documentation <https://python-poetry.org/docs/>`_.

**Install Python packages**

As a user, you need to download the ``pyproject.toml`` file from either the `master branch <https://github.com/IvS-KULeuven/PlatoSim3/blob/master/pyproject.toml>`_ or the `develop branch <https://github.com/IvS-KULeuven/PlatoSim3/blob/develop/pyproject.toml>`_, depending on which branch you have installed. We recommend to store this file in your ``$CONDA_PREFIX/python`` directory which accordingly points to your active Conda environment.
   	 
Move to the directory where the ``pyproject.toml`` file is located and install the packages with:

.. code-block:: shell

   poetry install --no-root

Verify the success of the installation by typing ``poetry show`` (any failed installations will show up in red). Also verify that Poetry successfully installed all packages into your Conda environment by typing ``conda list`` (see packages installed with PyPi).

You have now installed the default Python package distribution of PlatoSim, which is sufficient for running the `Jupyter noteboook tutorials <https://github.com/IvS-KULeuven/PlatoSim3/tree/master/docs/tutorials>`_ and to use most of PlatoSim's Python modules for scripting. Always remember to activate your Conda environment before using PlatoSim's Python interface.


.. warning::

   Some users may find that Poetry stalls during the installation procedure. If that happens it is typically caused by a bad *keyring* setting in Poetry. Simple run the following command and retry the installation:

   .. code-block:: shell
		   
      export PYTHON_KEYRING_BACKEND=keyring.backends.fail.Keyring



   
.. raw:: html

   <hr>

.. _install_conda_jenkins:

*A word about Jenkins*
----------------------

We are using `Jenkins <https://www.jenkins.io/>`_ to automatically build PlatoSim and make pre-built software available for a myriad of operating systems. The figure below summarises how we want to use it.

Each time code is pushed to the repository (either to the ``master`` or the ``develop`` branch) or a pull request is merged in GitHub, Jenkins will start building the new code, resolve the dependencies for you, and make it available via the Conda installation command. In case the build was successful, users can install this build on their system.

To monitor the status of the PlatoSim builds, check `here <https://jenkins.miricle.org/view/Platosim/>`_.

.. image:: ../figures/github_jenkins.png
   :align: center
   :width: 1000
