For Developers (from Source)
============================

The scheme below shows an overview of the procedure to download, update, and build the PlatoSim software on your system yourself. This scheme is valid for developers who may want to contribute to the code at some point. Also the procedure to make your own modifications to the code available for others is described.

.. image:: ../figures/flowchart_installPlatoSimFromSource.png
   :align: center
   :width: 600

.. admonition:: Overview of procedure to download, update, and build PlatoSim

   * :ref:`Prerequisites <install_conda_prerequisites>`: needed to download and update PlatoSim.
   * :ref:`Forking & cloning the repository <install_source_forking>`: for first-time installation on your local machine.
   * :ref:`Dependencies <install_source_dependencies>`: required for running PlatoSim.
   * :ref:`Building the code <install_source_build>`: related to the C++ code.
   * :ref:`Remote repositories <install_source_remote>`: to configure your GitHub setup.
   * :ref:`Update procedure <install_source_update>`: to get the latest version of the code (from upstream).
   * :ref:`Contributing to the code <install_source_contribute>`: by transfering your local changes to the upstream repo.
   * :ref:`Branching strategy <install_source_branching>`: which is used in PlatoSim.
   * :ref:`Python packages <install_source_python>`: to install if running PlatoSim with Python.

.. raw:: html

   <hr>




   
.. _install_source_prerequisites:

*Prerequisites*
---------------

**Git**

The PlatoSim code is under version control in GitHub. To be able to get the latest version of the code on your local machine and to share possible contributions with other people in the project, you need to install Git on your computer. Installation instructions can be found in the `Git reference documentation <https://git-scm.com/book/en/v2/Getting-Started-Installing-Git>`_.

.. image:: ../figures/logo_git.png
   :align: center
   :width: 200
   :alt: Alternative text 

**GitHub**

The PlatoSim code is located in a `private repository on GitHub <https://github.com/IvS-KULeuven/PlatoSim3>`_. To be able to access it, we have to grant you access explicitly. Please, send your GitHub username to the development team and you will be granted read access to the repository. You will then get an invitation by email to join.

If you have not done so already, you can make an `account on GitHub <https://github.com/>`_ for free.

.. image:: ../figures/logo_github.png
   :align: center
   :width: 300
   :alt: Alternative text 

**Installed software packages**

To be able to install the dependencies and build the code, the following software must be installed on your computer:

* `gcc <https://gcc.gnu.org/>`_ v5.1 or more recent, or `clang <https://clang.llvm.org/>`_ v3.3 or more recent
* `CMake <https://cmake.org/>`_: cross-platform open-source build system to control the software compilation process (using simple platform and compiler independent configuration files)
* `BLAS <http://www.openblas.net/>`_ and `LAPACK <https://netlib.org/lapack/>`_. Without these, the simulator will likely be slower. These libraries come pre-installed on Mac OS X (so Mac users do not have to do anything). Many Linux distributions also standardly have these libraries installed or offer a package manager to easily install them.

.. admonition:: Want to change the system default compiler?

   If you want to use a different compiler than the system default to build the code, you have to export two additional environment variables, CCX and CC, as follows:

   .. code-block:: shell
		   
      export CXX=g++-5
      export CC=gcc-5

   where you may want to adapt the right-hand side of these two lines to the compiler (version) of your choice.

.. raw:: html

   <hr>






   
.. _install_source_forking:

*Forking & Cloning*
-------------------

 The PlatoSim repository can be found on the `IvS-KULeuven GitHub pages <https://github.com/IvS-KULeuven/PlatoSim3>`_. This repository is referred to as *upstream*. This section describes the process shown in the diagram below.

.. image:: ../figures/github_forkAndClone.png
   :align: center
   :width: 900
   :alt: Alternative text 

**Fork**

To be able to not only use the code but also to contribute to it, you have to *fork* this repository. To do this, you have to go to the `upstream GitHub page <https://github.com/IvS-KULeuven/PlatoSim3>`_, shown in the figure below. Just press the ``Fork`` button at the top right (encircled in red in the screenshot above) and follow the instructions. Your personal copy of the PlatoSim repository will then show up on your personal GitHub pages. This copy is referred to as *origin*. 

.. image:: ../figures/github_forkPlatoSim.png
   :align: center
   :width: 900
   :alt: Alternative text 

**Clone**

From there you can *clone* it to a designated directory on your local machine. The following command will download the source code folder called ``PlatoSim3``:

.. code-block:: shell

   git clone https://github.com/<Your GitHub username>/PlatoSim3.git

It is possible that you will be asked for a usename and password. In that case, follow the next section *Credentials*.

After you have downloaded PlatoSim's source code, you first have to install a few packages (so-called dependencies) before you can actually build and run PlatoSim (:ref:`see next section <install_source_dependencies>`).

Note that it is also possible to clone the repository directly onto your local machine, without forking it first. Doing so you will be able to update the software but not to contribute to it, hence, we strongly discourage this approach. If you only want to use PlatoSim (without changing the code), you may want to follow the :doc:`user installation procedure <install_conda>` instead. 
	 
**GitHub credentials**

When cloning the repository to your local machine your GitHub username and `personal acces tokens <https://docs.github.com/en/authentication/keeping-your-account-and-data-secure/creating-a-personal-access-token>`_ will be needed. If you want your machine to remember these in the future, write the following command:

.. code-block:: shell
		
   git config credential.helper store

PlatoSim can then be cloned to your local machine using the command:

.. code-block:: shell

   git clone https://"<USERNAME>:<ACCES_TOKEN>"@github.com/<your-GitHub-username>/PlatoSim3.git .
   
.. raw:: html

   <hr>


   


   
.. _install_source_dependencies:

*Dependencies*
--------------

PlatoSim relies on a number software packages, which are all included in the PlatoSim distribution for your convenience. Everything concerning the dependencies can be found in the ``/dependencies`` directory. The ``/dependency/Downloads`` sub-directory contains the tarball or zipball file of all required packages. In the ``/installscripts`` sub-directory you can find Python scripts that help you to unzip or untar these files into the ``/dependency/Installs`` directory.

You can install all the required dependencies and build the code in one go by typing (from the PlatoSim3 directory): 

.. code-block:: shell

   ./install.sh

Running this script will create the ``platosim`` executable needed to run simulations from bash.
   
.. attention::

   If problems would arise when executing the ``./install.sh`` command, we recommend to install the dependencies one-by-one to pinpoint the problem. Hence, copy-and-paste and run each of the following commands (also from the PlatoSim3 directory):

   .. code-block:: shell

      python ./dependencies/installscripts/install_hdf5.py
      python ./dependencies/installscripts/install_yaml.py
      python ./dependencies/installscripts/install_armadillo.py
      python ./dependencies/installscripts/install_fftw.py
      python ./dependencies/installscripts/install_faddeeva.py
      python ./dependencies/installscripts/install_zeromq.py

   **Still experiencing problems?**

   If you are still experiencing problems with the installation above, :doc:`please raise an issue <basic_troubleshooting>`! It's convenient if you can send us the full error log, which you can get using the command:

   .. code-block:: shell
		   
      ./install.sh > output.txt 2> errors.txt

.. raw:: html

   <hr>




   
.. _install_source_build:

*Building the code*
-------------------

**Software changes**

In case of code changes (after retrieving the latest version from GitHub or after introducing changes to the code yourself), or in case you have updated the PlatoSim3 code but the dependencies remain unchanged, you only have to re-build the software but not resolve the dependencies again. This saves you a tremendous amount of time. From the directory ``/build`` type:

.. code-block:: shell

   cmake ..
   make -j 4

**Updated dependencies**

At some stage, we will want to update (some of) the dependencies. You will be notified by the developer team in case this happens. You will then have to:

* Clear the ``/dependencies/Installs`` directory
* Run the install script again (see `Denpendencies <install_source_dependencies>`_)

**Switching between branches**

As you can read :ref:`here <install_source_branching>`, we don't just use the master branch. If you switch to another branch and want to run simulations with the current branch, you will have to re-build the software.

.. raw:: html

   <hr>




   
.. _install_source_remote:

*Remote repositories*
---------------------

To *pull* or *fetch* the latest version of the upstream software into your local copy (we'll come back to this) we advise you to add *upstream* to your list of remote repositories. To check your list of remote repositories, execute the following command in the the installation folder (or one of its sub-folders):

.. code-block:: shell
		
   git remote -v

This should give output similar to this:

.. code-block::

   origin    https://github.com/<your GitHub username>/PlatoSim3.git (fetch)
   origin    https://github.com/<your GitHub username>/PlatoSim3.git (push)
   upstream  https://github.com/IvS-KULeuven/PlatoSim3.git (fetch)
   upstream  https://github;.com/IvS-KULeuven/PlatoSim3.git (push)

If there's no sign of the upstream (the last two lines), you can add it with the following command:

.. code-block::

   git remote add upstream https://github.com/IvS-KULeuven/PlatoSim3.git

In case you pointed upstream to the wrong location (i.e. you used the command from above with the wrong link), you can undo this by executing the following command in the the installation folder (or one of its sub-folders):

.. code-block::

   git remote rm upstream

.. raw:: html

   <hr>





   
.. _install_source_update:

*Update procedure*
------------------

In the course of the development process, the code in the upstream repository will be updated. The ``develop`` branch will always contain the latest version of the code, whereas the ``master`` branch is (supposed to be) stable and well-tested. The workflow is shown schematically in the figure below.

.. image:: ../figures/github_pullAndPush.png
   :align: center
   :width: 900

To get the latest version on your local machine, execute the following command for the ``develop`` and ``master`` branch, respectively:

.. code-block::
   
   git pull upstream/develop
   git pull upstream/master
   
Note that this will only work smoothly if you did not change any of the PlatoSim3 files or added files to the PlatoSim3 folders. The only exceptions are the ``/inputfiles`` and the ``/build`` folder, where you can add files.

Please note that you have to re-build the code each time you fetch software changes. How to do this is explained :ref:`here <install_source_build>`.






.. raw:: html

   <hr>

.. _install_source_contribute:

*Contribute to the code*
------------------------

The workflow to contribute to PlatoSim is shown schematically in the figure below. We assume you already went through the procedure to use the code. Do not forget to bring you local repository in line with the upstream.

.. image:: ../figures/github_pushPullRequest.png
   :align: center
   :width: 900

**Pushing your changes**

If you have made local changes to some files that are already under version control or you have added new files, and you want to transfer these changes to the ``origin`` repository, you must first *stage* these files and then *commit* them.

Staging a file can be done with the command:

.. code-block::
   
   git add <relative/path/to/file>

To inspect the status of the files in your local repository, execute the following command:

.. code-block::
   
   git status

To *commit* all staged files, execute the following command:

.. code-block::
   
   git commit -m "<commit-message>"

We advise you not to squeeze to much into one commit and to write clear commit messages. Before pushing your commit(s) you shall always fetch the upstream and merge potential changes (and perhaps resolve merge conflicts) before pushing your changes to your local GitHub account. This can be done using:

.. code-block::
   
   git fetch upstream
   git merge upstream/develop
   
You can now transfer the committed changes to the origin repository with the command:

.. code-block::
   
   git push origin develop

To further transfer these changes to the upstream repository, you must open a pull request (see below). Alternatively, you can use a *Graphical Git Clients* to perform the steps described above. An overview of the possibilities can be found `here <https://git-scm.com/downloads/guis>`_. 

**Pull requests**

Now that you have *pushed* your changes to the ``origin`` repository, you want your changes to be incorporated into the ``upstream`` repository, so other people can benefit from your efforts. To do this, go to the `upstream GitHub page <https://github.com/IvS-KULeuven/PlatoSim3>`_ and just select the ``Pull requests`` tab at the top (encircled in red in the screenshot below).

.. image:: ../figures/github_pullRequest.png
   :align: center
   :width: 900

From there you can open a new pull request by pressing the green ``New pull request`` button. Select ``compare across forks``, and then you can compare the upstream repository (on the left-hand side in the screenshot below) and the origin repository (on the right-hand side in the screenshot below). You will get an overview of the differences between the two..

.. image:: ../figures/github_compareAcrossForks.png
   :align: center
   :width: 900

To confirm you want to open a pull request, press the green ``Create pull request`` button and fill out the required information. The development team will review the changes, and accept or reject the request.

.. raw:: html

   <hr>




   
.. _install_source_branching:

*Branching strategy*
--------------------

We have adopted the branching strategy of Vincent Driessen, which means that the following two branches will be used permanently:

* The master branch, with the stable (i.e. tested) code, which is considered safe to use
* The develop branch, with the latest development, which might be (highly) experimental at times (use at your own risk).
  
**Branches**

To switch to a specific branch, use the command:

.. code-block::
   
   git checkout <branch-name>

To grab all branches and get an overview, use the following commands:

.. code-block::
   
   git fetch
   git branch

The first command pulls all remote branches (not only the one you are currently working on) and the second command gives you an overview of all branches you have on your system. 

**Release candidates and releases**

Release candidates and releases correspond to tagged versions of the master branch. To start using a specific release or release candidate, you have to check out the version of the master branch with a specific tag to a new branch, like this:

.. code-block::
   
   git checkout -b <new-branch-name> <tag-name>

We will send around the tag name for new release candidates and releases once they become available. To get an overview of the available tags, use

.. code-block::
   
   git tag -l
   
**Switching between branches**

If you switch to another branch and want to run simulations with the current branch, you will have to re-build the software, as described in :ref:`Building the Code <install_source_build>`.





.. raw:: html

   <hr>

.. _install_source_python:
   
*Python packages*
-----------------

.. important::
   
   To install PlatoSim's Python distribution we strongly recommend to install `Anaconda <https://docs.continuum.io/anaconda/install/>`_ and use Conda to manage your Python environments. Note that well tested package builds are currently with **Python version 3.8 and 3.9**.

**Create a Conda environment**

First, create a Conda environment (e.g. ``platosim``):

.. code-block:: shell
		
   conda create -n <environment name> python=<Python version>

Activate your new conda environment:

.. code-block:: shell

   conda activate <environment name>		
   
It is advisable to use multiple conda environments if you want to be able to switch between versions and/or branches in a smooth way (e.g. ``platosim_master`` and ``platosim_develop``). Please note that when you switch to a different version of Python it is advised to create a new Conda environment rather than trying to update your existing one. Find more information on how to use and manage `Conda environments <https://conda.io/projects/conda/en/latest/user-guide/tasks/manage-environments.html>`_.

**Install Poetry**

Minimally you need to install all packages listed under the section ``tool.poetry.dependencies`` in the file `pyproject.toml <https://github.com/IvS-KULeuven/PlatoSim3/blob/develop/pyproject.toml>`_. As an alternative to install each of the Python libraries one-by-one yourself using Conda/PyPi, we show how you can use `Poetry <https://python-poetry.org/>`_ to handle the package managing. Poetry manages your project installation easily over multiple platforms in a deterministic way, which has a clear advantage when you want to install an exact copy of your PlatoSim setup on a different machine (e.g. on a computing cluster).

The first step is to install Poetry by following `Poetry's documentation <https://python-poetry.org/docs/>`_.

**Install Python packages**

Move to the base of the PlatoSim repository and install the Python packages with:

.. code-block:: shell

   poetry install --with docs

This will install the minimal package distribution of PlatoSim, which is sufficient for running the `Jupyter noteboook tutorials <https://github.com/IvS-KULeuven/PlatoSim3/tree/master/docs/tutorials>`_ and to use most of PlatoSim's Python modules for scripting. Additionally it installs the Sphinx packages needed to modify this documentation (see ``sphinx`` branch).

.. warning::

   Some users may find that Poetry stalls during the installation procedure. If that happens it is typically caused by a bad *keyring* setting in Poetry. Simple run the following command and retry the installation:

   .. code-block:: shell
		   
      export PYTHON_KEYRING_BACKEND=keyring.backends.fail.Keyring

  
**Running the validation suite**

As a good practice in software building, PlatoSim has an automated test harnesses consisting of a collection of modules and test data configured to validate the software unit. Hence, as a developer, the test harnesses is a mandatory sanity check before your push code and open a pull-request.

In order to be able to run the test harnesses, you must first build the code and :doc:`export the required environment variables <run_prerequisites>`. From the PlatoSim3 directory, simply run the command:

.. code-block:: shell
		
   python tests/validations/scripts/run.py
   
