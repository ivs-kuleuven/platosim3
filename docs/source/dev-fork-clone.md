# Forking & Cloning the Repository {#dev-fork-clone}

The PlatoSim3 repository can be found on the [IvS-KULeuven GitHub pages](https://github.com/IvS-KULeuven/PlatoSim3).  This repository is referred to as **upstream**.

This section describes the process shown in the diagram below.

@image html /images/fork-clone.png ""

---

## Fork

To be able to not only use the code but also to contribute to it, you have to "fork" this repository.  To do this, you have to go to the [upstream GitHub page](https://github.com/IvS-KULeuven/PlatoSim3), shown below.

@image html /images/gitHub.png ""

Just press the `"Fork"` button at the top right (encircled in red in the screenshot above) and follow the instructions.  Your personal copy of the PlatoSim repository will then show up on your personal GitHub pages.  This copy is referred to as **origin**.

---

## Clone

From there you can "clone" it to a designated directory on your local machine, with the following command (mind the dot at the end of the command!):

    $ git clone https://github.com/<your GitHub username>/PlatoSim3.git .

It is possible that you will be asked for a usename and password. In that case, follow the section **Credentials** underneed this section.

After you have downloaded the PlatoSim3 code, you first have to install a few packages (so-called dependencies) before you can actually build and run the PLATO Simulator.  How to do this, is described @ref dev-dependencies "here".

Note that it is also possible to clone the repository directly onto your local machine, without forking it first.  You will be able to update the software but not to contribute to it.  We therefore strongly discourage this approach.  If you only want to use PlatoSim (without changing the code), you may want to follow the @ref user-overview "user installation procedure" instead.

---

## Credentials

When cloning a repository for the first time, GitHub might as for authentification. If this is the case, you will need to obtain a [Personal Acces Token](https://docs.github.com/en/authentication/keeping-your-account-and-data-secure/creating-a-personal-access-token) for you account. When cloning the repository to your local machine your GitHub username and personal acces tokens will be needed. If you want your machine to remember these in the future, write the following command 

    $ git config credential.helper store

PlatoSim can then be "cloned" to your local machine using the command

    $ git clone https://"<USERNAME>:<ACCES_TOKEN>"@github.com/<your GitHub username>/PlatoSim3.git .
	
