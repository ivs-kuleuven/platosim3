# Installing/Updating PlatoSim3 for Developers{#dev-overview}

The scheme below shows an overview of the procedure to download, update, and build the PlatoSim3 software on your system yourself, without using conda.  This scheme is valid for developers who may want to contribute to the code at some point (or for users who want to install the software without @ref InstallViaConda "using conda").  Also the procedure to make your own modifications to the code available for others is described.

@image html /images/download-update-build.png ""

---

## Overview

The following installation guides are available:

- [*Prerequisites*](#dev-reqs): Describes the prerequisites that must be fulfilled to be able to download, update, and build PlatoSim3.

- [*Forking & Cloning the Repository*](#dev-fork-clone): Describes the first-time installation of the software on your local machine.

- [*Installing Dependencies*](dev-dependencies.md): Describes how to install required dependencies.

- [*Building*](dev-building): Describes how to build the code.

- [*Remote Repositories*](dev-add-upstream.md): Describes how to configure your GitHub setup, so that you will be able to get hold of any changes to the code by others.

- [*Update procedure*](dev-pull.md): Describes how to get the latest version of the code (from the `upstream`) on your local machine.

- [*Contributing*](dev-push.md): Describes how to transfer your local changes to the `origin` repository and how to transfer these changes from the `origin` repository to the `upstream` repository.

- [*Branching Strategy*](dev-branching): Describes the branching strategy we will be following.

---

A wrap-up describing what you should know after you went through the instructions on how to contribute to the code, can be found [here](#dev-wrap-up).
