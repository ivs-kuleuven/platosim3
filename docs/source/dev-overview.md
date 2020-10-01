# Installing/Updating PlatoSim3 for Developers {#dev-overview}

The scheme below shows an overview of the procedure to download, update, and build the PlatoSim3 software on your system yourself, without using conda.  This scheme is valid for developers who may want to contribute to the code at some point (or for users who want to install the software without @ref InstallViaConda "using conda").  Also the procedure to make your own modifications to the code available for others is described.

@image html /images/download-update-build.png ""

---

## Overview

The following installation guides are available:

- @ref dev-prerequisites "**Prerequisites**": Describes the prerequisites that must be fulfilled to be able to download, update, and build PlatoSim3.

- @ref dev-fork-clone "**Forking & Cloning the Repository**": Describes the first-time installation of the software on your local machine.

- @ref dev-dependencies "**Installing Dependencies**": Describes how to install required dependencies.

- @ref dev-build "**Building**": Describes how to build the code.

- @ref dev-add-upstream "**Remote Repositories**": Describes how to configure your GitHub setup, so that you will be able to get hold of any changes to the code by others.

- @ref dev-pull "**Update procedure**": Describes how to get the latest version of the code (from the `upstream`) on your local machine.

- @ref dev-push "**Contributing**": Describes how to transfer your local changes to the `origin` repository and how to transfer these changes from the `origin` repository to the `upstream` repository.

- @ref dev-branching "**Branching Strategy**": Describes the branching strategy we will be following.


---

A wrap-up describing what you should know after you went through the instructions on how to contribute to the code, can be found @ref dev-wrap-up "here".
