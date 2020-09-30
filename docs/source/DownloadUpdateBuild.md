# Downloading, Updating, and Building PlatoSim3 for Developers{#DownloadUpdateBuild}

The scheme below shows an overview of the procedure to download, update, and build the PlatoSim3 software on your system yourself, without using conda.  This scheme is valid for developers who may want to contribute to the code at some point (or for users who want to install the software without @ref InstallViaConda "using conda").

We use Git as version control system for PlatoSim3, and the software can be forked and/or cloned from GitHub, depending on whether or not you plan to contribute to the software, as described @ref Downloading "here".  

Please, note there are a number of  @ref ReqsInstallUpdate "prerequisites" to be able to @ref Downloading "download", @ref Updating "update", and @ref Building "build" PlatoSim3.

The described procedure enables you to @ref Updating "update" your PlatoSim3 installation to the latest version on GitHub in a straightforward way, even without being (or having to become) a Git expert.

Once you have downloaded or updated the software, you have to @ref Building "build" it in order to be able to @ref Running "run" it.  The first time you download the software and on the rare occasions we update the dependencies (i.e. software packages we use in the code and that are distributed together with the PlatoSim3 code), you will have to resolve the dependencies (as described together with the build procedure).

@image html /images/download-update-build.png ""