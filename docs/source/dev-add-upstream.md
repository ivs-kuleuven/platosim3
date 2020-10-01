# Remote Repositories {#dev-add-upstream}

To be able to "pull" the latest version of the upstream software into your local copy (we'll come back to this), we advise you to add upstream to your list of remote repositories. To check your list of remote repositories, execute the following command in the the installation folder (or one of its sub-folders):

    $ git remote -v

This should give output, similar to this:

    $ origin    https://github.com/<your GitHub username>/PlatoSim3.git (fetch)
    $ origin    https://github.com/<your GitHub username>/PlatoSim3.git (push)
    $ upstream  https://github.com/IvS-KULeuven/PlatoSim3.git (fetch)
    $ upstream  https://github.com/IvS-KULeuven/PlatoSim3.git (push)
If there's no sign of the upstream (the last two lines), you can add it with the following command:

    $ git remote add upstream https://github.com/IvS-KULeuven/PlatoSim3.git

In case you pointed upstream to the wrong location (i.e. you used the command from above with the wrong link), you can undo this by executing the following command in the the installation folder (or one of its sub-folders):

    $ git remote rm upstream
