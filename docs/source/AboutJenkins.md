# A Word about Jenkins {#AboutJenkins}

We have started using <a href="https://jenkins.io/">Jenkins</a> to automatically build PlatoSim and make pre-built software available for a myriad of operating systems.  The figure below summarises how we want to use it.

Each time code is pushed to the repository (either to the <code>master</code> or the <code>develop</code> branch) or a pull request is merged in GitHub, Jenkins will start building the new code, resolves the dependencies for you, and makes it available via the conda installation command.  In case the build was successful, users can install this build on their system.

To monitor the status of the PlatoSim builds, check <a href="https://jenkins.miricle.org/view/Platosim/">here</a>.

@image html /images/jenkins.png "" 