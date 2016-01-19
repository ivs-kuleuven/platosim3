#
# install.sh runs automatically this script. It can be done manually running (from the same folder than install.sh):
# $ python ./installscripts/install_yaml.py
#


import os,shutil,subprocess


# Specify the dependency package name

packageName = "armadillo-6.400.3"

# Specify build and install folders

currentWorkingDir = os.getcwd()
buildDir   = currentWorkingDir + "/dependencies/Downloads/"
installDir = currentWorkingDir + "/dependencies/Installs/" + packageName 

# If there is an older version of the install dir, make sure to remove it first.

shutil.rmtree(installDir, ignore_errors=True)


# Build and nstall package

installProcedure = "cd {build};                                     \
					tar -xzvf {package}.tgz;                        \
					cd {package};                                   \
                    mkdir build;                                    \
                    cd build;                                       \
                    cmake ..;                                       \
                    make;                                           \
                    make install DESTDIR={install}".format(build=buildDir, package=packageName, install=installDir)

subprocess.call(installProcedure, shell=True)

# make installed Armadillo locally, but still with a /usr/local/ directory structure.
# Copy the include and lib folders in the install folder, and remove the /usr/local folder.

shutil.copytree(installDir+"/usr/local/include", installDir+"/include")
shutil.copytree(installDir+"/usr/local/lib", installDir+"/lib")
shutil.rmtree(installDir+"/usr", ignore_errors=True)

# After installation in the install folder, remove the decompressed package folder in 
# the Downloads dir so that only the .tgz file remains in the Downloads dir.

shutil.rmtree(buildDir+packageName, ignore_errors=True)


