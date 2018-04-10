#
# install.sh runs automatically this script. It can be done manually running (from the same folder than install.sh):
# $ python ./dependencies/installscripts/install_faddeeva.py
#


import os,sys,shutil,subprocess, glob


# Specify the dependency package name

packageName = "Faddeeva"

# Specify build and install folders

currentWorkingDir = os.getcwd()
buildDir   = currentWorkingDir + "/dependencies/Downloads/"
installDir = currentWorkingDir + "/dependencies/Installs/" + packageName 

# If there is an older version of the install dir, make sure to remove it first.

shutil.rmtree(installDir, ignore_errors=True)

# Print a banner

print("\n\n\n")
print("===================")
print("Installing FADDEEVA")
print("===================")
print("\n")

# Build and install package

installProcedure = "cd {build};                                     \
                    tar -xvf {package}.tgz;                         \
                    cd {package};                                   \
                    mkdir build;                                    \
                    cd build;                                       \
                    cmake ..;                                       \
                    make;                                           \
                    make install".format(build=buildDir, package=packageName)

subprocess.call(installProcedure, shell=True)


# For Mac systems, we still need to correct the relative path in the armadillo library

if sys.platform == "darwin":

    correctionProcedure = "cd {install}/lib;                                                            \
                           install_name_tool -id {install}/lib/libFaddeeva.dylib libFaddeeva.dylib    \
                          ".format(install=installDir)

    subprocess.call(correctionProcedure, shell=True)



# After installation in the install folder, remove the decompressed package folder in 
# the Downloads dir so that only the .tgz file remains in the Downloads dir.

shutil.rmtree(buildDir+packageName, ignore_errors=True)



