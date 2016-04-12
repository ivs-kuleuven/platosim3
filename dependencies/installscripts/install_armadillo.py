#
# install.sh runs automatically this script. It can be done manually running (from the same folder than install.sh):
# $ python ./dependencies/installscripts/install_armadillo.py
#


import os,sys,shutil,subprocess


# Specify the dependency package name

packageName = "armadillo-6.500.4"

# Specify build and install folders

currentWorkingDir = os.getcwd()
buildDir   = currentWorkingDir + "/dependencies/Downloads/"
installDir = currentWorkingDir + "/dependencies/Installs/" + packageName 

# If there is an older version of the install dir, make sure to remove it first.

shutil.rmtree(installDir, ignore_errors=True)

# Print a banner

print("\n\n\n")
print("====================")
print("Installing ARMADILLO")
print("====================")
print("\n")

# Build and install package

installProcedure = "cd {build};                                     \
                    tar -xzvf {package}.tgz;                        \
                    cd {package};                                   \
                    mkdir build;                                    \
                    cd build;                                       \
                    cmake ..;                                       \
                    make;                                           \
                    make install DESTDIR={install}".format(build=buildDir, package=packageName, install=installDir)

subprocess.call(installProcedure, shell=True)


# Armadillo installs the libraries and header files in a folder structure depending on the OS.
# This folder structure is different from what the other dependency packages use, so to harmonize
# structure, we copy the old folder structure into a new one.

# First try to copy the header files to the folder /include.

if os.path.isdir(installDir+"/usr/local/include"):
    shutil.copytree(installDir+"/usr/local/include", installDir+"/include")
elif os.path.isdir(installDir+"/usr/include"):
    shutil.copytree(installDir+"/usr/include", installDir+"/include")    
else:
    print("Armadillo installation failed: please email the folder structure in " + installDir + "to Joris De Ridder")
    exit(1)

# Then try to copy the library files to the folder /lib

if os.path.isdir(installDir+"/usr/local/lib"):
    shutil.copytree(installDir+"/usr/local/lib", installDir+"/lib")
elif os.path.isdir(installDir+"/usr/local/lib64"):
    shutil.copytree(installDir+"/usr/local/lib64", installDir+"/lib")
elif os.path.isdir(installDir+"/usr/lib64"):
    shutil.copytree(installDir+"/usr/lib64", installDir+"/lib")
elif os.path.isdir(installDir+"/usr/lib"):
    shutil.copytree(installDir+"/usr/lib", installDir+"/lib")
else:
    print("Armadillo installation failed: please email the folder structure in " + installDir + "to Joris De Ridder")
    exit(1)


# If we arrived here, the copying succeeded, and we need to remove the old folder

shutil.rmtree(installDir+"/usr", ignore_errors=True)


# For Mac systems, we still need to correct the relative path in the armadillo library

if sys.platform == "darwin":

    correctionProcedure = "cd {install}/lib;                                                            \
                           install_name_tool -id {install}/lib/libarmadillo.dylib libarmadillo.dylib    \
                          ".format(install=installDir)

    subprocess.call(correctionProcedure, shell=True)



# After installation in the install folder, remove the decompressed package folder in 
# the Downloads dir so that only the .tgz file remains in the Downloads dir.

shutil.rmtree(buildDir+packageName, ignore_errors=True)


