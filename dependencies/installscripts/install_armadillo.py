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


# Armadillo installs the libraries and header files differently depending on the operating system,
# which needs we need to distinguish between Mac OS X and Linux in what follows.

if sys.platform == "darwin":

    # make installed Armadillo locally, but still with a /usr/local/ directory structure.
    # Copy the include and lib folders in the install folder, and remove the /usr/local folder.

    shutil.copytree(installDir+"/usr/local/include", installDir+"/include")
    shutil.copytree(installDir+"/usr/local/lib", installDir+"/lib")
    shutil.rmtree(installDir+"/usr", ignore_errors=True)

    # Correct the relative path in the armadillo.6 library

    correctionProcedure = "cd {install}/lib;                                                            \
                           install_name_tool -id {install}/lib/libarmadillo.dylib libarmadillo.dylib    \
                          ".format(install=installDir)

    subprocess.call(correctionProcedure, shell=True)

elif sys.platform.startswith("linux"):

    # Note: in Python 2.x sys.platform = "linux2", in Python 3.x sys.platform = "linux",
    #       hence the startswith().

    # make installed Armadillo locally, but still with a /usr/local/ directory structure.
    # Copy the include and lib folders in the install folder, and remove the /usr/local folder.

    shutil.copytree(installDir+"/usr/local/include", installDir+"/include")
    shutil.copytree(installDir+"/usr/local/lib64", installDir+"/lib")
    shutil.rmtree(installDir+"/usr", ignore_errors=True)






# After installation in the install folder, remove the decompressed package folder in 
# the Downloads dir so that only the .tgz file remains in the Downloads dir.

shutil.rmtree(buildDir+packageName, ignore_errors=True)


