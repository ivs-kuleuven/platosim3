#
# install.sh runs automatically this script. It can be done manually running (from the same folder than install.sh):
# $ python ./dependencies/installscripts/install_fftw.py
#

import os,sys,shutil,subprocess


# Specify the dependency package name

packageName = "fftw-3.3.4"

# Specify build and install folders

currentWorkingDir = os.getcwd()
buildDir   = currentWorkingDir + "/dependencies/Downloads/"
installDir = currentWorkingDir + "/dependencies/Installs/" + packageName 

# If there is an older version of the install dir, make sure to remove it first.

shutil.rmtree(installDir, ignore_errors=True)

# Print a banner

print("\n\n\n")
print("===============")
print("Installing FFTW")
print("===============")
print("\n")

# Build and install package

installProcedure = "cd {build};                                                      \
                    tar -xzvf {package}.tgz;                                         \
                    cd {package};                                                    \
                    ./configure --prefix={install} --enable-threads --enable-float;  \
                    make;                                                            \
                    make install".format(build=buildDir, package=packageName, install=installDir)

subprocess.call(installProcedure, shell=True)


# After installation in the install folder, remove the decompressed package folder in 
# the Downloads dir so that only the .tgz file remains in the Downloads dir.

shutil.rmtree(buildDir+packageName, ignore_errors=True)
