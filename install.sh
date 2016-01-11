#!/bin/bash
#
# Just run this script:
# $ ./install.sh

#1. Install dependencies
#-----------------------

# There is one python script under installscripts/ folder for each dependency package.
# Each script Unzip/Untar packages in dependencyDownloads/ and
# installs dependencies packages in dependencyInstalls/

if [ ! -d ./dependencies/Installs ]; then
    mkdir ./dependencies/Installs
fi

python ./installscripts/install_hdf5.py
python ./installscripts/install_googletest.py


# If any of these packages is already installed, just comment the 
# corresponding line (please, be sure that the commented package is correctly installed in 
# your system, otherwise the installation will raise errors!!).


#2. Build and install the simulator
#----------------------------------

if [ ! -d ./build ]; then
    mkdir build
fi

cd build
cmake ..
make -j 4 

