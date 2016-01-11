#!/bin/bash
#
# Just run this script:
# $ ./install.sh

# Install dependencies
# There is one python script under dependencies/installscripts/ folder for each 
# dependency package. Each script Unzip/Untar packages in dependencies/Downloads/ 
# and installs dependencies packages in dependencies/Installs/

if [ ! -d ./dependencies/Installs ]; then
    mkdir ./dependencies/Installs
fi

python ./dependencies/installscripts/install_hdf5.py
python ./dependencies/installscripts/install_googletest.py


# Build the simulator

if [ ! -d ./build ]; then
    mkdir build
fi

cd build
cmake ..
make -j 4 

