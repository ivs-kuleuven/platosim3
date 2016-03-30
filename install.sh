#!/bin/bash
#
# Just run this script:
# $ ./install.sh


# Stop this script if we encounter an error with one of the packages

set -e

# Install dependencies
# There is one python script under dependencies/installscripts/ folder for each 
# dependency package. Each script Unzip/Untar packages in dependencies/Downloads/ 
# and installs dependencies packages in dependencies/Installs/

if [ ! -d ./dependencies/Installs ]; then
    mkdir ./dependencies/Installs
fi

python ./dependencies/installscripts/install_hdf5.py
python ./dependencies/installscripts/install_googletest.py
python ./dependencies/installscripts/install_yaml.py
python ./dependencies/installscripts/install_armadillo.py
python ./dependencies/installscripts/install_fftw.py


# Build the simulator

if [ ! -d ./build ]; then
    mkdir build
fi

cd build
cmake ..
make -j 4 

