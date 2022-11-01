#!/bin/bash
#
# Just run this script:
# $ ./install.sh


# Stop this script if we encounter an error with one of the packages

set -e


# Check if the "python" command exists, if not try "python3".

myPython=python
if ! [ -x "$(command -v python)" ]; then
  echo "Didn't find a python executable, trying python3." 
  myPython=python3
fi


# Install dependencies
# There is one python script under dependencies/installscripts/ folder for each 
# dependency package. Each script Unzip/Untar packages in dependencies/Downloads/ 
# and installs dependencies packages in dependencies/Installs/

$myPython ./dependencies/installscripts/install_hdf5.py
$myPython ./dependencies/installscripts/install_yaml.py
$myPython ./dependencies/installscripts/install_armadillo.py
$myPython ./dependencies/installscripts/install_fftw.py
$myPython ./dependencies/installscripts/install_faddeeva.py
$myPython ./dependencies/installscripts/install_zeromq.py


# Build the simulator

if [ ! -d ./build ]; then
    mkdir build
fi

cd build
cmake ..
make -j 4 

