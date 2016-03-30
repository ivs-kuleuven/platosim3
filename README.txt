INSTALLATION OF PLATOSIM
------------------------

The installation of PlatoSim assumes 3 prerequisites:

1) gcc v5.3 or more recent, or, clang v3.3 or more recent
2) cmake v2.8 or more recent (freely downloadable from https://cmake.org/download/)
3) BLAS and LAPACK. Without these, the simulator will likely be slower. These libraries
   come pre-installed on Mac OS X (so Mac users do not have to do anything). Many Linux 
   distributions also standardly have these libraries installed, or offer a package manager
   to easily install them. In case you do need to install them manually, we refer to
   the following websites where you can download them:
        http://www.openblas.net
        http://www.netlib.org/lapack/

PlatoSim relies on a number of other dependencies, which are all included in this package
for the convenience of the user. To build, and install them, you can run the install
bash script:

$ ./install.sh

This install script requires Python on your system, and runs the Python install 
scripts in dependencies/installscripts/. This may take a while. You may get some 
warnings, this is normal. Afterwards it also builds the PlatoSim simulator itself.

In the build/ folder you should find two executables: testplatosim and platosim. 
Executing 

$ ./testplatosim

runs all the unit tests.

If, for some reason, you ever want to recompile the simulator, simply:

$ cd build
$ rm *
$ cmake ..
# make



RUNNING PLATOSIM
----------------

Running PlatoSim is done with

$ cd build
$ ./platosim ../inputfiles/myInputfile.yaml myOutputfile.hdf5

The first argument is the simulation configuration input file, of which you can find an
example 'inputfile.yaml' in the folder inputfiles/. You can copy this configuration file
and adapt the line

ProjectLocation:         /Users/rik/Git/PlatoSim3

to your own location of PlatoSim.

The second argument is the name of the (non-existing!) HDF5 file to which all simulation 
output is written. Apart from this HDF5 file, the simulator also writes logging statements
to a log.txt file.




ACCESSING PLATOSIM'S OUTPUT
---------------------------

PlatoSim writes its output to an HDF5 file. HDF stands for Hierarchical Data Format, and is
a next generation file format that was specifically designed to store and organize large 
amounts of data.

HDF5 behaves much likes a unix-like folder structure, but where folders are called groups.
Each group can contain other groups, array datasets, and scalar attributes. For example,
the first subfield image is located in

/Images/image000000 

in the HDF5 file.

To quickly list the contents of the group structure of an HDF5 file on the command line, 
make sure that your PATH environment variable includes 

    dependencies/Installs/hdf5-1.8.16/bin/

so that you can execute

$ h5ls myOutputfile.hdf5

or e.g.

$ h5ls myOutputfile.hdf5/StarCatalog


For PYTHON USERS, we provided a simfile.py module in the python/ folder, with convenients 
tools to extract and plot the Simulator output. For example, one can plot a subfield image using

>>> from simfile import *
>>> myFile = SimFile("myOutputfile.hdf5")
>>> myFile.showImage(0)

The top of simfile.py contains documentation with several examples.


Still for PYTHON USERS, you can also access the HDF5 file using the pytables module. For 
example (using the latest version of PyTables):

>>> import tables as tbl
>>> myFile = tbl.open_file("myOutputfile.hdf5", "r")
>>> image = myFile.root.Images.image000000
>>> imshow(image, interpolation="nearest", origin="lower")
>>> myFile.root.InputParameters.CCD._v_attrs
[...]
>>> print(myFile.root.InputParameters.CCD._v_attrs.Gain)


For IDL USERS, you can access the HDF5 file using, for example, 

IDL> path = FILEPATH(“Simul01.hdf5")
IDL> file = H5F_OPEN(path)
IDL> contents = H5_PARSE(path)
IDL> help, contents, /STRUCTURE
[...]

IDL> help, contents.Images, /STRUCTURE
[...]

IDL> dataset = H5D_OPEN(file,'/Images/image000000') 
IDL> image = H5D_READ(dataset)
IDL> print, size(image)



