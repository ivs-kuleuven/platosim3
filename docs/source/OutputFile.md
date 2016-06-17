# Description of the Output File {#OutputFileDescription}

PlatoSim3 writes its output to an <a href="https://www.hdfgroup.org/HDF5/">HDF5</a> file. HDF stands for Hierarchical Data Format, and is a next generation file format that was specifically designed to store and organise large amounts of data.

HDF5 behaves much likes a Unix-like folder structure, but where folders are called groups.  Each group can contain other groups, array datasets, and scalar attributes. For example, the first subfield image is located in <code>/Images/image000000</code> in the HDF5 file.





<!-- ********* -->
<!-- Structure -->
<!-- ********* -->

## Structure

The general structure of the output file if shown in Fig. 1.  The colour code in this diagram is as follows:
- purple: groups (which contain other groups and/or datasets and/or attributes)
- green: datasets
- orange: attributes

@image html /images/outputStructure.png "Figure 1: Structure of the HDF5 output file."





<!-- ******************************************** -->
<!-- Inspecting the Structure & Content in Python -->
<!-- ******************************************** -->

## Inspecting the Structure & Content in Python

<!-- How-To -->
<!-- ****** -->

### How-To

Under <code>/docs/tutorials/InspectHDF5</code> you can find @ref Tutorials "tutorials" that demonstrate how you can inspect the structure and the content of the HDF5 output file.  This can be done with the <code>SimFile</code> class we offer in the <code>/python/simfile.py</code> module (as described in <code>Inspect HDF5 using simFile.ipynb</code>), or with one of the following packages:

- <a href="http://www.h5py.org/">h5py</a> (as decribed in <code>Inspect HDF5 using h5py.ipynb</code>)
- and <a href="http://www.pytables.org/">PyTables</a> (as described in <code>Inspect HDF5 using PyTables.ipynb</code>).

These packages have to be installed with <code>conda install</code>.  To use them, you must do the following imports and create the following objects: 

\code{.py}
# To use the SimFile class from PlatoSim3 itself

import simfile
from simfile import SimFile
simFile = SimFile(<full path to the HDF5 file>)

# To use the h5py module

import h5py
h5pyFile = h5py.File(<full path to the HDF5 file>)

# To use the PyTables module

import tables
pytablesFile = tables.openFile(<full path to the HDF5 file>)
\endcode

In the following sections, we will discuss which groups are present in the output HDF5 files and which information their datasets and attributes contain.  To access a dataset in a given group, the following commands can be used:

\code{.py}
simDataset = simFile.hdf5file["<group name>"]["<dataset name>"]

h5pyDataset = h5pyFile["/<group name>/<dataset name>"]

pytablesDataset = pytablesFile.root.<group name>.<dataset name>
\endcode 

and to access an attribute in a given group:

\code{.py}
simAttribute = simFile.hdf5file["<group name>"].attrs[attribute name>]

h5pyAttribute = h5pyFile["<group name>"].attrs["<attribute name>"]

pytablesAttribute = pytablesFile.root.<group name>._v_attrs.<attribute name>
\endcode

You can use <a href="http://matplotlib.org/users/image_tutorial.html">matplotlib.pyplot.imshow</a> to visualise 2D datasets (images, bias maps, smearing maps, etc.), and <a href="http://matplotlib.org/api/pyplot_api.html">matplotlib.pyplot.plot and matplotlib.pyplot.scatter</a> to plot the columns of the datasets against each other.  In some cases, <code>SimFile</code> offers alternative methods to access and visualise dataset.  They will be mentioned in the appropriate section below.

To get an overview of the groups and attributes, you can use the following commands:

<!-- \code{.py}
simFile.hdf5file.keys()								# Main groups
simFile.hdf5file.attrs.keys()						# Main attributes (None)
simFile.hdf5file["<group name>"].keys()				# Sub-groups in the given group
simFile.hdf5file["<group name>"].attrs.keys()		# Attributes in the given group

h5pyFile.hdf5file.keys()								# Main groups
h5pyFile.hdf5file.attrs.keys()						# Main attributes (None)
h5pyFile.hdf5file["<group name>"].keys()				# Sub-groups in the given group
h5pyFile.hdf5file["<group name>"].attrs.keys()		# Attributes in the given group

pytables.root										# Main groups
pytables.root._v_attrs								# Main attributes (None)
pytables.root.<group name>							# Sub-groups in the given group
pytables.<group name>._v_attrs						# Attributes in the given group
\endcode -->



<!-- Configuration Parameters -->
<!-- *********************** -->

### Configuration Parameters

To enable you to trace back which configuration parameters you have used to generate a certain HDF5 file, we have added the **InputFiles** group to the output HDF5 files, which contains a copy of the YAML file that was used as input.

To check the value of an individual configuration parameter (from the stored YAML file), the following commands can be used:

\code{.py}
simValue1 = simFile.getInputParameter("<block name in the YAML file>", "<parameter name in the YAML file>")
simValue2 = simFile.hdf5file["InputParameters"]["<block name in the YAML file>"].attrs["<parameter name in the YAML file>"]

h5pyValue = h5pyFile["InputParameters"]["<block name in the YAML file>"].attrs["<parameter name in the YAML file>"]

pytablesValue = pytablesFile.root.InputParameters.<block name in the YAML file>._v_attrs.<parameter name in the YAML file>
\endcode

Have a look @ref InputFileDescription "here" to find out the names of the blocks and parameters in the configuration file.



<!-- Images, Sub-Images, and Imagettes -->
<!-- ********************************* -->

### Images, Sub-Images, and Imagettes

The **Images** group contains one dataset, **images<exposure number with 6 digits>**, per exposure (counting starts at 000000).  An alternative way to get hold of such a dataset is

\code{.py}
simImage = simFile.getImage(<exposure number>)
\endcode

Visualisation can be done as follows:

\code{.py}
simFile.showImage(<exposure number>)
\endcode

In case you have stored the images at sub-pixel level too, you can find these in the **SubPixelImages** group, named **subPixelImage<exposure number with 6 digits>**.

To retrieve a small square imagette, centred about a star with the given ID, use the command

\code{.py}
simImagette = simFile.getImagette(<star ID>, <number of the exposure>, <radius>)
\endcode

where the output imagette has size 2 * radius + 1 pixels in both directions (expressed in pixels).



<!-- Smearing Maps -->
<!-- ************* -->

### Smearing Maps

The **SmearingMaps** group contains one dataset, **smearingMap<exposure number with 6 digits>, per exposure (counting starts at 000000).



<!-- Bias Maps -->
<!-- ********* -->

### Bias Maps

The **BiasMaps** group contains one dataset, **biasMap<exposure number with 6 digits>, per exposure (counting starts at 000000).



<!-- Star Catalogue -->
<!-- ************** -->

### Star Catalogue

The **StarCatalog** group contains the following datasets:

- **starIDs**: sequential number of those stars in the input star catalogue that were detected on the sub-field in at least one exposure;
- **RA**: right ascension of these stars, expressed in degrees;
- **Dec**: declination of these stars, expressed in degrees;
- **Vmag**: V magnitude of these stars;
- **xFPmm** and **yFPmm**: initial planar focal-plane coordinates of these stars, expressed in mm;
- **rowPix** and **colPix**: initial pixel coordinates of these stars (float values).

You can get hold of all these datasets in one go, with the following command:

\code{.py}
simStarIds, simRa, simDec, simVmag, simXFPmm, simYFPmm, simRowPix, simColPix = simFile.getStarCatalog()
\endcode



<!-- Star Positions -->
<!-- ************** -->

### Star Positions

The **StarPositions** group contains one group, **Exposure<exposure number with 6 digits>**, per exposure, which comprises the following datasets:

- **starID**: sequential number of those stars that are visible in the current sub-field;
- **rowPix** and **colPix**: pixel coordinates of those stars, averaged over the duration of the exposure (float values);
- **xFPmm** and **yFPmm**: planar focal-plane coordinates of those stars, averaged over the duration of the exposure and expressed in mm;
- **flux**: number of photons of those stars that got detected in the sub-field during the exposure.

To get hold of these dataset in one go, you can use the following command:

\code{.py}
simStarIds, simRowPix, simColPix, simXFPmm, simYFPmm = simFile.getStarCoordinates(<exposureNumber>, [minVmag = <minimum V magnitude in the selection>], [maxVmag = <maximum V magnitude in the selection>])
\endcode



<!-- Background -->
<!-- ********** -->

### Background

The **Background** contains one dataset, **skyBackground**, which contains the background (one entry per exposure), expressed in photons / pixel / exposure.



<!-- Spacecraft Attitude -->
<!-- ******************* -->

### Attitude

The **ACS** group contains the information concerning the spacecraft attitude, organised in the following groups:

- **Time**: time elapsed at the end of each exposure since the start of the simulation, expressed in seconds;
- **Yaw**: jitter yaw angle at the end of each exposure, expressed in arcsec;
- **Pitch**: jitter pitch angle at the end of each exposure, expressed in arcsec; 
- **Roll**: jitter roll angle at the end of each exposure, expressed in arcsec;
- **PlatformRA** and **PlatformDec**: right ascension and declination of the platform pointing at the end of each exposure, expressed in degrees; 

To get hold of the Euler angles and of the platform pointing, you can use the following commands:

\code{.py}
simJitterYaw, simJitterPitch, simJitterRoll = simFile.getYawPitchRoll()
simPlatformRa, simPlatformDec = simFile.getPlatformPointingCoordinates() 

\endcode



<!-- Telescope -->
<!-- ********* -->

### Telescope

The **Telescope** group contains the following datasets:

- **Time**: time elapsed at the end of each exposure since the start of the simulation, expressed in seconds;
- **TelescopeYaw**: drift yaw angle at the end of each exposure, expressed in arcsec;
- **TelescopePitch**: drift pitch angle at the end of each exposure, expressed in arcsec;
- **TelescopeRoll**: drift roll angle at the end of each exposure, expressed in arcsec;
- **TelescopeRA** and **TelescopeDec**: right ascension and declination of the telescope pointing at the end of each exposure, expressed in degrees;
- **Azimuth**: azimuth angle of the telescope for each exposure, expressed in degrees;
- **Tilt**: tilt angle of the telescope for each exposure, expressed in degrees;
- **FocalPlaneOrientation**: focal-plane orientation angle for each exposure, expressed in degrees.



<!-- PSF -->
<!-- *** -->

### PSF

The **PSF** group contains the following datasets:

- **rebinnedPSFpixel**: selected PSF rebinned to the image's number of pixels;
- **rebinnedPSFsubPixel**: selected PSF rebinned to the image's number of sub-pixels;
- **rotatedPSF**: 

and the following attributes:

- **rotationAngle**: angle over which the PSF was rotated, expressed in degrees;
- **selectedPSF**: information on which PSF was selected:
	- for a Gaussian PSF this reports the standard deviation of the Gaussian, expressed in pixels
	- for a pre-computed PSF this reports which one was selected from the PSF HDF5 file.

You can also get hold of the datasets as follows:

\code{.py}
simPsf = simFile.getPsf(<name of the dataset>)
\endcode

Visualisation can be done as follows:

\code{.py}
simFile.showPSF(<name of the dataset>)
\endcode



<!-- Vignetting -->
<!-- ********** -->

### Vignetting

The vignetting map is contained by the **vignettingMap** dataset in the **VignettingMap** group and can be accessed with the following commands:



<!-- Flatfield -->
<!-- ********* -->

### Flatfield

The **Flatfield** group contains two datasets:

- **PRNU**: pixel response non-uniformity map for the simulated sub-field;
- **IRNU**: intra-pixel response non-uniformity for the simulated sub-field.

These can be accessed with the following commands:

\code{.py}
simPrnu = simFile.getPRNU()
simIrnu = simFile.getIRNU()
\endcode



<!-- Version -->
<!-- ******* -->

### Version

The **Version** group contains two attributes:

- **Application**: the application that was used to generate the HDF5 file (i.c. PlatoSim3);
- **GitVersion**: version of the PlatoSim3 software that was used to generate the HDF5 file.





<!-- ********************** -->
<!-- Alternatives to Python -->
<!-- ********************** -->

## Alternatives to Python




<!-- IDL -->
<!-- *** -->

### IDL

IDL users can access the HDF5 file using, for example, 

\code{.idl}
path = FILEPATH(“Simul01.hdf5")
file = H5F_OPEN(path)
contents = H5_PARSE(path)
help, contents, /STRUCTURE
...
help, contents.Images, /STRUCTURE
...
dataset = H5D_OPEN(file,'/Images/image000000') 
image = H5D_READ(dataset)
print, size(image)
\endcode



<!-- Visualisation Tools -->
<!-- ******************* -->

### Visualisation Tools

There are two alternatives to Python that involve no coding we have been using ourselves (to quickly check the data):

- <a href="https://www.hdfgroup.org/products/java/hdfview/">HDFView</a>
- <a href="https://www.hdfgroup.org/projects/compass/">HDF Compass</a>