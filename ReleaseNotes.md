* IMPORTANT:

This new release has many changes that might be incompatble with previous versions. In order to keep using old input files make sure to make the following changes:

Delete:
"""
ObservingParameters:
    RApointing:
    DecPointing:
Sky:
    SkyBackground:
Platform:
    SolarPanelOrientation:
Telescope:
    UseDriftFromFile:
Camera:
    FieldDistortion:
        Type:
"""

Add:
"""
Sky:
    SkyBackground:
      UseConstantSkyBackground:
      BackgroundValue:
Platform:
    Orientation:
      Source:
      Angles:
         RAPointing:
        DecPointing:
        SolarPanelOrientation:
      Quaternion:
        Components:
Telescope:
    DriftSource:
ControlHDF5Content:
    GroupByExposure:
    WriteFlatfieldMap:
    WriteTransmissionEfficiency:
    WriteBackgroundMap:
    WriteCTI:
    WriteDiffusedPSF:
    WriteHighResolutionPSF:
        WriteTelescopeACS:
    WriteStarPositions:
    WriteGhostPositions:
    WriteCosmics:
"""

Also the PSF files have been changed. There was a bug where the old mapped PSF files were not consistent with the analytic PSF files. If such a thing is needed for your
simulation make sure you download the most recent PSF files. (see website)




* Changelog PlatoSim 3.6.0

** Improvements

*** Change of YAML inputfile `Telescope` block named `UseDriftFromFile` to `DriftSource` (Issue #766)

*** Changed input structure of input files:
    - `ObservingParameters/RApointing` -> `Platform/Orientation/Angles/RAPointing`
    - `ObservingParameters/DecPointing` -> `latform/Orientation/Angles/DecPointing`

*** Changed validation tests to deal with new input structure.

*** Changed the naming convention SC (spacecraft) into PLM (Payload Module) in python code.

*** Made absolute aberation test easier to read.

*** Removed ``doc/Validation`` directory. The content of this directory is mostly containd the ``tests/validationTests`` directory.

*** Interpolated the ``skyBackground`` map so that the entire sky is filled.

*** When the ``WriteBackgroundMap`` option is true, we will either save a time series of the background value (if we use a constant background) or the background map at the beginning of the simulation (if we use a variable background map).

*** Changed `simfile.py` functions to remove and cleanup many of the repeating routines used.

*** Changed default jitter timescale to 250s instead of 3600s.

*** Updated the old Doxygen documentation to a more user-friendly Sphinx documentation.

*** Updated the old Jupyter notebook tutorials. Each notebook follows a cronological order and is inline with the changes made to the YAML file.

*** Updated all Python function docstrings (to a appropiate Pythonic version).

*** Changed directory structure of `/python/platosim`. All command line scripts are now placed within the folder `script`.

*** Change of HDF5 structure for `Cosmics`, now using consistent upper case letter for `Exposure` (Issue #765)





** Bug fixes

*** Changed Mapped distortion and  Mapped inverse distortion routines so that they now are each others inverse. Also changed the corresponding python scripts.

*** Fixed that `simfile.showImage()` was not scaling the correct flux in ADU.

*** Fixed log color scaling bug in `simfile.showImage()`.

*** Fixed bug in `mappedGaussianPSF` validation test.

*** Fixed bug where generating throughputmap, automatically assumed distortion should be used. For analytic PSF without distortion this gave a error.





** New features/functionality

*** Added option to save high resolution, analytic non-Gaussian PSF to HDF5.

*** Added variable background (Issue #729)

*** Added the option to use Quaternions in the input file. (Issue #709)

*** Added option to change the output format of the HDF5 output file. This speeds up long simulations.

*** Added option to use differently formatted HDF5-PSF files. (Issue #811)

*** Added PLATOnium toolkit (see the new documentation page).

*** Added Poetry installation for developers (see the new documentation).
