"""
Auxiliary scripts to facilitate the validation and verification of PlatoSim:
    - switch off all effects for the PlatoSim simulations
    - conversion between equatorial and galactic coordinates
    - absolute stellar aberration
"""

from astropy.coordinates import SkyCoord
from math import sin, cos, degrees

def switchOffAllEffects(sim):

    """
    PURPOSE: Switch off all effects for the given simulation:
                - variability;
                - cosmics;
                - jitter;
                - thermo-elastic drift;
                - aberration correction (absolute + differential);
                - field distortion;
                - charge diffusion;
                - jitter smoothing;
                - flatfield;
                - dark signal;
                - brighter-fatter effect (BFE);
                - photon noise;
                - readout noise;
                - charge transfer inefficiency (CTI);
                - open-shutter smearing;
                - quantum efficiency;
                - natural vignetting;
                - mechanical vignetting;
                - polarisation;
                - particulate contamination;
                - molecular contamination;
                - convolution with the PSF;
                - full-well saturation (blooming);
                - digital saturation;
                - quantisation.

    INPUT:
        - sim: Simulation for which to switch off all effects.

    OUTPUT:
        - Simulation with all effects switched off.
    """

    # Sky parameters

    sim["Sky/IncludeVariableSources"] = "no"
    sim["Sky/IncludeCosmicsInSubField"] = "no"
    sim["Sky/IncludeCosmicsInSmearingMap"] = "no"
    sim["Sky/IncludeCosmicsInBiasMap"] = "no"

    # Platform parameters
    
    sim["Platform/UseJitter"] = "no"
    
    # Telescope parameters

    sim["Telescope/UseDrift"] = "no"
    
    # Camera parameters

    sim["Camera/IncludeAberrationCorrection"] = "no"
    sim["Camera/IncludeFieldDistortion"] = "no"
    
    # PSF parameters

    sim["PSF/MappedGaussian/IncludeChargeDiffusion"] = "no"
    sim["PSF/MappedFromFile/IncludeChargeDiffusion"] = "no"
    sim["PSF/AnalyticNonGaussian/IncludeChargeDiffusion"] = "no"
    sim["PSF/MappedGaussian/IncludeJitterSmoothing"] = "no"
    sim["PSF/MappedFromFile/IncludeJitterSmoothing"] = "no"
    
    # CCD parameters

    sim["CCD/IncludeFlatfield"] = "no"
    sim["CCD/IncludeDarkSignal"] = "no"
    sim["CCD/IncludeBFE"] = "no"
    sim["CCD/IncludePhotonNoise"] = "no"
    sim["CCD/IncludeReadoutNoise"] = "no"
    sim["CCD/IncludeCTIeffects"] = "no"
    sim["CCD/IncludeOpenShutterSmearing"] = "no"
    sim["CCD/IncludeQuantumEfficiency"] = "no"
    sim["CCD/IncludeNaturalVignetting"] = "no"
    sim["CCD/IncludeMechanicalVignetting"] = "no"
    sim["CCD/IncludePolarization"] = "no"
    sim["CCD/IncludeParticulateContamination"] = "no"
    sim["CCD/IncludeMolecularContamination"] = "no"
    sim["CCD/IncludeConvolution"] = "no"
    sim["CCD/IncludeFullWellSaturation"] = "no"
    sim["CCD/IncludeDigitalSaturation"] = "no"
    sim["CCD/IncludeQuantisation"] = "no"
    
    return sim





def equatorial2galactic(ra, dec):

    """
    PURPOSE: Convert the given equatorial coordinates to galactic coordinates.

    INPUT:
        - ra: Right ascension [degrees]
        - dec: Declination [degrees]
    
    OUTPUT:
        - galactic longitude [radians]
        - galacti latitude [radians]
    """
    
    coordinates = SkyCoord(ra, dec, frame = 'icrs', unit = 'deg').transform_to('geocentrictrueecliptic')

    return (coordinates.lon.radian, coordinates.lat.radian)





def galactic2equatorial(lon, lat):

    """
    PURPOSE: Convert the given galactic coordinates to equatorial coordinates.

    INPUT:
        - lon: Galactic longitude [radians]
        - dec: Galactic latitude [radians]
    
    OUTPUT:
        - right ascension [degrees]
        - declination [degrees]
    """

    coordinates = SkyCoord(lon, lat, frame = 'geocentrictrueecliptic', unit = 'rad').transform_to('icrs')
    
    return (coordinates.ra.deg, coordinates.dec.deg)





def aberration(lon, lat, lonSpacecraft):

    """
    PURPOSE: Aberrate the given galactic coordinates for the given galactic longitude 
             of the spacecraft in its (circular) orbit around the Sun.
    
    INPUT:
        - lon: Galactic longitude [radians]
        - lat: Galactic latitude [radians]
        - lonSpacecraft: Galactic longitude of the spacecraft in its (circular) orbit
                         around the Sun [radians]
    
    OUTPUT:
        - aberrated galactic longitude [radians]
        - aberrated galactic latitude [radians]
    """

    amplitude = degrees(20.496 / 3600.)

    deltaLon = - amplitude * cos(lon - lonSpacecraft) / cos(lat)
    deltaLat   = - amplitude * sin(lon - lonSpacecraft) * sin(lat)

    return (lon + deltaLon, lat + deltaLat)
