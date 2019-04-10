def switchAllEffectsOff(sim):

    """
    PURPOSE: Switch off all effects for the given simulation.

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
    