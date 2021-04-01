
import numpy as np





def stellarFlux(Vmag, exposureTime, fluxm0=1.00238e8, throughputBandwidth=400, transmissionEfficiency=0.76, 
                lightCollectingArea=0.01131, quantumEfficiency=0.87):

    """
    PURPOSE: compute the stellar flux (electrons / exposure) given the instrumental characteristics

    INPUT: Vmag:                   Johnson V magnitude
           exposureTime:           Exposure time (without the readout) [s]
           fluxm0:                 Photon flux of a V=0 star (default SpT=G2V) [phot/s/m^2/nm]
           throughputBandwidth:    FWHM [nm]
           transmissionEfficiency: In [0,1]
           lightCollectingArea:    Of the telescope [m^2]
           quantumEfficiency:      In [0,1]

    OUTPUT: flux: [e-/exposure]
    """

    photonFlux = fluxm0 * throughputBandwidth * transmissionEfficiency * lightCollectingArea * pow(10.0, -0.4 * Vmag) * exposureTime
    electronFlux = photonFlux * quantumEfficiency

    return electronFlux


