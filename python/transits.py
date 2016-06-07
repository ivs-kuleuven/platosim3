"""
This file contains tools to generate artificial stellar flux time series containing exoplanet transits

"""


import sys
import numpy as np
from numpy.random import normal







def simpleTransit(time, t0, flatPartDuration, transitDuration, orbitalPeriod, relativeDepth):

    """
    PURPOSE: generate a time series containing exoplanet transits. The transits are modelled 
             with a linear ingress, a flat part, and a linear egress.

    INPUT: time:             Time points for which the signal needs to be constructed [Ms]             
           t0:               Reference time [Ms]            
           flatPartDuration: Duration of the time when a  [Ms] 
           transitDuration:  Duration of the whole transit = ingress + flat part + egress [Ms]
           orbitalPeriod:    Time span between two transits [Ms]
           relativeDepth:    Relative depth of the transit in [0,1]         

    OUTPUT: relativeFlux:    delta F / F0   lightcurve

    """

    phase = np.fmod(time - t0, orbitalPeriod)   # in [0, orbitalPeriod]

    ingressDuration = (transitDuration-flatPartDuration)/2.0
    egressDuration = (transitDuration-flatPartDuration)/2.0
    ingressStart = 0.0
    flatPartStart = ingressDuration
    egressStart = flatPartStart + flatPartDuration
    
    # Start with everything out of transit

    relativeFlux = np.ones_like(time)

    # Ingress part

    ingress = (phase >= ingressStart) & (phase < ingressStart + ingressDuration)
    relativeFlux[ingress] = 1.0 - relativeDepth / ingressDuration * phase[ingress]

    # Flat part

    flat = (phase >= flatPartStart) & (phase < flatPartStart + flatPartDuration)
    relativeFlux[flat] = 1.0-relativeDepth

    # Egress part

    egress = (phase >= egressStart) & (phase < egressStart + egressDuration)
    relativeFlux[egress] = 1.0-relativeDepth + relativeDepth / egressDuration * (phase[egress] - egressStart)

    # That's it!

    return phase, relativeFlux
