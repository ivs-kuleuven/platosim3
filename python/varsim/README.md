# README: Stellar Oscillation & Plnaetary Eclipse Simulator

## Overview


Model Convection-Driven Variability in the passbands

In this notebook we use asteroseismically derived scaling relations to
simulate the granulation and pulsation signal of stars to-be observed
with PLATO (PLAnetary Transits and Oscillations of stars) mission. This
notebook was developed for the ARIEL (Atmospheric Remote-sensing Exoplanet
Large-survey) mission.

The granulation, pulsations and bolmetric coefficient comes from this [paper](https://academic.oup.com/mnras/article/481/3/2871/5092616?login=true)

## Input parameters

List of spectral types or star selections
  Spectral type  = ["M5V"]
  Selected stars = ["Sun", "GJ1214", "HD209458"]

List of scaling relations for the convective turbulence (granulation, pulsations)
and active regions. A consice review for the pulsational scaling relations,
developed from 1995 to 2013, is presented in: Corsaro et al. 2013 "A Bayesian
approach to scaling relations for amplitudes of solar-like oscillations in Kepler stars"

  ScalingRelation_gran = ['None',      'KjeldsenBedding2011', 'Kallinger2014']
  ScalingRelation_puls = ['None',      'KjeldsenBedding1995', 'KB1995Brown1991', 'Mosser2010',
                          'Huber2011', 'KjeldsenBedding2011', 'Corsaro2013', 'Kallinger2014']
  ScalingMethod_actr   = ['None', 'Sun']
  Correction           = ['None', 'Kepler Bandpass']


    - In case solarlike oscillations and granulation was included:
        - Plot the noiseless input time series of the oscgran time series
        - Plot the PSD of the noiseless oscgran time series
        - Plot an echelle diagram of the frequencies (like Fig 10 top panel in article of Reza)
        - Plot the damping time as a function of the frequency,
        - Plot the amplitudes of the oscillation modes as a function of frequency
    - For other (non-stochastic) variabiliy types (e.g. binarity)
        - Plot the noiseless time series (or a zoom to see better what's going on)
        - Plot a power spectrum
