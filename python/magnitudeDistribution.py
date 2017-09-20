"""
Extracting the magnitude values from a star catalogue (containing right ascension, 
declination, and magnitude), plotting a histogram, fitting a function
a * EXP(-b * magnitude) + c to it, and showing the resulting fit.

- readMagnitudes: Reading the magnitude values from a catalogue
- expFunction: Function f(x) = a * EXP(-b * x) + c
- inverseExpFunction: Inverse of function f(x) = a * EXP(-b * x + c)
- getMagnitudeDistribution: Plotting the magnitude histogram and the fit.
"""





#########
# Imports
#########

import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit





####################
# Auxialiary methods
####################

def readMagnitudes(filename):
    
    """
    Reads the magnitudes from the star catalogues with the given filename.  Assumed is
    that this file comprises three (space-separated) columns:
        (1) right ascension of the stars, expressed in degrees
        (2) declination of the stars, expressed in degrees
        (3) magnitude of the stars
    
    INPUT: filename: Filename (full path) of the star catalogue from which to read the magnitude.
    """

    # Open the star catalogue

    catalog = open(filename, 'r')
    
    # Used to store the magnitude values from the catalogue
    
    magnitudes = []
    
    # Read the file line-by-line and store the 3rd column in the array of magnitudes
    
    for star in catalog:
        
        magnitude = star.split()[2]     # 3rd column
        magnitudes.append(magnitude)
        
    magnitudes = np.asarray(magnitudes, dtype = np.double)  # Conversion to numpy array
    
    return magnitudes
    
    
    


def expFunction(x, a, b, c):
    
    """
    Evaluation of the function a * EXP(-b * x) + c, in variable x.
    
    INPUT: x: Variable for which the evaluate the exponential function.
    INPUT: a: Scale factor for the exponential component.
    INPUT: b: Scale factor for the variable in the argument of the exponential function.
    INPUT: c: Constant to add to the exponential function.
    
    OUTPUT: Result of evaluating the function y = a * EXP(-b * x) + c, for x.
    """
    
    return a * np.exp(-b * x) + c




    
def inverseExpFunction(y, a, b, c):
    
    """
    Evaluation of the inverse of the function y = a * EXP(-b * x) + c, in variable x,
    which is x = -LN((y - c) / a) / b.
    
    INPUT: y: Result of the evaluation of the function, for variable x.
    INPUT: x: Variable for which the evaluate the exponential function.
    INPUT: a: Scale factor for the exponential component.
    INPUT: b: Scale factor for the variable in the argument of the exponential function.
    INPUT: c: Constant to add to the exponential function.
    
    OUTPUT: Result of evaluating the inverse of function y = a * EXP(-b * x) + c, for y.
    """
    
    return -np.log((y - c) / a) / b
    
        



def getMagnitudeDistribution(filename):
    
    """
    Reads the magnitudes from a given star catalogue, plots a histogram and fits the
    magnitude distribution with a function a * EXP(-b * magnitude) + c.

    INPUT: filename: Filename (full path) of the star catalogue from which to read the magnitudes.
    
    OUTPUT: Optimised parameters (a, b, c) of fitting function a * EXP(-b * magnitude) + c 
            to the magnitude distribution in the given star catalogue.
    """
    
    # Read the magnitudes from a real catalogue
    
    magnitudes = readMagnitudes(filename)

    # Make a histogram

    hist, bins, patches = plt.hist(magnitudes, bins = 100, normed = 'True')

    # Fit a function a * EXP(-b * magnitude) + c to the magnitude distribution
    
    x = (bins[1:] + bins[:-1]) / 2                                                  # Centre of the bins
    magnitudeFit, covariance = curve_fit(expFunction, x, hist, p0 = (1, 1e-6, 1))

    xx = np.linspace(2, 15, 100)
    yy = expFunction(xx, *magnitudeFit)
    plt.plot(xx, yy, 'r')

    plt.xlabel("Magnitude")
    plt.show()
    
    return magnitudeFit
    