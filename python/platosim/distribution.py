import random
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit


class Distribution(object):

    """
    This class provides the Python interface to different distribution models.
    """


    def __init__(self, filename):
        """
        PURPOSE:

        INPUT: 
        """
        
        self.filename = filename




    def __del__(self):
        """
        PURPOSE: Destructor.
        """

        pass





    def readMagnitudes(self):
        """
        Reads the magnitudes from the star catalogues with the given filename.  Assumed is
        that this file comprises three (space-separated) columns:
            (1) right ascension of the stars, expressed in degrees
            (2) declination of the stars, expressed in degrees
            (3) magnitude of the stars

        INPUT: filename: Filename (full path) of the star catalogue from which to read the magnitude
        """

        # Open the star catalogue

        catalog = open(self.filename, 'r')

        # Used to store the magnitude values from the catalogue

        magnitudes = []

        # Read the file line-by-line and store the 3rd column in the array of magnitudes

        for star in catalog:

            magnitude = star.split()[2]     # 3rd column
            magnitudes.append(magnitude)

        magnitudes = np.asarray(magnitudes, dtype = np.double)  # Conversion to numpy array

        return magnitudes





    def expFunction(self, x, a, b, c):

        """
        Evaluation of the function a * EXP(-b * x) + c, in variable x.

        INPUT: x: Variable for which the evaluate the exponential function.
        INPUT: a: Scale factor for the exponential component.
        INPUT: b: Scale factor for the variable in the argument of the exponential function.
        INPUT: c: Constant to add to the exponential function.

        OUTPUT: Result of evaluating the function y = a * EXP(-b * x) + c, for x.
        """

        return a * np.exp(-b * x) + c





    def inverseExpFunction(self, y, a, b, c):

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







    def getRandomPosition(self, numRows, numColumns):

        """
        Generates random position (row, column) in the sub-field with the given dimensions, 
        using a uniform distribution for row and column.  This does not take the position of the
        sub-field on the CCD into account.

        INPUT: numRows: Number of rows in the sub-field, expressed in pixels.
        INPUT: numColumns: Number of columns in the sub-field, expressed in pixels.

        OUTPUT: Randomly generated position (row, column) in the sub-field.  Note that these shoud NOT
                be integer values (i.e. stars don't necessarily fall in the centre of a pixel).
        """

        randomRow    = random.uniform(0, numRows    - 1)  # Not necessarily an integer value
        randomColumn = random.uniform(0, numColumns - 1)  # Not necessarily an integer value

        return randomRow, randomColumn


    




    def getRandomMagnitude(self, minMagnitude, maxMagnitude, a, b, c):
        """
        Generates random magnitude in the given magnitude range, using the (fit to the)
        distribution of the magnitudes in a real catalogue (a * EXP(-b * magnitude) + c).

        INPUT: minMagnitude: Lowest magnitude (i.e. of the brightest star).
        INPUT: maxMagnitude: Highest magnitude (i.e. of the faintest star).
        INPUT: a: Scale factor for the exponential component.
        INPUT: b: Scale factor for the variable in the argument of the exponential function.
        INPUT: c: Constant to add to the exponential function.
        """

        # Inverse transform sampling method

        y = 1e3 + 1e4*random.random()
        magnitude = self.inverseExpFunction(y, a, b, c)

        # Make sure the generated magnitude is withing the allowed range
        # (keep on generating new values until it is)
        
        while (magnitude < minMagnitude) or (magnitude > maxMagnitude):     
            magnitude = self.getRandomMagnitude(minMagnitude, maxMagnitude, a, b, c)

        return magnitude

    


    

    

    def getMagnitudeDistribution(self, minMagnitude=2, maxMagnitude=15, plot=False):
        """
        Extracting the magnitude values from a star catalogue (containing right ascension, 
        declination, and magnitude), plotting a histogram, fitting a function
        a * EXP(-b * magnitude) + c to it, and showing the resulting fit.

        - readMagnitudes: Reading the magnitude values from a catalogue
        - expFunction: Function f(x) = a * EXP(-b * x) + c
        - inverseExpFunction: Inverse of function f(x) = a * EXP(-b * x + c)
        - getMagnitudeDistribution: Plotting the magnitude histogram and the fit.
        Reads the magnitudes from a given star catalogue, plots a histogram and fits the
        magnitude distribution with a function a * EXP(-b * magnitude) + c.

        INPUT: filename: Filename (full path) of the star catalogue from which to read the magnitudes.

        OUTPUT: Optimised parameters (a, b, c) of fitting function a * EXP(-b * magnitude) + c 
                to the magnitude distribution in the given star catalogue.
        """

        # Read the magnitudes from a real catalogue

        magnitudes = self.readMagnitudes()

        # Make a histogram
        
        hist, bins, patches = plt.hist(magnitudes, bins = 100)

        # Fit a function a * EXP(-b * magnitude) + c to the magnitude distribution

        x = (bins[1:] + bins[:-1]) / 2  # Centre of the bins
        magnitudeFit, covariance = curve_fit(self.expFunction, x, hist, p0 = (1, 1e-6, 1))
        
        # Plot if requested
        
        if plot:
            xx = np.linspace(minMagnitude, maxMagnitude, 100)
            yy = self.expFunction(xx, *magnitudeFit)
            plt.plot(xx, yy, 'r')
            plt.xlabel("Magnitude")
            plt.ylabel("Number of stars")
            plt.show()

        return magnitudeFit






    def getStarCatalog(self, numStars, minMagnitude, maxMagnitude, subFieldDimensions, plot=False):

        # Fetch the SPF starcatalog used by default
        a, b, c = self.getMagnitudeDistribution(minMagnitude, maxMagnitude, plot=plot)

        rows = []
        columns = []
        magnitudes = []

        for starIndex in range(numStars - 1):

            # Random position in the sub-field (uniform distribution)
            randomRow, randomColumn = self.getRandomPosition(subFieldDimensions, subFieldDimensions)

            # Random magnitude between choices
            randomMagnitude = self.getRandomMagnitude(minMagnitude, maxMagnitude, a, b, c)

            rows.append(randomRow)
            columns.append(randomColumn)
            magnitudes.append(randomMagnitude)

        # Add star which shows blooming

        # Random position in the sub-field (uniform distribution)
        # bloomingRow, bloomingColumn = 25, 25
        # #bloomingRow, bloomingColumn = getRandomPosition(subFieldDimensions, subFieldDimensions)
        # rows.append(bloomingRow)
        # columns.append(bloomingColumn)
        # magnitudes.append(7.0)

        return np.array(rows), np.array(columns), np.array(magnitudes)

    
