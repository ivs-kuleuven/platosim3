#!/usr/bin/env python

import os
import scipy
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm

#from platosim.utilities import powerDensityFFT




# Function to for axes

def axes_yminmax(y):
    ymin = np.min(y) - (np.max(y)-np.min(y))*pt
    ymax = np.max(y) + (np.max(y)-np.min(y))*pt
    return ymin, ymax

#-------------------




