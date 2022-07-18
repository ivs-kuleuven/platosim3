#!/usr/bin/env python

"""
Placeholder for small utility scripts.
"""

import os
import gzip
import pickle
import logging
import numpy as np
from scipy.integrate import cumtrapz
from colorama import Fore, Style

#==============================================================#
#                           UTILITIES                          #
#==============================================================#

def convertQuarterRange(dQ):
    """
    Small function that takes a string of numbers (here quarters)
    and split it up into readable float values used as real number
    ranges. If a single number is given, an quarter integer is 
    returned.
    """
    quarters = []
    for part in dQ.split(','):
        if '-' in part:
            # If a range in mag is provided
            q1, q2 = part.split('-')
            q1, q2 = int(q1), int(q2)
            quarters.append(q1)
            quarters.append(q2)
        else:
            # If only one mag-value is given select 1 mag around it
            q1 = int(part)
            quarters.append(q1)
    return quarters

        
def find_nearest(array, value):
    a = (np.abs(array - value))
    index = np.argmin(a)
    return index


def diff(new, old):
    result = (new - old) / old
    return result


def superLorentzian(nu, b, sigma):
    xi = 2. * np.sqrt(2) / np.pi
    return (xi * sigma**2. / b) / (1+(nu/b)**4.)


def rebin3(x, xp, fp):
    if np.diff(xp).min() < np.diff(x).min():
        # Binning
        x_cum = xp[1:]
        c =  cumtrapz(fp,xp)
        x_diff =  np.diff(x)
        b = x[:-1] + x_diff/2.
        # Deal with edge points - estimate x diff in the outer directions
        b = np.hstack((x[0] - x_diff[0]/2. , b, x[-1] + x_diff[-1]/2. ))
        c_new = np.interp(b, x_cum, c)
        d = 0.5*(x_diff[:-1] + x_diff[1:])
        # Deal with edge points - estimate x diff in the outer directions
        d = np.hstack((x_diff[0] , d, x_diff[-1]))
        new_f = (c_new[1:] - c_new[:-1] ) / d
    else:
        # Interpolate!
        new_f = np.interp(x, xp, fp, left=0.0, right=0.0)
    return x, new_f


# save and load pickle objects
def save_obj(obj,name,path=''):
    with open(path + name + '.pkl', 'wb') as f:
        pickle.dump(obj, f, pickle.HIGHEST_PROTOCOL)


def load_obj(name,path=''):
    with open(path + name + '.pkl', 'rb') as f:
        return pickle.load(f)
