#!/usr/bin/env python3

"""
This script is used to configure the matplotlibrc file in order to consistantly
use the same settings for all plots.
"""

import matplotlib.pyplot as plt


def setup():

    # Control tickers
    plt.rcParams['xtick.top']           = True
    plt.rcParams['xtick.bottom']        = True
    plt.rcParams['xtick.top']           = True
    plt.rcParams['xtick.labelbottom']   = True
    plt.rcParams['xtick.direction']     = 'out'
    plt.rcParams['xtick.minor.visible'] = True
    plt.rcParams['xtick.major.top']     = True
    plt.rcParams['xtick.minor.top']     = True
    plt.rcParams['xtick.minor.bottom']  = True
    plt.rcParams['xtick.alignment']     = 'center'
    plt.rcParams['ytick.left']          = True
    plt.rcParams['ytick.right']         = True
    plt.rcParams['ytick.labelleft']     = True
    plt.rcParams['ytick.minor.visible'] = True
    plt.rcParams['ytick.major.left']    = True
    plt.rcParams['ytick.major.right']   = True
    plt.rcParams['ytick.minor.left']    = True
    plt.rcParams['ytick.minor.right']   = True

    # Legends
    plt.rcParams['legend.loc']        = 'best'
    plt.rcParams['legend.frameon']    = True
    plt.rcParams['legend.framealpha'] = 0.8
    plt.rcParams['legend.fancybox']   = True

    # Font
    plt.rcParams['font.family'] = 'serif'
    plt.rcParams['font.size']   = 17

    
    
    
def latex():

    # Set Latex font
    plt.rcParams['text.usetex'] = True
    setup()
    

    
    
def setup_notebook():

    # Set default
    plt.rcParams['font.size']   = 20
    setup()
    latex()
