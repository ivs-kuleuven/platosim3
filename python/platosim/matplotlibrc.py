#!/usr/bin/env python3

"""
This script is used to configure the matplotlibrc file in order to 
consistantly use the same settings for all plots.
"""

# Built-in
import os

# PlatoSim standard
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
    plt.rcParams['legend.fancybox']   = True
    plt.rcParams['legend.framealpha'] = 0.8
    plt.rcParams['legend.fontsize']   = 15

    # Font
    plt.rcParams['font.family'] = 'serif'
    plt.rcParams['font.size']   = 17

    
        
def latex():
    setup()
    plt.rcParams['text.usetex'] = True

    
    
def setup_notebook():
    setup()
    latex()

    
    
def setup_paper():
    setup()
    latex()
    plt.rcParams['font.size']       = 20
    plt.rcParams['legend.fontsize'] = 17
