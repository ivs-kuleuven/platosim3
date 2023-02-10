#!/usr/bin/env python3

"""
"""

import os
import datetime
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt 
from pathlib import Path

import platosim.utilities       as ut
import platosim.referenceFrames as rf
from platosim.matplotlibrc import latex
latex()

day2sec = 86400.

#==============================================================#
#                          FUNCTIONS                           #
#==============================================================#


def getPRE(ra, dec, kappa, quarter, sigma=3, outfile=False, show_table=False):
    """
    TODO under development!

    PURPOSE: 
             
    INPUT:

    OUTPUT:
    """

    # Sort input quarters
    n = len(quarter)
    
    # Coordinates
    ICRS = np.array([ra, dec, kappa])
    
    # Pointing Reproducibility Error (PRE) in P/L reference frame (yaw, pitch, roll)
    # Here t stands for transverse direction and [deg]
    t = 3.0/3600
    b = 6.0/3600

    # Find distribution within 3 sigma of req.
    tt = np.array([np.random.normal(0, t/sigma) for i in range(n)])
    bb = np.array([np.random.normal(0, b/sigma) for i in range(n)])

    # Corresponding yaw, pitch, roll
    y = tt
    z = 3 * y
    x = bb - z

    # ICRS pointing angles
    phi   = np.deg2rad(ra)
    theta = np.deg2rad(dec)

    # Find change to pointing for quarters
    PRE = np.zeros((n , 4))
    for i in range(n):
        data = rf.changeOfPointing(x[i], y[i], z[i], phi, theta)[0]
        PRE[i,:] = np.append(quarter[i], data)

    df = pd.DataFrame(PRE, columns=["quarter", "yaw", "pitch", "roll"])
    df = df.astype({"quarter":int, "yaw":np.float64, "pitch":np.float64, "roll":np.float64})
    
    # Print generated values
    if show_table:         

        print('\nChange of coordinates [arcsec]')
        df0 = df.copy()
        df0.iloc[:,1:] = df0.iloc[:,1:] * 3600
        print(df0)

        print('\nNew ICRS coordinates [deg]')
        df1 = df.copy()
        df1.rename(columns={"yaw":"RA", "pitch":"Dec", "roll":"kappa"}, inplace=True)
        df1.iloc[:,1] = df1.iloc[:,1] + ra
        df1.iloc[:,2] = df1.iloc[:,2] + dec
        print(df1)

    # Save file with relative pointing errors [deg]
    if outfile:
        df.to_csv(outfile, sep=" ", header=False, index=False)

    return PRE





def getAPE(ra, dec, kappa, sigma=3, outfile=False, show_table=False):

    # Pointing Reproducibility Error (PRE) in P/L reference frame (yaw, pitch, roll)
    t = 4.5/60  # [deg]
    b = 9.0/60  # [deg]

    # Find distribution within 3 sigma of req.
    tt = np.array([np.random.normal(0, t/sigma) for i in range(24)])
    bb = np.array([np.random.normal(0, b/sigma) for i in range(24)])

    # Corresponding yaw, pitch, roll
    dy = tt
    dz = 3 * dy
    dx = bb - dz

    # Mean and standard deviation
    mu, sigma = 0, sigma
    s = np.random.normal(mu, sigma, 1000)

    # Store APE
    APE = np.transpose([tt, bb])
    df  = pd.DataFrame(APE, columns=["tilt", "azimuth"])

    # Plot histogram and data
    if show_table:
        # Print generate values
        print('\nCamera alignment errors for all 24 N-CAMs [arcmin]')
        APE0 = np.transpose([tt, bb, dx, dy, dz]) * 60
        df0  = pd.DataFrame(APE0, columns=["Alt", "Az", "Yaw", "Pitch", "Roll"])
        print(df0)
        
        # Plot
        #count, bins, ignored = plt.hist(s, 30, density=True)
        # plt.plot(bins, 1/(sigma * np.sqrt(2 * np.pi)) * 
        #          np.exp( - (bins - mu)**2 / (2 * sigma**2) ),
        #          linewidth=2, color='r')
        # plt.show()

    # Save APE camera misalignments
    if outfile:
        df.to_csv(outfile, sep=" ", header=False, index=False)
                       
    return APE





def getTED(quarter, model="poly", outfile=False, plot=False):
        """
        Function to create a Themo-Elastic Distortion (TED) file.
        """
        
        # Constants
        time0 = np.arange(0, 90*day2sec, 25)
        cols  = ["yaw", "pitch", "roll"]
        N     = len(quarter)

        # Create data frame and store default time0 for fit
        df  = pd.DataFrame()
        df1 = pd.DataFrame()
            
        # Loop over each quarter
        
        for Q in range(quarter[0]-1, quarter[-1]):

            # Time column
            t0 = round(90. *  Q    * day2sec)
            t1 = round(90. * (Q+1) * day2sec)
            df1["time"] = np.arange(t0, t1, 25)

            # Generate linear model
            
            # Generate a random 2nd order polynomial
            
            for col in cols:

                if model == 'linear':
                    a = 1.3 * 15      
                    if col == "roll":
                        df1[col] = np.zeros(len(df1.time))
                    else:
                        df1[col] = np.linspace(0, a, len(df1.time))

                else:
                    # NOTE these parameters has been compared to Prime TED
                    a = np.random.uniform(-10, 10) * 1e-14
                    b = np.random.uniform(-15, 15) * 1e-7
                    # Secure that c (the y offset) is always zero
                    c = 0
                    # Make sure that a and b always has opposite signs
                    if np.sign(a) == np.sign(b): b *= -1
                    # Get model fit 
                    poly = np.array([a, b, c])
                    df1[col] = np.polyval(poly, time0)

            # File to save
            df = pd.concat([df, df1])

        # Plot model
        if plot:
            fig, ax = plt.subplots(3,1,figsize=(9, 6))
            # Plots
            for i, col in zip(range(3), cols):
                ax[i].plot(df["time"]/day2sec, df[col], 'k-')
                ax[i].axhline(y=0, linestyle=':', color='k')
                for k in range(N-1):
                    ax[i].axvline(x=quarter[k]*90, linestyle='--', color='b')
            # Settings
            ax[2].set_xlabel("Time [days]")
            ax[0].set_ylabel("Yaw [arcsec]")
            ax[1].set_ylabel("Pitch [arcsec]")
            ax[2].set_ylabel("Roll [arcsec]")
            for i in range(3):
                ax[i].set_xlim(df.time.min()/day2sec, df.time.max()/day2sec)
            # Layout
            ax[0].set_xticklabels([])
            ax[1].set_xticklabels([])
            plt.tight_layout(h_pad=0.2, w_pad=0)
            plt.show()

        # Save data in one big drift text file for PlatoSim
        if outfile:
            fig.savefig(f"{outfile[:-4]}.png", bbox_inches='tight', dpi=200)
            df.to_csv(outfile, sep=" ", header=False, index=False)
