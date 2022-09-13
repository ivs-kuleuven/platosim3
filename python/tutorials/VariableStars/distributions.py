#!/usr/bin/env python3

"""
This script..
"""

import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker

#==============================================================#
#                    EXOPLANET DISTRIBUTION                    #
#==============================================================#

# Load NASA exoplanet file
# https://exoplanetarchive.ipac.caltech.edu/

filename = os.getcwd() + '/NASA_archive/nasa_exoplanets.csv'
planetID = np.loadtxt(filename, delimiter=',', usecols=[0], dtype=str)
data     = np.genfromtxt(filename, delimiter=',')

# Sort out planets that do not contain {P, a, R, i}={2, 6, 10, 22}:

cols = [2, 6, 10, 22]
dex0 = np.ones_like(data[:,0], dtype=bool)

for row in range(len(data)):
    for col in cols:
        dex0[row] *= ~np.isnan(data[row,col])

# Check if either {a, Mp*sini}={6, 14}

cols = [6, 14]
dex1 = np.zeros_like(data[:,0], dtype=bool)
for row in range(len(data)):
    for col in cols:
        dex1[row] += ~np.isnan(data[row,col])

dex = dex0 * dex1

# Ommit targets tabulated several times

d = {'ID': planetID[dex],  'P': data[:,2][dex], 'a': data[:,6][dex],
     'R': data[:,10][dex], 'M': data[:,14][dex],
     'e': data[:,18][dex], 'i': data[:,22][dex]}

# Create a data frame for easy handling

df = pd.DataFrame(d, columns=['ID', 'R', 'M', 'P', 'a', 'i', 'e'])

# Fetch columns of interest

df = df.drop_duplicates(subset=['ID'])

# Convert upper mass limit to actual mass

df['M'] /= np.sin(np.deg2rad(df['i']))

# Save data frame into ascii

#df.to_csv('cat.txt')

# PLOT DISTRIBUTION OF PLANETS

fig, ax = plt.subplots(1, 3, figsize=(12,4))

sc = ax[0].scatter(df['M'], df['R'], c=df['e'], edgecolor='w', cmap='coolwarm')
ax[0].set_xscale('log')
ax[0].set_title('Radius vs. Mass')
ax[0].set_xlabel(r'Mass [$M_{\oplus}$]')
ax[0].set_ylabel(r'Radius [$R_{\oplus}$]')

ax[1].scatter(df['P'], df['R'], c=df['e'], edgecolor='w', cmap='coolwarm')
ax[1].set_xscale('log')
ax[1].set_yscale('log')
ax[1].set_title('Radius vs. Period')
ax[1].set_xlabel('Period [days]')
ax[1].set_ylabel('Radius [$R_{\oplus}$]')

ax[2].scatter(df['a'], df['M'], c=df['e'], edgecolor='w', cmap='coolwarm')
ax[2].set_xscale('log')
ax[2].set_yscale('log')
ax[2].set_title('Mass vs. Semimajor axis')
ax[2].set_xlabel('Semimajor axis [AU]')
ax[2].set_ylabel(r'Mass [$M_{\oplus}$]')

ax[0].xaxis.set_major_formatter(mticker.ScalarFormatter())
ax[1].xaxis.set_major_formatter(mticker.ScalarFormatter())
ax[1].yaxis.set_major_formatter(mticker.ScalarFormatter())
ax[2].xaxis.set_major_formatter(mticker.ScalarFormatter())
ax[2].yaxis.set_major_formatter(mticker.ScalarFormatter())

plt.tight_layout(h_pad=1.0)
cbar_ax = fig.add_axes([0.07, 0.58, 0.015, 0.3])
cbar = fig.colorbar(sc, cax=cbar_ax)
cbar.set_label('Eccentricity')
plt.show()
