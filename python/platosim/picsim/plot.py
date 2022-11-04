#!/usr/bin/env python3

import os
import h5py
import numpy as np
import ligo.skymap.plot
from matplotlib import pyplot as plt
from matplotlib.ticker import ScalarFormatter

import astropy.units as u
from astropy.coordinates import SkyCoord

import platosim.referenceFrames as rf
import platosim.utilities       as ut
from platosim.matplotlibrc import setup
setup()

#==============================================================#
#                      GRAPHICAL FUNCTIONS                     #
#==============================================================#


def plotPlatoFOV(pointingField, raStars, decStars, magStars=None,
                 nCamVis=None, skymap=None, title=None, fs=20):
    """
    Funtion to plot 

    Parameters
    ----------

    Return
    ------
    fig : object
        Axes matplotlib.pyplot handle object to be modified by the user
        Use fig..savefig('<plot.png>', bbox_inches='tight', dpi=200)
    """

    # Select field
    indir = os.getenv('PLATO_PROJECT_HOME') + '/python/platosim/picsim/'
    if pointingField == 'NPF': PF_gal = [65.0, 30.0]
    if pointingField == 'SPF': PF_gal = [253.0, -30.0]

    PF_gal  = SkyCoord(PF_gal[0], PF_gal[1], frame='galactic', unit='deg')  # [deg]
    PF_icrs = PF_gal.icrs  # [deg]
    
    # Load PIC stars for each N-CAM visibility
    PF06 = np.load(indir + f'{pointingField}-NCAM06.npy')
    PF12 = np.load(indir + f'{pointingField}-NCAM12.npy')
    PF18 = np.load(indir + f'{pointingField}-NCAM18.npy')
    PF24 = np.load(indir + f'{pointingField}-NCAM24.npy')
    starPF06 = SkyCoord(PF06[:,0]*u.deg, PF06[:,1]*u.deg, frame='icrs', unit='deg')
    starPF12 = SkyCoord(PF12[:,0]*u.deg, PF12[:,1]*u.deg, frame='icrs', unit='deg')
    starPF18 = SkyCoord(PF18[:,0]*u.deg, PF18[:,1]*u.deg, frame='icrs', unit='deg')
    starPF24 = SkyCoord(PF24[:,0]*u.deg, PF24[:,1]*u.deg, frame='icrs', unit='deg')
    
    # Load brightest stars
    starPF = SkyCoord(raStars*u.deg, decStars*u.deg, frame='icrs', unit='deg')

    # MAKE PLOTS
    
    fig = plt.figure(figsize=(9,9))
    ax = plt.axes(projection='astro zoom', center=PF_icrs, radius='30 deg', rotate='184 deg')

    # Plot PIC1.1.0 stars after N-CAM visibility
    ax.plot(starPF06.ra.deg, starPF06.dec.deg, '.', c='skyblue',
            transform=ax.get_transform('world'), markersize=1, zorder=1)
    ax.plot(starPF12.ra.deg, starPF12.dec.deg, '.', c='deepskyblue',
            transform=ax.get_transform('world'), markersize=1, zorder=2)
    ax.plot(starPF18.ra.deg, starPF18.dec.deg, '.', c='dodgerblue',
            transform=ax.get_transform('world'), markersize=1, zorder=3)
    ax.plot(starPF24.ra.deg, starPF24.dec.deg, '.', c='royalblue',
            transform=ax.get_transform('world'), markersize=1, zorder=4)

    # Plot stars and add legend scaled to the stellar magnitudes
    if magStars is not None and len(magStars) > 0:
        maxMarkerSize = 30
        dm = (max(magStars) - magStars) * maxMarkerSize
        mag_range = np.arange(min(magStars), max(magStars)).astype(int)
        dm_range  = (max(magStars) - mag_range) * maxMarkerSize/10
        mark, color = 'o', 'gold'
        handle = [plt.plot([],[], "o", c='gray', ms=dm_range[i], ls="")[0] for i in range(len(dm_range))]
        ax.legend(handles=handle, labels=mag_range.tolist(), loc='upper right', title=r"P [mag]", fontsize=16, title_fontsize=16)
    else:
        dm, mark, color = 20, '*', 'none'
    # Plot all stars
    scatter = ax.scatter(starPF.ra.deg, starPF.dec.deg, transform=ax.get_transform('world'), 
                         s=dm, marker=mark, c=color, ec='k', lw=1, zorder=5)

    # Plot pointing of each camera group
    #raGroups, decGroups = rf.getCameraGroupCoordinates(PF_icrs.ra.deg, PF_icrs.dec.deg, -8)
    #camPointing = SkyCoord(raGroups*u.deg, decGroups*u.deg, frame='icrs', unit='deg')  
    #ax.plot(camPointing.ra.deg, camPointing.dec.deg, 'r.', transform=ax.get_transform('world'), markersize=10, zorder=6)

    # Plot pointing of PIC1.1.0 and PIC2.0.0
    ax.plot(PF_icrs.ra.deg, PF_icrs.dec.deg, '*', transform=ax.get_transform('world'), ms=20, c='k', mfc='r', zorder=6)
    #ax.plot(277.18, 52.85, '*', transform=ax.get_transform('world'), ms=20, c='k', mfc='b', zorder=7)

    # Add-on's
    ax.scalebar((0.05, 0.05), 10 * u.deg).label()
    ax.compass(0.95, 0.05, 0.1)
    ax.grid(color='gray')

    # Settings
    if title is not None:
        ax.set_title(title, fontsize=fs+2)
    ax.set_xlabel('RA',  fontsize=fs)
    ax.set_ylabel('Dec', fontsize=fs)
    plt.xticks(fontsize=fs)
    plt.yticks(fontsize=fs)
    ax.tick_params(axis='both', labelsize=fs)
    
    # Return figure
    return fig











def plotTeffvsRadius(fig0, starSample, title,
                     ds, df_dK, df_dG, df_dF,
                     sg, df_sgK, df_sgG, df_sgF,
                     df, ms_limit):

    # Fontsize
    ms = 3
    if starSample == 'P5': da = 0.2
    else: da = 0.0
    
    # Plots sub-giants
    plt.plot(df_sgK['Teff'], df_sgK['R'], 'o', alpha=0.5-da, markersize=ms, color='orange',      label=r'K$\,$IV')
    plt.plot(df_sgG['Teff'], df_sgG['R'], 'o', alpha=0.7-da, markersize=ms, color='greenyellow', label=r'G$\,$IV')
    plt.plot(df_sgF['Teff'], df_sgF['R'], 'o', alpha=0.7-da, markersize=ms, color='skyblue',     label=r'F$\,$IV')

    # Plot dwarfs
    plt.plot(df_dK['Teff'], df_dK['R'], 'o', alpha=0.4-da, markersize=ms, color='orangered', label=r'K$\,$V')
    plt.plot(df_dG['Teff'], df_dG['R'], 'o', alpha=0.4-da, markersize=ms, color='limegreen', label=r'G$\,$V')
    plt.plot(df_dF['Teff'], df_dF['R'], 'o', alpha=0.4-da, markersize=ms, color='royalblue', label=r'F$\,$V')

    # Plot selected targets
    plt.plot(df['Teff'], df['R'], 'o', mfc='none', mec='k', alpha=0.5, markersize=ms)

    # Compute main sequence devision
    dt = np.arange(np.min(ds['Teff']), np.max(ds['Teff']), 10)
    plt.plot(dt, ms_limit(dt), 'k-')

    # Settings
    plt.title(title)
    plt.xlabel(r'Effective temperature, $T_{\mathrm{eff}}$ [K]', fontsize=16)
    plt.ylabel(r'Stellar radius, $R$ [$R_{\odot}$]', fontsize=16)

    # Legend
    order = [3, 4, 5, 0, 1, 2]
    handles, labels = plt.gca().get_legend_handles_labels()
    h = [handles[idx] for idx in order]
    l = [labels[idx] for idx in order]
    plt.legend(h, l, ncol=2, loc='upper left', prop={'size':12},
               columnspacing=0.5, handletextpad=0)

    # Finito!
    return

    







    

def plotStellarSampleDistributions(fig, mag, magCon, magRange, numConPerTar, distCon):
    """
    This function plots 4 different stellar sample distribution plots
    for an PLATO Input Catalogue (PIC)
    1) Magnitude distribution of PIC targets
    2) Magnitude distribution of PIC contaminants
    3) Number distribution of contaminants per target
    4) Distance distribution of contaminants

    Parameters
    ----------
    mag : list, array
        The stellar target magnitudes
    magCon : list, array
        The stellar contaminant magnitudes
    magRange : list, array
        Upper and lower magnitude limit for input sample
    numConPerTar : list, array
        The number of contaminants (integer) per target star
    distCon : list, array
        Distances of each contaminant star w.r.t. their target

    Return
    ------
    axes : object
        Axes matplotlib.pyplot handle object to be modified by the user
    """

    # Copy figure object not to overwrite it

    fig1, axes = fig

    # Prepare bins and plot magnitude distribution of targets

    magbinTar  = 0.1
    if magRange[1]-magRange[0] > 10: magbinTar = 0.2
    binsizeTar = int((magRange[1] - magRange[0]) / magbinTar) + 1
    binlistTar = np.linspace(magRange[0], magRange[1], binsizeTar)

    axes[0,0].hist(mag, binlistTar, facecolor='b', edgecolor='b', fill=True, alpha=0.3)
    axes[0,0].set_title('Magnitude distribution of PIC targets')
    axes[0,0].set_xlabel(r'$P$ passband')
    #axes[0,0].set_xlabel(r'$V$ Johnson-Cousin')
    axes[0,0].set_ylabel('Number of stars')
    axes[0,0].locator_params(axis='y', integer=True)
    axes[0,0].tick_params(axis='x', which='minor', bottom=True, top=False)
    axes[0,0].tick_params(axis='x', which='major', bottom=True, top=False)
    axes[0,0].tick_params(axis='y', which='minor', left=False, right=False)
    axes[0,0].tick_params(axis='y', which='major', left=True, right=False)
    axes[0,0].grid(axis='y', color='gray', alpha=0.3)

    # Prepare bins and plot magnitude distribution of contaminants

    magbinCon  = 0.2
    binsizeCon = int((np.max(magCon) - np.min(magCon)) / magbinCon) + 1
    binlistCon = np.linspace(round(np.min(magCon)), round(np.max(magCon)), binsizeCon)

    axes[0,1].hist(magCon, binlistCon, facecolor='m', edgecolor='m', fill=True, alpha=0.3)
    axes[0,1].set_title('Magnitude distribution of PIC contaminants')
    #axes[0,1].set_xlabel(r'$V$ Johnson-Cousin')
    axes[0,1].set_xlabel(r'$P$ passband')
    axes[0,1].set_ylabel('Number of stars')
    axes[0,1].tick_params(axis='x', which='minor', bottom=True, top=False)
    axes[0,1].tick_params(axis='x', which='major', bottom=True, top=False)
    axes[0,1].tick_params(axis='y', which='minor', left=False, right=False)
    axes[0,1].tick_params(axis='y', which='major', left=True, right=False)
    axes[0,1].grid(axis='y', color='gray', alpha=0.3)

    # Prepare bins and plot number distribution of contaminants per target

    numbinCon  = 1 + int(np.max(numConPerTar)/50)
    binsizeNum = int((np.max(numConPerTar) - 0) / numbinCon) + 2
    binlistNum = np.linspace(-0.5, np.max(numConPerTar)+0.5, binsizeNum)  # -0.5 because num x-axis

    axes[1,0].hist(numConPerTar, binlistNum, facecolor='g', edgecolor='g', fill=True, log=True, alpha=0.3)
    axes[1,0].yaxis.set_major_formatter(ScalarFormatter())
    axes[1,0].set_title('Number distribution of contaminants per target')
    axes[1,0].set_xlabel('Number of contaminants')
    axes[1,0].set_ylabel('Number of targets')
    axes[1,0].tick_params(axis='x', which='minor', bottom=False, top=False)
    axes[1,0].tick_params(axis='x', which='major', bottom=True, top=False)
    axes[1,0].tick_params(axis='y', which='minor', left=False, right=False)
    axes[1,0].tick_params(axis='y', which='major', left=True, right=False)
    axes[1,0].grid(axis='y', color='gray', alpha=0.3)

    # Prepare bins and plot distance distribution of contaminants in respect to their target star

    distbinCon  = 1.0
    binsizeDist = int((np.max(distCon) - np.min(distCon)) / distbinCon) + 2  # +1 extra because zero is rare
    binlistDist = np.linspace(round(np.min(distCon)), round(np.max(distCon)), binsizeDist)

    axes[1,1].hist(distCon, binlistDist, facecolor='orange', edgecolor='orange', fill=True, alpha=0.4)
    axes[1,1].set_title('Distance distribution of contaminants')
    axes[1,1].set_xlabel('Distances [arcsec]')
    axes[1,1].set_ylabel('Number of stars')
    axes[1,1].locator_params(axis='y', integer=True)
    axes[1,1].tick_params(axis='x', which='minor', bottom=True, top=False)
    axes[1,1].tick_params(axis='x', which='major', bottom=True, top=False)
    axes[1,1].tick_params(axis='y', which='minor', left=False, right=False)
    axes[1,1].tick_params(axis='y', which='major', left=True, right=False)
    axes[1,1].grid(axis='y', color='gray', alpha=0.3)

    # Layout

    plt.tight_layout()

    # Finito!

    return axes


