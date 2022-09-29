#!/usr/bin/env python

import os
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm
from scipy.interpolate import make_interp_spline
from scipy.ndimage import median_filter
from PyAstronomy import pyasl
#from platosim.utilities import powerDensityFFT

# Hard-code values

fs = 14   # Font-size
ms = 0.2  # Marker-size
lw = 0.5  # Line-width
pt = 0.1  # Percentage 
pp = 0.01

# Define some nice colors

colors_sea = ['royalblue', 'lightseagreen', 'limegreen']
colors_hot = ['tomato', 'darkorange', 'gold']
colors_new = ['royalblue', 'limegreen', 'darkorange', 'tomato', 'gold']

# Function to for axes

def axes_yminmax(y):
    ymin = np.min(y) - (np.max(y)-np.min(y))*pt
    ymax = np.max(y) + (np.max(y)-np.min(y))*pt
    return ymin, ymax

#-------------------




def plot_phoenix_sed(wvl, wvl1_in, wvl2_in, wvl_equi,
                     flux, bb_flux, flux1_in, flux2_in, flux_equi,
                     Teff, Teff_upper, Teff_lower):
    """
    Plot PHOENIXS SED and blackbody model
    """

    fig, ax = plt.subplots(2, 1, figsize=(12,7))

    # Plot PHOENIXS SED and blackbody model
    #ax[0].plot(wvl, flux, c='k', label=r'PHOENIX model $T_{\mathrm{eff}}$'+' = {0} K'.format(int(Teff)), lw=lw)
    #ax[0].plot(wvl*1000, bb_flux, label='Blackbody model')
    #ax[0].set_xlim(0, 20000)
    #ax[0].legend(fontsize=12)

    # Plot the interpolation of the grid
    ax[0].plot(wvl2_in, flux2_in, c='blue', lw=lw, alpha=0.8, label=r'$T_{\mathrm{eff}}$'+' = {0} K'.format(int(Teff_upper)))
    ax[0].plot(wvl,     flux,     c='k', lw=lw, alpha=0.8, label=r'$T_{\mathrm{eff}}$'+' = {0} K'.format(int(Teff)))
    ax[0].plot(wvl1_in, flux1_in, c='green', lw=lw, alpha=0.8, label=r'$T_{\mathrm{eff}}$'+' = {0} K'.format(int(Teff_lower)))
    ax[0].set_xlim(2000, 10000)
    ax[0].legend(fontsize=12)

    # Plot the final equidistant grid used for further calculations
    ax[1].plot(wvl, flux,  'k', lw=lw,  alpha=0.8, label='Zoom-in on original grid')
    ax[1].plot(wvl_equi, flux_equi, 'r', lw=0.8, alpha=1.0, label=r'Equidistant grid: by Subhajit Sarkar')
    ax[1].set_xlim(wvl_equi[0]-500, wvl_equi[-1]+500)
    ax[1].set_xlabel('$\lambda$ [AA]')
    ax[1].legend(fontsize=12)
    plt.tight_layout()

    fig.text(0.001, 0.5, r'Flux [ergs sec$^{-1}$ cm$^{-2}$ AA$^{-1}$ sr$^{-1}$]', va='center', rotation='vertical')
    
    # Finito!
    plt.show()
    #fig.savefig('/home/nicholas/phoenixSED.png', bbox_inches='tight', dpi=300)







    

def plot_amplitude_time_series(time, signal_gran, signal_puls, signal_total, star):
    """
    Plot bolometric luminosity amplitude timeseries.

    """

    # Correct time points from Ms to days

    time = time * 1e6 / 86400.

    # Plot

    fig, (ax1, ax2, ax3) = plt.subplots(3, figsize=(12, 12), sharex=True)
    ax1.plot(time, signal_gran,  colors_hot[0], linewidth=lw, label = 'Granulation')
    ax2.plot(time, signal_puls,  colors_hot[1], linewidth=lw, label = 'Pulsations')
    ax3.plot(time, signal_total, 'k',           linewidth=lw, label = 'Total Aemplitude')

    # Limits

    ax1.set_xlim(0, time[-1])
    ax1.set_ylim(signal_gran.min()  - signal_gran.std(),  signal_gran.max()  + signal_gran.std())
    ax2.set_ylim(signal_puls.min()  - signal_puls.std(),  signal_puls.max()  + signal_puls.std())
    ax3.set_ylim(signal_total.min() - signal_total.std(), signal_total.max() + signal_total.std())

    # Labels

    ax1.set_title('Bolometric luminosity amplitude time series of ' + star, fontsize = fs)
    ax3.set_xlabel('Time [days]',           fontsize = fs-2)
    ax1.set_ylabel('Granulation [ppm]',     fontsize = fs-2)
    ax2.set_ylabel('Pulsation [ppm]',       fontsize = fs-2)
    ax3.set_ylabel('Total Amplitude [ppm]', fontsize = fs-2)

    # Extra settings

    fig.subplots_adjust(hspace=0)
    plt.setp([a.get_xticklabels() for a in fig.axes[:-1]], visible=False)
    plt.tight_layout()

    # Finito!

    plt.show()








def plot_amplitude_spectrum(time, signals, sampling, freqlim=1e-2, title=False, save=False):
    """
    POWER DENSITY SPECTRA

    PARAMETERS
    ----------
    time : narray
        Time points
    datasets : narray, list-narray
        Either single signal array or a list of signal arrays
    title : str (optional)
        Title for plot
    labels : list-str (optinal)
        List of string labels where the first is the xlabel and the rest is ylabels

    OUTPUT:
    Plot or/and saved plot to PNG.
    """

    # Compute frequencies uptil the Nyquist frequency

    medfilt = 144  # [hour for N-Cams]
    Nfreq   = int(len(time)/2.+1)

    PSD  = np.zeros((3, Nfreq))
    med  = np.zeros((3, Nfreq))

    for i in range(3):
        freq, PSD[i,:] = powerDensityFFT(signals[i], sampling)
        med[i,:] = median_filter(PSD[i,:], medfilt)

    # PLOT SEPERATE

    fig, ax = plt.subplots(3,1, figsize=(12,12))

    # Plot subplots

    ax[0].plot(freq, PSD[0], '-', color='gold',   linewidth=lw)
    ax[1].plot(freq, PSD[1], '-', color='tomato', linewidth=lw)
    ax[2].plot(freq, PSD[2], '-', color='gray',   linewidth=lw)

    ax[0].plot(freq, med[0], '-', color='darkorange', linewidth=lw+2)
    ax[1].plot(freq, med[1], '-', color='r',          linewidth=lw+2)
    ax[2].plot(freq, med[2], '-', color='k',          linewidth=lw+2)

    # Limits

    ax[0].set_ylim(PSD[0].min(), PSD[0].max())
    ax[1].set_ylim(PSD[1].min(), PSD[1].max())
    ax[2].set_ylim(PSD[2].min(), PSD[2].max())

    # Common settings

    for plot in range(3):
        ax[plot].set_xlim(100, max(freq)+100)
        ax[plot].set_xscale("log")
        ax[plot].set_yscale("log")

    # Labels

    if title is False: ax[0].set_title('Amplitude spectrum - log scale', fontsize=fs)
    else: ax[0].set_title(title, fontsize=fs)
    ax[2].set_xlabel(r'Frequency [$\mu$Hz] ',  fontsize=fs-2)
    ax[0].set_ylabel(r'Granulation [ppm$^2$ $\mu$Hz$^{-1}$]', fontsize=fs-2)
    ax[1].set_ylabel(r'Pulsation [ppm$^2$ $\mu$Hz$^{-1}$]',   fontsize=fs-2)
    ax[2].set_ylabel(r'Total Power [ppm$^2$ $\mu$Hz$^{-1}$]', fontsize=fs-2)

    # Settings

    plt.setp([a.get_xticklabels() for a in fig.axes[:-1]], visible=False)
    fig.subplots_adjust(hspace=0)
    plt.tight_layout()
    plt.show()

    # PLOT TOGETHER

    # plt.figure(figsize=(12,10))

    # # Plotting

    # plt.plot(freq, PSD[2], '-', color='gray',   linewidth=lw, label='Total amplitude')
    # plt.plot(freq, PSD[0], '-', color='gold',   linewidth=lw, label='Granulation')
    # plt.plot(freq, PSD[1], '-', color='tomato', linewidth=lw, label='Pulsations')
    # plt.plot(freq, med[2], '-', color='k', linewidth=lw+2)

    # # Limits

    # plt.ylim(1e-1, 1e7)
    # plt.xscale("log")
    # plt.yscale("log")

    # # Lables

    # if title is False: ax[0].set_title('Amplitude spectrum - log scale', fontsize=fs)
    # else: ax[0].set_title(title, fontsize=fs)
    # plt.xlabel(r'Frequency [$\mu$Hz]', fontsize=fs-2)
    # plt.ylabel(r'Amplitude [ppm$^2$]', fontsize=fs-2)

    # # Settings

    # plt.tight_layout()
    # plt.show()







def plot_passband_ldc(wvl_int_plato, tran_int_plato, grid_no,
                      mu_trunc, intensity_VTA_trunc, LD_values, ldc):
    """
    Plot... TODO
    """

    # Import TESS passband

    wvl_tess = np.loadtxt(os.getcwd() + '/data/Passbands/response_tess.txt')[:,0]*10.  # [Å]
    tra_tess = np.loadtxt(os.getcwd() + '/data/Passbands/response_tess.txt')[:,1]      # Norm.
    wvl_int_tess  = np.linspace(wvl_tess[0], wvl_tess[-1], grid_no)
    passband_tess = make_interp_spline(wvl_tess, tra_tess, k=3)
    tran_int_tess = passband_tess(wvl_int_tess)

    # Import Kepler passband

    wvl_kepler  = np.loadtxt(os.getcwd() + '/data/Passbands/response_kepler.txt')[:,0]*10.
    tran_kepler = np.loadtxt(os.getcwd() + '/data/Passbands/response_kepler.txt')[:,1]
    wvl_int_kepler  = np.linspace(wvl_kepler[0], wvl_kepler[-1], grid_no)
    passband_kepler = make_interp_spline(wvl_kepler, tran_kepler, k=3)
    tran_int_kepler = passband_kepler(wvl_int_kepler)

    # Create the plot

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 6))#13,6))

    # Response functions:

    ax1.plot(wvl_int_plato,  tran_int_plato,  'r-', label='PLATO')
    ax1.plot(wvl_int_tess,   tran_int_tess,   'g--', label='TESS')
    ax1.plot(wvl_int_kepler, tran_int_kepler, 'b:', label='Kepler')
    ax1.set_xlabel(r'Wavelength, $\lambda$ [Å]')
    ax1.set_ylabel(r'Norm. Spectral Response, $S_{\lambda}$')
    ax1.set_title('Bandpass')
    ax1.legend(fontsize=12)

    # Quadratic LD coeffients fitting:

    ax2.plot(mu_trunc, intensity_VTA_trunc, 'ko', alpha=0.2, label='Data')
    lab = 'Model:'+'\n'+r'$u_1$ = %.5s'%ldc[0]+'\n'+r'$u_2$ = %.5s'%ldc[1]
    ax2.plot(mu_trunc, LD_values.model, 'r-', label=lab)
    ax2.set_xlabel(r'Norm. Wavelength, $\lambda$')
    ax2.set_ylabel(r'Norm. Intensity, $I_{\lambda}$')
    ax2.set_title('Quadratic LD coeffients fitting')
    ax2.legend(fontsize=12)
    plt.tight_layout()
    plt.show()
    #fig.savefig('/home/nicholas/Nextcloud/presentations/presentation_PW12/plotPassbandPLATO.png', bbox_inches='tight', dpi=300)






    
def plot_orbital_phase_curve(fig, time, lc_tra, lc_occ, lc_beam, lc_elli, lc_final,
                             t0, P, dt_c, t_tra_cen, t_tra_tot, t_occ_cen, t_occ_tot,
                             A_beam, A_elli, colors=None):

    """
    Module to plot any given TODO
    """

    # Input parameters
    if colors is None: colors = colors_new

    # Phase-fold light curve
    phase = pyasl.foldAt(time, P, T0=t0)

    # Locate first transit
    t_tra_min = t_tra_cen - t_tra_tot
    t_tra_max = t_tra_cen + t_tra_tot
    dex_tra = (time >= t_tra_min) * (time < t_tra_max)

    # Locate first occultation
    t_occ_min = t_occ_cen - t_occ_tot
    t_occ_max = t_occ_cen + t_occ_tot
    dex_occ = (time >= t_occ_min) * (time < t_occ_max)

    # Setup
    plt.subplots_adjust(wspace=0.15, hspace=0.20)
    fig.text(0.5, 0.92, 'Time [days]', ha='center', fontsize=fs)
    fig.text(0.5, 0.06, 'Phase', ha='center', fontsize=fs)
    fig.text(0.05, 0.4, 'Relative Flux [ppm]', va='center', rotation='vertical', fontsize=fs)

    # Final light curve in time

    ax0 = fig.add_subplot(4,2,(1,2))
    # Plot
    ax0.axvline(t0,      color='gray', linestyle='--')
    ax0.axvline(t0+dt_c, color='gray', linestyle=':')
    ax0.plot(time, lc_final/1e6+1, 'k-')
    # Axes
    ax0.xaxis.set_label_position('top')
    ax0.xaxis.tick_top()
    ymin, ymax = axes_yminmax(lc_final/1e6+1)
    ax0.set_ylim(ymin, ymax)
    ax0.set_xlim(time[0], time[-1])
    # Color fill areas of interest
    ax0.fill_between((t_tra_min, t_tra_max), ymin, ymax, facecolor=colors[0], alpha=0.3)
    ax0.fill_between((t_occ_min, t_occ_max), ymin, ymax, facecolor=colors[1], alpha=0.2)
    ax0.ticklabel_format(style='plain', useOffset=False)
    # Labels
    ax0.set_ylabel('Relative Flux')

    # Transit

    ax1 = fig.add_subplot(4,2,3)
    # Plot
    ax1.plot(time[dex_tra], lc_tra[dex_tra], '-', c=colors[0], label='Transit')
    ax1.legend(loc='upper center')
    # Text
    x_pos     = t_tra_max - t_tra_tot*1.3
    delta_tra = np.max(lc_tra[dex_tra]) - np.min(lc_tra[dex_tra])
    y_pos     = np.max(lc_tra[dex_tra]) - delta_tra/2.
    ax1.text(x_pos, y_pos, r'$\delta_{\mathrm{tra}}=%.1f$ ppm' % delta_tra, fontsize=fs-4)
    # Axes
    ax1.set_xlim(time[dex_tra][0], time[dex_tra][-1])
    ax1.set_ylim(axes_yminmax(lc_tra))
    ax1.xaxis.set_label_position('top')
    ax1.xaxis.tick_top()

    # Occultation and phase curve

    ax2 = fig.add_subplot(4,2,4)
    # Plot
    ax2.plot(time[dex_occ], lc_occ[dex_occ], '-', c=colors[1], label='Occultation')
    ax2.legend(loc='upper center')
    # Text
    x_pos     = t_occ_max - t_occ_tot*1.3
    delta_occ = np.max(lc_occ[dex_occ]) - np.min(lc_occ[dex_occ])
    y_pos     = np.max(lc_occ[dex_occ]) - delta_occ/2.
    ax2.text(x_pos, y_pos, r'$\delta_{\mathrm{occ}}=%.1f$ ppm' % delta_occ, fontsize=fs-4)
    # Axes
    ax2.set_xlim(time[dex_occ][0], time[dex_occ][-1])
    ax2.set_ylim(axes_yminmax(lc_occ))
    ax2.xaxis.set_label_position('top')
    ax2.xaxis.tick_top()

    # Beaming and Ellipsoidal

    sort = np.argsort(phase)

    ax3 = fig.add_subplot(4,2,(5,6))
    # Lines
    ax3.axhline(0.00, color='gray', linestyle='--')
    ax3.axvline(0.00, color='gray', linestyle='--')
    ax3.axvline(0.25, color='gray', linestyle=':')
    ax3.axvline(0.50, color='gray', linestyle='-.')
    ax3.axvline(0.75, color='gray', linestyle=':')
    ax3.axvline(1.00, color='gray', linestyle='--')
    # Plots
    ax3.plot(phase[sort], lc_beam[sort], '-', c=colors[2], ms=ms, label='Beaming')
    ax3.plot(phase[sort], lc_elli[sort], '-', c=colors[3], ms=ms, label='Ellipsoidal')
    ypos_max   = np.max([lc_beam, lc_elli])
    ypos_text0 = ypos_max - ypos_max*2.0*pt
    ypos_text1 = ypos_max - ypos_max*5.0*pt
    ax3.text(0.347, ypos_text0, r'$A_{\mathrm{beam}}=%.3f$ ppm' % A_beam, fontsize=fs-4)
    ax3.text(0.360, ypos_text1, r'$A_{\mathrm{elli}}=%.3f$ ppm' % A_elli, fontsize=fs-4)
    ax3.legend(loc='upper left')
    # Axes
    ax3.xaxis.set_major_formatter(plt.NullFormatter())
    ax3.set_xlim(0-pp, 1+pp)

    # Combined model

    ax4 = fig.add_subplot(4,2,(7,8))
    # Lines
    ax4.axvline(0.00, color='gray', linestyle='--', zorder=0)
    ax4.axvline(0.25, color='gray', linestyle=':',  zorder=1)
    ax4.axvline(0.50, color='gray', linestyle='-.', zorder=2)
    ax4.axvline(0.75, color='gray', linestyle=':',  zorder=3)
    ax4.axvline(1.00, color='gray', linestyle='--', zorder=4)
    # Text labels
    ymin, ymax = axes_yminmax(lc_final - lc_tra)
    ydif = (ymax-ymin)*pp
    ypos_text = ymax + ydif + ymax*pt
    ax4.text(0.00-0.02, ypos_text, 'Transit',     fontsize=fs-5)
    ax4.text(1.00-0.02, ypos_text, 'Transit',     fontsize=fs-5)
    ax4.text(0.50-0.03, ypos_text, 'Occultation', fontsize=fs-5)
    ax4.text(0.25-0.03, ypos_text, 'Quadrature',  fontsize=fs-5)
    ax4.text(0.75-0.03, ypos_text, 'Quadrature',  fontsize=fs-5)
    # Plot
    yy = lc_final/(np.max(lc_final)+np.max(lc_final*2*pt))
    ax4.plot(phase[sort], lc_final[sort], 'k-', zorder=5, label='Combined Model')
    ax4.scatter(phase, lc_final, marker='o', s=5, c=cm.hot(yy), ec='None', zorder=6)
    # Axes
    ax4.set_ylim(ymin-ydif, ymax+ydif)
    ax4.set_xlim(0-pp, 1+pp)
    ax4.legend(loc='upper left')

    # Show plot
    plt.show()
    #fig.savefig('/home/nicholas/Nextcloud/presentations/presentation_PW12/plotPhaseCurve.png', bbox_inches='tight', dpi=300)



    

    
def plot_phasefold_lightcurve(fig, time, flux, t0, P, tdur, tdepth, odepth, mark='o'):

    """
    Module to plot any given TODO
    """

    # Initial calculations
    dt_c = P/2 * (1 + 4*e*np.cos(w)/np.pi)

    # Setup
    fs, ms = 15, 1
    fig.text(0.5, 0.04, 'Time [d]', ha='center', fontsize=fs)
    fig.text(0.04, 0.5, 'Relative Flu--x', va='center', rotation='vertical', fontsize=fs)
    plt.subplots_adjust(wspace=0.25, hspace=0.25)

    # Find phases
    phase = np.mod(time, per)

    # # Full phase-folded light curve
    ax0 = fig.add_subplot(3,2,(1,2))
    ax0.plot(time, flux, 'bo', markersize=ms)
    ax0.plot(time, flux, 'k-')

    # Full phase-folded light curve
    ax1 = fig.add_subplot(3,2,(3,4))
    ax1.plot(phase, flux, 'bo', markersize=ms)

    # Transit zoom-in
    ax2 = fig.add_subplot(3,2,5)
    ax2.set_xlim(phi-tdur, phi+tdur)
    #ax2.set_ylim(tdepth, ])
    # xpos  = phi-dT+((phi+dT)-(phi-dT))*0.70
    # ypos1 = depth+0.0000
    # ypos2 = depth+0.0002
    # plt.text(xpos, ypos2, '$P={:.4f}$ days'.format(P), fontsize=18)
    # plt.text(xpos, ypos1, '$\phi={:.4f}$ days'.format(phi), fontsize=18)
    # plt.title('Phase  folded', fontsize=18)
    # plot_settings('$t$ $mod$ $P$ [days]', 'Flux')


    # ax2.plot(phase, flux, 'go')
    # ax2.plot(phase, flux, 'k-')
    # # Position text
    # xpos = min(phase)+(max(phase)-min(phase))*0.70
    # ypos = max(flux)-(max(flux)-min(flux))*0.98
    # ax1.text(xpos, ypos,'$P={:.4f}$ d'.format(per), fontsize=fs)
    # # Labels
    # ax2.set_title('Phase folded', fontsize=15)
    #ax2.set_xlim(phase-tdur, phase+tdur)
    #ax1.set_ylim(min(flux), max(flux))

    # Occultation zoom-in
    ax4 = fig.add_subplot(3,2,6)
    #sub3.plot(x, y)
    
    # Finito!
    plt.show()





    
def plot_transit_zoom(): # TODO in progress

    plt.figure(figsize=(4,3))
    # Plot
    plt.plot(time[dex_tra], lc_tra[dex_tra], '-', c='orange', linewidth=2, label='Transit')
    # Text
    x_pos     = t_tra_max - t_tra_tot*1.5
    delta_tra = np.max(lc_tra[dex_tra])-min(lc_tra[dex_tra])
    y_pos     = np.max(lc_tra[dex_tra]) + delta_tra/10.
    plt.text(x_pos, y_pos, r'$\delta_{\mathrm{tra}}=%.1f$ ppm' % delta_tra, fontsize=fs-4)
    # Axes
    plt.xlim(time[dex_tra][0], time[dex_tra][-1])
    plt.ylim(-14000, 3000)
    plt.xlabel('Time [days]')
    plt.ylabel('Relative flux [ppm]')
    plt.tight_layout()
    plt.show()





    

def plot_final_lc(lc):

    zeros = np.zeros(len(lc['time']))
    try: lc['spot']
    except: lc['spot'] = zeros.tolist()
    try: lc['gran']
    except: lc['gran'] = zeros.tolist()
    try: lc['puls']
    except: lc['puls'] = zeros.tolist()
    try: lc['tran']
    except: lc['tran'] = zeros.tolist()
    
    time = lc['time']/86400.
    lc_med = median_filter(lc['sum'], 144)

    fig, ax = plt.subplots(4, 1, figsize=(12, 10), sharex=True)

    ax[0].plot(time, lc['gran'] + lc['puls'], 'g-', label='Gran + Puls')
    ax[1].plot(time, lc['spot'], 'b-', label='Spots')
    ax[2].plot(time, lc['tran'], 'r-', label='Transits')
    ax[3].plot(time, lc['sum'],  'k-', label='Combined')
    ax[3].plot(time, lc_med,     'm-', label='1h median')            

    plt.xlabel('Time [days]')
    fig.text(0.01, 0.5, 'Relative flux [ppm]', va='center', rotation='vertical')

    for i in range(4):
        ax[i].legend(loc='lower right')
        ax[i].set_xlim(time.iloc[0], time.iloc[-1])

    plt.subplots_adjust(wspace=-2, hspace=0)
    
    plt.tight_layout()
    plt.show()


    # plt.figure(figsize=(12,3.5))
    # plt.plot(time, lc['sum'],  'k-', label='Combined')
    # plt.plot(time, lc_med,     'm-', label='1h median')            
    # plt.xlabel('Time [days]')
    # plt.ylabel('Relative flux [ppm]')
    # plt.legend(loc='lower left')
    # plt.xlim(time.iloc[0], time.iloc[-1])
    # plt.tight_layout()
    # plt.show()
