#!/usr/bin/env python3

"""
This file contains tools to generate artificial noise time series
"""

# Python standard
import os
import sys
import datetime

# PlatoSim standard
import scipy
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt 
from scipy.signal import periodogram
from pathlib import Path

# PlatoSim imports
import platosim.plot            as pt
import platosim.utilities       as ut
import platosim.referenceFrames as rf
from platosim.utilities    import errorcode
from platosim.matplotlibrc import latex
latex()

# Global parameters
day2sec = 86400.


#--------------------------------------------------------------#
#                   TIME AND FREQUENCY DOMAIN                  #
#--------------------------------------------------------------#


def numpyFFT(signal, timestep):
    """
    Computes power spectrum of an equidistant time series 'signal'
    using the FFT algorithm. The length of the time series need not
    be a power of 2 (zero padding is done automatically).
    Normalisation is such that a signal A*sin(2*pi*nu_0*t)
    gives power A^2 at nu=nu_0  (IF nu_0 is in the 'freq' array)
    @param signal: the time series [0..Ntime-1]
    @type signal: ndarray
    @param timestep: time step fo the equidistant time series
    @type timestep: float
    @return: frequencies and the power spectrum
    @rtype: array,array
    """

    # Compute the FFT of a real-valued signal. If N is the number
    # of points of the original signal, 'Nfreq' is (N/2+1).

    fourier = np.fft.rfft(signal)
    Ntime = len(signal)
    Nfreq = len(fourier)

    # Compute the power

    power = np.abs(fourier)**2 * 4.0 / Ntime**2

    # Compute the frequency array.
    # First compute an equidistant array that goes from 0 to 1 (included),
    # with in total as many points as in the 'fourier' array.
    # Then rescale the array that it goes from 0 to the Nyquist frequency
    # which is 0.5/timestep

    freq = np.arange(float(Nfreq)) / (Nfreq-1) * 0.5 / timestep

    # That's it!

    return freq, power





def timeSeriesFromFourier(time, freq, ampl, phase, power=1,
                          plot=False, title=False):

    """Generate light curve from Fourier info.

    Paramters
    ---------
    time : ndarray, pdframe
        Time points of which light curve will be generated [s]
    freq : ndarray, pdframe
        Frequencies of sinusoids [d]
    ampl : ndarray, pdframe
        Amplitudes of sinusoids [mmag]
    phase : ndarray, pdframe
        Phases of sinusoids [rad]

    Returns
    -------
    mag : ndarray
        Signal for each time point [mag]
    """

    # Number of pulsation modes
    nmodes = len(freq)

    # Loop over the number of modes and sum of every mode
    mag = np.zeros_like(time)
    for i in range(nmodes):
        mag += ampl[i] * np.sin((2*np.pi * freq[i]) * time + phase[i])

    # Normalize the magnitude so its values is in [-1, 1] (so roots are not undefined)
    # Then add 1, raise the power and substract 1
    A = np.max(np.abs(mag))
    mag = A * (((1 + mag/A)**power) - 1)

    # If requested, plot model 
    if plot:

        fig, ax = plt.subplots(2, 1, figsize=(12, 7))

        mag0 = mag / 1e3  # [mag -> mmag] 
        
        # Plot time series
        ax[0].plot(time, mag0, 'k-', lw=0.3)
        ax[0].set_xlabel(r'Time [days]')
        ax[0].set_ylabel(r'$\delta m$ [mmag]')
        ax[0].set_xlim(time.min(), time.max())
        #ax[0].set_ylim(self.mag.min(),  self.mag.max())
        if title:
            ax[0].set_title(str(title))
        
        # Plot power spectrum
        cadence = np.diff(time)[0]
        freq, power = numpyFFT(mag0, cadence)
        ymin = np.mean(power)
        ymax = np.max(power) + 0.1*np.max(power)
        carray = np.linspace(0, 1, 1000)
        colors = plt.cm.rainbow(carray)
        #freq = 1/freq
        # for i in range(N):
        #     if self.starname is not 'g-Dor':
        #         snr_norm = self.snr[i]/self.snr.max()
        #         dex = ut.findNearestIndex(carray, snr_norm)
        #         ax[1].plot([self.freq[i], self.freq[i]], [ymin, ymax],
        #                    c=colors[dex], lw=snr_norm*3)
        #     else:
        #         ax[1].plot([self.freq[i], self.freq[i]], [ymin, ymax], c='royalblue')
        ax[1].plot(freq, power, '-', c='deeppink', lw=1)
        # Settings
        ax[1].set_ylabel(r'Amplitude [mmag]')
        ax[1].set_xlabel(r'Period, $P$ [day]')
        #ax[1].set_yscale('log')
        ax[1].set_xlim(0, np.max(freq)+1)
        ax[1].set_ylim(ymin, ymax)
        plt.tight_layout()
        plt.show()

    # Return signal
    
    return mag





def timeSeriesFromMeanPSD(freq, psd):

    """
    PURPOSE: Given the average noise profile in the power spectral density, compute the 
             corresponding noisy time series. 

    INPUT: freq: Array containing {n / N / deltat } with n in [0,..,Ntime/2],
                 where N is the number of time points, and deltat the time step.
                 Unit: [Hz | microHz | mHz]
           psd:  Power spectral density: abs(fourier)**2 / Ntime * timestep
                 Unit: [ppm^2/Hz | ppm^2/microHz | ppm^2/mHz]

    OUTPUT: time:    time points [s | Ms | KHz]
            signal:  [ppm]

    EXAMPLE: 
        >>> timeStep = 25.0e-6    # in Ms
        >>> Nfreq = 15001
        >>> freq = arange(float(Nfreq)) / (Nfreq-1) * 0.5 / timeStep
        >>> omega = 2*np.pi*freq
        >>> omegaMin  = 2. * np.pi *  20.0
        >>> omegaKnee = 2. * np.pi * 200.0
        >>> meanPSD = (omega**2 + omegaKnee**2)/(omega**2 + omegaMin**2) *  timeStep
        >>> time, signal = timeSeriesFromMeanPSD(freq, meanPSD)
        >>> nu, noisyPSD = FFTpowerdensity(signal, timeStep)
        >>> plt.loglog(nu, noisyPSD, c="b")
        >>> plt.loglog(freq, meanPSD, c="r")

    NOTE: - N points in the time domain correspond to N/2+1 points in the frequency 
            domain (not N/2), where the division is an _integer_ division (rounded down).

    """

    # Determine the number of time points. 

    Nfreq = len(freq)
    Ntime = 2*(Nfreq-1)

    # Determine the frequency resolution
    # Note: N points in the time domain corresponds to N/2+1 points in the frequency domain (not N/2).

    freqStep = freq[1] - freq[0]
    timeStep = 1./(freqStep*Ntime)

    # Generate the time points of the signal

    time = np.arange(Ntime) * timeStep

    # Generate the real and imaginary parts of the fourier transform with gaussian noise

    realPart = normal(0., 1., Nfreq)
    imagPart = normal(0., 1., Nfreq)

    # Construct the full fourier spectrum, using the proper scale factor
    
    fourier = np.sqrt(psd * Nfreq / timeStep) * (realPart + imagPart * 1j)

    # Inverse-FFT the fourier transform to generate the time series [ppm]

    signal = np.real(np.fft.irfft(fourier)) 

    # That's it!

    return time, signal






#--------------------------------------------------------------#
#                    STOCASHIC NOISE SOURCES                   #
#--------------------------------------------------------------#


#@njit
def getRedNoise(time, currenttime, kicktimestep, Ntime, timescale, varscale, noise, mu, sigma):

    
    signal = np.zeros(Ntime)
    
    for i in range(Ntime):

        # Compute the contribution of each component separately.
        # First advance the time series right *before* the time point i,

        while( (currenttime + kicktimestep) < time[i]):
            noise = noise * (1.0 - kicktimestep/timescale) + np.random.normal(mu[0], sigma[0])
            currenttime = currenttime + kicktimestep

        # Then advance the time series with a small time step right *on* time[i]

        delta  = time[i] - currenttime
        #sigma1 = np.sqrt(delta/timescale)*varscale
        sigma1 = varscale * 0.66  # Correction factor to have varscale in RMS arcsec
        noise  = noise * (1.0 - delta/timescale) + np.random.normal(mu[0], sigma1)
        currenttime = time[i]

        # Add the different components to the signal. 

        signal[i] = np.sum(noise)

    return signal





def modelRedNoise(time, timescale, varscale):

    """Function to generate a red noise time series.
    
    Parameters
    ----------
    time : ndarray
        Time points: time[0..Ntime-1]
    timescale : ndarray
        Time scale tau of each red noise component: timescale[0..Ncomp-1]
    varscale : ndarray
        Variation scale of each red noise component: varscale[0..Ncomp-1]
            
    Returns
    -------
    signal : ndarray
        Signal containing all red noise components: signal[0..Ntime-1]
    """

    Ntime = len(time)
    Ncomp = len(timescale)

    # Set the kick (= excitation) timestep to be one 100th of the
    # shortest noise time scale (i.e. kick often enough).

    kicktimestep = min(timescale) / 100.0
    currenttime  = time[0] - kicktimestep
    
    # Predefine some arrays

    delta = 0.0
    noise = np.zeros(Ncomp)
    mu    = np.zeros(Ncomp)
    sigma = np.sqrt(kicktimestep/timescale)*varscale
    
    # Warm up the first-order autoregressive process

    for i in range(2000):
        noise = noise * (1.0 - kicktimestep / timescale) + np.random.normal(mu, sigma)

    # Start simulating the granulation time series
    
    return getRedNoise(time, currenttime, kicktimestep, Ntime,
                       timescale, varscale, noise, mu, sigma)
    







def modelRedNoisePSD(freq, timescale, varscale):

    """Function to generate a red noise model from the PSD.

    Compute the mean power spectral density (PSD) corresponding to the 
    red noise time series that is generated by modelRedNoise().

    INPUT: freq:       frequency points of the PSD  [microHz | mHz | Hz]
           timescale:  see the function rednoise()
           varscale:   see the function rednoise()

    OUTPUT: psd:  power spectral density   [ppm^2/microHz | ppm^2/microHz | ppm^2/microHz]
    """

    psd = np.zeros_like(freq)

    for n in range(len(timescale)):
        sigma = varscale[n]
        tau = timescale[n]
        psd += sigma * sigma * tau / (1.0 + (2.0*np.pi*freq*tau)**2)
    
    return psd





#--------------------------------------------------------------#
#                   SYSTEMATIC NOISE SOURCES                   #
#--------------------------------------------------------------#


def getPRE(alpha, delta, kappa, quarter, sigma=3,
           ofile=False, table=False, plot=False):

    """Pointing Reproducibility Error (PRE) in PLM reference frame.
    
    Paramters
    ---------
    alpha : float
        Equtorial right ascension [rad]
    delta : float
        Equatorial declination [rad]
    kappa : float
        Orientation of the solar panel [rad]
        This corresponds to (Q1, Q2, Q3, Q4) = (0, pi/2, pi, 3pi/2)
    quarter : range(int, int)
        Mission quarter range, e.g. [Q1, Q8] = range(1, 9)
    sigma : float
        Standard deviation to draw within
    ofile : str
        Output filename pointing to storage directory
    table : bool
        Option to print the table produced
    plot : bool
        Option to show the plot of the PRE data

    Return
    ------
    - A pandas dataframe with PRE data
    - Optionally a feather output file
    """

    # Sort input quarters
    
    n = len(quarter)
        
    # PRE in the PLM reference frame (yaw, pitch, roll)
    # Here t stands for transverse direction and [deg]
    # NOTE: Performance values "as required"
    
    t = 3.0/3600 
    b = 6.0/3600

    # Find distribution within 3 sigma of req.
    
    tt = np.array([np.random.normal(0, t/sigma) for i in range(n)])
    bb = np.array([np.random.normal(0, b/sigma) for i in range(n)])

    # Corresponding yaw, pitch, roll

    # y = tt
    # z = 3 * y
    # x = bb - z
    y = tt
    z = bb - y
    x = tt

    # ICRS pointing angles
    
    phi   = np.deg2rad(alpha)
    theta = np.deg2rad(delta)

    # Find change of pointing for each quarter
    
    PRE = np.zeros((n, 4))
    for i in range(n):
        data = rf.perturbPlatformPointing(x[i], y[i], z[i], phi, theta)[0]
        PRE[i,:] = np.append(quarter[i], data)

    df = pd.DataFrame(PRE, columns=["quarter", "yaw", "pitch", "roll"])
    df = df.astype({"quarter":int, "yaw":np.float64, "pitch":np.float64, "roll":np.float64})
    
    df0 = df.copy()
    df0.iloc[:,1:] = df0.iloc[:,1:] * 3600
        
    df1 = df.copy()
    df1.rename(columns={"yaw":"alpha", "pitch":"delta", "roll":"kappa"}, inplace=True)
    df1.iloc[:,1] = df1.iloc[:,1] + alpha
    df1.iloc[:,2] = df1.iloc[:,2] + delta
    df1.iloc[:,3] = df1.iloc[:,3] + kappa

    # Print generated values
           
    if table:
        print('\nChange of coordinates [arcsec]')
        print(df0)
        print('\nNew ICRS coordinates [deg]')
        print(df1)

    # Plot distributions
    
    t *= 3600
    b *= 3600
    y = t/sigma
    z = 3 * y
    x = np.abs(b/sigma - z)
    xx = np.linspace(-10*x, 10*x, 1000)
    alpha0 = (df1.alpha - alpha) * 3600 / 18
    delta0 = (df1.delta - delta) * 3600 / 18

    fig, ax = plt.subplots(1, 2, figsize=(9,4))

    # Plot PDF
    ax[0].set_title(f'PRE distributions at {sigma}$\sigma$')        
    ax[0].plot(xx, scipy.stats.norm.pdf(xx, 0, x)*100, '-', c='b', label='Trans.')
    ax[0].plot(xx, scipy.stats.norm.pdf(xx, 0, z)*100, '-', c='m', label='Rot.')
    ax[0].set_xlabel('Platform pointing errors in FPA [pixel]')
    ax[0].set_ylabel('Probability (PDF) [\%]')
    ax[0].set_xlim(xx[0], xx[-1])
    ax[0].legend()

    # Show distribution on sky 
    ax[1].grid(zorder=0)
    ax[1].plot(0, 0, 'k*', ms=10, zorder=2)
    for i in range(len(alpha0)):
        ax[1].scatter(alpha0[i], delta0[i], marker=f'${i+1}$', s=50, alpha=0.8, zorder=2)
    ax[1].set_title('Distribution on Sky')
    ax[1].set_xlabel('RA [pixel]')
    ax[1].set_ylabel('Dec [pixel]')
    ax[1].set_aspect('equal', adjustable='box')

    # Settings
    lim = np.max([np.max(np.abs(alpha0)), np.max(np.abs(delta0))])
    lim += lim/10.
    ax[1].set_xlim(-lim, +lim)
    ax[1].set_ylim(-lim, +lim)    
    plt.tight_layout()
        
    # Plot figure above

    if plot: plt.show()
    
    # Save file with relative pointing errors [deg]
    
    if ofile:
        df.to_csv(ofile, sep=" ", header=False, index=False)
        fig.savefig(f"{ofile[:-4]}.png", bbox_inches='tight', dpi=200)

    # That's it!
    
    return PRE





def getAPE(alpha, delta, kappa, sigma=3,
           ofile=False, table=False, plot=False):

    """Pointing Reproducibility Error (PRE) in P/L reference frame.

    The mission reuquirements are 4.5 arcmin in X and Y, and 9 arcmin
    for the rotation. This gives a RMS value of 1.5 arcsec (6 pixel)
    
    Paramters
    ---------
    alpha : float
        Equtorial right ascension [rad]
    delta : float
        Equatorial declination [rad]
    kappa : float
        Orientation of the solar panel [rad]
        This corresponds to (Q1, Q2, Q3, Q4) = (0, pi/2, pi, 3pi/2)
    sigma : float
        Standard deviation to draw within
    ofile : str
        Output filename pointing to storage directory
    table : bool
        Option to print the table produced
    plot : bool
        Option to show the plot of the APE data

    Return
    ------
    - A pandas dataframe with APE data
    - Optionally a feather output file
    """

    # APE in the PLM reference frame (yaw, pitch, roll)
    # Here t stands for transverse direction and [deg]    
    # NOTE: Performance values "as required"
    
    t = 4.5/60  # [deg]
    b = 9.0/60  # [deg]
        
    # Find distribution within 3 sigma of req.

    tt = np.array([np.random.normal(0, t/sigma) for i in range(26)])
    bb = np.array([np.random.normal(0, t/sigma) for i in range(26)])

    # Corresponding yaw, pitch, roll

    dy = tt
    dz = 3 * dy
    dx = bb - dz

    # Store APE
    
    APE = np.transpose([tt, bb])
    df  = pd.DataFrame(APE, columns=["tilt", "azimuth"])

    # Print distributions to bash
    
    if table:
        
        print(f'\nCamera alignment errors for all 26 cameras [pixel]')
        APE0 = np.transpose([tt, bb, dx, dy, dz]) * 3600 / 15
        df0  = pd.DataFrame(APE0, columns=["Alt", "Az", "Yaw", "Pitch", "Roll"])
        print(df0)
        
    # Create figure object
    
    t *= 3600 / ( 15 * sigma * (sigma-1))
    b *= 3600 / ( 15 * sigma * (sigma-1))
    xx = np.linspace(-10*t, 10*t, 1000)

    fig, ax = plt.subplots(1, 2, figsize=(10,5))

    # Plot PDF
    ax[0].plot(xx, scipy.stats.norm.pdf(xx, 0, t)*100, '-', c='b', label='Trans.')
    ax[0].plot(xx, scipy.stats.norm.pdf(xx, 0, b)*100, '-', c='m', label='Rot.')
    ax[0].set_title(f'APE distributions at {sigma}$\sigma$')
    ax[0].set_xlabel('Camera misalignment in FPA [pixel]')
    ax[0].set_ylabel('Probability (PDF) [\%]')
    ax[0].set_xlim(xx[0], xx[-1])
    ax[0].legend()

    # Plot distribution on sky
    azim = df.azimuth*3600/15
    tilt = df.tilt*3600/15
    ax[1].grid(zorder=0)
    ax[1].plot(0, 0, 'k*', ms=10, zorder=2)
    for i in range(24):
        if i < 6:
            ax[1].scatter(azim[i], tilt[i], marker=rf'${i+1}$', s=50, c='g', alpha=0.8,zorder=2)
        else:
            ax[1].scatter(azim[i], tilt[i], marker=rf'${i+1}$', s=70, c='g', alpha=0.8,zorder=2)
    ax[1].scatter(azim[24], tilt[24], marker=r'$1$', s=50, c='b', alpha=0.8, zorder=3)
    ax[1].scatter(azim[25], tilt[25], marker=r'$2$', s=50, c='r', alpha=0.8, zorder=3)
    ax[1].set_title('Relative offset')
    ax[1].set_xlabel('Azimuth [pixel]')
    ax[1].set_ylabel('Tilt [pixel]')
    ax[1].set_aspect('equal', adjustable='box')

    # Settings
    lim = np.max([np.max(np.abs(azim)), np.max(np.abs(tilt))])
    lim += lim/10.
    ax[1].set_xlim(-lim, lim)
    ax[1].set_ylim(-lim, lim)
    plt.tight_layout()

    # Plot figure above 
    
    if plot: plt.show()
    
    # Save APE camera misalignments
    
    if ofile:
        df.to_csv(ofile, sep=" ", header=False, index=False)
        fig.savefig(f"{ofile[:-4]}.png", bbox_inches='tight', dpi=200)

    # That's it!
    
    return APE





def getTED(quarter, model="poly", ofile=False, table=False, plot=False):

    """Generate a Themo-Elastic Drift (TED) file.
   
    This function generates a complete TED model returned in euler angles.
    
    Paramters
    ---------
    quarter : range
        Range of quarters (e.g. range(1,8) for Q1-Q8).
    model : str
        Model to produce TED: 'linear' or 'poly'.
    ofile : bool, str
        Parse string to save the model to a ascii file.
    plot : bool
        True will make a plot of the models.
        
    Return
    ------
    Output file if requested.
    """
        
    # Constants
    time0 = np.arange(0, ut.year()/4, 25)
    cols  = ["yaw", "pitch", "roll"]
    N     = len(quarter)
    n     = len(time0)

    # Create data frame and store default time0 for fit
    df  = pd.DataFrame()
    df1 = pd.DataFrame()
    A = np.zeros((len(quarter), 4))
    # Loop over each quarter

    for Q in range(quarter[0]-1, quarter[-1]):

        # Time column
        t0 = round(ut.year()/4 * Q)
        df1["time"] = t0 + np.arange(0, n) * 25

        # Generate linear model

        # Generate a random 2nd order polynomial

        for col in cols:

            if model == 'linear':
                a = 1.3 * 15      
                if col == "roll":
                    df1[col] = np.zeros(n)
                else:
                    df1[col] = np.linspace(0, a, n)

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
        # Array with all amplitudes
        A[Q-quarter[0],0] = Q+1
        A[Q-quarter[0],1] = df1.yaw.max()   - df1.yaw.min()
        A[Q-quarter[0],2] = df1.pitch.max() - df1.pitch.min()
        A[Q-quarter[0],3] = df1.roll.max()  - df1.roll.min()

    # Show amplitudes
    
    if table:
        print('TED model amplitudes [arcsec]')
        names = ['Quarter', 'A_yaw', 'A_pitch', 'A_roll']
        da = pd.DataFrame(A, columns=names)
        da = da.sort_values(['Quarter'])
        da = da.astype({'Quarter':np.int})
        da = da.reset_index(drop=True)
        print(da)

    # Plot model
    
    fig, ax = plt.subplots(3,1,figsize=(9, 6))

    # Plots
    for i, col in zip(range(3), cols):
        ax[i].plot(df["time"]/day2sec, df[col], 'k-')
        ax[i].axhline(y=0, linestyle=':', color='k')
        Qday = ut.year()/86400/4
        for k in range(N-1):
            ax[i].axvline(x=quarter[k]*Qday, linestyle='--', color='b')

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
        
    # Plot figure above 
    
    if plot: plt.show()
        
    # Save data in one big drift text file for PlatoSim
    if ofile:
        df.to_csv(ofile, sep=" ", header=False, index=False)
        fig.savefig(f"{ofile[:-4]}.png", bbox_inches='tight', dpi=200)




            
def getACS(time, rms=[0.038, 0.038, 0.040], ofile=False, plot=False):

    """Generate a Attitude and orbit Control System (ACS) jitter file.
   
    This function generates a complete ACS model returned in euler angles.
    
    Paramters
    ---------
    quarter : range
        Range of quarters (e.g. range(1,8) for Q1-Q8).
    ofile : bool, str
        Parse string to save the model to a ascii file.
    plot : bool
        True will make a plot of the models.
        
    Return
    ------
    Output file if requested.
    """

    # Generate red noise components

    tau = time[1] - time[0]
    df = pd.DataFrame()
    df['t'] = time
    for i, n, in zip(range(3), ['x', 'y', 'z']):
        print(f'Generating {n} red noise component using {rms[i]}')
        df[n] = modelRedNoise(time, timescale=np.array([tau]), varscale=np.array(rms[i]))

    # Plot time series and PSD

    fig1, ax = pt.plotYawPitchRollTimeSeries(df.t/3600, np.array([df.x, df.y, df.z]),
                                             units=["hours", "arcsec"])
    fig2, ax = pt.plotYawPitchRollPSD(df.t, np.array([df.x, df.y, df.z]))
    if plot: plt.show()
        
    # Save data in one jitter file
    
    if ofile:
        df.to_csv(ofile, sep=" ", header=False, index=False)
        fig1.savefig(f"{ofile[:-4]}_timeseries.png", bbox_inches='tight', dpi=200)
        fig2.savefig(f"{ofile[:-4]}_psd.png", bbox_inches='tight', dpi=200)
    
    



def getDataGaps(time, quarter=range(1,9), ofile=False, plot=False):

    """Function to create a job script to be used on the VSC.

    All time gaps are based on knowledge from the Kepler mission.
    The following time gaps are considered:

      - Quarterly rolls
      - Monthly data donlinks
      - Loss of fine guidance
      - Save mode events

    The Quarterly rolls introduce a time gap of 1-2 days happening every
    three months (hence four times in a year). These are are needed in 
    order to align the spacecraft's solar panels towards the Sun. Note that
    the PES called PRE takes into account for the error in every pointing
    from quarter to quarter. 

    Montly data downlinks is needed during the mission. During these events
    the spacecraft changed its orientation once a month to point at Earth
    and downlink the last month of data. This caused a gap in the time series
    data, as the spacecraft can not collect while it was downlinking. When 
    the spacecraft returns to the regular pointing, its motion induced what
    is known as a “thermal transient.” This means that the telescope components
    and detector electronics underwent a temperature change, and the electron
    count reading was temporarily increased. A change in temperature can slightly
    change the camera focus. This manifests in the data as a downward slope 
    caused by “reheating,” while the flux returns to its previous level. 

    Loss of fine guidance and coarse pointing events which results in
    a lower photometric precision. Typically these cadences are not
    suggested for use in photometry. Since these event do not introduce
    thermal transients in the light curve, they can simply be removed 
    after the simulations have been generated by the user. 

    Safe modes (Kepler Data Characteristics Handbook, Section 4.2)
    are another type of thermal transient that appears in Kepler data. 
    A safe mode occured when the telescope temporarily shut off operation
    due to an unexpected event, usually caused by an issue with the 
    detector electronics. Typically a positive flux jump happens each event
    (probably due to an increased gain) an exponentially decrease to the
    count level before the event over a period of a few days.
    """

    # Initialise random generator
    rng = np.random.default_rng()

    # Storage arrays
    roll   = np.zeros_like(time, dtype=bool)
    #link   = np.zeros_like(time, dtype=bool)
    jitter = np.zeros_like(time, dtype=bool)
    safe   = np.zeros_like(time, dtype=bool)

    # QUARTERLY ROLLS

    roll_period   = ut.quarter()
    roll_duration = 2.0
    roll_anomaly  = 0.5

    n_roll = len(quarter)
    roll_events = quarter
    roll_event0 = np.zeros(n_roll)
    roll_event1 = np.zeros(n_roll)
    
    for i, Q in zip(range(n_roll), roll_events):
        roll_gap = roll_duration + np.random.uniform(-roll_anomaly, roll_anomaly) # [d]
        roll_event0[i] = (roll_period * Q - roll_gap/2) * day2sec                 # [s]
        roll_event1[i] = (roll_period * Q + roll_gap/2) * day2sec                 # [s]
        roll_dex       = np.where((time>=roll_event0[i]) & (time<=roll_event1[i]))[0]
        roll[roll_dex] = True

    # DOWNLINK GAPS
    # NOTE: Not applicable for PLATO!
    
    # link_period   = 365.25/4/3
    # link_duration = 5/24.           # [d] i.e. 5 hours
    # link_anomaly  = 0.5/24.         # [d] i.e. 0.5 hour

    # # Remove overlaps with quarters
    # array0 = np.arange((quarter[0]-1)*roll_period, quarter[-1]*roll_period, link_period)[1:]
    # array1 = array0[2::3]
    # link_events = [i for i in array0 if i not in array1]

    # n_link = len(link_events)
    # link_event0 = np.zeros(n_link)
    # link_event1 = np.zeros(n_link)

    # for i, L in zip(range(n_link), link_events):
    #     link_gap = link_duration + np.random.uniform(-link_anomaly, link_anomaly)
    #     link_event0[i] = (L - link_gap/2) * day2sec
    #     link_event1[i] = (L + link_gap/2) * day2sec
    #     link_dex       = np.where((time>=link_event0[i]) & (time<=link_event1[i]))[0]
    #     link[link_dex] = True
        
    # LOSS OF FINE GUIDANCE

    jitter_freq = 120
    jitter_offset = np.random.uniform(0, jitter_freq)
    t = (quarter[0] - 1) * roll_period - jitter_offset
    jitter_duration = 0.5/24.
    jitter_anomaly  = 0.1/24.
    jitter_events   = []

    while t < quarter[-1] * roll_period:            
        jitter_event = np.random.poisson(lam=jitter_freq)
        t += jitter_event
        jitter_events.append(t)

    n_jitter = len(jitter_events)
    jitter_event0 = np.zeros(n_jitter)
    jitter_event1 = np.zeros(n_jitter)
    event0 = roll_event0 # np.concatenate((roll_event0, link_event0))
    event1 = roll_event1 # np.concatenate((roll_event1, link_event1))

    for i, J in zip(range(n_jitter), jitter_events):
        jitter_gap = jitter_duration + np.random.uniform(-jitter_anomaly, jitter_anomaly)
        jitter_event0[i] = (J - jitter_gap/2) * day2sec
        jitter_event1[i] = (J + jitter_gap/2) * day2sec
        jitter_dex       = np.where((time>=jitter_event0[i]) & (time<=jitter_event1[i]))[0]
        jitter[jitter_dex] = True

        for j in range(len(event0)):
            if ((jitter_event0[i] >= event0[j]) and (jitter_event1[i] <= event1[j]) or
                (jitter_event0[i] <= event0[j]) and (jitter_event1[i] >= event1[j]) or
                (jitter_event0[i] <= event0[j]) and (jitter_event1[i] >= event0[j]) or
                (jitter_event0[i] <= event1[j]) and (jitter_event1[i] >= event1[j])):
                jitter_event0[i] += event1[j] - event0[j]
                jitter_event1[i] += event1[j] - event0[j]

    # SAFE MODE EVENTS

    safe_freq = 270
    safe_offset = np.random.uniform(0, safe_freq)
    t = (quarter[0] - 1) * roll_period - safe_offset
    safe_duration = 1
    safe_anomaly  = 12/24.
    safe_events   = []

    while t < quarter[-1] * roll_period:            
        safe_event = np.random.poisson(lam=safe_freq)
        t += safe_event
        safe_events.append(t)
        
    n_safe = len(safe_events)
    safe_event0 = np.zeros(n_safe)
    safe_event1 = np.zeros(n_safe)
    event0 = np.concatenate((roll_event0, jitter_event0))
    event1 = np.concatenate((roll_event1, jitter_event1))        

    for i, S in zip(range(n_safe), safe_events):
        safe_gap = safe_duration + np.random.uniform(-safe_anomaly, safe_anomaly)
        safe_event0[i] = (S - safe_gap/2) * day2sec
        safe_event1[i] = (S + safe_gap/2) * day2sec
        safe_dex       = np.where((time>=safe_event0[i]) & (time<=safe_event1[i]))[0]
        safe[safe_dex] = True
        
        for j in range(len(event0)):
            if ((safe_event0[i] >= event0[j]) and (safe_event1[i] <= event1[j]) or
                (safe_event0[i] <= event0[j]) and (safe_event1[i] >= event1[j]) or
                (safe_event0[i] <= event0[j]) and (safe_event1[i] >= event0[j]) or
                (safe_event0[i] <= event1[j]) and (safe_event1[i] >= event1[j])):
                safe_event0[i] += 2*(event1[j] - event0[j])
                safe_event1[i] += 2*(event1[j] - event0[j])
        
    # Show figure

    if n_roll != 0:

        fig, ax = plt.subplots(figsize=(9, 3))
        ax.axhline(y=0, linestyle=':', color='k')

        for i in range(n_roll):
            ax_roll = ax.axvspan(roll_event0[i]/day2sec, roll_event1[i]/day2sec,
                                 color='b', alpha=0.5)
            if i == n_roll-1: ax_roll.set_label('Quarter rolls')

        # for i in range(n_link):
        #     ax_link = ax.axvspan(link_event0[i]/day2sec, link_event1[i]/day2sec,
        #                          color='m', alpha=0.5)
        #     if i == n_link-1: ax_link.set_label('Downlinks')

        for i in range(n_jitter):
            ax_jitter = ax.axvspan(jitter_event0[i]/day2sec, jitter_event1[i]/day2sec,
                                   color='orange', alpha=0.5)
            if i == n_jitter-1: ax_jitter.set_label('Loss of FGS')

        for i in range(n_safe):
            ax_safe = ax.axvspan(safe_event0[i]/day2sec, safe_event1[i]/day2sec,
                                 color='r', alpha=0.5)
            if i == n_safe-1: ax_safe.set_label('Safe modes')

        # Labels
        ax.set_xlabel("Time [days]")
        ax.set_ylabel("Arb.")
        ax.legend(ncol=4, loc="center", bbox_to_anchor=(0.5, 1.2))

        # Layout
        ax.set_yticklabels([])
        ax.set_xlim(time[0]/day2sec, time[-1]/day2sec+5)
        ax.set_ylim(-1,1)
        plt.tight_layout()

    # Show plot
    if plot: plt.show()

    # Create pandas data frame for different flags    

    df = pd.DataFrame()
    df["time"]   = time
    df["roll"]   = roll
    #df["link"]   = link
    df["jitter"] = jitter
    df["safe"]   = safe
    df["all"]    = roll + jitter + safe
        
    # Compute event times

    t0 = np.concatenate((roll_event0, safe_event0))
    t1 = np.concatenate((roll_event1, safe_event1))
    dt = pd.DataFrame({'t0':t0, 'td':t1-t0})
    
    # Save file (and plot) if requested
    
    if ofile:
        df.reset_index(drop=True, inplace=True)
        df.to_feather(ofile)
        dt.to_feather(f"{ofile[:-4]}.tab")
        fig.savefig(f"{ofile[:-4]}.png", bbox_inches='tight', dpi=200)

    # That's it!
    
    return df, t0, t1-t0





def temperatureTransients(time, t0, td, tempCCD=203.15, tempConst=10, gapSize=0.1, timeSpan=30,
                          timeScale=False, amplitude=False, ofile=False, plot=False):

    """Function to model detector temperature transients.

    Parameters
    ----------
    time : ndarray
        Time points with gaps [s]
    tempCCD : float
        Average CCD temperature [K]
    gapSize : float, ndarray
        Minimum duration to detect the gaps [days]
    timeSpan : float, ndarray
        Time duration to model each transient [days]
        Parameter used to minimize run time.
    timeScale : float, ndarray
        Time scale duration of the flare(s) [days]
        Linear model is used if False.
    amplitude : float, ndarray
        Temperature amplitude of the transient(s) [K]
        Linear model is used if False.
    plot : bool
        Flag to plot the generated model.

    NOTE: The following code snippet removes the data gaps from the df:
    >> df0 = df.drop(df[(df.roll == True)   | (df.link == True) |
                        (df.jitter == True) | (df.safe == True)].index)
    """
    
    # Secure numpy syntax
    
    if isinstance(time, pd.Series):
        time = time.to_numpy()
    if isinstance(t0, list):
        t0 = np.array(t0)
    if isinstance(td, list):
        td = np.array(td)
        
    # Create temperature array to write to
    
    time = time / 86400.
    temp = np.zeros_like(time)

    # Convert to days
    
    timeGap0 = t0   / 86400.
    tdurGap  = td   / 86400.
    timeGap1 = timeGap0 + tdurGap
    
    # Unit parameters
    
    cadence = np.diff(time)[0]
    ndex    = int(timeSpan / cadence)
    n       = len(t0)

    # Indices for start and end of each event
    
    timeDex0 = [np.argmin(np.absolute(time - timeGap0[i])) for i in range(n)]
    timeDex1 = [np.argmin(np.absolute(time - timeGap1[i])) for i in range(n)]
    
    # Linear CCD temperature dependece with gap size
    
    if not timeScale:
        timeScale = 1/tempConst * tdurGap
    
    if not amplitude:
        amplitude = tempConst * tdurGap

    # Secure that a single event works
    
    try: len(timeGap1)
    except: timeGap1 = [timeGap1]
    try: len(timeScale)
    except: timeScale = [timeScale]
    try: len(amplitude)
    except: amplitude = [amplitude]

    # Loop over each transient event

    for i in range(n):

        # Time array during transient event
        
        tn = time/timeScale[i]

        # Model parameters of transients
        
        a0 = 0.689
        a1 = -1.6
        b0 = 1 - a0
        b1 = -0.2783
        
        # Secure last event is within time series
        
        if timeDex1[i]+ndex > len(time):
            timeDex2 = len(time)-1
        else:
            timeDex2 = timeDex1[i] + ndex

        # Loop over every time-step in the transient time interval

        for j,k in zip(range(timeDex1[i], timeDex2), range(len(tn))):
            temp[j] += (a0 * np.exp(a1 * tn[k]) +
                        b0 * np.exp(b1 * tn[k]) * amplitude[i])

            # Add amplitude for overlapping events
            
            if (temp[timeDex0[i]] > tempCCD+0.001) and (j==timeDex1[i]):
               amplitude[i] += temp[timeDex0[i]]

    # Add CCD zero point temperature
    
    temp += np.ones_like(temp) * tempCCD
               
    # Plot if requested
    
    if ofile:
        
        fig, ax = plt.subplots(figsize=(9,3))
        for i in range(n):
            ax.axvspan(timeGap0[i], timeGap1[i], color='b', alpha=0.2)
        ax.plot(time, temp, '-', lw=1, c='deeppink')
        ax.set_xlabel('Time [d]')
        ax.set_ylabel('CCD temperature [K]')
        ax.set_xlim(np.min(time), np.max(time))
        plt.tight_layout()

    # Plot figure above
        
    if plot: plt.show()        

    # Save data if requested

    if ofile:
        np.savetxt(ofile, np.transpose([time*day2sec, temp]), fmt=['%.1f', '%.6f'])
        fig.savefig(f"{ofile[:-4]}.png", bbox_inches='tight', dpi=200)
        
    # That's it!
        
    return temp






def getGain(sigma=3, gain0CCD=False, ofile=False, plot=False):

    """ointing Reproducibility Error (PRE) in PLM reference frame.
    
    TODO under construction!

    Paramters
    ---------

    Return
    ------
    """

    # Total number of CCD x 2 halves

    nCCD = 104
    
    # Find distribution within 3 sigma of req.

    if not gain0CCD: gain0CCD = 1.8
    gainRef = gain0CCD / sigma
    deltaGainF = np.array([np.random.normal(0, gainRef) for i in range(nCCD)])
    deltaGainE = np.array([np.random.normal(0, gainRef) for i in range(nCCD)])
    
    # Plot distributions
    
    t *= 3600
    b *= 3600
    y = t/sigma
    z = 3 * y
    x = np.abs(b/sigma - z)
    xx = np.linspace(-10*x, 10*x, 1000)

    ra0  = (df1.RA-ra)   * 3600 / 18
    dec0 = (df1.Dec-dec) * 3600 / 18

    fig, ax = plt.subplots(1, 2, figsize=(9,4))

    ax[0].set_title(f'PRE distributions at {sigma}$\sigma$')        
    ax[0].plot(xx, scipy.stats.norm.pdf(xx, 0, x)*100, '-', c='b', label='Trans.')
    ax[0].plot(xx, scipy.stats.norm.pdf(xx, 0, z)*100, '-', c='m', label='Rot.')
    ax[0].set_xlabel('Platform pointing errors in FPA [pixel]')
    ax[0].set_ylabel('Probability (PDF) [\%]')
    ax[0].set_xlim(xx[0], xx[-1])
    ax[0].legend()

    ax[1].grid(zorder=0)
    ax[1].plot(0, 0, 'k*', ms=10, zorder=2)
    for i in range(len(ra0)):
        ax[1].scatter(ra0[i], dec0[i], marker=f'${i+1}$', s=50, alpha=0.8, zorder=2)
    ax[1].set_title('Distribution on Sky')
    ax[1].set_xlabel('RA [pixel]')
    ax[1].set_ylabel('Dec [pixel]')
    ax[1].set_aspect('equal', adjustable='box')

    lim = np.max([np.max(np.abs(ra0)), np.max(np.abs(dec0))])
    lim += lim/10.
    ax[1].set_xlim(-lim, +lim)
    ax[1].set_ylim(-lim, +lim)        
    plt.tight_layout()

    # Plot figure above 
    
    if plot: plt.show()
        
    # Save file with relative pointing errors [deg]
    
    if ofile:
        df.to_csv(ofile, sep=" ", header=False, index=False)
        fig.savefig(f"{ofile[:-4]}.png", bbox_inches='tight', dpi=200)

    # That's it!
    
    return PRE
