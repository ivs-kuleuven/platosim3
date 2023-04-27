
import h5py
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.pyplot import cm
from scipy.interpolate import make_interp_spline

from platosim.plot import plotPSD
from platosim.noise import redNoise, rednoiseModel, powerDensityFFT


# Change matplotlib defautls

plt.rc('xtick', labelsize=13)
plt.rc('ytick', labelsize=13)
plt.rcParams.update({'font.size': 13})

# Constants

rad2arcsec = 206265.
sec2day    = 1/86400.
#----------------------------

# Prime TED

filename = '/home/nicholas/software/python/platonium/models/aocs_prime_2021-01/04_rev349_ARC3_EOL_90d_3600.csv'
data = np.loadtxt(filename, delimiter=';')
time_ted = data[:,0]  # [s]
# Extract mean dirft time series
ted_dir = np.zeros((24,len(data)))
ted_rot = np.zeros((24,len(data)))
for i in range(data.shape[1]-1):
    if i % 2 == 0:
        ted_dir[int(i/2),:] = data[:,i+1]
    if i % 2 != 0:
        ted_rot[int(i/2),:] = data[:,i+1]
ted_dir_mean = np.mean(ted_dir, axis=0)
ted_rot_mean = np.mean(ted_rot, axis=0)

# Prime Jitter + TED

# filename = '/home/nicholas/software/python/platonium/models/aocs_prime_2021-01/03_PLATO_PDR_AOCSandTED.csv'
# data = np.loadtxt(filename, delimiter=';')
# time_ted_jitter = data[:,0]  # [s]
# x    = data[:,1]  # [arcsec]
# y    = data[:,2]  # [arcsec]
# z    = data[:,3]  # [arcsec]
# data_ted_jitter = [x, y, z]

# Prime jitter

filename_jitter_prime = '/home/nicholas/software/python/platonium/models/aocs_prime_2021-01/01_PLATO_PDR_FPM_02_longrun_APE.csv'
data = np.loadtxt(filename_jitter_prime, delimiter=';')
time_jitter = data[:,0] - data[:,0][0]  # [s] Measure from 1000s needs to be subtracted
x = data[:,1] * rad2arcsec              # [arcsec]
y = data[:,2] * rad2arcsec              # [arcsec]
z = data[:,3] * rad2arcsec              # [arcsec]
data_jitter = [x, y, z]

#----------------------------

test = 'yes'




if test == 'jitterPrime':

    # Find frequencies

    sampling = np.diff(time)[0]/1e3    # [Ms]
    freq, PSD0 = powerDensityFFT(data_jitter_prime[0], sampling)
    freq, PSD1 = powerDensityFFT(data_jitter_prime[1], sampling)
    freq, PSD2 = powerDensityFFT(data_jitter_prime[2], sampling)
    PSD = [PSD0, PSD1, PSD2]

    # Plot PSD

    fig = plt.figure(figsize=(8,4))
    plotPSD(fig, freq, PSD[0], title='PSD jitter from Prime', carbox=False)
    plt.show()


    timescale = np.array([1e5])    # [Ms]
    varscale  = np.array([1e-3])   # [arcsec/microHz] 

    #signal    = redNoise(time, timescale, varscale)  # [arcsec]

    #plt.figure(figsize=(10,6))
    #plt.plot(time/3600., signal,'b-')
    #plt.show()

    # Find frequencies

    sampling = np.diff(time)[0]/1e6    # [Ms]
    fourier  = np.fft.rfft(signal)
    Nfreq    = len(fourier)
    freq     = np.arange(float(Nfreq)) / (Nfreq-1) * 0.5 / sampling

    #freq = 1e-6/time  # [microHz]

    # Plot with IvS code

    #freq0, psd0 = powerDensityFFT(signal, timescale[0]/1e6) # [arcsec/microHz]
    #freq1, psd1 = powerDensityFFT(signal, timescale[1]/1e6)

    # Plot with noise code

    psd = rednoiseModel(freq, timescale, varscale) # [arcsec/microHz]




if test == 'driftPrime':

    cams = ['CAM 1.1', 'CAM 1.2', 'CAM 1.3', 'CAM 1.4', 'CAM 1.5', 'CAM 1.6',
            'CAM 2.1', 'CAM 2.2', 'CAM 2.3', 'CAM 2.4', 'CAM 2.5', 'CAM 2.6',
            'CAM 3.1', 'CAM 3.2', 'CAM 3.3', 'CAM 3.4', 'CAM 3.5', 'CAM 3.6',
            'CAM 4.1', 'CAM 4.2', 'CAM 4.3', 'CAM 4.4', 'CAM 4.5', 'CAM 4.6',]
    color = cm.rainbow(np.linspace(0,1,49))


    # Plot each camera and mean value

    fig = plt.figure(figsize=(15,10))
    for i,c in zip(range(data_ted_dir.shape[1]-1), color):
        #Plot even number i.e. directional in file
        if i % 2 == 0:
            plt.plot(time_ted, data_ted_dir[:,i+1], '-', c=c, label='{}'.format(cams[int(i/2)]))
    plt.plot(time_ted, data_dir_mean, 'k-', lw=5, label='Mean')
    #Plot settings
    plt.title('TED Cartesian (yaw + pitch) drift from Prime')
    plt.xlabel('Time [days]')
    plt.ylabel('Amplitude [arcsec]')
    plt.legend(bbox_to_anchor=(1.01, 1))
    plt.show()


if test == 'yes':

    # Find freq and PSD

    unit = 1e-6
    samp_ted  = np.diff(time_ted)[0] * unit
    freq_ted, PSD_ted = powerDensityFFT(ted_dir_mean/unit, samp_ted)

    samp_jitter = np.diff(time_jitter)[0] * unit
    freq_jitter, PSD_jitter = powerDensityFFT(data_jitter[0]/unit, samp_jitter)

    freq = [freq_ted, freq_jitter]
    PSD  = [PSD_ted, PSD_jitter]

    # Plot drift PSD

    #fig = plt.figure(figsize=(10,6))
    #plotPSD(fig, freq, PSD, units=['$\mu$Hz', 'ppm'], title='PSD drift for Prime')
    #plt.show()

    # Interpolate (piecewise cubic) into higher resolution grid

    grid_no  = int(time_ted[-1]/25.)
    time_int = np.linspace(time_ted[0], time_ted[-1], grid_no)
    tedgrid  = make_interp_spline(time_ted, ted_dir_mean, k=3)
    ted_int  = tedgrid(time_int)

    plt.figure(figsize=(15,10))
    plt.plot(time_int, ted_int, 'r,')
    plt.plot(time_ted, ted_dir_mean, 'k.')
    plt.show()



