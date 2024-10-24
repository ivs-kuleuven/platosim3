#!/usr/bin/env python

"""
This script uses the PIC targets and their contaminants to simulate realistic imagettes
"""

import os
import sys
import glob
import h5py
import argparse
import datetime
import numpy as np
import pandas as pd
from colorama import Back, Fore, Style
from scipy.ndimage import median_filter

from platosim.utilities import errorcode, normalize
from platosim.photometryfile import PhotometricFile
from platosim.plot import axes_minmax, axes_maskupdates

# Import repo matplotlibrc
import matplotlib.pyplot as plt
import matplotlibrc

# Turn off warnings
import warnings
warnings.filterwarnings("ignore")

# Constants
ppmh = 144   # For a cadence of 25s


#===========================================================================#
#                             PLATOnium's QUICK TOOLS                       #
#===========================================================================#

parser = argparse.ArgumentParser(epilog=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter,
                                 description=errorcode('software', '\nPLATOnium Quick Tools'))

parser.add_argument('--combineFeathers',     action='append', type=str, nargs='*', metavar='PICKLE',  help='')
parser.add_argument('--convertCSV2TXT',     action='append', type=str, nargs=3,   help='')
parser.add_argument('--generateVarSourceIndex', action='store_true',  help='')
parser.add_argument('--generateDriftFile',  action='append', type=str, nargs=3,   metavar='FILE', help='')
parser.add_argument('--plotLCfromPipeline', action='append', type=str, nargs=1,   metavar='FILE', help='')
parser.add_argument('--plotLCfromAscii',    action='append', type=str, nargs=1,   metavar='FILE', help='')
parser.add_argument('--plotLCfromHDF5',     action='append', type=str, nargs=1,   metavar='FILE', help='')
args = parser.parse_args()




def addTimeColToTimeSeries(args):
    """
    """

    # Load file and data
    filename = str(sys.argv[1])
    data     = np.loadtxt(filename)
    numData  = len(data)
    time     = np.arange(numData) * 25.0

    # Write new file with col1: time, col2: signal
    with open(filename, "w") as f:
        for i, x_i, y_i in zip(range(numData), time, data):
            line = "{:.1f} ".format(x_i) + '{:.20e}'.format(y_i) + "\n"
            f.write(line)
    print('A regular time-grid has been added to the file {0}'.format(filename))





def adjustMicroscanFile(args):


    cadence  = 25.  # [seconds]
    startQ   = 24
    filename_input  = os.getenv('PLATO_PROJECT_HOME') + '/inputfiles/spiral_8Hz_airbus_3h.txt' 

    filename_output = f'spiral_8Hz_airbus_3h_Q{startQ}.txt'

    data = np.loadtxt(filename_input)
    x, y, z = data[:,1], data[:,2], data[:,3]
    time_points = len(data)

    day = 86400.
    timeStart = round(90. * startQ * day)
    timeEnd   = round(timeStart + time_points * cadence)
    t = np.arange(timeStart, timeEnd, cadence)

    # Plot data
    plt.plot(t/day, x, 'b-', label='Yaw')
    plt.plot(t/day, y, 'm-', label='Pitch')
    plt.plot(t/day, z, 'g-', label='Roll')
    plt.xlabel('Time [days]')
    plt.ylabel('Jitter [arcsec]')
    plt.legend(loc='best')
    plt.tight_layout()
    plt.show()

    # Save data
    np.savetxt(filename_output, np.transpose([t,x,y,z]))



    
    
def combineFeathers(args):
    """
    Small utility to combine Pandas data frame tables in "pickle" (.pkl) format.
    This function takes only 1 input argument which is the file name of all the
    files that should be combined. Hence the files needs to have the same start
    prefix -> prefix**.pkl
    """

    # Fetch all data frames
    files = args.combineFeathers[0]

    # Load all tables
    df = []
    for i in range(len(files)):
        df_i = pd.read_feather(files[i])
        df_i['ID'] = f'{i+1}'.zfill(9)
        df.append(df_i)

    # Combine tables
    table = pd.concat(df)

    # Rearange ID to first col
    cols = table.columns.tolist()
    cols = cols[-1:] + cols[:-1]
    table = table[cols]

    # Feather needs proper indicing to be able to save
    indices = pd.Series(range(len(files)))
    table.set_index([indices], inplace=True, drop=True)

    # Save and print table
    table.to_feather('table.ftr')
    #table.to_csv('table.txt', float_format='%.3f')
    print(table)









def convertCSV2TXT(fileCSV, fileTXT, extension):

    # Paramters
    fileCSV = args.fileCSV
    fileTXT = args.fileTXT
    extension = args.extension
    
    # Open CSV file
    with open(fileCSV, newline='') as csv_file:

        # Skip comments
        csv_file = pandas.read_csv(csv_file, sep=';', comment='#')

        # Load data with pandas
        dataset = pandas.DataFrame(csv_file)
        data    = dataset.values

        # Check the number of columns
        if numCols is None:
            numCols = data.shape[1]

        # Select number columns to save
        data = data[:, :numCols]

    # Check filename
    if fileASCII is None:
        fileASCII = fileCSV[:-4]

    # Check file extension
    if extension is None:
        extension = 'txt'

    # Save ascii file
    fileNameOut = os.getcwd() + '/' + fileASCII + '.' + extension
    print('Saving ascii file {0}'.format(fileNameOut))
    np.savetxt(fileNameOut, data)





def generateDriftFile(args):
    """
    """

    # Parsed parameters
    filename  = args.generateDriftFile
    drifttype = args.generateDriftFile
    driftamp  = args.generateDriftFile  # [arcsec]
    cadence   = args.generateDriftFile  # [seconds]
    tdur      = args.generateDriftFile  # [days] 
    startQ    = args.generateDriftFile  # [quarter nr.]

    # Timings
    day = 86400.
    timeStart = round(90. * (startQ - 1) * day)
    timeEnd   = round(timeStart + tdur * day)
    t = np.arange(timeStart, timeEnd, cadence)

    # Model of yaw, pitch, roll
    #if drifttype == 'linear':
    x = np.linspace(0, drift, len(t))
    y = x
    z = np.zeros(len(t))

    # Plot result
    plt.plot(t/day, x, 'b-')
    plt.plot(t/day, y, 'r-')
    plt.show()

    # Save file
    np.savetxt(filename+'.txt', np.transpose([t,x,y,z]), fmt=['%i', '%f', '%f', '%f'])






def generateVarSourceIndex(args):

    # Load file
    ifile = '../sims/L1-pipeline/input/data_hpc_P5.txt'
    ofile = '../sims/L1-pipeline/input/data_starcat_hpc_P5.txt'
    a = np.zeros(7)
    b = np.ones(100)
    c = np.ones(445) * 2
    d = np.ones(445) * 3
    e = np.ones(3) * 4
    seed = 1234

    # Make array and shuffle it
    var = np.concatenate((a,b,c,d,e))
    rng = np.random.default_rng(seed)
    rng.shuffle(var)
    var = var.reshape(len(var), 1)
    
    data = np.loadtxt(ifile, skiprows=1, delimiter=',')
    data0 = np.append(data[:,:-1], var, axis=1)
    header = 'PIC, M, R, Teff, xsource'
    np.savetxt(ofile, data0, header=header, comments='', delimiter=',',
               fmt=['%i', '%0.3f', '%0.3f', '%i', '%i'])

    print(data0)
    print(np.shape(data0))







def plotLCfromAscii(args):
    
    # Command line arguments
    filename = args.plotLCfromAscii[0][0]

    # Load photometry class
    data = np.loadtxt(filename)
    time = data[:,0] / 86400           # [days]
    dmag = data[:,1]                   # [mag]
    flux = (10**(-dmag/2.5) - 1)* 1e6  # [ppm]
    
    # Plot the input variable source
    fig, ax = plt.subplots(1,1,figsize=(12,6))
    plt.plot(time, flux, 'k-', markersize=1)
    plt.xlabel('Time [days]')
    plt.ylabel(r'Relative flux, $f_P$')

    # Finito!

    plt.show()






def generateAOCSRedNoise(args):



    ifile = '/lhome/nicholas/software/workdir/L1-pipeline/output/pic13448892_Ncam1.1_Q23.hdf5'
    ofile = 'AOCS_4mHz_4mas_180d_Q23-24.txt'
    f = h5py.File(ifile, 'r')

    data_AOCS = np.transpose([f['ACS/Time'], f['ACS/Yaw'], f['ACS/Pitch'], f['ACS/Roll']])
    
    print(f['TransmissionEfficiency'])
    exit()

    np.savetxt(ofile, data_AOCS, fmt=['%.1f', '%.9f', '%.9f', '%.9f'])







    
def plotLCfromPipeline(args):
    
    # Command line arguments
    filename = args.plotLCfromPipeline[0][0]

    # Load photometry class
    data  = np.loadtxt(filename)
    flux0 = normalize(data[:,0])  # [ppm]
    flux1 = median_filter(flux0, ppmh)
    time = np.arange(len(flux0)) * 25 / 86400  # [day]

    # Input light curve
    #lc_input = filename[:-3]
    #data_input = np.loadtxt(filename)
    
    
    # Plot the input variable source
    fig, ax = plt.subplots(1,1,figsize=(15,6))
    plt.plot(time, flux0, 'b.', markersize=1)
    plt.plot(time, flux1, 'k-', label='Median per hour')
    plt.xlabel('Time [days]')
    plt.ylabel(r'Flux [e-]')
    plt.show()

    



def plotLCfromHDF5(args):
    
    filename = args.plotLCfromHDF5[0][0]
    # Fetch photometry
    photometryClass = PhotometricFile(filename)
    data = photometryClass.getPhotometricTimeSeries(1)
    time = data[0] / day2sec                  # [days]
    flux = normalize(data[2])                 # [ppm]
    maskupdates = (data[3] * 25.) / day2sec   # [days]

    # Stellar pixel coordinates
    #if COB is True:
    #    COB = photometryClass.getStellarPixelCoordinates(1)

    # Create figure and subplots
    #fig = plt.subplots(figsize=(5,9))
    #plt.subplots_adjust(wspace=0.15, hspace=0.20)

    # Plot photometric time series
    fig, ax = plt.subplots(1,1,figsize=(8,4))
    plt.plot(time, flux, '.', color='orange', markersize=1, alpha=0.2, label='Raw flux')
    plt.plot(time, median_filter(flux, ppmh), 'k-', label='Median per hour')
    axes_maskupdates(ax, time, maskupdates)
    # Labels
    plt.legend(loc='lower center', fancybox=True, ncol=2)
    plt.ylabel('Relative flux [ppm]')
    # Limits
    flux_min = np.min(flux)
    #plt.ylim(np.min(flux), np.abs(flux_min))
    plt.xlim(axes_minmax(x=time))
    plt.ylim(np.min(flux), 20000)
    # Finito!
    plt.show()


# Run functions

if args.generateVarSourceIndex: generateVarSourceIndex(args)
if args.combineFeathers: combineFeathers(args)
if args.plotLCfromPipeline: plotLCfromPipeline(args)
if args.plotLCfromAscii: plotLCfromAscii(args)
if args.plotLCfromHDF5: plotLCfromHDF5(args)
