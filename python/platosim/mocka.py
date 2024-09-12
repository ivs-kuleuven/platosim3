#!/usr/bin/env python3

"""
This python module contains plot utilities used in the minimal 
PlatoSim installation and in the extra PLATOnium installation.
"""

# Built-in
import os
import glob

# PlatoSim standard
import natsort
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from tqdm import tqdm
from pathlib import Path

# PlatoSim imports
import platosim.noise           as ns
import platosim.utilities       as ut
import platosim.referenceFrames as rf                                


#--------------------------------------------------------------#
#                         GRAPHICAL TOOLS                      #
#--------------------------------------------------------------#


def create_timeseries(ids, idir, odir, power=2.2):
    
    """Function create time series from the pulsation modes.
    """

    # Time array of 2 years
    duration = ut.quarter() * 8
    time_sec = np.arange(0, duration * 86400, 25)
    time_day = np.arange(0, duration, 25 / 86400)

    for i in tqdm(ids, bar_format=ut.tqdmBar()):

        # Fetch simulation table
        starID = f'{i}'.zfill(9)

        # Create varsource from pulsations
        dx = pd.read_feather(f'{idir}/pulsations_{starID}_001.ftr')
        dv = pd.DataFrame()
        dv['time'] = time_sec
        dv['dmag'] = ns.timeSeriesFromFourier(time_day, dx.freq, dx.ampl, dx.phase, power=power)

        # Save light curve
        output_dir = Path(f'{odir}/{starID}')
        output_dir.mkdir(parents=True, exist_ok=True)
        ofile = f'{output_dir}/varsource_001.txt'
        data = np.transpose([dv.time, dv.dmag])
        np.savetxt(ofile, data, fmt=['%.1f', '%.8f'])
        os.system(f'chmod 755 {ofile}')




        
def fetch_amplitude_correction(path):
    
    """Function to fetch the amplitudes before and after the passband correction
    """
    
    folders_parameters = natsort.natsorted(glob.glob(f'{path}/parameters/*'))
    folders_pulsations = natsort.natsorted(glob.glob(f'{path}/pulsations/*'))

    # Load all pulsation modes
    N = len(folders_parameters)
    N_modes = np.zeros(N)
    PC = np.zeros(N)
    df = pd.DataFrame()
    for m,i,j in zip(range(N), folders_parameters, tqdm(folders_pulsations, bar_format=ut.tqdmBar())):    

        # Parameter file
        df_parameters = pd.read_feather(i)
        PC[m] = df_parameters.PC_kepler

        # Pulsation file
        df_pulsations = pd.read_feather(j)
        n = df_pulsations.shape[0]
        N_modes[m] = n
        df_pulsations['A_PC'] = df_pulsations.ampl / PC[m]

        # Save only pulsation modes
        df = pd.concat([df, df_pulsations[:n]])

    # Convert to mmag
    A_sim = df.ampl * 1e3
    A_PC  = df.A_PC * 1e3
    
    return df, PC, A_sim, A_PC




def match_modes(file_parameters, file_pulsations, file_modes, file_table):

    """Function to match input with output modes detected above BIC/SNR threshold.
    """

    # Load files from for varsource
    dp = pd.read_feather(file_parameters)
    do = pd.read_feather(file_pulsations)
    
    # Load results from simulations
    dt = pd.read_feather(file_table)
    dm = pd.read_feather(file_modes)
    
    # Correct for gamma factor
    do.ampl /= 2.2

    # Convert amplitudes [dmag -> ppm]
    do.ampl = (1 - ut.fromMagToFlux(do.ampl)) * 1e6

    # All modes passed SNR
    dm_all_snr = dm[dm.passed_snr].reset_index(drop=True)
    
    # Fetch input frequencies in pettern [c/d]
    f_i = 1 / np.array([dp.DeltaP0_day * ((1 + dp.slope)**i - 1)/dp.slope + dp.P0_day 
                        for i in range(dp.N_modes[0])])

    # Find indices of matching frequencies
    dex_do = np.array([ut.findNearestIndex(do.freq, f_i[i]) for i in range(dp.N_modes[0])])
    dex_dm = np.array([ut.findNearestIndex(dm.freq, f_i[i]) for i in range(dp.N_modes[0])])

    # print(file_parameters)
    # print(file_pulsations)
    # print(file_table)
    # print(file_modes)    
    # print(dex_do)
    # print(dex_dm)
    # exit()
    
    # Get pattern passed BIC criterion
    do_bic = do.loc[dex_do].reset_index(drop=True)
    dm_bic = dm.loc[dex_dm].reset_index(drop=True)
    
    # Get pattern passed SNR criterion
    do_snr = do_bic[dm_bic.passed_snr].reset_index(drop=True)
    dm_snr = dm_bic[dm_bic.passed_snr].reset_index(drop=True)

    # Compute O-C values [ppm]
    f_oc_bic = (dm_bic.freq.to_numpy() - do_bic.freq.to_numpy())
    f_oc_snr = (dm_snr.freq.to_numpy() - do_snr.freq.to_numpy())
    A_oc_bic = (dm_bic.ampl.to_numpy() - do_bic.ampl.to_numpy())
    A_oc_snr = (dm_snr.ampl.to_numpy() - do_snr.ampl.to_numpy())

    # Remove matches above O-C threshold
    x = 0.0005
    dex_bic = np.where((np.abs(f_oc_bic) > x))[0]
    dex_snr = np.where((np.abs(f_oc_snr) > x))[0]
    dm_bic = dm_bic.drop(index=dex_bic)
    dm_snr = dm_snr.drop(index=dex_snr)
    f_oc_bic = np.delete(f_oc_bic, dex_bic) * 1e6
    f_oc_snr = np.delete(f_oc_snr, dex_snr) * 1e6
    A_oc_bic = np.delete(A_oc_bic, dex_bic)
    A_oc_snr = np.delete(A_oc_snr, dex_snr)

    # Store values into data frame
    df_oc_bic = pd.DataFrame({'freq':dm_bic.freq, 'ampl':dm_bic.ampl,
                              'freq_oc':f_oc_bic, 'ampl_oc':A_oc_bic})
    df_oc_snr = pd.DataFrame({'freq':dm_snr.freq, 'ampl':dm_snr.ampl,
                              'freq_oc':f_oc_snr, 'ampl_oc':A_oc_snr})
    
    # Store in data frame
    dx = pd.DataFrame({'Pmag': dp.Pmag,
                       'ncam': int(dt.shape[0]/8),
                       'rOA': dt.rOA.mean(),
                       'SPR': dt.SPR.mean(),
                       'N_bic': dm_bic.shape[0],
                       'N_snr': dm_snr.shape[0],
                       'N_input': dp.N_modes[0],
                       'A_max': dm.ampl.max(),
                       'A_limit_bic': dm.ampl.min(), 
                       'A_limit_snr': dm_all_snr.ampl.min(),
                       'A_ampl_bic': dm_bic.ampl.min(),
                       'A_ampl_snr': dm_snr.ampl.min(),
                       'f_rms_bic': ut.rootMeanSquare(f_oc_bic),
                       'f_rms_snr': ut.rootMeanSquare(f_oc_snr),
                       'A_rms_bic': ut.rootMeanSquare(A_oc_bic),
                       'A_rms_snr': ut.rootMeanSquare(A_oc_snr),
                      })
    
    return dx, df_oc_bic, df_oc_snr 





def fetch_all_modes(path, star='GDOR', batch='finals_affogato', old=False):
    
    # files_parameters = natsort.natsorted(glob.glob(f'{path}/{star}/varsource/parameters/*'))
    # files_pulsations = natsort.natsorted(glob.glob(f'{path}/{star}/varsource/pulsations/*'))
    
    # files_table = natsort.natsorted(glob.glob(f'{path}/{star}/{batch}/table/*'))

    files_modes = natsort.natsorted(glob.glob(f'{path}/{star}/{batch}/modes/*'))
    names_stars = [Path(files_modes[i]).stem[6:] for i in range(len(files_modes))]

    path_parameters = f'{path}/{star}/varsource/parameters'
    path_pulsations = f'{path}/{star}/varsource/pulsations'
    path_modes = f'{path}/{star}/{batch}/modes'
    path_table = f'{path}/{star}/{batch}/table'
    
    dx = pd.DataFrame()
    df_bic = pd.DataFrame()
    df_snr = pd.DataFrame()

    for i in tqdm(names_stars, bar_format=ut.tqdmBar()):

        if old:
            file_parameters = f'{path_parameters}/{i}_parameters.ftr'
            file_pulsations = f'{path_pulsations}/{i}_pulsations.ftr'
        else:
            file_parameters = f'{path_parameters}/parameters_{i}_001.ftr'
            file_pulsations = f'{path_pulsations}/pulsations_{i}_001.ftr'
            
        try:
            dx0, df0_bic, df0_snr = match_modes(file_parameters,
                                                file_pulsations,
                                                f'{path_modes}/modes_{i}.ftr',
                                                f'{path_table}/table_{i}.ftr')
        except:
            print(i)
            pass
        else:
            dx = pd.concat([dx, dx0])
            df_bic = pd.concat([df_bic, df0_bic])
            df_snr = pd.concat([df_snr, df0_snr])
    
    # Sort rows
    try:
        dx = dx.sort_values(by=['ncam', 'Pmag'])
    except: pass

    return dx, df_bic, df_snr






def plot_percentage_heatmap(dx, mag_min=7.5, mag_max=16.5, col='limit', cmap='rainbow',
                            title=None, clab=None, clim=None, cext=None):

    fig, ax = plt.subplots(1, 1, figsize=(9,5))

    # Plot a detection matrix
    P_bins = np.arange(mag_min, mag_max+1, 1) 
    ncams = [6, 12, 18, 24]
    xlab = P_bins[:-1] + 0.5
    ylab = ncams 
    
    # Create data to be plotted
    if isinstance(dx, pd.DataFrame):
        data = np.zeros((4, len(P_bins)-1))
        for i,n in zip(range(4), ncams):
            dx0 = dx[dx.ncam == n]
            data0 = [dx0[dx0["Pmag"].between(P_bins[i], P_bins[i+1])] for i in range(len(P_bins)-1)]
            if col == 'limit':
                Ebins = [data0[i]["N_snr"].mean() / data0[i]["N_input"].mean() * 100 for i in range(len(data0))]
            elif col == 'N_bic-N_snr':
                Ebins = [data0[i]["N_bic"].mean() - data0[i]["N_snr"].mean() for i in range(len(data0))]
            elif col == 'N_bic-N_snr/N_input':
                Ebins = [(np.mean(data0[i]["N_bic"]/data0[i]["N_input"]) - np.mean(data0[i]["N_snr"]/data0[i]["N_input"])) * 100 for i in range(len(data0))]
            else:
                Ebins = [data0[i][col].mean() for i in range(len(data0))]
            data[i,:] = np.array([Ebins])
    else:
        data = dx

    # Plot data
    if clim:
        im = ax.imshow(data, cmap=cmap, vmin=clim[0], vmax=clim[1])
    else:
        im = ax.imshow(data, cmap=cmap)
    for i in range(len(ylab)):
        for j in range(len(xlab)):
            if col in ['limit', 'phi', 'A_rms_snr', 'f_rms_snr']:
                text = ax.text(j, i, f'{data[i, j]:.1f}', ha="center", va="center", color="k")
            else:
                text = ax.text(j, i, f'{data[i, j]:.2f}', ha="center", va="center", color="k")

    ax.set_xticks(np.arange(len(xlab)), labels=xlab)
    ax.set_yticks(np.arange(len(ylab)), labels=ylab)

    ax.set_xticks(np.arange(data.shape[1] + 1) - 0.5, minor = True)
    ax.set_yticks(np.arange(data.shape[0] + 1) - 0.5, minor = True)
    ax.grid(which='minor', color='k', linestyle='-', linewidth=1) 
    ax.tick_params(which =  "minor", bottom = False, left = False)

    if clim == [0, 100]:
        ticks = [0, 25, 50, 75, 100]
    elif clim == [0, 90]:
        ticks = [0, 30, 60, 90]
    elif clim == [0, 55]:
        ticks = [0, 10, 20, 30, 40, 50]        
    elif clim == [0, 30]:
        ticks = [0, 10, 20, 30]
    elif clim == [40, 95]:
        ticks = [40, 50, 60, 70, 80, 90]        
    else:
        ticks = None

    cbar = plt.colorbar(im, extend=cext, pad=0.02, shrink=0.65, orientation='vertical', ticks=ticks)
    
    if clab:
        cbar.set_label(clab)
        
    ax.set_xlabel('PLATO magnitude, $\mathcal{P}$')
    ax.set_ylabel(r'Camera visibility, $n_{\rm CAM}$')
    
    if title: 
        ax.set_title(title, pad=15)
                
    plt.tight_layout()
    
    return fig, ax, data
