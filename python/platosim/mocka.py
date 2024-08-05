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

# PlatoSim imports
import platosim.noise           as ns
import platosim.utilities       as ut
import platosim.referenceFrames as rf                                

#--------------------------------------------------------------#
#                         GRAPHICAL TOOLS                      #
#--------------------------------------------------------------#

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

    # Fetch input frequencies in pettern
    f_i = 1 / np.array([dp.DeltaP0_day * ((1 + dp.slope)**i - 1)/dp.slope + dp.P0_day 
                        for i in range(dp.N_modes[0])])

    # Find indices of matching frequencies
    dex_do = np.array([ut.findNearestIndex(do.freq, f_i[i]) for i in range(dp.N_modes[0])])
    dex_dm = np.array([ut.findNearestIndex(dm.freq, f_i[i]) for i in range(dp.N_modes[0])])
    
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
                       'A_max': dm.ampl.max(),
                       'A_limit_bic': dm_bic.ampl.min(),
                       'A_limit_snr': dm_snr.ampl.min(),
                       'N_input': dp.N_modes[0],
                       'N_bic': dm_bic.shape[0],
                       'N_snr': dm_snr.shape[0],
                       'f_rms_bic': ut.rootMeanSquare(f_oc_bic),
                       'f_rms_snr': ut.rootMeanSquare(f_oc_snr),
                       'A_rms_bic': ut.rootMeanSquare(A_oc_bic),
                       'A_rms_snr': ut.rootMeanSquare(A_oc_snr),
                      })
    
    return dx, df_oc_bic, df_oc_snr 





def fetch_all_modes(path, star='GDOR', batch='finals_affogato'):
    
    files_parameters = natsort.natsorted(glob.glob(f'{path}/{star}/varsource/parameters/*'))
    files_pulsations = natsort.natsorted(glob.glob(f'{path}/{star}/varsource/pulsations/*'))
    
    files_modes = natsort.natsorted(glob.glob(f'{path}/{star}/{batch}/modes/*'))
    files_table = natsort.natsorted(glob.glob(f'{path}/{star}/{batch}/table/*'))
    
    dx = pd.DataFrame()
    df_bic = pd.DataFrame()
    df_snr = pd.DataFrame()
    for i in tqdm(range(len(files_modes)), bar_format=ut.tqdmBar()):    
        dx0, df0_bic, df0_snr = match_modes(files_parameters[i], 
                                             files_pulsations[i], 
                                             files_modes[i], 
                                             files_table[i])
        dx = pd.concat([dx, dx0])
        df_bic = pd.concat([df_bic, df0_bic])
        df_snr = pd.concat([df_snr, df0_snr])
    
    # Sort rows
    dx = dx.sort_values(by=['ncam', 'Pmag'])

    return dx, df_bic, df_snr
