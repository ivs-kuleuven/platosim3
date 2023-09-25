#!/usr/bin/env python3

"""
Run the PlatoSim functions in parallel.
"""

# Python standard
import os
import argparse
import datetime
from datetime import datetime

# PlatoSim standard
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

# PlatoSim functions
from platosim.utilities import errorcode


#==============================================================#
#                        SLURM UTILITIES                       #
#==============================================================#


def convertWorkerLog(workerLog):

    """Get a pandas data frame of the worker log.

    This function takes the worker log file in ascii format and convert
    it to a pandas data frame for easier data handling.
    """
    
    month_convert = {'Jan': 1, 'Feb': 2, 'Mar': 3, 'Apr': 4, 'May': 5, 'Jun': 6, 
                     'Jul': 7, 'Aug': 8, 'Sep': 9, 'Oct': 10, 'Nov': 11, 'Dec': 12}

    df = pd.read_csv(workerLog, sep=' ', 
                     names=['state', 'a', 'core', 'b', 'c',
                            'month', 'date', 'time', 'year'])
    df = df.drop(['a', 'b', 'c'], axis=1)
    df = df.reset_index()
    df = df.rename(columns={'index':'sim'})
    N = len(df.sim)
    df['datetime'] = [datetime(df.year.iloc[i],
                               month_convert[df.month.iloc[i]],
                               df.date.iloc[i],
                               int(df.time.iloc[i][0:2]),
                               int(df.time.iloc[i][3:5]),
                               int(df.time.iloc[i][6:8])) for i in range(N)]

    return df.drop(['month', 'date', 'time', 'year'], axis=1)




    
def workerOverview(workerLog, paramFile, ofile=False, plot=False):

    """Overview of the run time for worker using SLURM.
    """

    # Load data frame and manipulate worker log 

    df = convertWorkerLog(workerLog)

    # If requested, plot an overview of start and end times
    
    if plot:
        fig, ax = plt.subplots(1, 1, figsize=(8, 10))
        dex = []

        for i in df.sim.unique():

            df0 = df[df.sim == i]

            if df0.shape[0] == 2:
                ax.hlines(df0.sim.iloc[0], df0.datetime.iloc[0], df0.datetime.iloc[1], color='b', alpha=0.2)
                ax.plot(df0.datetime.iloc[0], df0.sim.iloc[0], 'blue', marker='>', mec='k', ms=6)
                ax.plot(df0.datetime.iloc[1], df0.sim.iloc[1], 'lime', marker='s', mec='k', ms=6)
            elif df0.shape[0] == 1:
                ax.plot(df0.datetime.iloc[0], df0.sim.iloc[0], 'r', marker='>', mec='k', ms=6)
                dex.append(df0.sim.iloc[0]-1)
                print(f'Simulation {df0.sim.iloc[0]} did not finish!')

        plt.title('Overview of worker')
        plt.xlabel('Date-time')
        plt.ylabel('Simulation no.')
        plt.xticks(rotation=30)
        plt.tight_layout()
        plt.show()

    # Check how many simulations that did not finish in time
    # Also save a parameterisation file in that case
        
    if dex:
        paramFile = Path(paramFile)
        df1 = pd.read_csv(paramFile)
        df1 = df1.loc[dex]
        errorcode('warning', f'{df1.shape[0]} simulations did not finish within the walltime!')
        if ofile:
            df1.to_csv(ofile)
    else:
        errorcode('module', 'All simulations finished successfully!')
