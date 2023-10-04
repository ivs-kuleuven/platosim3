#!/usr/bin/env python3

"""
Run the PlatoSim functions in parallel.
"""

# Python standard
import os
import shutil
import random
import argparse
import datetime
from joblib import Parallel, delayed, parallel_config

# PlatoSim standard
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

# Disable warnings
#warnings.simplefilter("ignore")



#==============================================================#
#                        PLATOnium CLASS                       #
#==============================================================#

class HPC(object):

    """Class for parallisation of PLATOnium functions.
    """
    
    def __init__(self, cpus=8, backend='threading'):

         # PATHS AND INPUT

        # Global variables
        self.cpus    = cpus
        self.backend = backend

        # Absolute pwd path
        self.path = Path(__file__).parent.resolve()
        self.platoInputDir   = self.path.joinpath(os.getenv('PLATO_PROJECT_HOME'), 'inputfiles')





    def run(self, script, project, paramFile):

        """Function to run the parallelisation.
        """

        # Fetch parameterisation file
        projectDir = self.path.joinpath(os.getenv('PLATO_WORKDIR'), project)
        params     = pd.read_csv(paramFile).to_numpy()

        # Configure parallel computing
        
        with parallel_config(backend=self.backend, n_jobs=self.cpus):

            if script == 'platonium':
                S = params[:,0]
                G = params[:,1]
                C = params[:,2]
                Q = params[:,3]
                N = len(S)
                Parallel()(delayed(self.run_platonium)
                           (i, N, S[i], G[i], C[i], Q[i], project)
                           for i in range(N))

            if script == 'varsim':
                S = params[:,0].astype(int)
                N = len(S)
                odir = f'{projectDir}/varsource'
                Parallel()(delayed(self.run_varsim)
                           (i, N, S[i], odir)
                           for i in range(N))

            # Print when script is done
            errorcode('message', 'Parallel simulations are done!')


                
                
            
    def run_platonium(self, i, N, S, G, C, Q, project):

        """Function to run the parallelisation.
        """
        
        print(f'Simulation {i+1}/{N} -> ({round(i/N*100,2)}%)', end='\r')
        PLATOnium = os.getenv('PLATO_PROJECT_HOME') + '/python/platosim/platonium/platonium.py'

        #if mode == 'normal':
        #    os.system(f'{PLATOnium} {S} {G} {C} {Q} --project {project} --nexp 1 -w -v 0')
        #elif mode == 'statistics':
        os.system(f'{PLATOnium} {S} {G} {C} {Q} --project {project} --statistics')





    def run_varsim(self, i, N, S, odir):

        """Function to run the parallelisation.
        """

        print(f'Simulation {i+1}/{N} -> ({round(i/N*100,2)}%)', end='\r')
        VARSIM = os.getenv('PLATO_PROJECT_HOME') + '/python/platosim/platonium/varsim.py'
        starID = f'{i+1}'.zfill(9)
        ofile  = f'{odir}/varsource_{starID}.txt'
        os.system(f'{VARSIM} --star SMBHB --quarter 1-8 -o {ofile} -v 0')
        
