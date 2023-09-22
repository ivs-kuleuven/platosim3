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
        self.platoInputDir  = self.path.joinpath(os.getenv('PLATO_PROJECT_HOME'), 'inputfiles')
        self.platoPythonDir = self.path
        

    def run(self, script, project, paramFile):

        """Function to run the parallelisation.
        """

        # Fetch parameterisation file
        params = pd.read_csv(paramFile).to_numpy()

        # Configure parallel computing
        
        with parallel_config(backend=self.backend, n_jobs=self.cpus):

            if script == 'platonium':
                S = params[:,0]
                G = params[:,1]
                C = params[:,2]
                Q = params[:,3]
                N = len(S)
                Parallel()(delayed(self.run_platonium)(i,N,S[i],G[i],C[i],Q[i],project)
                           for i in range(N))

            if script == 'varsim':
                Parallel()(delayed(function)(i,S[i],G[i],C[i],Q[i],project) for i in range(N))


                
            
    def run_platonium(self, i, N, S, G, C, Q, project):

        """Function to run the parallelisation.
        """

        print(f'Simulation {i}/{N} -> ({round(i/N*100,2)}%)', end='\r')
        PLATOnium = os.getenv('PLATO_PROJECT_HOME') + '/python/platosim/platonium/platonium.py'

        #if mode == 'normal':
        #    os.system(f'{PLATOnium} {S} {G} {C} {Q} --project {project} --nexp 1 -w -v 0')
        #elif mode == 'statistics':
        os.system(f'{PLATOnium} {S} {G} {C} {Q} --project {project} --statistics')

