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

# PlatoSim imports
from platosim.utilities import errorcode

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
                N, M = params.shape
                Parallel()(delayed(self.run_platonium)
                           (i, N, M, params, project)
                           for i in range(N))

            if script == 'varsim':
                S = params[:,0].astype(int)
                N = len(S)
                odir = f'{projectDir}/varsource'
                Parallel()(delayed(self.run_varsim)
                           (i, N, S[i], odir)
                           for i in range(N))

            # Print when script is done
            errorcode('message', 'All simulations finished successfully!')



            
    #def set_compress(self, compress=False): return compress
    #def set_nexp(self,    nexp    = ''):  self.nexp     = nexp
    def set_seed(self,    seed    = ''):  self.seed     = seed
    def set_cadence(self, cadence = ''):  self.cadence  = cadence
    def set_vardir(self,  vardir  = ''):  self.vardir   = vardir

        
    def run_platonium(self, i, N, M, params, project):

        """Function to run the PLATOnium in parallel.
        """
        
        print(f'Simulation {i+1}/{N} -> ({round(i/N*100,2)}%)', end='\r')
        platonium = os.getenv('PLATO_PROJECT_HOME') + '/python/platosim/platonium/platonium.py'
        
        S = int(params[i,0])
        G = int(params[i,1])
        C = int(params[i,2])
        Q = int(params[i,3])

        # Nine digit star ID 
        starID = f'{S}'.zfill(9)
        odir   = os.getenv('PLATO_WORKDIR') + f'/{project}/output/{starID}'
        # Parse flags
        #compress = self.set_compress()

        # Parse arguments
        #if self.set_nexp    != '': nexp     = f'--nexp {self.nexp}'
        if self.set_seed    != '': seed     = f'--seed {self.seed}'
        if self.set_cadence != '': cadence  = f'--cadence {self.cadence}'
        if self.set_vardir  != '': varfile  = f'--varfile {self.vardir}/varsource_{starID}.txt'

        # Run PlatoSim simulation
        os.system(f'{platonium} {S} {G} {C} {Q} --project {project} -o {odir} ' +
                  f'{seed} {cadence} {varfile} ' +
                  '--compress -v 0 -w')




        
    def run_varsim(self, i, N, S, odir):

        """Function to run the parallelisation.
        """

        print(f'Simulation {i+1}/{N} -> ({round(i/N*100,2)}%)', end='\r')
        VARSIM = os.getenv('PLATO_PROJECT_HOME') + '/python/platosim/platonium/varsim.py'
        starID = f'{i+1}'.zfill(9)
        ofile  = f'{odir}/varsource_{starID}.txt'
        os.system(f'{VARSIM} --star SMBHB --quarter 1-8 -o {ofile} -v 0')
        
