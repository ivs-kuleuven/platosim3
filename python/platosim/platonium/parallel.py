#!/usr/bin/env python3

# Python standard
import os
import shutil
import random
import argparse
import datetime

# PlatoSim standard
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

# Disable warnings
#warnings.simplefilter("ignore")

from joblib import Parallel, delayed, parallel_config

#------------------------
cpus    = 8
backend = 'threading'
project = 'test'
#------------------------

df = pd.read_csv('/lhome/nicholas/software/workdir/test/input/cluster.data')
params = df.to_numpy()
S = params[:,0]
G = params[:,1]
C = params[:,2]
Q = params[:,3]
N = len(S)

def run_platonium(i, S,G,C,Q):
    print(f'Simulation {i}', end='\r')
    PLATOnium = os.getenv('PLATO_PROJECT_HOME') + '/python/platosim/platonium/platonium'
    os.system(f'{PLATOnium} {S} {G} {C} {Q} --project {project} --nexp 1 -w -v 0')

with parallel_config(backend=backend, n_jobs=cpus):
    Parallel()(delayed(run_platonium)(i,S[i],G[i],C[i],Q[i]) for i in range(N))
