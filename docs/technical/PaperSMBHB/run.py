import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from astropy.time import Time
from astropy import units as u
from astropy import constants as c
from ipywidgets import *
from pathlib import Path

# PlatoSim libraries
import platosim.utilities as ut
import platosim.plot      as pt
import platosim.noise     as ns
from platosim.varsource    import SMBHB
from platosim.lightcurve   import LightCurve
from platosim.slurm        import workerOverview
from platosim.matplotlibrc import setup_notebook
setup_notebook()

# Paths to where data is stored
path = Path(os.getenv('PLATO_WORKDIR')) / 'smbhb'
fdir = path / 'figures'

# Initialise model
dt   = 3600
tdur = 3 * ut.year()
time = np.arange(0, tdur, dt) * u.s
model = SMBHB(time, seed=123456789)

# Parameters
t0  = u.yr * 1.05
z  = 0.962
P  = u.yr * 1.144
M1 = 10**7.4 * u.M_sun                                                           
M2 = 10**6.7 * u.M_sun
i  = np.arccos(0.140) * u.rad
I  = np.pi/2 * u.rad - i
J  = np.pi / 4 * u.rad
e  = 0.524
w  = 1.477 * u.rad
alpha = 2.09
f_lum = 0.89

model.doppler_boosting(t0, z, P, M1, M2, i, e, w, alpha, f_lum, plot=True)
