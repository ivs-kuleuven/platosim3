#!/usr/bin/env Python3

import ligo.skymap.plot
from matplotlib import pyplot as plt


from astropy.coordinates import SkyCoord
from astropy.io import fits
from astropy import units as u

import ligo.skymap.plot
from matplotlib import pyplot as plt


raPointing = 86.79
decPointing = -46.39

pointing = SkyCoord(raPointing, decPointing, unit='deg')  # Defaults to ICRS

# Plot

ax = plt.axes(projection='astro zoom', center=pointing, radius='10 deg', rotate='0 deg')
ax.grid()

ax.plot(pointing, 'r*')


ax.set_xlabel('x')
ax.set_ylabel('y')

plt.show()
