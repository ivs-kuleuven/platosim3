import os
import warnings
import numpy as np
import matplotlib.pyplot as plt

from scipy import constants

from platosim.simulation import Simulation
from platosim.plot import drawCCDsInCameraFocalPlane

# Change default settings

#warnings.simplefilter("ignore")

plt.rc('xtick', labelsize=13) 
plt.rc('ytick', labelsize=13) 
plt.rcParams.update({'font.size': 13})
plt.rcParams['text.usetex'] = True

# LOAD AND CREATE USER DEFINED STAR CATALOGUE

# Loading PIC targets
pic_tar = np.loadtxt(os.getcwd() + '/starcat-SPF-P2-Targets.txt')
ID  = pic_tar[:,0].astype(int)
ra  = pic_tar[:,1]
dec = pic_tar[:,2]
mag = pic_tar[:,3]

# Loading PIC contaminants
pic_con = np.loadtxt(os.getcwd() + '/starcat-SPF-P2-Contaminants.txt')
ID_con   = pic_con[:,0].astype(int)
ra_con   = pic_con[:,1]
dec_con  = pic_con[:,2]
mag_con  = pic_con[:,3]

# Instrumental parameters

pixelSize   = 18        # [µm]
plateScale  = 15        # [arcsec]
fovGhost    = 8.0       # [deg]
focalLength = 247.52    # [mm]

# Calculate the maximum FOV for a point-like ghost creation

fovGhost = focalLength * np.tan(np.radians(fovGhost))  # [mm]

# START PLOT

fig = plt.figure(figsize=(14,10))

# Plotting

ax = drawCCDsInCameraFocalPlane(fig)
c = plt.Circle((0, 0), radius=fovGhost, color="none", linewidth=2, label=r"Ghost FOV $<8^{\circ}$", zorder=4)
c.set_edgecolor("m")
ax.add_patch(c)
ax.plot(-15, 15, '*', c='m', ms=15, label='Star position', zorder=5)
ax.plot(15, -15, 'o', c='m', ms=13, fillstyle="none", label='Ghost position', zorder=6)
ax.plot([-15, 15], [15, -15], '--', c='m', zorder=7)

# Settings

plt.legend(prop={'size':20}, bbox_to_anchor=(1.0,1.0))
ax.set_aspect('equal', 'box')
plt.tight_layout()

plt.show()

# Save figrue
fig.savefig(os.getcwd() + '/plotGhostInFocalPlane.png', dpi=300)

