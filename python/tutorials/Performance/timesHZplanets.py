
import numpy as np
import matplotlib.pyplot as plt
from astropy import constants as c
from astropy import units as u
from prettytable import PrettyTable
plt.rcParams['text.usetex'] = True
plt.rcParams['text.latex.preamble'] = [r'\usepackage{amsmath}']

# Select range of temperatures

w = 0.
e = 0.
i = (90. * u.deg).to('rad').value

# Value for F0-M5 main-sequence dwarf stars
s = ['F0', 'F5', 'G0', 'G2', 'G5', 'K0', 'K5', 'M0', 'M5']
Teff = np.array([7240, 6420, 5920, 5780, 5610, 5240, 4410, 3800, 3120]) * u.K
L = np.array([6.0, 2.5, 1.26, 1.0, 0.79, 0.40, 0.16, 0.072, 0.0027]) * u.L_sun
M = np.array([1.7, 1.3, 1.10, 1.0, 0.93, 0.78, 0.69, 0.60,  0.15]) * u.M_sun
R = np.array([1.3, 1.2, 1.05, 1.0, 0.93, 0.85, 0.74, 0.51,  0.18]) * u.R_sun

def timeHZ(Rp):

    # Find effective stellar flux
    Teff_sun = 5780. * u.K
    T = (Teff - Teff_sun).value

    c_inn = [8.1774e-5, 1.7063e-9, -4.3241e-12, -6.6462e-16]
    c_out = [5.8942e-5, 1.6558e-9, -3.0045e-12, -5.2983e-16]

    Seff_inn = 1.0140 + c_inn[0]*T + c_inn[1]*T**2 + c_inn[2]*T**3 + c_inn[3]*T**4
    Seff_out = 0.3438 + c_out[0]*T + c_out[1]*T**2 + c_out[2]*T**3 + c_out[3]*T**4

    # Semi-major axis

    a_inn = ((L.value/Seff_inn)**0.5 * u.AU).to('m')   # [AU] => [m]
    a_out = ((L.value/Seff_out)**0.5 * u.AU).to('m')   # [AU] => [m]

    # Orbital period (K3)

    P_inn = np.sqrt(4 * np.pi**2 * a_inn**3 / (c.G * M.to('kg')))  # [s]
    P_out = np.sqrt(4 * np.pi**2 * a_out**3 / (c.G * M.to('kg')))  # [s]

    
    # fig = plt.figure(figsize=(12,5))
    # # Plot Teff vs. Seff
    # ax1 = fig.add_subplot(131)
    # ax1.plot(Seff_inn, Teff, 'r-')
    # ax1.plot(Seff_out, Teff, 'b-')
    # ax1.set_xlabel(r'$S_{eff}$ ($S_{\odot}$)')
    # ax1.set_ylabel(r'$T_{eff}$ (K)')
    # ax1.invert_xaxis()
    # # Plot Mass vs. Distance
    # ax2 = fig.add_subplot(132)
    # ax2.plot(a_inn.to('AU'), M, 'r-')
    # ax2.plot(a_out.to('AU'), M, 'b-')
    # ax2.set_xscale('log')
    # ax2.set_yscale('log')
    # ax2.set_xlabel(r'Distance, $D$ (AU)')
    # ax2.set_ylabel(r'Stellar mass, $M$ ($M_{\odot}$)')
    # # Plot Mass vs. Distance
    # ax3 = fig.add_subplot(133)
    # ax3.plot(P_inn.to('d'), Teff, 'r-')
    # ax3.plot(P_out.to('d'), Teff, 'b-')
    # ax3.set_xlabel(r'Orbital Period, $P$ (days)')
    # ax3.set_ylabel(r'$T_{eff}$ (K)')
    # # Plot
    # plt.tight_layout()
    # plt.show()
    # exit()

    # Impact parameter: Winn (2014) Eq. 7 & 8

    b_tra_inn = 0#a_inn.to('m')*np.cos(i)/R.to('m') * (1 - e**2)/(1 + e*np.sin(w))
    b_tra_out = 0#a_inn.to('m')*np.cos(i)/R.to('m') * (1 - e**2)/(1 + e*np.sin(w))

    #b_occ_inn = a_out*np.cos(i)/R * (1 - e**2)/(1 - e*np.sin(w))
    #b_occ_out = a_out*np.cos(i)/R * (1 - e**2)/(1 - e*np.sin(w))

    # Transit depth first approximation

    k = Rp.to('R_sun')/R.to('R_sun')

    # Transit times: Winn (2014) Eq. 14, 15 & 16
    # NOTE on circular orbits the transit and occultation times are equal

    x  = np.sqrt(1 - e**2)   # Optimization constant
    e_tra = x/(1 + e*np.sin(w))
    e_occ = x/(1 - e*np.sin(w))

    t_tra_tot_inn = P_inn.to('s')/np.pi * np.arcsin( R.to('m')/a_inn.to('m') * np.sqrt((1 + k)**2 - b_tra_inn**2)/np.sin(i) ) * e_tra / u.rad
    t_tra_tot_out = P_out.to('s')/np.pi * np.arcsin( R.to('m')/a_out.to('m') * np.sqrt((1 + k)**2 - b_tra_out**2)/np.sin(i) ) * e_tra / u.rad

    t_tra_ful_inn = P_inn.to('s')/np.pi * np.arcsin( R.to('m')/a_inn.to('m') * np.sqrt((1 - k)**2 - b_tra_inn**2)/np.sin(i) ) * e_tra / u.rad
    t_tra_ful_out = P_out.to('s')/np.pi * np.arcsin( R.to('m')/a_out.to('m') * np.sqrt((1 - k)**2 - b_tra_out**2)/np.sin(i) ) * e_tra / u.rad

    tau_tra_inn = (t_tra_tot_inn - t_tra_ful_inn)/2.
    tau_tra_out = (t_tra_tot_out - t_tra_ful_out)/2.

    # Occultation times

    t_occ_tot_inn = P_inn.to('s')/np.pi * np.arcsin( R.to('m')/a_inn.to('m') * np.sqrt((1 + k)**2 - b_tra_inn**2)/np.sin(i) ) * e_occ / u.rad
    t_occ_tot_out = P_out.to('s')/np.pi * np.arcsin( R.to('m')/a_out.to('m') * np.sqrt((1 + k)**2 - b_tra_out**2)/np.sin(i) ) * e_occ / u.rad

    t_occ_ful_inn = P_inn.to('s')/np.pi * np.arcsin( R.to('m')/a_inn.to('m') * np.sqrt((1 - k)**2 - b_tra_inn**2)/np.sin(i) ) * e_occ / u.rad
    t_occ_ful_out = P_out.to('s')/np.pi * np.arcsin( R.to('m')/a_out.to('m') * np.sqrt((1 - k)**2 - b_tra_out**2)/np.sin(i) ) * e_occ / u.rad

    tau_occ_inn = (t_occ_tot_inn - t_occ_ful_inn)/2.
    tau_occ_out = (t_occ_tot_out - t_occ_ful_out)/2.

    # Finito!

    return ([t_tra_tot_inn, t_tra_tot_out, tau_tra_inn, tau_tra_out],
            [t_occ_tot_inn, t_occ_tot_out, tau_occ_inn, tau_occ_out],
            [P_inn, P_out])



# Compute HZ for different planet radii

#mars = timeHZ(0.5*u.R_earth)
earth = timeHZ(1.0*u.R_earth)
neptune = timeHZ(5.0*u.R_earth)
jupiter = timeHZ(10.0*u.R_earth)

# Make a pretty table for an Earth analog of transit times

p = jupiter
t = PrettyTable(['Type', 'T_tra1 [h]', 'T_tra2 [h]', 'tau_tra1 [min]', 'tau_tra2 [min]', 'P_inn [days]', 'P_out [days]'])
for i in range(len(s)):
    t.add_row([str(s[i]),
               '%.2f' % p[0][0].to('h').value[i], '%.2f' % p[0][1].to('h').value[i],
               '%.2f' % p[0][2].to('min').value[i], '%.2f' % p[0][3].to('min').value[i],
               '%.2f' % p[2][0].to('d').value[i], '%.2f' % p[2][1].to('d').value[i]])

print(t)
exit()
# Plot figure showing results

fig = plt.figure(figsize=(4.3,5))
ax1 = fig.add_subplot(111)

for i in range(len(Teff)):
    ax1.hlines(y=Teff.value[i], xmin=0, xmax=27, color='k', linestyle=':', lw=1)
    ax1.text(27.2, Teff.value[i]-50, s[i], fontsize=12)

#ax1.plot(mars[0][0].to('h'), Teff, 'r-', label=r'$0.5 R_{\oplus}$')
#ax1.plot(mars[0][1].to('h'), Teff, 'r-')

ax1.plot(earth[0][0].to('h'), Teff, 'b-', label=r'$1 R_{\oplus}$')
ax1.plot(earth[0][1].to('h'), Teff, 'b-')

ax1.plot(neptune[0][0].to('h'), Teff, 'c-', label=r'$5 R_{\oplus}$')
ax1.plot(neptune[0][1].to('h'), Teff, 'c-')

ax1.plot(jupiter[0][0].to('h'), Teff, '-', c='orange', label=r'$10 R_{\oplus}$')
ax1.plot(jupiter[0][1].to('h'), Teff, '-', c='orange')

ax1.set_xlabel('Eclipse times (h)')
ax1.set_ylabel(r'$T_{\text{eff}}$ (K)')

ax1.set_xlim(0,30)
ax1.legend()
plt.tight_layout()
plt.show()
