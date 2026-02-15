# Beware that this does not use the common cpp header 
# file, so if there is a change in it, this code has to be 
# updated manually.

import numpy as np
import matplotlib.pyplot as plt

kT_center = np.array([0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5, 0.55, 0.6, 0.65]);
energies = np.array([3.0, 3.2, 3.5, 3.9, 4.5, 7.7, 9.2, 11.5, 14.5, 19.6, 27.0]);

rmaxdef = 72.0

def rhofitmax(kT, energy, syst=0):
    ktparam = (kT-0.2) / 0.05
    systpm = 0.0
    if syst == 1:
        systpm = -1.0
    elif syst == 2:
        systpm = 1.0
    return (rmaxdef - ktparam * 5.0) * (1.0 + systpm * 0.6) * np.sqrt(energy/11.5)

plt.figure(figsize=(8,6))
plt.title(r'$\rho_{fit}^{max}$ vs $k_T$ for different energies')
plt.xlabel(r'$k_T$ (GeV)')
plt.ylabel(r'$\rho_{fit}^{max}$')
for energy in energies:
    plt.plot(kT_center, rhofitmax(kT_center, energy), 'o-', label=f'{energy} GeV')
plt.legend()
fig1 = plt.gcf()
fig1.savefig("figs/rhofitmax_vs_kt_snn.png")
plt.show()
