# Healpy test file

import matplotlib.pyplot as plt
import numpy as np
import healpy as hp

planck = np.loadtxt("COM_PowerSpect_CMB-base-plikHM-TT-lowTEB-minimum-theory_R2.02.txt")

# ls = []
# TT = []

# ls = 0
# TT = 0

ls = np.array(planck[:,0])
TT = np.array(planck[:,1])/(ls*(ls+1)/(2*np.pi))
TE = np.array(planck[:,2])/(ls*(ls+1)/(2*np.pi))
EE = np.array(planck[:,3])/(ls*(ls+1)/(2*np.pi))
BB = np.array(planck[:,4])/(ls*(ls+1)/(2*np.pi))

TTmap = hp.sphtfunc.synfast(TT,2048,lmax=500)
hp.mollview(TTmap)
plt.show()

cltest = hp.sphtfunc.anafast(TTmap,lmax=500)
plt.plot(cltest)
plt.plot(TT)
plt.xmax=500
plt.show()

# pols = [TT,EE,BB,TE]

# polmap = hp.sphtfunc.synfast(pols,2048,new=True)

# hp.mollview(polmap[0])
# plt.show()

# hp.mollview(polmap[1])
# plt.show()

# hp.mollview(polmap[2])
# plt.show()