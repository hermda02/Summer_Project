# loads most recent Planck release data for map test

import matplotlib.pyplot as plt
import numpy as np
import healpy as hp

planck = np.loadtxt("COM_PowerSpect_CMB-base-plikHM-TT-lowTEB-minimum-theory_R2.02.txt")

ls = np.array(planck[:,0])
TT = np.array(planck[:,1])/(ls*(ls+1)/(2*np.pi))
TE = np.array(planck[:,2])/(ls*(ls+1)/(2*np.pi))
EE = np.array(planck[:,3])/(ls*(ls+1)/(2*np.pi))
BB = np.array(planck[:,4])/(ls*(ls+1)/(2*np.pi))
