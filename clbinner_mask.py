# For use after running polspicetest.py to initialize arrays (such as computed C_ls)

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import healpy as hp

polspicetest = np.loadtxt('polspicetest_mask.txt')
planck = np.loadtxt("COM_PowerSpect_CMB-base-plikHM-TT-lowTEB-minimum-theory_R2.02.txt")
fsky = hp.read_map('mask_latlon_512.fits')

ls = np.array(planck[:,0])
TT = np.array(planck[:,1])
TE = np.array(planck[:,2])
EE = np.array(planck[:,3])
BB = np.array(planck[:,4])

newls = polspicetest[2:,0]
newTT = polspicetest[2:,1]
newEE = polspicetest[2:,2]
newBB = polspicetest[2:,3]
newTE = polspicetest[2:,4]

lmax = max(newls)

ls    = ls[:int(lmax)-1]
TT    = TT[:int(lmax)-1]
TE    = TE[:int(lmax)-1]
BB    = BB[:int(lmax)-1]
EE    = EE[:int(lmax)-1]

for i in range(int(lmax-1)):
	newTT[i] = newTT[i]*(newls[i]*(newls[i]+1))/(2*np.pi)
	newEE[i] = newEE[i]*(newls[i]*(newls[i]+1))/(2*np.pi)
	newBB[i] = newBB[i]*(newls[i]*(newls[i]+1))/(2*np.pi)
	newTE[i] = newTE[i]*(newls[i]*(newls[i]+1))/(2*np.pi)

binl   = np.empty(int(lmax/25))
binsTT = np.empty(int(lmax/25))
binTT  = np.empty(int(lmax/25))
declTT = np.empty(int(lmax/25))
binsEE = np.empty(int(lmax/25))
binEE  = np.empty(int(lmax/25))
declEE = np.empty(int(lmax/25))
binsBB = np.empty(int(lmax/25))
binBB  = np.empty(int(lmax/25))
declBB = np.empty(int(lmax/25))
binsTE = np.empty(int(lmax/25))
binTE  = np.empty(int(lmax/25))
declTE = np.empty(int(lmax/25))

yes = 0.
no  = 0.

for i in range(len(fsky)):
	if fsky[i] == 0.:
		no = no + 1.
	else:
		yes = yes + 1.

frac = yes / (no + yes)

for i in range(int(lmax/25)):
	binl[i]   = 12 + i*25
	binTT[i]  = np.sum(TT[i*25:24+i*25])/25
	binsTT[i] = np.sum(newTT[i*25:24+i*25])/25
	declTT[i] = np.sqrt(2/((2*binl[i]+1)*frac))*binsTT[i]

for i in range(int(lmax/25)):
	binEE[i]  = np.sum(EE[i*25:24+i*25])/25
	binsEE[i] = np.sum(newEE[i*25:24+i*25])/25
	declEE[i] = np.sqrt(2/((2*binl[i]+1)*frac))*binsEE[i]

for i in range(int(lmax/25)):
	binBB[i]  = np.sum(BB[i*25:24+i*25])/25
	binsBB[i] = np.sum(newBB[i*25:24+i*25])/25
	declBB[i] = np.sqrt(2/((2*binl[i]+1)*frac))*binsBB[i]

for i in range(int(lmax/25)):
	binTE[i]  = np.sum(TE[i*25:24+i*25])/25
	binsTE[i] = np.sum(newTE[i*25:24+i*25])/25
	declTE[i] = np.sqrt(2/((2*binl[i]+1)*frac))*binsTE[i]

plt.errorbar(binl,binsTT,yerr=declTT,color='grey')
plt.plot(binl,binTT)
plt.xlim(0,lmax)
plt.xlabel(r'$l$')
plt.ylabel(r'$l(l+1)C_l/2\pi$')
plt.title('PolSpice TT Power Spectrum with error')
plt.legend(loc='best')
plt.show()

plt.errorbar(binl,binsEE,yerr=declEE,color='grey')
plt.plot(binl,binEE)
plt.xlim(0,lmax)
plt.xlabel(r'$l$')
plt.ylabel(r'$l(l+1)C_l/2\pi$')
plt.title('PolSpice EE Power Spectrum with error')
plt.legend(loc='best')
plt.show()

plt.errorbar(binl,binsBB,yerr=declBB,color='grey')
plt.plot(binl,binBB)
plt.xlim(0,lmax)
plt.xlabel(r'$l$')
plt.ylabel(r'$l(l+1)C_l/2\pi$')
plt.title('PolSpice BB Power Spectrum with error')
plt.legend(loc='best')
plt.show()

plt.errorbar(binl,binsTE,yerr=declTE,color='grey')
plt.plot(binl,binTE)
plt.xlim(0,lmax)
plt.xlabel(r'$l$')
plt.ylabel(r'$l(l+1)C_l/2\pi$')
plt.title('PolSpice TE Power Spectrum with error')
plt.legend(loc='best')
plt.show()