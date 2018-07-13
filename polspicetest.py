import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import healpy as hp

from scipy.stats import chisquare

mpl.rcParams['text.usetex'] = True

polspicetest = np.loadtxt('polspicetest_mask.txt')
planck = np.loadtxt("COM_PowerSpect_CMB-base-plikHM-TT-lowTEB-minimum-theory_R2.02.txt")

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

print('lmax = ', lmax)

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

chiTT = chisquare(newTT,TT)
chiEE = chisquare(newEE,EE)
chiBB = chisquare(newBB,BB)
chiTE = chisquare(newTE,TE)

print('TT Chi square results:', chiTT)

plt.plot(ls,TT,label='Planck TT')
plt.plot(newls,newTT, label='PolSpice TT')
plt.xlim(0,lmax)
plt.xlabel(r'$l$')
plt.ylabel(r'$l(l+1)C_l/2\pi$')
plt.title('TT Power Spectrum')
plt.legend(loc='best')
plt.show()

print('EE Chi square results:', chiEE)

plt.plot(ls,EE,label='Planck EE')
plt.plot(newls,newEE, label='PolSpice EE')
plt.xlim(0,lmax)
plt.xlabel(r'$l$')
plt.ylabel(r'$l(l+1)C_l/2\pi$')
plt.title('EE Power Spectrum')
plt.legend(loc='best')
plt.show()

print('BB Chi square results:', chiBB)

plt.plot(ls,BB,label='Planck BB')
plt.plot(newls,newBB, label='PolSpice BB')
plt.xlim(0,lmax)
plt.xlabel(r'$l$')
plt.ylabel(r'$l(l+1)C_l/2\pi$')
plt.title('BB Power Spectrum')
plt.legend(loc='best')
plt.show()

print('TE Chi square results:', chiTE)

plt.plot(ls,TE,label='Planck TE')
plt.plot(newls,newTE, label='PolSpice TE')
plt.xlim(0,lmax)
plt.xlabel(r'$l$')
plt.ylabel(r'$l(l+1)C_l/2\pi$')
plt.title('TE Power Spectrum')
plt.legend(loc='best')
plt.show()