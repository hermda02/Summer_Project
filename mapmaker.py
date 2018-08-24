import matplotlib.pyplot as plt
import numpy as np
import healpy as hp
import random 

simulated = np.loadtxt('./data/r_lensedtotCls_sim.dat')
planck = np.loadtxt("./data/COM_PowerSpect_CMB-base-plikHM-TTTEEE-lowl-lowE-lensing-minimum-theory_R3.01.txt")

filename = './fits/'

ell       = planck[:,0]
planckTT  = planck[:,1]/(ell*(ell+1)/(2*np.pi))
planckTE  = planck[:,2]/(ell*(ell+1)/(2*np.pi))
planckEE  = planck[:,3]/(ell*(ell+1)/(2*np.pi))
planckBB  = planck[:,4]/(ell*(ell+1)/(2*np.pi))

simell    = simulated[:,0]
simTT     = simulated[:,1]/(simell*(simell+1)/(2*np.pi))
simEE     = simulated[:,2]/(simell*(simell+1)/(2*np.pi))
simBB     = simulated[:,3]/(simell*(simell+1)/(2*np.pi))
simTE     = simulated[:,4]/(simell*(simell+1)/(2*np.pi))

planckpol = [planckTT,planckEE,planckBB,planckTE]
simupol   = [simTT,simEE,simBB,simTE]

maptype   = input('Map for Planck (p) or Simulated (s)?: ')
maptype   = str(maptype)

nside     = input('N_side for map: ')
nside     = int(nside)
pixel     = nside*nside*12

arcmin    = input('Beam width (in arcmin): ')
beam      = np.radians(arcmin/60)

noise     = input('How much noise should be added? (microK^2): ')

if noise == '':
	noise = 0.0
else:
	noise = float(noise)
	TT_noise  = np.empty(pixel)
	EE_noise  = np.empty(pixel)
	BB_noise  = np.empty(pixel)


if maptype == 'p':
	filename = filename + 'planckmap_' + str(nside)
	syn      = hp.synfast(planckpol,nside,fwhm=beam,new=True)

if maptype == 's':
	filename = filename + 'simulatedmap_' + str(nside)
	syn   = hp.synfast(simupol,nside,fwhm=beam,new=True)

if beam != 0.0:
	filename = filename+'_beam_'+str(arcmin)

if noise != 0.0:
	filename = filename+'_w_'+str(noise)+'_noise.fits'
	synTT = syn[0]
	synEE = syn[1]
	synBB = syn[2]
	for i in range(pixel):
		TT_noise[i] = synTT[i] + random.gauss(0,noise)
		EE_noise[i] = synEE[i] + random.gauss(0,noise)
		BB_noise[i] = synBB[i] + random.gauss(0,noise)
	newmap = [TT_noise,EE_noise,BB_noise]
else:
	newmap = [syn[0],syn[1],syn[2]]
	filename = filename+'.fits'

hp.write_map(filename,newmap)
print 'Map written to ', filename