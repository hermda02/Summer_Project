import matplotlib.pyplot as plt
import numpy as np
import healpy as hp
from scipy.optimize import leastsq
from scipy.optimize import minimize
from scipy.stats import chisquare

plt.rc('text', usetex=True)

# Import text files
lens_temp  = np.loadtxt('./data/lens_template.txt')
r_temp     = np.loadtxt('./data/r0.1_template.txt')
simulated  = np.loadtxt('./data/r_lensedtotCls_sim.dat')
polspice   = np.loadtxt('./data/sim_map_512.dat')
spiders    = np.loadtxt('./data/bb_data_spec.txt')
foreground = np.loadtxt('./data/dust_template.dat')

polls = polspice[2:,0]
polTT = polspice[2:,1]
polEE = polspice[2:,2]
polBB = polspice[2:,3]
polTE = polspice[2:,4]

ells  = spiders[:19,0]
dls   = spiders[:19,1]
err   = spiders[:19,2]

lmax  = int(max(polls))

lens  = lens_temp[:lmax-1,1]
r     = r_temp[:lmax-1,1]

simTT = simulated[:lmax-1,1]
simEE = simulated[:lmax-1,2]
simBB = simulated[:lmax-1,3]
simTE = simulated[:lmax-1,4]

forBB = foreground[:lmax-1,3]

# Weight Power spectra for plotting
for i in range(int(lmax-1)):
	polTT[i] = polTT[i]*(polls[i]*(polls[i]+1))/(2*np.pi)
	polEE[i] = polEE[i]*(polls[i]*(polls[i]+1))/(2*np.pi)
	polBB[i] = polBB[i]*(polls[i]*(polls[i]+1))/(2*np.pi)
	polTE[i] = polTE[i]*(polls[i]*(polls[i]+1))/(2*np.pi)
	forBB[i] = forBB[i]*(polls[i]*(polls[i]+1))/(2*np.pi)

inp = input('Create fit for simulated (sim) or SPIDER (spi) data?: ')
inp = str(inp)

if inp == 'sim':
	# Binning simulated spectra
	binl   = np.empty(int(lmax/25))
	binsBB = np.empty(int(lmax/25))
	binBB  = np.empty(int(lmax/25))
	declBB = np.empty(int(lmax/25))
	binlen = np.empty(int(lmax/25))
	binr   = np.empty(int(lmax/25))

	for i in range(int(lmax/25)):
		binl[i]   = 12 + i*25
		binBB[i]  = np.sum(simBB[i*25:24+i*25])/24
		binsBB[i] = np.sum(polBB[i*25:24+i*25])/24
		declBB[i] = np.sqrt(2/((2*binl[i]+1)))*binsBB[i]
		binlen[i] = np.sum(lens[i*25:24+i*25])/24
		binr[i]   = np.sum(r[i*25:24+i*25])/24

	x0 = np.empty(2,dtype=float)
	x = input("Lensing weight estimate: ")
	x0[0] = float(x)
	y = input("r component estimate: ")
	x0[1] = float(y)

	# Minimize the simulated fit!
	def sim_fit(x):
		model = x[0]*binlen + x[1]*binr
		chisq = np.sum((binsBB - model)**2/declBB)
		return chisq

	bnds = ((0,None),(0,None))

	result = minimize(sim_fit,x0,bounds=bnds)

	z         = np.empty(2,dtype=float)
	z[0],z[1] = result.x
	fit       = z[0]*binlen + z[1]*binr 

	# chi = chisquare(fit,binsBB)

	print "Optimized values: "
	print "Lensing weight = ", z[0]
	print "r weight = ", z[1]
	print "Fit chi-squared: ", sim_fit(z)

	plt.errorbar(binl,binsBB,yerr=declBB,color='grey',label='Binned PolSpice w/ error')
	# plt.plot(polBB,label='PolSpice Data')
	plt.plot(binl,fit,label='Binned Best fit')
	plt.xlim(0,lmax)
	plt.xlabel(r'$l$')
	plt.ylabel(r'$l(l+1)C_l/2\pi$')
	plt.title('BB Power Spectrum')
	plt.legend(loc='best')
	plt.show()
#-------------------------------------------------------------------------
#-------------------------------------------------------------------------
#-------------------------------------------------------------------------
# Binning stuff for SPIDER

if inp == 'spi':

	binlen = np.empty(int(lmax/25)-1)
	binr   = np.empty(int(lmax/25)-1)
	binfor = np.empty(int(lmax/25)-1)
	declBB = np.empty(int(lmax/25)-1)
	signal = np.empty(int(lmax/25)-1)
	vari   = np.empty(int(lmax/25)-1)
	total  = np.empty(int(lmax/25)-1)

	for i in range(int(lmax/25-1)):
		binlen[i] = np.sum(lens[7+i*25:32+i*25])/24
		binr[i]   = np.sum(r[7+i*25:32+i*25])/24
		vari[i]   = np.sqrt(2/((2*ells[i]+1)))*dls[i]
		declBB[i] = np.sqrt(vari[i]**2+err[i]**2)
		signal[i] = np.sum(forBB[7+i*25:32+i*25])/24

	x0 = np.empty(2,dtype=float)
	x = input("Lensing weight estimate: ")
	x0[0] = float(x)
	y = input("r component estimate: ")
	x0[1] = float(y)

	# Minimize the SPIDER data fit!
	def spider_fit(x):
		model = x[0]*binlen + x[1]*binr + signal
		chisq = np.sum(((dls - model)/declBB)**2)
		return abs(chisq)

	bnds = ((0,None),(0,None))

	result = minimize(spider_fit,x0,bounds=bnds)

	z         = np.empty(2,dtype=float)
	z[0],z[1] = result.x
	fit       = z[0]*binlen + z[1]*binr + signal

	chi = str(spider_fit(z))

	print "Optimized values: "
	print "Lensing weight = ", z[0]
	print "r weight = ", z[1]
	print "Fit chi-squared: ", chi

	plt.errorbar(ells,dls,yerr=declBB,color='grey',marker='o',label='SPIDER data w/ error')
	# plt.plot(ells,z[0]*binlen,label='lensing')
	# plt.plot(ells,z[1]*binr,label='r')
	plt.plot(ells,fit,label='Binned Best fit')
	plt.text(250,2.7,r'$\chi^2 =$ %s' % (chi))
	plt.text(250,1.7,'lensing weight = %s \n r = %s' % (z[0],0.1*z[1]))
	plt.xlim(0,lmax)
	plt.xlabel(r'$l$')
	plt.ylabel(r'$l(l+1)C_l/2\pi$')
	plt.title('BB Power Spectrum')
	plt.legend(loc='best')
	plt.show()