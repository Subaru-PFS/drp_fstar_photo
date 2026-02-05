import numpy as np
import emcee
import matplotlib.pyplot as plt
import pandas as pd
from scipy import interpolate
from scipy.stats import sigmaclip
from scipy.integrate import quad

import sys, os, glob, time,datetime, resource
import corner
from astropy import units as u
from astropy.coordinates import SkyCoord

import h5py
	

from memory_profiler import profile

# See https://docs.astropy.org/en/stable/api/astropy.coordinates.galactocentric_frame_defaults.html
from astropy.coordinates import galactocentric_frame_defaults
_ = galactocentric_frame_defaults.set('v4.0')

state = galactocentric_frame_defaults.get_from_registry("v4.0")
r_sun = (state["parameters"]["galcen_distance"] * 1000.).value
z_sun = (state["parameters"]["z_sun"]).value

from multiprocessing import Pool, RawArray
import multiprocessing as mp



from brutus import filters
from brutus import seds
from brutus import utils as butils

import brutus




#state = galactocentric_frame_defaults.get_from_registry("v4.0")



def log_prob_test(x, mu, cov):

	diff = x - mu

	return -0.5 * np.dot(diff, np.linalg.solve(cov, diff))



def log_like(theta, magobs, cov, sedmaker):
	
	# The extinction vector R for PS1 grizy from Table 1 of Green et al. 2019, ApJ, 887, 27
	extvect = np.array([3.518, 2.617, 1.971, 1.549, 1.263])

	# 50% completeness limit from https://outerspace.stsci.edu/display/PANSTARRS/PS1+Photometric+Depth#PS1PhotometricDepth-DepthofGalaxiesinComparisontoPointSources
	#           => should be revised to the selection function 
	maglim = np.array([23.2, 23.2, 23.1, 22.3, 21.2]) 

	# Obtain model magnitudes:
	eep, feh, mini, ebv, dmod = theta
	rv = 3.1
	av = ebv * rv
	dist = 10**((dmod + 5.)/5.) 

	#modelabsmags, teff, logg = get_magmodel_from_isochrone(lage, feh, mass)	

	modelmags, params, _ = sedmaker.get_sed(mini=mini, feh=feh, eep=eep, av=av, rv=rv, dist=dist)


	if np.count_nonzero(np.isnan(modelmags))!=0:
		return -1.*float('inf')


	# Incorpolating survey selection function 
	for i, modelmag in enumerate(modelmags):
		if modelmag > maglim[i]:
			return -1.*float('inf') 


	modelplx = 1000./dist 
	model = np.hstack((modelmags, modelplx))

	diff = magobs - model 

	return -0.5 * np.dot(diff, np.linalg.solve(cov, diff))
	

def log_prior(params, galcoord, ebvmax):

	eep, feh, mass, ebv, dmod = params

	if ebv>ebvmax:
		return -1.*float('inf')

	# Calculate coordinates in the Galactocentric Cylinderical frame 
	l, b = galcoord
	d = 10**((dmod + 5.)/5.)

	c = SkyCoord(l = l*u.degree, b = b*u.degree, distance = d*u.pc, frame = 'galactic')


	x_star = c.galactocentric.x.value
	y_star = c.galactocentric.y.value
	z_star = c.galactocentric.z.value
	r_star = np.sqrt(x_star**2 + y_star**2)

	return lp_eep(eep, r_star, z_star) + lp_feh(feh, r_star, z_star) + lp_mass(mass) + lp_ebv(ebv, ebvmax) + lp_dmod(dmod, r_star, z_star) 


def lp_eep(eep, r_star, z_star):

	if eep < 0 or eep>500:
		return -1.*float('inf')
	else:
		return 0.

    

def lp_age(lage, r_star, z_star):


	age = 10.**lage

	mean_thin, std_thin = 5.0e+9, 4.0e+9
	mean_thick, std_thick = 10.5e+9, 2.0e+9
	mean_halo, std_halo = 12.5e+9, 1.0e+9



	if lage > 10.13 or lage < 6.6:
		return -1.*float('inf')


	p_thin = nthin(r_star, z_star) /(nthin(r_star, z_star) + nthick(r_star, z_star) + nhalo(r_star, z_star))
	p_thick = nthick(r_star, z_star) /(nthin(r_star, z_star) + nthick(r_star, z_star) + nhalo(r_star, z_star))
	p_halo = 1.-p_thin - p_thick

	
	p_age_thin = p_thin * (1./(np.sqrt(2*np.pi)*std_thin)) * np.exp(-(age-mean_thin)**2/(2*std_thin**2))
	p_age_thick = p_thick * (1./(np.sqrt(2*np.pi)*std_thick)) * np.exp(-(age-mean_thick)**2/(2*std_thick**2))
	p_age_halo = p_halo  * (1./(np.sqrt(2*np.pi)*std_halo)) * np.exp(-(age-mean_halo)**2/(2*std_halo**2))
	
	
	return np.log(p_age_thin + p_age_thick + p_age_halo)

	

def lp_mabs(mabs_r):

	if mabs_r < -1. or mabs_r>16:
		return -1.*float('inf')

	else:
		return 0.

	return

def lp_mass(mass):

	# Kroupa 2001:
	#  https://ui.adsabs.harvard.edu/abs/2001MNRAS.322..231K/abstract

	if mass <= 0.01 or mass > 50.:
		return -1.*float('inf')
	elif mass <= 0.08:
		return np.log(mass**(-0.3))
	elif mass <= 0.50:
		return np.log(mass**(-1.3))
	elif mass <= 1.00:
		return np.log(mass**(-2.3))
	else:
		return np.log(mass**(-2.3))


def lp_ebv(ebv, ebvmax):

	if ebv < 0. or ebv > ebvmax:
		return -1.*float('inf')
	else:
		return 0.


	return




def lp_feh(feh, r_star, z_star):


	# [Fe/H] parameters from Table 2 of Green+14
	a_D = -0.89
	sig_D = 0.20
	c = 0.63
	Delta_a = 0.14
	Delta_mu = 0.55
	H_mu = 0.5 * 1000.
	a_H = -1.46
	sig_H = 0.30




	p_disk = (nthin(r_star, z_star) + nthick(r_star, z_star))/(nthin(r_star, z_star) + nthick(r_star, z_star) + nhalo(r_star, z_star))
	p_halo = 1.-p_disk


	a_Z = a_D + Delta_mu * np.exp(-1.*np.abs(z_star)/H_mu)

	a = 1./(np.sqrt(2.*np.pi)*sig_D)
	p_feh_disk = c *  a * np.exp(-1.*(feh - a_Z)**2/(2.*sig_D**2)) + \
	(1.-c) * a * np.exp(-1.*(feh - (a_Z+Delta_a))**2/(2.*sig_D**2))

	a = 1./(np.sqrt(2.*np.pi)*sig_H)
	p_feh_halo = a * np.exp(-1.*(feh - a_H)**2/(2.*sig_H**2))
	
	p = p_feh_disk * p_disk + p_feh_halo * p_halo

	lp = np.log(p)
	return lp 



def lp_dmod(dmod, r_star, z_star):

	p = (nthin(r_star, z_star) + nthick(r_star, z_star) + nhalo(r_star, z_star)) * 10**(3.*dmod/5.)

	return np.log(p)


def nthin(r_star, z_star):
	## thin disk
	rthin = 2150.
	zthin = 245.
	return np.exp(-1.*((r_star-r_sun)/rthin + (np.abs(z_star)-np.abs(z_sun))/zthin))


def nthick(r_star, z_star):
	## thick disk
	rthick = 3261.
	zthick = 743.
	fthick = 0.13
	return fthick * np.exp(-1.*((r_star-r_sun)/rthick + (np.abs(z_star)-np.abs(z_sun))/zthick))


def nhalo(r_star, z_star):

	## halo
	qh = 0.70
	fhalo = 0.003
	eota=2.62
	re = 500.
	reff = np.sqrt(r_star**2 + (z_star/qh)**2 + re**2)
	return fhalo * (reff/r_sun)**(-1*eota)
	

def log_prob(theta, galcoord, ebvmax, magobs, cov, sedmaker):
	
	lp = log_prior(theta, galcoord, ebvmax)

	return log_like(theta, magobs, cov, sedmaker) + lp





def get_magmodel(lage, feh, mass):


	modelmags = np.full(5, np.nan)


	if feh<=-2.1:
		isofile = "PS1_MH-2.5_-2.0_logAge6.6_10.20.dat"
	elif feh>-2.1 and feh<=-1.1:
		isofile = "PS1_MH-2.0_-1.0_logAge6.6_10.20.dat"
	elif feh>-1.1 and feh<=-0.1:
		isofile = "PS1_MH-1.0_-0.0_logAge6.6_10.20.dat"
	elif feh>-0.1:
		isofile = "PS1_MH0.0_1.0_logAge6.6_10.20.dat"


	df0 = pd.read_table(isopath + isofile, comment="#", delimiter="\s+", usecols = ["MH", "logAge", "Mini", "logTe", "logg", "label", "gP1mag", "rP1mag", "iP1mag", "zP1mag", "yP1mag"])


	filt = ((df0["MH"]).round(1) == np.round(feh, decimals=1)) & ((df0["logAge"]).round(1) == np.round(lage, decimals=1)) & (df0["label"] < 4)
	
	df = df0[filt]


	if len(df)==0:
		return(modelmags, np.nan)

	elif mass < np.min(df["Mini"]) or mass > np.max(df["Mini"]):
		return(modelmags, np.nan)
	
	



	#plt.plot(10**df["logTe"], df["logg"])
	#plt.savefig("../figs/isotest.png")
	

	for i, band in enumerate(["gP1mag", "rP1mag", "iP1mag", "zP1mag", "yP1mag", "logTe"]):

		intpfunc = interpolate.interp1d(df["Mini"], df[band], kind = "linear")
		if band == "logTe":	
			modelteff = 10**intpfunc(mass)
		else:
			modelmags[i] = intpfunc(mass)

	return(modelmags, modelteff)



def get_magmodel_cached(lage, feh, mass, df0):


    modelmags = np.full(5, np.nan)




    filt = ((df0["MH"]).round(1) == np.round(feh, decimals=1)) & ((df0["logAge"]).round(1) == np.round(lage, decimals=1)) & (df0["label"] < 4)


    df = df0[filt]


    if len(df)==0:
        return(modelmags, np.nan)

    elif mass < np.min(df["Mini"]) or mass > np.max(df["Mini"]):
        return(modelmags, np.nan)

    #plt.plot(10**df["logTe"], df["logg"])
    #plt.savefig("../figs/isotest.png")


    for i, band in enumerate(["gP1mag", "rP1mag", "iP1mag", "zP1mag", "yP1mag", "logTe"]):

        intpfunc = interpolate.interp1d(df["Mini"], df[band], kind = "linear")
        if band == "logTe":
            modelteff = 10**intpfunc(mass)
        else:
            modelmags[i] = intpfunc(mass)

    return(modelmags, modelteff)



def get_magmodel_arr(lage, feh, mass, isoarr0):


	modelmags = np.full(5, np.nan)

	filt = (np.round(isoarr0[:, 0], 1) == np.round(feh, 1)) & (np.round(isoarr0[:, 1], 1) == np.round(lage, 1)) & (isoarr0[:, 5]<4)
	isoarr = isoarr0[filt]

	if len(isoarr[:,0])==0:
		return(modelmags, np.nan)
	elif mass < np.min(isoarr[:, 2]) or mass > np.max(isoarr[:, 2]):
		return(modelmags, np.nan)

	#"MH", "logAge", "Mini", "logTe", "logg", "label", "gP1mag", "rP1mag", "iP1mag", "zP1mag", "yP1mag"

	for i in [6, 7, 8, 9, 10, 11, 3]:

		intpfunc = interpolate.interp1d(isoarr[:, 2], isoarr[:, i], kind = "linear")
		if i == 3:
			modelteff = 10**intpfunc(mass)
		else:
			modelmags[i] = intpfunc(mass)

	return(modelmags, modelteff)



def get_magmodel_from_isochrone(lage, feh, mass):

	# Output absolute magnidudes, Teff, logg given age and mass

	modelmags = np.full(5, np.nan)


	filt = (np.round(isoarr0[:, 0], 1) == np.round(feh, 1)) & (np.round(isoarr0[:, 1], 1) == np.round(lage, 1)) & (isoarr0[:, 5] < 4 )
	isoarr = isoarr0[filt]

	if len(isoarr[:,0])==0:
		return(modelmags, np.nan, np.nan)
	elif mass < np.min(isoarr[:, 2]) or mass > np.max(isoarr[:, 2]):
		return(modelmags, np.nan, np.nan)
	
	#"MH", "logAge", "Mini", "logTe", "logg", "label", "gP1mag", "rP1mag", "iP1mag", "zP1mag", "yP1mag"

	for i in [ 6, 7, 8, 9, 10, 3, 4 ]:

		intpfunc = interpolate.interp1d(isoarr[:, 2], isoarr[:, i], kind = "linear")
		if i == 3:
			modelteff = 10**intpfunc(mass)
		elif i == 4:
			modellogg = intpfunc(mass)
		else:
			j = i-6
			modelmags[j] = intpfunc(mass)

	return(modelmags, modelteff, modellogg)




def get_ini(param0, nwalkers, stepfrac):


	p0 = param0 + np.abs(param0)*stepfrac * np.random.rand(nwalkers, len(param0))

	return(p0)


def init_worker(X, X_shape):

	var_dict['X'] = X
	var_dict['X_shape'] = X_shape



def sample_posterior_segue():
	# Read stellar data 
	path = "../catalogs/MAST_PS1_crossmatch/"
	catalog = "sppParams_PS1DR1_crossmatch_GaiaEDR3_magerr.csv"
	df_segue = pd.read_csv(path + catalog)
	


	bands = ["gmag","rmag","imag","zmag","ymag"]
	banderrs = ["e_gmag","e_rmag","e_imag","e_zmag","e_ymag"] 

	ct = 0	
	for index, row in df_segue.iterrows():
		starname = row['specobjid']
		obsmags = [ row[band] for band in bands]
		obsmagerrs = [ row[banderr] for banderr in banderrs]  
		plx = row["parallax"]
		plxerr = row["parallax_error"]
		if index==0:
			continue
		if plx<0 or plxerr/plx > 0.2 or row["teffadopunc"] > 100.:
			continue 

		l = row["l"]
		b = row["b"]	

		sample_posterior(starname, l, b, obsmags, obsmagerrs, plx, plxerr)

		break
		ct = ct + 1

	return


@profile()
def sample_posterior(starname, l, b, obsmags, obsmagerrs, plx, plxerr, sedmaker):
	

	# Turn off NumPy's parallelizaion operation 
	os.environ["OMP_NUM_THREADS"] = "1"


	# Read stellar data 
	#path = "../catalogs/MAST_PS1_crossmatch/"
	#catalog = "sspParams_PS1DR1_crossmatch_GaiaEDR3_ExtinctionCorrected_selected.csv"
	#df0 = pd.read_csv(path + catalog)
	#df = df0.query('specobjid == 3228029599975761920')

	#print(df["raMean"], df["decMean"])
	#bands = ["gMeanPSFMag","rMeanPSFMag","iMeanPSFMag","zMeanPSFMag","yMeanPSFMag"]

	
	#obsmags = [ (df[band].values)[0] for band in bands]
	#obsmagerrs = [0.00168, 0.00348, 0.004935,0.005148, 0.000842]

	#plx = (df["parallax"].values)[0]
	#plx_err = (df["parallax_error"].values)[0]

	cov = np.diag(np.hstack((obsmagerrs, plxerr)))

	#l = (df["l"].values)[0]
	#b = (df["b"].values)[0]
	galcoord = [l, b]

	dist = 1000./plx
	coords = SkyCoord(l*u.deg, b*u.deg, distance = dist*u.pc, \
                          frame = 'galactic')

	sfd = SFDQuery()
	ebvmax = sfd(coords)


	# A grid of SED model
  
	#sedmaker = create_sedmodel()

	# Initial guess of the free parameters
	nwalkers = 32

	eep0 = 400
	feh0 = -1.0
	mass0 = 1.0
	ebv0 = 0.01
	dmod0 = 11.



	param0 = [eep0, feh0, mass0, ebv0, dmod0]
	ndim = len(param0)


	# fraction of each parameter for one step of the walkers
	stepfrac = [1, 0.1, 0.1, 0.01, 0.1]


	# Setup the backend 
	filename = "../outputs/MCMC/mcmc_brutus_" + str(starname) + ".h5"
	backend = emcee.backends.HDFBackend(filename)	
	backend.reset(nwalkers, ndim)

	theta0 = get_ini(param0, nwalkers, stepfrac)


	#max_n = 10000
	#index = 0
	#autocorr = np.empty(max_n)

	#old_tau = np.inf

	# Get isochrone data
	#df = pd.read_pickle(isopath + "PS1.pkl")
	#isoarr = df.to_numpy()
	#X_shape = isoarr.shape
	#X = RawArray('d', X_shape[0] * X_shape[1])
	#X_np = np.frombuffer(X).reshape(X_shape)
	#np.copyto(X_np, isoarr)


	with Pool(12) as pool:
	#with Pool(processes = 50, initializer = init_worker, initargs = (X, X_shape)) as pool:

		sampler = emcee.EnsembleSampler(nwalkers, ndim, log_prob, args = [galcoord, ebvmax, np.hstack((obsmags, plx)), cov, sedmaker], pool = pool, backend = backend)

		#sampler = emcee.EnsembleSampler(nwalkers, ndim, log_prob, args = [galcoord, np.hstack((obsmags, plx)), cov], backend = backend)

		start = time.time()


		# Start sampling

		print("Start sampling!")

		#for sample in sampler.sample(p0, iterations = max_n, progress=True):
		#if sampler.iteration % 100:
			#	continue
			#tau = sampler.get_autocorr_time(tol=0)
			#autocorr[index] = np.mean(tau)
			#index += 1
			
		# Check convergence
		#converged = np.all(tau * 100 < sampler.iteration)
		#converged &= np.all(np.abs(old_tau-tau)/tau < 0.01)

		#if converged:
		#	break
		#old_tau = tau
	
		sampler.run_mcmc(theta0, 10)

		end = time.time()
		multi_time = end - start
		print("Processing took {0:.1f} seconds".format(multi_time))

	print(
    	"Mean acceptance fraction: {0:.3f}".format(
        	np.mean(sampler.acceptance_fraction)
    	)
	)


	#print(
    # 	"Mean autocorrelation time: {0:.3f} steps".format(
    #     	np.mean(sampler.get_autocorr_time())
    # 	)
	#)



	# Print status of the chane

	#tau = sampler.get_autocorr_time()
	#burnin = int(2 * np.max(tau))
	#thin = int(0.5 * np.min(tau))
	#samples = sampler.get_chain(discard=burnin, flat=True, thin=thin)
	#log_prob_samples = sampler.get_log_prob(discard=burnin, flat=True, thin=thin)
	#log_prior_samples = sampler.get_blobs(discard=burnin, flat=True, thin=thin)

	#print("burn-in: {0}".format(burnin))
	#print("thin: {0}".format(thin))
	#print("flat chain shape: {0}".format(samples.shape))
	#print("flat log prob shape: {0}".format(log_prob_samples.shape))
	#print("flat log prior shape: {0}".format(log_prior_samples.shape))



	#all_samples = np.concatenate((samples, log_prob_samples[:, None], log_prior_samples[:, None]), axis = 1)

	#labels = list(map(r"$\theta_{{{0}}}".format, range(1, ndim + 1)))
	#labels += ["log prob", "log prior"]

	#labels = ["log(Age)", "[Fe/H]", "mass", "E(B-V)", "dmod" ]
	#corner.corner(samples, labels=labels)
	#plt.savefig("../figs/emcee_" + str(starname) + "_corner.png")
	
	#samples = sampler.get_chain(flat=True)
	#plt.hist(samples[:, 4], 100, color="k", histtype="step")
	#plt.savefig("../figs/emcee_startype.png")

	return



def calc_posterior_segue():

	path = "../catalogs/MAST_PS1_crossmatch/"
	catalog = "sppParams_PS1DR1_crossmatch_GaiaEDR3_magerr.csv"
	df_segue = pd.read_csv(path + catalog)


	bands = ["gmag","rmag","imag","zmag","ymag"]
	banderrs = ["e_gmag","e_rmag","e_imag","e_zmag","e_ymag"]



	ct = 0
	for index, row in df_segue.iterrows():
		starname = row['specobjid']
		obsmags = [ row[band] for band in bands]
		obsmagerrs = [ row[banderr] for banderr in banderrs]
		plx = row["parallax"]
		plxerr = row["parallax_error"]



		if index==0:
			continue
		if plx<0 or plxerr/plx > 0.2 or row["teffadopunc"] > 100.:
			continue

		l = row["l"]
		b = row["b"]


		start = time.time()

		print("Start calculating an unnormalized posterior for %s: Teff = %f, logg = %f, [Fe/H] = %f"%(starname, row["teffadop"], row["loggadop"], row["fehadop"]))

		calc_posterior(starname, l, b, obsmags, obsmagerrs, plx, plxerr)

		end = time.time()
		multi_time = end - start
		print("Processing took {0:.1f} seconds".format(multi_time))

		
		break
		ct = ct + 1

	return

def calc_posterior(starname, l, b, obsmags, obsmagerrs, plx, plxerr):
	
	# Calculate a posterior probability distribution for a fixed grid of parameters

	## Stellar parameter range and step size
	logages = np.arange(8.0, 10.2, 0.1)
	fehs = np.arange(-2.5, 1.0, 0.2)
	logmasses = np.arange(-1.1, 1.7, 0.05)
	masses = 10**logmasses

	#logages = np.arange(10.0, 10.2, 0.1)
	#fehs = np.arange(-1.0, -0.8, 0.2)
	#logmasses = np.arange(-0.1, 0.1, 0.1)




	## Define a covariance matrix
	cov = np.diag(np.hstack((obsmagerrs, plxerr)))


	## Galactic corrdinates to be used for the 3D dust map 
	galcoord = [l, b]


	## Consider +-2 sigma for parallax
	if plx - 2.*plxerr > 0:
		dists = [1000. / (plx + 2.*plxerr), 1000. / (plx - 2.*plxerr)]
	else:
		dists = [1000. / (plx + 2.*plxerr), 100000.]

	coords = SkyCoord([l, l]*u.deg, [b, b]*u.deg, distance = dists*u.pc, \
                          frame = 'galactic')

	# 3D map
	bayestar = BayestarQuery(max_samples=3)
	ebv_baye, flags = bayestar(coords, mode = 'percentile', pct = [16., 50., 84.], return_flags = True)

	
	if flags['reliable_dist'][0] == True:
		ebv_min = (ebv_baye[0])[0]
	else:
		ebv_min = 0.0

	if flags['reliable_dist'][1] == True:
		ebv_max = (ebv_baye[1])[2]
	else:
		sfd = SFDQuery()
		ebv_max = sfd(SkyCoord(l*u.deg, b*u.deg, distance = dists[1]*u.pc, frame='galactic')) 


	ebvstep_min = 0.01
	ebv_nsteps = 10.
	ebvstep = (ebv_max - ebv_min) / ebv_nsteps
	
	if ebvstep < ebvstep_min:
		ebvstep = ebvstep_min

	ebvs = np.arange(ebv_min, ebv_max, ebvstep)

	if len(ebvs)==0:
		ebvs = list[ebv_min]

	dmodstep = 0.05
	dmodminmax = 5.*np.log10(dists) - 5.
	dmods = np.arange(dmodminmax[0], dmodminmax[1], dmodstep)


	paramall = []

	for logage in logages:
		for feh in fehs:
			for mass in masses:
				for ebv in ebvs:
					for dmod in dmods:
						paramall.append((logage, feh, mass, ebv, dmod))

	
	ncpus = 10
	print("Number of CPUs used: ", ncpus)
	pool = mp.Pool(ncpus)

	args = [[params, galcoord, ebv_max, np.hstack((obsmags, plx)), cov] for params in paramall ]


	log_posterior = pool.map(wrap_log_prob, args)

	pool.close()


	# Get Teff and logg 
	pool = mp.Pool(ncpus)

	args = [ params[:3] for params in paramall ]

	results = pool.map(wrap_get_magmodel, args)

	pool.close()

	teff = [ line[1] for line in results ]
	logg = [ line[2] for line in results ]


	data = {'params': paramall, 'log_posterior': log_posterior, 'Teff': teff, 'logg': logg }
	dfout = pd.DataFrame(data)
	dfout.to_pickle('test3.pkl')

	return


def plot_hist_teff_logg():

	df = pd.read_pickle('test3.pkl')

	tstep = 50.
	t_bin_mins = np.arange(3000., 10000., tstep)
	p_m = np.zeros_like(t_bin_mins)

	for log_p, teff in zip(df['log_posterior'], df['Teff']):

		p = 10**log_p
		
		indx = np.where((t_bin_mins <= teff) & (t_bin_mins + tstep > teff))
		p_m[indx] = p_m[indx] + p

	print(t_bin_mins)
	print(p_m)

	fig, ax = plt.subplots(1, 1)
	ax.plot(t_bin_mins + 0.5*tstep, np.log(p_m))
	ax.set_xlim(3500., 6000.)

	plt.savefig("../figs/test_posterior.png")



def wrap_log_prob(args):
	return log_prob(*args)




def wrap_get_magmodel(args):
	return get_magmodel_from_isochrone(*args)


def wrap_get_magmodel_brutus(args):

	eep, feh, mini = args
	av = 0.0
	rv = 3.1
	dist = 1000.
	_, params, _ = sedmaker.get_sed(mini=mini, feh=feh, eep=eep, av=av, rv=rv, dist=dist)
	return(params)

def read_hdf(hdffile):

	starname = ((hdffile.split("_"))[2]).strip(".h5")

	reader = emcee.backends.HDFBackend(hdffile)
	samples = reader.get_chain(flat = True)
	
	labels = ["EEP", "[Fe/H]", "mass", "E(B-V)", "dmod" ]
	corner.corner(samples, labels = labels)
	plt.savefig("../figs/emcee_" + starname + "_corner.png")

	

	ncpus = 30
	print("Number of CPUs used: ", ncpus)
	pool = mp.Pool(ncpus)

	args = [sample[:3] for sample in samples[5000:, :] ]

	results = pool.map(wrap_get_magmodel_brutus, args)

	pool.close()

	teffs = [ 10**line['logt'] for line in results ]

	print(teffs)

	fig, ax = plt.subplots(1, 1)	
	ax.hist(teffs, 100, histtype="step")
	plt.savefig("../figs/emcee_" + starname + "_teff_hist.png")



	return




def test_emcee():
	ndim = 5

	np.random.seed(42)
	means = np.random.rand(ndim)
	cov = 0.5 - np.random.rand(ndim ** 2).reshape((ndim, ndim))

	cov = np.triu(cov)


	cov += cov.T - np.diag(cov.diagonal())

	cov = np.dot(cov, cov)



	nwalkers = 32
	p0 = np.random.rand(nwalkers, ndim)

	print(p0[0])
	sampler = emcee.EnsembleSampler(nwalkers, ndim, log_prob, args = [means, cov])


	print(log_prob(p0[0], means, cov))

	state = sampler.run_mcmc(p0, 100)
	sampler.reset()


	sampler.run_mcmc(state, 10000)


	samples = sampler.get_chain(flat=True)
	plt.hist(samples[:, 0], 100, color="k", histtype="step")
	plt.savefig("../figs/emceetest.png")


	print("Mean acceptance fraction: {0:.3f}".format(
	        np.mean(sampler.acceptance_fraction)
	    ))

	print(
    	"Mean autocorrelation time: {0:.3f} steps".format(
        	np.mean(sampler.get_autocorr_time())
    	)
	)


def test_get_magmodel():


	time_start = time.perf_counter()
	lage=10.0
	feh = -1.0
	mass = 1.0
	results = get_magmodel(lage, feh, mass)

	time_elapsed = (time.perf_counter() - time_start)
	memMb = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss/1024.0/1024.0
	print("%5.5f secs %5.1f MByte" % (time_elapsed,memMb))


	return

def test_get_magmodel2():


	df = pd.read_pickle(isopath + "PS1.pkl")
	lage=10.0
	feh = -1.0
	mass = 1.0


	time_start = time.perf_counter()
	results = get_magmodel_cached(lage, feh, mass, df)

	time_elapsed = (time.perf_counter() - time_start)
	memMb = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss/1024.0/1024.0
	print("%5.5f secs %5.1f MByte" % (time_elapsed,memMb))
	return


def test_get_magmodel3():

	df = pd.read_pickle(isopath + "PS1.pkl")
	lage=10.0
	feh = -1.0
	mass = 1.0

	isoarr0 = df.to_numpy()


	time_start = time.perf_counter()
	results = get_magmodel_arr(lage, feh, mass, isoarr0)



	time_elapsed = (time.perf_counter() - time_start)
	memMb = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss/1024.0/1024.0
	print("%5.5f secs %5.1f MByte" % (time_elapsed,memMb))
	return


def pickle_isofiles():


	isofiles = ["PS1_MH-2.5_-2.0_logAge6.6_10.20.dat", "PS1_MH-2.0_-1.0_logAge6.6_10.20.dat", "PS1_MH-1.0_-0.0_logAge6.6_10.20.dat", "PS1_MH0.0_1.0_logAge6.6_10.20.dat"]

	for i in range(0, 4):
		df0 = pd.read_table(isopath + isofiles[i], comment="#", delimiter="\s+", usecols = ["MH", "logAge", "Mini", "logTe", "logg", "label", "gP1mag", "rP1mag", "iP1mag", "zP1mag", "yP1mag"])
		print("Data size = %i is read from "%(len(df0)), isofiles[i])
		if i==0:
			df = df0
		else:
			df = pd.concat([df, df0], ignore_index = True)
	df.to_pickle(isopath + "PS1.pkl")

	return




def view_models():

	mistfile = '../../brutus/data/DATAFILES/MIST_1.2_EEPtrk.h5'
	
	f = h5py.File(mistfile, 'r')

	print(f.keys())
	grp = f['-1.00']

	dset = grp['0.00']
	#grp_name = grp.name
	#dset = f[grp_name]
	print(dset.keys())


	return




def create_sedmodel():

	from brutus import filters
	from brutus import seds
	from brutus import utils as butils

	import brutus


	# --- Importing filter names --------------------------------------------------------------------------------

	filt = filters.ps[:-2]

	# --- Generating a grid of models --------------------------------------------------------------------------------

	nnfile = '../../brutus/data/DATAFILES/nn_c3k.h5'
	mistfile = '../../brutus/data/DATAFILES/MIST_1.2_EEPtrk.h5'
	sedmaker = seds.SEDmaker(filters=filt, nnfile=nnfile, mistfile=mistfile)

	return(sedmaker)




def run_brutus_segue(outfile, include_plx=True):

	from brutus import filters
    #from brutus import seds
	from brutus import utils as butils
	from brutus.fitting import BruteForce
	import brutus

	filt = filters.ps[:-2]

	gridfile = "../../brutus/data/DATAFILES/grid_mist_v9.h5"
	(models_mist, labels_mist, lmask_mist) = butils.load_models(gridfile, filters=filt)

	# Initialize BruteForce class
	BF_mist = BruteForce(models_mist, labels_mist, lmask_mist)
	
	# Load data

	path = "../catalogs/MAST_PS1_crossmatch/"
	catalog = "sppParams_PS1DR1_crossmatch_GaiaEDR3_magerr.csv"

	cols = ['specobjid', 'gmag', 'e_gmag', 'rmag', 'e_rmag', 'imag', 'e_imag', 'zmag', 'e_zmag', 'ymag', 'e_ymag', 'parallax', 'parallax_error', 'l', 'b' ]
	df0 = pd.read_csv(path + catalog, usecols = cols, nrows=500)   # add nrows=X to limit the number of entries.
	
	# Remove entries with NaN
	#filt = (df0.isnull().any(axis=1)==False) & (df0['e_gmag']>0.)  & (df0['e_rmag']>0.)  & (df0['e_imag']>0.)  & (df0['e_zmag']>0.)  & (df0['e_ymag']>0.)
	df = df0[df0.isnull().any(axis=1)==False]
	#filt = (df1['e_gmag']>0.)  & (df1['e_rmag']>0.)  & (df1['e_imag']>0.)  & (df1['e_zmag']>0.)  & (df1['e_ymag']>0.)
	#df = df1[filt]
	#df = df1
	print("Analyzing %i/%i stars"%(len(df), len(df0)))


	objid = df['specobjid']
	mag = [ [line[0], line[1], line[2], line[3], line[4]] for line in zip(df["gmag"], df["rmag"], df["imag"], df["zmag"], df["ymag"])]
	magerr = [ [line[0], line[1], line[2], line[3], line[4]] for line in zip(df["e_gmag"], df["e_rmag"], df["e_imag"], df["e_zmag"], df["e_ymag"]) ]
	mask = np.isfinite(magerr)
	phot, err = butils.inv_magnitude(np.array(mag), np.array(magerr))
      

	print(len(phot[phot<0.]))
	print(len(err[err<=0.0]))
	
	plx = df["parallax"].values
	plxerr =df["parallax_error"].values


	l = df["l"].values
	b = df["b"].values
	coords = np.c_[l, b]


	# Photometric offset 
	off_file_mist = "../../brutus/data/DATAFILES/offsets_mist_v8.txt"
	off_mist = butils.load_offsets(off_file_mist, filters=filt)

	# Fitting the data
	np.random.seed(862)

	dustfile = "../../brutus/data/DATAFILES/bayestar2019_v1.h5"

	
	if os.path.isfile(outfile + ".h5"):
		os.remove(outfile + ".h5")


	if include_plx==True:
		
		parallax = plx
		parallax_err = plxerr

	else:
		parallax = np.zeros_like(plx)
		parallax[:]=np.nan
		parallax_err = np.zeros_like(plxerr)
		parallax_err[:] = np.nan
	
	BF_mist.fit(phot, err, mask, objid, outfile, data_coords=coords, \
	parallax=parallax, parallax_err=parallax_err, phot_offsets=off_mist, mag_max=50., merr_max=0.5, \
	dustfile=dustfile, Ndraws=2500, Nmc_prior=50, logl_dim_prior=True, save_dar_draws=True, running_io=True, verbose=True)

	return


def run_brutus(inputcatalog, include_plx=False):

	from brutus import filters
    #from brutus import seds
	from brutus import utils as butils
	from brutus.fitting import BruteForce
	import brutus
	

	outfile = "../outputs/PS1_Gaia_mist_woplx" + ((inputcatalog.split("GaiaEDR3"))[-1]).replace(".csv", "")

	filt = filters.ps[:-2]

	gridfile = "../../brutus/data/DATAFILES/grid_mist_v9.h5"
	(models_mist, labels_mist, lmask_mist) = butils.load_models(gridfile, filters=filt)

	# Initialize BruteForce class
	BF_mist = BruteForce(models_mist, labels_mist, lmask_mist)
	
	# Load data


	cols = ['objID', 'source_id','gmag', 'e_gmag', 'rmag', 'e_rmag', 'imag', 'e_imag', 'zmag', 'e_zmag', 'ymag', 'e_ymag', 'parallax', 'parallax_error', 'l', 'b' ]
	df0 = pd.read_csv(inputcatalog, usecols = cols)   # add nrows=X to limit the number of entries.
	
	# Remove entries with NaN
	#filt = (df0.isnull().any(axis=1)==False) & (df0['e_gmag']>0.)  & (df0['e_rmag']>0.)  & (df0['e_imag']>0.)  & (df0['e_zmag']>0.)  & (df0['e_ymag']>0.)
	df = df0[df0.isnull().any(axis=1)==False]
	#filt = (df1['e_gmag']>0.)  & (df1['e_rmag']>0.)  & (df1['e_imag']>0.)  & (df1['e_zmag']>0.)  & (df1['e_ymag']>0.)
	#df = df1[filt]
	#df = df1
	print("Analyzing %i/%i stars"%(len(df), len(df0)))


	objid = df['source_id']
	mag = [ [line[0], line[1], line[2], line[3], line[4]] for line in zip(df["gmag"], df["rmag"], df["imag"], df["zmag"], df["ymag"])]
	magerr = [ [line[0], line[1], line[2], line[3], line[4]] for line in zip(df["e_gmag"], df["e_rmag"], df["e_imag"], df["e_zmag"], df["e_ymag"]) ]
	mask = np.isfinite(magerr)
	phot, err = butils.inv_magnitude(np.array(mag), np.array(magerr))
      

	
	plx = df["parallax"].values
	plxerr =df["parallax_error"].values


	l = df["l"].values
	b = df["b"].values
	coords = np.c_[l, b]


	# Photometric offset 
	off_file_mist = "../../brutus/data/DATAFILES/offsets_mist_v8.txt"
	off_mist = butils.load_offsets(off_file_mist, filters=filt)

	# Fitting the data
	np.random.seed(862)

	dustfile = "../../brutus/data/DATAFILES/bayestar2019_v1.h5"

	
	if os.path.isfile(outfile + ".h5"):
		os.remove(outfile + ".h5")


	if include_plx==True:
		
		parallax = plx
		parallax_err = plxerr

	else:
		parallax = np.zeros_like(plx)
		parallax[:]=np.nan
		parallax_err = np.zeros_like(plxerr)
		parallax_err[:] = np.nan
	
	BF_mist.fit(phot, err, mask, objid, outfile, data_coords=coords, \
	parallax=parallax, parallax_err=parallax_err, phot_offsets=off_mist, mag_max=50., merr_max=0.5, \
	dustfile=dustfile, Ndraws=2000, Nmc_prior=50, logl_dim_prior=True, save_dar_draws=True, running_io=True, verbose=True)

	return



def run_brutus_multiprocess():

	#start_time = time.process_time()

	start_time = time.perf_counter()

	catalogs = glob.glob("../comissioning/Nov2022/PS1DR2_GaiaEDR3_ra?.*csv")

	#catalogs = ["../comissioning/Nov2022/PS1DR2_GaiaEDR3_ra0.0_1.0_dec-40.0_90.0.csv", "../comissioning/Nov2022/PS1DR2_GaiaEDR3_ra1.0_2.0_dec-40.0_90.0.csv"]


	

	ncpus = mp.cpu_count()
	print("Number of CPUs used: ", ncpus)
	pool = mp.Pool(ncpus)

	results = pool.map(run_brutus, [catalog for catalog in catalogs])

	pool.close()

	#end_time = time.process_time()
	end_time = time.perf_counter()

	elapse_time = end_time - start_time
	print("CPU time = " ,elapse_time)

	return()





def read_brutus_results_segue(filename, paramfile, plot=False):
	
	from brutus import utils as butils
	from brutus import filters

	from brutus import plotting as bplot
	
	filt = filters.ps[:-2]

	# import MIST model grid
	gridfile = '../../brutus/data/DATAFILES/grid_mist_v9.h5'
	(models_mist, labels_mist, lmask_mist) = butils.load_models(gridfile, filters=filt)


	# Photometric offset 
	off_file_mist = "../../brutus/data/DATAFILES/offsets_mist_v8.txt"
	off_mist = butils.load_offsets(off_file_mist, filters=filt)

	print(filename +  '.h5')
	f = h5py.File(filename + '.h5', 'r')



	
	idxs_mist = f['model_idx'][:]  # model indices
	chi2_mist = f['obj_chi2min'][:]  # best-fit chi2
	nbands_mist = f['obj_Nbands'][:]  # number of bands in fit
	dists_mist = f['samps_dist'][:]  # distance samples
	reds_mist = f['samps_red'][:]  # A(V) samples
	dreds_mist = f['samps_dred'][:]  # R(V) samples
	lnps_mist = f['samps_logp'][:]  # log-posterior of samples



	# Load SEGUE catalog data 
	path = "../catalogs/MAST_PS1_crossmatch/"
	catalog = "sppParams_PS1DR1_crossmatch_GaiaEDR3_magerr.csv"
	df0 = pd.read_csv(path + catalog)

	teff = () 
	teff_low = ()
	teff_high = ()

	logg = ()
	logg_low = ()
	logg_high = ()

	
	oteff = ()
	oteff_sig = ()	

	ologg = ()
	ologg_sig = ()

	ofeh = ()
	ofeh_sig = ()





	for i, lab in enumerate(f['labels'][:]):
		df = df0[df0['specobjid']==lab]	

		mag = [ [line[0], line[1], line[2], line[3], line[4]] for line in zip(df["gmag"], df["rmag"], df["imag"], df["zmag"], df["ymag"])]
		magerr = [ [line[0], line[1], line[2], line[3], line[4]] for line in zip(df["e_gmag"], df["e_rmag"], df["e_imag"], df["e_zmag"], df["e_ymag"]) ]
		mask = np.isfinite(magerr)
		phot, err = butils.inv_magnitude(np.array(mag), np.array(magerr))
		plx = df["parallax"].values
		plxerr =df["parallax_error"].values

		oteff = np.hstack((oteff, df["teffadop"]))
		oteff_sig = np.hstack((oteff_sig, df["teffadopunc"]))
		ologg = np.hstack((ologg, df["loggadop"]))
		ologg_sig = np.hstack((ologg_sig, df["loggadopunc"]))
		ofeh = np.hstack((ofeh, df["fehadop"]))
		ofeh_sig = np.hstack((ofeh_sig, df["fehadopunc"]))
			

		#l = df["l"].values
		#b = df["b"].values
		#coords = np.c_[l, b]
		if plot==True:
			fig, ax, parts = bplot.posterior_predictive(models_mist,  # stellar model grid
                                            idxs_mist[i],  # model indices
                                            reds_mist[i],  # A(V) draws
                                            dreds_mist[i],  # R(V) draws
                                            dists_mist[i],  # distance draws
                                            data=phot[0], data_err=err[0],  # data
                                            data_mask=mask[0],  # band mask
                                            offset=off_mist,  # photometric offsets
                                            psig=2.,  # plot 2-sigma errors
                                            labels=filt,  # filters 
                                            vcolor='blue',  # "violin plot" colors for the posteriors
                                            pcolor='black')  # photometry colors for the data
			figname = "../figs/sed_" + str(lab) + ".png"

			fig.savefig(figname)

		print('Best-fit chi2 (MIST):', chi2_mist[i])


		# Get Teff
		#print(labels_mist.dtype.names)
		params = labels_mist[idxs_mist[i]]
		lteff = [ line[5] for line in params ]
		q1, q2, q3 = butils.quantile(lteff, [0.159, 0.5, 0.841])  # 1sigma quantiles
		teff = np.hstack((teff,10**q2))
		teff_low = np.hstack((teff_low, 10**q1))
		teff_high = np.hstack((teff_high, 10**q3))
			
		loggs = [ line[6] for line in params ]
		q1, q2, q3 = butils.quantile(loggs, [0.159, 0.5, 0.841])
		logg = np.hstack((logg, q2))
		logg_low = np.hstack((logg_low, q1))
		logg_high = np.hstack((logg_high, q3))

		if plot==True:
			fig2, axes = bplot.cornerplot(idxs_mist[i], 
                             (dists_mist[i], reds_mist[i], dreds_mist[i]),
                             labels_mist,  # model labels
                             parallax=plx[0], parallax_err=plxerr[0],  # overplot if included
                             show_titles=True, color='blue', pcolor='orange',
                             fig=plt.subplots(11, 11, figsize=(55, 55)))  #  custom figure
			plt.xticks(fontsize=12)
			figname = "../figs/corner_" + str(lab) + ".png"
			fig2.savefig(figname)

		if i==300:
			break

	outdata = { "teff":teff, "teff_low":teff_low, "teff_high":teff_high, "logg":logg, "logg_low":logg_low, "logg_high":logg_high, \
	"oteff":oteff, "oteff_sig":oteff_sig, "ologg":ologg, "ologg_sig":ologg_sig, "ofeh":ofeh, "ofeh_sig":ofeh_sig }

	dfout = pd.DataFrame(outdata)
	dfout.to_csv(paramfile, index = False)

	return


def asymmetric_gaussian(x, mu, sd_low, sd_high):
	
	if x <= mu:
		sd = sd_low
	else:
		sd = sd_high
	return 1 / (sd*(2.*np.pi)**0.5) * np.exp((x-mu)**2/(-2.*sd**2))

def fstar_prob(teff, teff_low, teff_high):

	# Fstar range
	tmin = 6000.
	tmax = 7500.

	mu = teff
	sd_low = teff-teff_low
	sd_high = teff_high-teff

	fstar_prob = quad(asymmetric_gaussian, tmin, tmax, args=(mu, sd_low, sd_high))
	#fstar_prob = quad(lambda x: 1 / (sd*(2.*np.pi)**0.5) * np.exp((x-mu)**2/(-2.*sd**2)), tmin, tmax)
	#print(quad(lambda x: 1 / (sd*(2.*np.pi)**0.5) * np.exp((x-mu)**2/(-2.*sd**2)), -np.inf, np.inf), "  should be 1.")


	return fstar_prob

def read_brutus_results(filename, catalog, paramfile, plot=False):
	
	from brutus import utils as butils
	from brutus import filters

	from brutus import plotting as bplot
	
	filt = filters.ps[:-2]

	# import MIST model grid
	gridfile = '../../brutus/data/DATAFILES/grid_mist_v9.h5'
	(models_mist, labels_mist, lmask_mist) = butils.load_models(gridfile, filters=filt)


	# Photometric offset 
	off_file_mist = "../../brutus/data/DATAFILES/offsets_mist_v8.txt"
	off_mist = butils.load_offsets(off_file_mist, filters=filt)

	f = h5py.File(filename, 'r')



	
	idxs_mist = f['model_idx'][:]  # model indices
	chi2_mist = f['obj_chi2min'][:]  # best-fit chi2
	nbands_mist = f['obj_Nbands'][:]  # number of bands in fit
	dists_mist = f['samps_dist'][:]  # distance samples
	reds_mist = f['samps_red'][:]  # A(V) samples
	dreds_mist = f['samps_dred'][:]  # R(V) samples
	lnps_mist = f['samps_logp'][:]  # log-posterior of samples



	# Load SEGUE catalog data 
	df0 = pd.read_csv(catalog)


	# Information from the catalog

	obj_ids = np.array((),dtype=int)
	catalogs = ()
	ras = ()
	decs = ()
	epochs = ()
	plxs =()
	plx_errs = ()
	pmras = ()
	pmra_errs = ()
	pmdecs = ()
	pmdec_errs = ()




	teff = () 
	teff_low = ()
	teff_high = ()

	logg = ()
	logg_low = ()
	logg_high = ()

	
	fstar_probs = ()




	photos = None 

	for i, lab in enumerate(f['labels'][:]):

		

		df = df0[df0['source_id']==lab]	
		


		mag = [ [line[0], line[1], line[2], line[3], line[4]] for line in zip(df["gmag"], df["rmag"], df["imag"], df["zmag"], df["ymag"])]
		magerr = [ [line[0], line[1], line[2], line[3], line[4]] for line in zip(df["e_gmag"], df["e_rmag"], df["e_imag"], df["e_zmag"], df["e_ymag"]) ]
		mask = np.isfinite(magerr)
		phot, err = butils.inv_magnitude(np.array(mag), np.array(magerr), zeropoints = 3631.*1.e+9)

		if photos is None:
			photos = np.round(phot, decimals = 2)
			photo_errs = np.round(err, decimals = 2)

		else:
			photos = np.vstack((photos, np.round(phot[0,:], decimals = 2)))
			photo_errs = np.vstack((photo_errs, np.round(err[0,:], decimals = 2)))


		obj_ids = np.hstack((obj_ids, lab))
		catalogs = np.hstack((catalogs, "GaiaDR3"))
		plxs = np.hstack((plxs, df["parallax"].values))
		plx_errs = np.hstack((plx_errs, df["parallax_error"].values))
		ras = np.hstack((ras, df["ra"].values))
		decs = np.hstack((decs, df["dec"].values))
		epochs = np.hstack((epochs, df["ref_epoch"].values))
		pmras = np.hstack((pmras, df["pmra"].values))
		pmra_errs = np.hstack((pmra_errs, df["pmra_error"].values))
		pmdecs = np.hstack((pmdecs, df["pmdec"].values))
		pmdec_errs = np.hstack((pmdec_errs, df["pmdec_error"].values))




		if plot==True:
			fig, ax, parts = bplot.posterior_predictive(models_mist,  # stellar model grid
                                            idxs_mist[i],  # model indices
                                            reds_mist[i],  # A(V) draws
                                            dreds_mist[i],  # R(V) draws
                                            dists_mist[i],  # distance draws
                                            data=phot[0], data_err=err[0],  # data
                                            data_mask=mask[0],  # band mask
                                            offset=off_mist,  # photometric offsets
                                            psig=2.,  # plot 2-sigma errors
                                            labels=filt,  # filters 
                                            vcolor='blue',  # "violin plot" colors for the posteriors
                                            pcolor='black')  # photometry colors for the data
			figname = "../figs/sed_" + str(lab) + ".png"

			fig.savefig(figname)

		print('Best-fit chi2 (MIST):', chi2_mist[i])


		# Get Teff
		#print(labels_mist.dtype.names)
		params = labels_mist[idxs_mist[i]]
		lteff = [ line[5] for line in params ]
		q1, q2, q3 = butils.quantile(lteff, [0.159, 0.5, 0.841])  # 1sigma quantiles
		teff = np.hstack((teff,np.round(10**q2, decimals = 1)))
		teff_low = np.hstack((teff_low, np.round(10**q1, decimals = 1)))
		teff_high = np.hstack((teff_high, np.round(10**q3, decimals = 1)))
			

		# Probability of being F-star
		P_fstar, P_fstar_err = fstar_prob(10**q2, 10**q1, 10**q3)

		fstar_probs = np.hstack((fstar_probs, P_fstar))

		if P_fstar<0 or P_fstar>1.:
			print("Probability of Fstar for %i is out of range:  %f"%(lab, P_fstar))

		if P_fstar>0.5:
			print("Fstar!")
		

		loggs = [ line[6] for line in params ]
		q1, q2, q3 = butils.quantile(loggs, [0.159, 0.5, 0.841])
		logg = np.hstack((logg, np.round(q2, decimals = 2)))
		logg_low = np.hstack((logg_low, np.round(q1, decimals = 2)))
		logg_high = np.hstack((logg_high, np.round(q3, decimals = 2)))

		if plot==True:
			fig2, axes = bplot.cornerplot(idxs_mist[i], 
                             (dists_mist[i], reds_mist[i], dreds_mist[i]),
                             labels_mist,  # model labels
                             parallax=plx[0], parallax_err=plxerr[0],  # overplot if included
                             show_titles=True, color='blue', pcolor='orange',
                             fig=plt.subplots(11, 11, figsize=(55, 55)))  #  custom figure
			plt.xticks(fontsize=12)
			figname = "../figs/corner_" + str(lab) + ".png"
			fig2.savefig(figname)

		if i==300:
			break
	
	#print(photos)
	#print(np.shape(photos[:,0]))	
	outdata = { "obj_id": obj_ids,  "catalog": catalogs, "ra": ras, "dec": decs, \
	"epoch": epochs, "parallax":plxs, "parallax_error":plx_errs, "pmra":pmras, \
	"pmra_error":pmra_errs, "pmdec":pmdecs, "pmdec_error":pmdec_errs, "gFluxJy":photos[:, 0], \
	"rFluxJy":photos[:, 1], "iFluxJy":photos[:, 2], "zFluxJy":photos[:, 3], \
	"yFluxJy":photos[:, 4], "gFluxJy_err":photo_errs[:, 0], "rFluxJy_err":photo_errs[:, 1], \
	"iFluxJy_err":photo_errs[:, 2], "zFluxJy_err":photo_errs[:, 3], \
	"yFluxJy_err":photo_errs[:, 4], "probfstar": fstar_probs, \
	"teff":teff, "teff_low":teff_low, "teff_high":teff_high, "logg":logg, "logg_low":logg_low, "logg_high":logg_high}

	dfout = pd.DataFrame(outdata)
	dfout.to_csv(paramfile, index = False)

	return



	









def plot_param(paramfiles):

	cmap = plt.get_cmap("Oranges")
	nc = 5


	for paramfile in paramfiles:

		df = pd.read_csv(paramfile)

		t_p = df["teff"]
		t_p_err = 0.5 * (np.abs(df["teff_low"] - df["teff"]) + np.abs(df["teff_high"] - df["teff"]))	

		logg_p = df["logg"]
		logg_p_err = 0.5 * (np.abs(df["logg_low"] - df["logg"]) + np.abs(df["logg_high"] - df["logg"]))	




		t_p_errbar = np.median(t_p_err)
		logg_p_errbar = np.median(logg_p_err)


		t_o = df["oteff"]
		t_o_err = df["oteff_sig"]
		t_o_errbar = np.median(t_o_err)

		logg_o = df["ologg"]
		logg_o_err = df["ologg_sig"]
		logg_o_errbar = np.median(logg_o_err)


		feh_o = df["ofeh"]


		t_diff = t_p-t_o
		logg_diff = logg_p - logg_o

		fact = 3.

		# Temperature 
		c, low, upp = sigmaclip(t_diff, fact, fact)
		diff_mean = np.mean(c)
		diff_std = np.std(c) 

		plt.rcParams["font.size"] = 20
		fig, ax = plt.subplots(1, 1, figsize=(15, 7))

		for i, femax in enumerate([-2., -1.5, -1.0, -0.5, 0.0]):
			if i==0:
				femin = -5.0
			else:
				femin = femax-0.5
			filt = (feh_o>=femin) & (feh_o<femax)
			ax.errorbar(t_o[filt], t_diff[filt], linestyle = "", marker = "o", ms = 10, mec = cmap(i/nc), mfc = cmap(i/nc), label="%.1f<=[Fe/H]<%.1f"%(femin, femax))
		ax.text(3200., 500., r"Mean$=%.0f$, $\sigma=%.0f$ (3$\sigma$-clipped)"%(diff_mean, diff_std))
		x = np.arange(3000, 9000, 10)
		y = np.zeros_like(x)
		ax.plot(x, y, linewidth=2, linestyle = ":", color='gray' )
		ax.errorbar([3200], [-400.], xerr = t_o_errbar, yerr = t_p_errbar, linewidth=2, mew = 2, capsize=2, ecolor='gray')
		ax.set_xlabel(r"$T_{eff}$(SEGUE)[K]")
		ax.set_ylabel(r"$T_{eff}$(Predicted)$-T_{eff}$(SEGUE)[K]")	
		ax.set_xlim(3000., 9000)
		ax.set_ylim(-600., 600)
		ax.legend(loc=4, prop={"size":12})
		figname = ((paramfile.split("/"))[-1])[:-4] + "_teff.png"
		plt.savefig("../figs/" + figname)


		# logg
		c, low, upp = sigmaclip(logg_diff, fact, fact)
		diff_mean = np.mean(c)
		diff_std = np.std(c) 

		plt.rcParams["font.size"] = 20
		fig, ax = plt.subplots(1, 1, figsize=(15, 7))

		for i, femax in enumerate([-2., -1.5, -1.0, -0.5, 0.0]):
			if i==0:
				femin = -5.0
			else:
				femin = femax-0.5
			filt = (feh_o>=femin) & (feh_o<femax)
			ax.errorbar(logg_o[filt], logg_diff[filt], linestyle = "", marker = "o", ms = 10, mec = cmap(i/nc), mfc = cmap(i/nc), label="%.1f<=[Fe/H]<%.1f"%(femin, femax))

		ax.text(0.2, 2.0, r"Mean$=%.1f$, $\sigma=%.1f$ (3$\sigma$-clipped)"%(diff_mean, diff_std))
		x = np.arange(0.0, 5.0, 0.1)
		y = np.zeros_like(x)
		ax.plot(x, y, linewidth=2, linestyle = ":", color='gray' )
		ax.errorbar([0.2], [-2.0], xerr = logg_o_errbar, yerr = logg_p_errbar, linewidth=2, mew = 2, capsize=2, ecolor='gray')
		ax.set_xlabel(r"$\log g$(SEGUE)[K]")
		ax.set_ylabel(r"$\log g$(Predicted)$-\log g$(SEGUE)[K]")
		ax.set_xlim(0.0, 5.5)
		ax.set_ylim(-3., 3.)
		ax.legend(loc=1, prop={"size":12})

		figname = ((paramfile.split("/"))[-1])[:-4] + "_logg.png"
		plt.savefig("../figs/" + figname)




	return



def produce_readme_brutus():

	samplefile = "../outputs/params/PS1_Gaia_mist_woplx_ra0.0_1.0_dec-40.0_90.0_params.csv"
	df = pd.read_csv(samplefile)
	colnames = df.columns
	
	line = "="*100
	empty = " "*100 
	empty5 = "\n"*5
	empty10 = "\n"*10
	empty3 = "\n"*3

	line2 = " --------------------------------------"
	title = "        Candidate F-type stars"
	virsion = "            Virsion: 2.0"
	dt_now = datetime.datetime.now()
	datestring = "     Table created on %i/%i/%i"%(dt_now.day, dt_now.month, dt_now.year)	
	author = "      Authors: PFS obsproc team"
	section = " File description"


	f=open("../outputs/params/Readme", "w")
	for string in empty, line2, empty, title, empty, virsion, empty, datestring, empty, author, empty, line2, empty3, section, empty:
		f.write(string + "\n")

	for i, column in enumerate(colnames):
		if column == "obj_id":
			unit = " - "
			description = "Gaia source ID"
		elif column == "catalog":
			unit = " - "
			description = "Catalog name"
		elif column == "ra":
			unit = "deg" 
			description = "Right ascension " 
		elif column == "dec":
			unit = "deg" 
			description = "Declination"
		elif column == "epoch":
			unit = "year"
			description = "The reference epoch for ra, dec, parallax, pmra, and pmdec"
		elif column == "parallax":
			unit = "mas"
			description = "Parallax"
		elif column == "parallax_error": 
			unit = "mas"
			description = "Standard error of parallax"
		elif column == "pmra":
			unit = "mas/yr"
			description = "Proper motion in right ascension direction"
		elif column == "pmra_error":
			unit = "mas/yr"
			description = "Standard error of pmra"
		elif column == "pmdec":
			unit = "mas/yr"
			description = "Proper motion in declination direction"
		elif column == "pmdec_error":
			unit = "mas/yr"
			description = "Standard error of pmdec"
		elif column == "gFluxJy":
			unit = "nJy"
			description = "g-band flux (PS1DR2)"
		elif column == "rFluxJy":
			unit = "nJy"
			description = "r-band flux (PS1DR2)"
		elif column == "iFluxJy":
			unit = "nJy"
			description = "i-band flux (PS1DR2)"
		elif column == "zFluxJy":
			unit = "nJy"
			description = "z-band flux (PS1DR2)"
		elif column == "yFluxJy":
			unit = "nJy"
			description = "y-band flux (PS1DR2)"
		elif column == "gFluxJy_err":
			unit = "nJy"
			description = "g-band flux error (PS1DR2)"
		elif column == "rFluxJy_err":
			unit = "nJy"
			description = "r-band flux error (PS1DR2)"
		elif column == "iFluxJy_err":
			unit = "nJy"
			description = "i-band flux error (PS1DR2)"
		elif column == "zFluxJy_err":
			unit = "nJy"
			description = "z-band flux error (PS1DR2)"
		elif column == "yFluxJy_err":
			unit = "nJy"
			description = "y-band flux error (PS1DR2)"
		elif column == "probfstar":
			unit = " - "
			description = "Probability of being an F-type star"
		elif column == "teff":
			unit = "K"
			description = "Effective temperature"
		elif column == "teff_low":
			unit = "K"
			description = "16% quantile of the teff posterior distribution"
		elif column == "teff_high":
			unit = "K"
			description = "84% quantile of the teff posterior distribution"
		elif column == "logg":
			unit = "dex"
			description = "Surface gravity"
		elif column == "logg_low":
			unit = "dex"
			description = "16% quantile of the logg posterior distribution"
		elif column == "logg_high":
			unit = "dex"
			description = "84% quantile of the logg posterior distribution"
		elif column == "probfstar":
			unit = " - "
			description = "Probability of being an F-type star"
		else:
			continue

		if i+1 < 10:
			margin1 = " "*(19-len(column))
			margin2 = " "*(10-len(unit))
		else:
			margin1 = " "*(18-len(column))
			margin2 = " "*(10-len(unit))

		f.write(" %i "%(i+1) + column + margin1 + unit + margin2 + description + "\n")


	section = " Notes"
	note1 = "  - For ra, dec, epoch, parallax, pmra, pmdec and their associated uncertainties."
	note2 = "    see the Gaia archive website and the documentation."

	for string in empty3, section, empty, note1, note2, empty3:

		f.write(string + "\n")

	f.close()


def read_brutus_results_all(filelist):

	for filename in filelist:
		
		coord_range = ((((filename.split("/"))[-1]).split("woplx_"))[-1]).replace(".h5", "")
		catalog = "../comissioning/Nov2022/PS1DR2_GaiaEDR3_" + coord_range + ".csv"
		paramfile = filename.replace(".h5", ".csv")
		
		read_brutus_results(filename, catalog, paramfile)


		print(catalog)
		sys.exit()


		




if __name__ == "__main__":
	
	#calc_posterior_segue()

	#plot_hist_teff_logg()

	#sample_posterior()
	#sample_posterior_segue()
	#hdffile = "../outputs/MCMC/mcmc_brutus_2123448927880505344.h5"
	#read_hdf(hdffile)
	#test_get_magmodel3()
	#pickle_isofiles()
	#lookup_magmodel(lage, feh, mabs_r)

	#sample_posterior_segue()

	#mini, feh, eep = 1.2, -0.7, 450.
	#av, rv, dist = 0.5, 3.1, 2000.
	#mag, params, c = sedmaker.get_sed(mini=mini, feh=feh, eep=eep, av=av, rv=rv, dist=dist)
	#inputcatalog = "../comissioning/Nov2022/PS1DR1_GaiaEDR3_ra0.0_1.0_dec-40.0_90.0.csv"
	#outfile="../outputs/ps1_ra0-1_mist_woplx"
	#run_brutus(inputcatalog)
	#run_brutus_multiprocess()

	#filename="../outputs/segue_mist_woplx"
	#filename="../outputs/segue_mist"
	#paramfile = "../outputs/param_compare_woplx.csv"
	#paramfile = "../outputs/param_compare.csv"
	#filename = "../outputs/PS1_Gaia_mist_woplx_ra0.0_1.0_dec-40.0_90.0.h5"
	filelist = glob.glob("../outputs/PS1_Gaia_mist_woplx_ra?.*.h5")
	read_brutus_results_all(filelist)	

	#catalog = "../comissioning/Nov2022/PS1DR2_GaiaEDR3_ra0.0_1.0_dec-40.0_90.0.csv"
	#paramfile = "../outputs/params/PS1_Gaia_mist_woplx_ra0.0_1.0_dec-40.0_90.0_params.csv"
	#read_brutus_results(filename, catalog, paramfile)
	

	#produce_readme_brutus()

	#paramfiles = ["../outputs/param_compare.csv", "../outputs/param_compare_woplx.csv"]
	#plot_param(paramfiles)
	#view_models()
