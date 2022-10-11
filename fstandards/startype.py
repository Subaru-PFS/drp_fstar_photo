import numpy as np
import emcee
import matplotlib.pyplot as plt
import pandas as pd
from scipy import interpolate
import sys, os, time, resource
import corner
from astropy import units as u
from astropy.coordinates import SkyCoord


	
from dustmaps.bayestar import BayestarQuery
from dustmaps.sfd import SFDQuery



# See https://docs.astropy.org/en/stable/api/astropy.coordinates.galactocentric_frame_defaults.html
from astropy.coordinates import galactocentric_frame_defaults
_ = galactocentric_frame_defaults.set('v4.0')

state = galactocentric_frame_defaults.get_from_registry("v4.0")
r_sun = (state["parameters"]["galcen_distance"] * 1000.).value
z_sun = (state["parameters"]["z_sun"]).value

from multiprocessing import Pool, RawArray
import multiprocessing as mp


#state = galactocentric_frame_defaults.get_from_registry("v4.0")


# --- Should be removed ---- 
# Path to the isochrone files
#isopath = "../isochrones/PARSEC/"
#df_isochrone = pd.read_pickle(isopath + "PS1.pkl")
#isoarr0 = df_isochrone.to_numpy()


import brutus





# A globa dictionary storing the variables passed from the initialize.
#var_dict = {}


def log_prob_test(x, mu, cov):

	diff = x - mu

	return -0.5 * np.dot(diff, np.linalg.solve(cov, diff))



def log_like(params, magobs, cov):
	
	# The extinction vector R for PS1 grizy from Table 1 of Green et al. 2019, ApJ, 887, 27
	extvect = np.array([3.518, 2.617, 1.971, 1.549, 1.263])

	# 50% completeness limit from https://outerspace.stsci.edu/display/PANSTARRS/PS1+Photometric+Depth#PS1PhotometricDepth-DepthofGalaxiesinComparisontoPointSources
	#           => should be revised to the selection function 
	maglim = np.array([23.2, 23.2, 23.1, 22.3, 21.2]) 

	# Obtain model magnitudes:
	lage, feh, mass, ebv, dmod = params

	modelabsmags, teff, logg = get_magmodel_from_isochrone(lage, feh, mass)	

	if np.count_nonzero(np.isnan(modelabsmags))!=0:
		return -1.*float('inf')

	modelmags = modelabsmags + ebv * extvect + dmod


	# Incorpolating survey selection function 
	for i, modelmag in enumerate(modelmags):
		if modelmag > maglim[i]:
			return -1.*float('inf') 


	modelplx = 10**(-0.2 * (dmod + 5.)) 
	model = np.hstack((modelmags, modelplx))

	diff = magobs - model 

	return -0.5 * np.dot(diff, np.linalg.solve(cov, diff))
	

def log_prior(params, galcoord, ebvmax):

	lage, feh, mass, ebv, dmod = params

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

	return lp_age(lage, r_star, z_star) + lp_feh(feh, r_star, z_star) + lp_mass(mass) + lp_ebv(ebv, ebvmax) + lp_dmod(dmod, r_star, z_star) 


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
	

def log_prob(params, galcoord, ebvmax, magobs, cov):
	
	lp = log_prior(params, galcoord, ebvmax)

	return log_like(params, magobs, cov) + lp

	




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

		if ct  > 10:
			break
		ct = ct + 1

	return



def sample_posterior(starname, l, b, obsmags, obsmagerrs, plx, plxerr):
	

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



	# Initial guess of the free parameters
	nwalkers = 32

	lage0 = 8.95
	feh0 = -1.0
	mass0 = 1.0
	ebv0 = 0.01
	dmod0 = 11.



	param0 = [lage0, feh0, mass0, ebv0, dmod0]
	ndim = len(param0)


	# fraction of each parameter for one step of the walkers
	stepfrac = [0.1, 0.1, 0.1, 0.01, 0.1]


	# Setup the backend 
	filename = "../MCMC/mcmc_" + str(starname) + ".h5"
	backend = emcee.backends.HDFBackend(filename)	
	backend.reset(nwalkers, ndim)


	p0 = get_ini(param0, nwalkers, stepfrac)


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


	with Pool(64) as pool:
	#with Pool(processes = 50, initializer = init_worker, initargs = (X, X_shape)) as pool:

		sampler = emcee.EnsembleSampler(nwalkers, ndim, log_prob, args = [galcoord, ebvmax, np.hstack((obsmags, plx)), cov], pool = pool, backend = backend)

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
	
		sampler.run_mcmc(p0, 200000)

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

def read_hdf(hdffile):

	starname = ((hdffile.split("_"))[1]).strip(".h5")

	reader = emcee.backends.HDFBackend(hdffile)
	samples = reader.get_chain(flat = True)
	
	labels = ["log(Age)", "[Fe/H]", "mass", "E(B-V)", "dmod" ]
	corner.corner(samples, labels = labels)
	plt.savefig("../figs/emcee_" + starname + "_corner.png")

	

	ncpus = 10
	print("Number of CPUs used: ", ncpus)
	pool = mp.Pool(ncpus)

	args = [sample[:3] for sample in samples[5000:, :] ]

	results = pool.map(wrap_get_magmodel, args)

	pool.close()

	teffs = [ line[1] for line in results ]

	fig, ax = plt.subplots(1, 1)	
	ax.hist(teffs, 100, color="k", histtype="step")
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








if __name__ == "__main__":
	
	#calc_posterior_segue()

	#plot_hist_teff_logg()

	#sample_posterior()
	#sample_posterior_segue()
	#hdffile = "../MCMC/mcmc_2122426107071326208.h5"
	#read_hdf(hdffile)
	#test_get_magmodel3()
	#pickle_isofiles()
	#lookup_magmodel(lage, feh, mabs_r)

