from __future__ import print_function

import pandas as pd



from astropy.coordinates import SkyCoord
import astropy.units as u

from dustmaps.bayestar import BayestarQuery
from dustmaps.sfd import SFDQuery


import numpy as np
import matplotlib.pyplot as plt
import sys, glob, os, gc


import multiprocessing as mp
import time





code_dir = os.path.abspath(os.getcwd()) + "/../"
sys.path.append(code_dir)

import fstandards as fs

from fstandards.database import *



#from dustmaps.bayestar import BayestarWebQuery


def evec(band):

    # From Table 1 of
    #    https://iopscience.iop.org/article/10.3847/1538-4357/ab5362/pdf

    if band == "g":
        R = 3.518
    elif band == "r":
        R = 2.617
    elif band == "i":
        R = 1.971
    elif band == "z":
        R = 1.549
    elif band == "y":
        R = 1.263
    return(R)


def correct_extinction(plot = True):

    # Load bayestar
    bayestar = BayestarQuery(max_samples = 2)


    
    path = "../catalogs/"
    filelist = ["GEDR3_PS1DR1_RA86_4_DEC28_9_r0_65", \
                "GEDR3_PS1DR1_RA96_7_DEC33_8_r0_65", \
                "GEDR3_PS1DR1_RA102_2_DEC35_8_r0_65", \
                "GEDR3_PS1DR1_RA108_0_DEC37_6_r0_65", \
                "GEDR3_PS1DR1_RA133_6_DEC41_6_r0_65", \
                "GEDR3_PS1DR1_RA159_9_DEC39_6_r0_65", \
                "GEDR3_PS1DR1_RA182_9_DEC32_2_r0_65"]

    b = [0.0, 10., 15., 20., 40., 60., 80.]
    l = 180.

    for i, bb in enumerate(b):


        #coord = SkyCoord(l = l, b = bb, unit = 'deg', frame = 'galactic')
        #coord_icrs = coord.icrs
        #print(coord_icrs)


        df = pd.read_csv(path + filelist[i])

        print("Data size: %i"%(len(df)))

        #filt = (df0['parallax'] > 0.)
        
        #df = df0[filt]

        #print("Filtered data size: %i"%(len(df)))

        
        
        
        lgal = df['l'].values
        bgal = df['b'].values
        dist0 = 1./df['parallax'].values
        dist = np.where(dist0<0., 20., dist0)
        
        if plot == True:
            
            figname = "../Figs/Distances_l%.0f_b%.0f.png"%(l, bb)
            plothist(dist, 'Distance', "D[kpc]", figname)

        
        coords = SkyCoord(l = lgal * u.deg, b = bgal * u.deg, \
                          distance = dist * u.kpc, frame = 'galactic')
        

        ebv = bayestar(coords, mode = 'median')

        df['ebv'] = ebv
        df['g0'] = df['gmeanpsfmag'] + ebv * evec('g')
        df['r0'] = df['rmeanpsfmag'] + ebv * evec('r')
        df['i0'] = df['imeanpsfmag'] + ebv * evec('i')
        df['z0'] = df['zmeanpsfmag'] + ebv * evec('z')
        df['y0'] = df['ymeanpsfmag'] + ebv * evec('y')


        if plot == True:
            figname = "../Figs/Extinction_g_l%.0f_b%.0f.png"%(l, bb)
            plothist(ebv*evec('g'), 'Extinction', "A(g)", figname)
            #plothist(ebv*evec('i'), 'Extinction', "A(i)")
            #plothist(ebv*evec('g') - ebv*evec('i'), 'Extinction', "E(g-i)")


        df.to_csv(path + filelist[i] + "_ExtinctionCorrected.csv")
        
        
    return()


def correct_extinction_PS1(catalog, dustmap = "bayestar", exclude_lowGalLat = True):


    path = "../PS1_ExtinctionCorrected/"
    
    outcatname = path + ((catalog.split("/"))[-1])[:-4] + "_ExtinctionCorrected.csv"


    if os.path.isfile(outcatname):
        print("Output file exists. STOP.")
        return

    start_time = time.process_time()

    print("Processing ", catalog, " STARTED")

    df0 = pd.read_csv(catalog)

	
    if exclude_lowGalLat==True:
        filt = ((df0['b'] > 10.) | (df0['b'] < -10.)) & (df0['parallax']>0.)
        df = df0[filt].copy()
    else:
        filt = df0['parallax']>0.
        df = df0[filt].copy()

    del df0
    gc.collect()	


    print("Data size: %i"%(len(df)))
    
    l = df['l'].values
    b = df['b'].values
    dist = 1./(df['parallax'].values)

    coords = SkyCoord(l*u.deg, b*u.deg, distance = dist*u.kpc, \
                          frame = 'galactic')
    


    if dustmap=="SFD":
        sfd = SFDQuery()
        ebv = sfd(coords)
        ebv_source = ["SFD"] * len(ebv)

    elif dustmap=="bayestar":
        # 3D map
        bayestar = BayestarQuery(max_samples=2)
#,Unnamed: 0,ra,dec,specobjid,u,g,r,i,z,teffadop,teffadopunc,loggadop,loggadopunc,fehadop,fehadopunc,flag,snr,qa,elodiervfinal,elodiervfinalerr,spectypehammer,objID,raMean,decMean,l,b,nDetections,ng,nr,ni,nz,ny,gMeanPSFMag,rMeanPSFMag,iMeanPSFMag,zMeanPSFMag,yMeanPSFMag,rmag_difference,ra_epoch2000,dec_epoch2000,errHalfMaj,errHalfMin,errPosAng,source_id,ra_x,ra_error,dec_x,dec_error,parallax,parallax_error,parallax_over_error,pm,pmra,pmra_error,pmdec,pmdec_error,astrometric_n_good_obs_al,astrometric_gof_al,astrometric_chi2_al,astrometric_excess_noise,astrometric_excess_noise_sig,astrometric_params_solved,pseudocolour,pseudocolour_error,visibility_periods_used,ruwe,duplicated_source,phot_g_mean_flux,phot_g_mean_flux_error,phot_g_mean_mag,phot_bp_mean_flux,phot_bp_mean_flux_error,phot_bp_mean_mag,phot_rp_mean_flux,phot_rp_mean_mag,phot_bp_rp_excess_factor,bp_rp,dr2_radial_velocity,dr2_radial_velocity_error,dr2_rv_nb_transits,dr2_rv_template_teff,dr2_rv_template_logg,panstarrs1,sdssdr13,skymapper2,urat1,phot_g_mean_mag_error,phot_bp_mean_mag_error,phot_rp_mean_mag_error,phot_g_mean_mag_corrected,phot_g_mean_mag_error_corrected,phot_g_mean_flux_corrected,phot_bp_rp_excess_factor_corrected,ra_epoch2000_error,dec_epoch2000_error,ra_dec_epoch2000_corr,angDist,ebv,g0,r0,i0,z0,y0,fstars
        ebv_baye, flags = bayestar(coords, mode='median', return_flags = True)
        print(flags['reliable_dist'])
        # 2D map
        sfd = SFDQuery()
        ebv_sfd = sfd(coords)
        print(ebv_sfd)

        # Evaluation of the ebv values
        ebv = ()
        ebv_source = ()
        for j, val in enumerate(ebv_baye):
             if flags['reliable_dist'][j] == True:
                 ebv = np.hstack((ebv, val))
                 ebv_source = np.hstack((ebv_source,'Bayestar'))
             else:
                 ebv = np.hstack((ebv, ebv_sfd[j]))
                 ebv_source = np.hstack((ebv_source, 'SFD'))


    df.loc[:,'ebv'] = ebv
    df.loc[:,'g0'] = df['gmag'] + ebv * evec('g')
    df.loc[:,'r0'] = df['rmag'] + ebv * evec('r')
    df.loc[:,'i0'] = df['imag'] + ebv * evec('i')
    df.loc[:,'z0'] = df['zmag'] + ebv * evec('z')
    df.loc[:,'y0'] = df['ymag'] + ebv * evec('y')

    df.loc[:, 'ebv_source'] = ebv_source

    df.to_csv(outcatname)

    print(catalog, " COMPLETED") 

   
    return()



def correct_extinction_PS1_region(region, dustmap = "bayestar", exclude_lowGalLat = True):

	outfile = "../PS1Gaia_ExtinctionCorrected/PS1DR1_GaiaEDR3_ra%.1f_%.1f_dec%.1f_%.1f.csv"%(region)
	

	if os.path.isfile(outfile):
		print("Output file exists. STOP.")
		return

	start_time = time.process_time()

	print("Processing ra=%.1f-%.1f, dec=%.1f-%.1f"%(region), " STARTED")

	df0 = query_ps1_region(region[0], region[1], region[2], region[3])

	size_org = len(df0)


	if exclude_lowGalLat==True:
		filt = (df0['b'] > 10.) | (df0['b'] < -10.)  #) & (df0['parallax']>0.)
		df = df0[filt].copy()
	else:
		df = df0.copy()

	del df0
	gc.collect()

	print("|b|>0: original/selected = %i / %i"%(size_org, len(df)))


	#print(df["parallax"].isnull().values.any())


	for k in range(0,2):
		
		if k==0:
			plx_filt = df['parallax']>0.
		elif k==1:
			plx_filt = (df['parallax']<=0. ) | (df['parallax'].isnull())


		df_plx = df[plx_filt]

		print("Parallax iteration %i: Nstar=%i"%(k, len(df_plx)))

		l = df_plx['l'].values
		b = df_plx['b'].values
		dist = 1./(df_plx['parallax'].values)
		flags_dist = (df_plx['parallax_error']/np.abs(df_plx['parallax']) > 0.2 ) | (df_plx['parallax']<=0.)


		if k==0:
			coords = SkyCoord(l*u.deg, b*u.deg, distance = dist*u.kpc, frame = 'galactic')

			## Call the 3D dust map
			bayestar = BayestarQuery(max_samples=3)
			ebv_baye, flags = bayestar(coords, mode = 'percentile', pct = [16., 50., 84.], return_flags = True)

			## Evaluation of the ebv values
			ebv = ()
			ebv_source = ()
			flags_ebv = ()
			for j, val in enumerate(ebv_baye):
				ebv = np.hstack((ebv, val[1]))
				ebv_source = np.hstack((ebv_source,'Bayestar'))

				if flags['reliable_dist'][j] == True and (np.abs(val[0]-val[1])/val[1] < 0.2) and (np.abs(val[2]-val[1])/val[1] < 0.2) :
					flags_ebv = np.hstack((flags_ebv, False))
				else:
					flags_ebv = np.hstack((flags_ebv, True))


		elif k==1:

			coords = SkyCoord(l*u.deg, b*u.deg, frame = 'galactic')

			## Call the 2D dust map 
			sfd = SFDQuery()
			ebv = sfd(coords)
			ebv_source = ["SFD"] * len(ebv)
			flags_ebv = np.array([False] * len(ebv))

		df_plx.loc[:,'ebv'] = ebv
		df_plx.loc[:,'g0'] = df_plx['gmag'] + ebv * evec('g')
		df_plx.loc[:,'r0'] = df_plx['rmag'] + ebv * evec('r')
		df_plx.loc[:,'i0'] = df_plx['imag'] + ebv * evec('i')
		df_plx.loc[:,'z0'] = df_plx['zmag'] + ebv * evec('z')
		df_plx.loc[:,'y0'] = df_plx['ymag'] + ebv * evec('y')

		df_plx.loc[:, 'ebv_source'] = ebv_source
		df_plx.loc[:, 'flags_dist'] = flags_dist

		flags_ebv = flags_ebv.astype(bool)
		df_plx.loc[:, 'flags_ebv'] = flags_ebv
    
		if k==0:
			df_out1 = df_plx
		elif k==1:
			df_out2 = df_out1.append(df_plx)


	df_out2.to_csv(outfile)

	print("COMPLETED")

	return()


def correct_extinction_PS1_multiprocess():

	start_time = time.process_time()

	catalogs = glob.glob("../PS1_Gaia_crossmatch_stilts/PS1DR1_GaiaEDR3_ra*.csv")
	ncpus = 30
	print("Number of CPUs used: ", ncpus)
	pool = mp.Pool(ncpus)

	results = pool.map(correct_extinction_PS1, [catalog for catalog in catalogs])

	pool.close()

	end_time = time.process_time()

	elapse_time = end_time - start_time
	print("CPU time = " ,elapse_time)


	return()


def correct_extinction_PS1_query():


	start_time = time.process_time()

	ncpus = 60
	print("Number of CPUs used: ", ncpus)
	pool = mp.Pool(ncpus)

	ramins = np.arange(0., 360., 1.) 
	#ramins = np.arange(0., 2., 1.)
	regions = [(ramin, ramin+1., -40., 90.) for ramin in ramins ]


	results = pool.map(correct_extinction_PS1_region, [region for region in regions])

	pool.close()

	end_time = time.process_time()

	elapse_time = end_time - start_time
	print("CPU time = " ,elapse_time)
	
	return()


def query_ps1_region(ramin, ramax, decmin, decmax):

    with get_connection() as conn:
        with conn.cursor() as cur:
            #cur.execute('CREATE TABLE test_match AS SELECT gp.*, ps.* FROM gaia_ps1_bestneighbour gp JOIN ps1dr1 ps ON gp.original_ext_source_id=ps."objID"')
            #cur.execute('SELECT gp.* FROM gaia_ps1_bestneighbour gp')
            cur.execute('SET max_parallel_workers_per_gather=4')
            df = pd.read_sql(sql='SELECT gp."objID", source_id, ref_epoch, ra, dec, l, b, parallax, parallax_error, pmra, pmra_error, pmdec, pmdec_error, ruwe, phot_g_mean_mag, phot_bp_mean_mag, phot_rp_mean_mag, gmag, e_gmag, rmag, e_rmag, imag, e_imag, zmag,e_zmag, ymag, e_ymag FROM gaia_ps1_crossmatch gp WHERE ra>=%f AND ra<%f AND dec>=%f AND dec<%f;'%(ramin, ramax, decmin, decmax), con=conn)
            
			#colnames = [col.name for col in cur.description]
            print("FINISHED.")
    return(df)



def correct_extinction_segue(plot = False):

    # Load bayestar
    bayestar = BayestarQuery(max_samples = 2)


    
    path = "../../../catalogs/MAST_PS1_crossmatch/"
    filename = "sppParams_PS1DR1_crossmatch_GaiaEDR3.csv"


    df = pd.read_csv(path + filename)

    print("Data size: %i"%(len(df)))

    #filt = (df0['parallax'] > 0.)
    
    #df = df0[filt]

    #print("Filtered data size: %i"%(len(df)))


    #coord = SkyCoord(l = l, b = bb, unit = 'deg', frame = 'galactic')
    #coord_icrs = coord.icrs
    #print(coord_icrs)
    
    #coords = SkyCoord(ra = df['ra'].values, dec = df['dec'].values, \
    #                  unit = 'deg', frame = 'icrs')

        
    #lgal = coords.galactic.l
    #bgal = coords.galactic.b
    dist0 = 1./df['parallax'].values
    dist = np.where(dist0<0., 20., dist0)
        
    if plot == True:
            
        figname = "../Figs/Distances_Segue.png"
        plothist(dist, 'Distance', "D[kpc]", figname)

        
    coords = SkyCoord(l = df['l'].values * u.deg, b = df['b'].values * u.deg, \
                      distance = dist * u.kpc, frame = 'galactic')

    ebv = bayestar(coords, mode = 'median')
    
    df['ebv'] = ebv
    df['g0'] = df['gMeanPSFMag'] + ebv * evec('g')
    df['r0'] = df['rMeanPSFMag'] + ebv * evec('r')
    df['i0'] = df['iMeanPSFMag'] + ebv * evec('i')
    df['z0'] = df['zMeanPSFMag'] + ebv * evec('z')
    df['y0'] = df['yMeanPSFMag'] + ebv * evec('y')
    #df['lgal'] = lgal
    #df['bgal'] = bgal


    if plot == True:
        figname = "../Figs/Extinction_g_Segue.png"
        plothist(ebv*evec('g'), 'Extinction', "A(g)", figname)
        #plothist(ebv*evec('i'), 'Extinction', "A(i)")
        #plothist(ebv*evec('g') - ebv*evec('i'), 'Extinction', "E(g-i)")


    df.to_csv(path + filename[:-4] + "_ExtinctionCorrected.csv")
        
        
    return()










def plothist(a, dataclass, xlab, figname):

    fig, ax = plt.subplots(1, 1)
    if dataclass == "Extinction":
        ax.hist(a)
        ax.set_xlabel(xlab)
    elif dataclass == "Distance":
        bins = np.arange(0., 25., 0.2)
        ax.hist(a, bins = bins)
        ax.set_xlabel(xlab)

    plt.savefig(figname)

    return()
    


if __name__ == "__main__":
	correct_extinction_PS1_query()

#correct_extinction_PS1_multiprocess()
        
#correct_extinction_segue(plot = False)

#catalogs = glob.glob("../PS1/*.csv")
#for catalog in catalogs:
#    path = "../PS1_ExtinctionCorrected/"

#    outcatname = path + ((catalog.split("/"))[-1])[:-4] + "_ExtinctionCorrected.csv"
#    if os.path.isfile(outcatname):
#        continue

#    correct_extinction_PS1(catalog)

