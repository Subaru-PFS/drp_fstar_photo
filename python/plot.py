import pandas as pd
import matplotlib.pyplot as plt
import matplotlib
from pylab import meshgrid 
import astropy.coordinates as coord
import astropy.units as u
from astropy.io import ascii
from astropy.coordinates import SkyCoord, FK5
import numpy as np
import sys, glob, os
import matplotlib.cm as cm


def func_r(bprp, G):

    a = [-0.12879, 0.24662, -0.027464, -0.049465]

    func = a[0] + a[1] * bprp + a[2] + bprp**2 + a[3] * bprp**3
    
    r = G - func

    return(r)




def func_g(bprp, G):

    #https://gea.esac.esa.int/archive/documentation/GDR2/Data_processing/chap_cu5pho/sec_cu5pho_calibr/ssec_cu5pho_PhotTransf.html

    a = [0.13518, -0.46245, -0.25171, 0.021349]
    
    func = a[0] + a[1] * bprp + a[2] * bprp**2 + a[3] * bprp**3
    
    g = G - func

    return(g)



def func_i(bprp, G):

    #https://gea.esac.esa.int/archive/documentation/GDR2/Data_processing/chap_cu5pho/sec_cu5pho_calibr/ssec_cu5pho_PhotTransf.html
    
    a = [-0.29676, 0.64728, -0.10141]
    
    func = a[0] + a[1] * bprp + a[2] * bprp**2 
    
    i = G - func

    return(i)





def plot_gum(glim, gilim):


    path = "../catalogs/"
    filelist = ["GUM_RA86.4_DEC28.9_r0.65.csv", \
                "GUM_RA108.0_DEC37.6_r0.65.csv", \
                "GUM_RA133.6_DEC41.6_r0.65.csv", \
                "GUM_RA159.9_DEC39.6_r0.65.csv", \
                "GUM_RA182.9_DEC32.2_r0.65.csv"]

    b = [0.0, 20., 40., 60., 80.]
    l = 180.

    fig, ax = plt.subplots(1, len(filelist), figsize = (20, 3))
    
    for i, catalog in enumerate(filelist):

        ra_cen = np.float(((catalog.split("_"))[1])[2:])
        dec_cen = np.float(((catalog.split("_"))[2])[3:])
        
        
        df0 = pd.read_csv(path + catalog)

        
        

        bprp = df0['mag_bp']-df0['mag_rp']

        gmag = func_g(bprp, df0['mag_g'])

        imag = func_i(bprp, df0['mag_g'])


        
        df0['gmag'] = gmag
        df0['imag'] = imag
        df0['gi'] = gmag - imag


        
        df0.to_csv(catalog.replace(".csv", "_MockPS1Photo.csv"))

        
        # Get SDSS12 relation ship (temporarily!! this should be PS1)

        

        
        filt = (df0['gi'] > gilim[0]) & (df0['gi'] < gilim[1]) & \
            (df0['gmag'] < glim[1]) & \
            (df0['gmag'] > glim[0])

        df = df0[filt]


        # get the number of unique entry
        unique_id = ()
        for id in df["source_extended_id"]:
            unique_id = np.append(unique_id, id[:17])

        
        df["source_extended_id_Unique"] = unique_id

        N_selected = len(df)
        

        N_unique = df["source_extended_id_Unique"].nunique()

        print("Number of stars %i  => %i "%(len(df0), N_unique))

        
        print("Number of UNIQUE object %i"%(N_unique) )
       
        #print(" (Teff not available: %i)"%(df['teff'].isnull().sum()))


        
        filt = (df['teff'] > 6000.) & (df['teff'] < 7500.)
        df_F = df[filt]
        #print(" (Teff not available: %i)"%(df_F['teff'].isnull().sum()))
        print("Fraction of genuine Fstars: %.1f"%(len(df_F)/N_unique))


        filt = (df['teff'] > 6000.) & (df['teff'] < 7500.) & (df['feh'] < -1.0)
        df_F_MP = df[filt]

        #print(" (Teff not available: %i)"%(df_F['teff'].isnull().sum()))
        print("Fraction of genuine F MP stars: %.1f"%(len(df_F_MP)/N_unique))


        

        ax[i].plot(df['ra'], df['dec'], linestyle = "", marker = ".", ms = 3, \
                   label = "PS1-selected")
        ax[i].plot(df_F['ra'], df_F['dec'], linestyle = "", marker = "o", \
                   mfc = "none", ms = 5, label = "F-stars")
        ax[i].plot(df_F_MP['ra'], df_F_MP['dec'], linestyle = "", \
                   marker = "s", mfc = "none", ms = 10, \
                   label = "F, MP stars")
            
        w = 0.9
        ax[i].set_xlim(ra_cen - w, ra_cen + w)
        ax[i].set_ylim(dec_cen - w, dec_cen + w)

        ax[i].text(ra_cen-0.85, dec_cen-0.85, "$N(%.1f<g-i<%.1f)=%i$"%(gilim[0], gilim[1], N_unique))
        ax[i].text(ra_cen-0.85, dec_cen+0.78, "$%.0f<g<%.0f$"%(glim[0], glim[1]))
        
        ax[i].set_xlabel("RA")
        if i==0:
            ax[i].set_ylabel("DEC")
        ax[i].set_title(r"$(l, b)=(%.0f^{\circ}, %.0f^{\circ})$"%(l, b[i]))

        #ax[i].axis('equal')
                
    ax[0].legend()

    plt.tight_layout()
    plt.savefig("mock_gaia.png")
    plt.show()

    return()



def plot_segue():

	import matplotlib.cm as cm


	catalog = "../catalogs/MAST_PS1_crossmatch/sspParams_PS1DR1_crossmatch_GaiaEDR3_ExtinctionCorrected_selected.csv"
	
	df = pd.read_csv(catalog)


	fig, ax = plt.subplots(1, 3, figsize=(15, 5))


	snrmins = [ 0., 20., 50., 100., 150.]
	snrmaxs = [20., 50., 100., 150., 500.]
	labs = [r'$0<$SNR$\leq 20$', r'$20<$SNR$\leq 50$', r'$50<$SNR$\leq 100$', r'$100<$SNR$\leq 150$', r'$150<$SNR']


	nsnr = len(snrmins)

	for i, snrmin in enumerate(snrmins):
		
			filt = (df['snr']>snrmin) & (df['snr']<=snrmaxs[i])
			ax[0].plot(df['g0'][filt]-df['r0'][filt], df['gMeanPSFMag'][filt], label = labs[i], color = cm.magma(i/nsnr), marker = "o", linestyle = "", ms = 1, alpha =0.7)
	#ax[0].hist2d(df['bp_rp'], df['phot_g_mean_mag'], bins = [np.arange(0.1, 1.6, 0.01), np.arange(13.0, 20.6, 0.05)], cmap = cm.Blues)
	ax[0].set_xlim(0.0, 1.0)
	ax[0].set_ylim(20.5, 13.0)
	ax[0].set_xlabel(r"$(g-r)_0$")
	ax[0].set_ylabel("$g$")

	ax[0].legend()

	fehmins = [ -6., -2., -1., 0.]
	fehmaxs = [-2., -1., 0., 1.]
	labs = [r'[Fe/H]$<-2$', r'$-2\leq$[Fe/H]$<-1$', r'$-1\leq$[Fe/H]$<0$', r'$0\leq$[Fe/H]$<1$']	

	nfehs = len(fehmins)

	for i, fehmin in enumerate(fehmins):
			fehmax = fehmaxs[i]
			filt = (df['fehadop']>fehmin) & (df['fehadop']<=fehmax)
	
			ax[1].plot(df['teffadop'][filt], df['loggadop'][filt], label = labs[i], color = cm.magma(i/nfehs), marker = "o", linestyle = "", ms = 1, alpha = 0.7)
			#ax[1].hist2d(df['teffadop'], df['loggadop'], bins = [np.arange(3000.,9000., 100.), np.arange(-0.5, 5.5, 0.1)], cmap = cm.Blues)
		
	ax[1].legend()
	ax[1].set_xlim(9000., 4000.)
	ax[1].set_ylim(5.0, -0.5)
	ax[1].set_xlabel(r"$T_{eff}$ [K]")
	ax[1].set_ylabel(r"$\log g$")

	filt = (df['parallax'] > 0. ) & (df['parallax'] < 10.)
	ax[2].hist(df['parallax'][filt], bins = np.arange(0.0, 5., 0.05))
	ax[2].set_xlabel("Parallax [mas]")

	plt.savefig("../figs/segue.png")
	return


def plot_gaia_ps1_cmatch(glim, gilim, extinctioncorr = True):

    path = "../catalogs/"


    filelist = \
        ["GEDR3_PS1DR1_RA86_4_DEC28_9_r0_65_ExtinctionCorrected.csv", \
         "GEDR3_PS1DR1_RA96_7_DEC33_8_r0_65_ExtinctionCorrected.csv", \
         "GEDR3_PS1DR1_RA102_2_DEC35_8_r0_65_ExtinctionCorrected.csv", \
         "GEDR3_PS1DR1_RA108_0_DEC37_6_r0_65_ExtinctionCorrected.csv", \
         "GEDR3_PS1DR1_RA133_6_DEC41_6_r0_65_ExtinctionCorrected.csv", \
         "GEDR3_PS1DR1_RA159_9_DEC39_6_r0_65_ExtinctionCorrected.csv", \
         "GEDR3_PS1DR1_RA182_9_DEC32_2_r0_65_ExtinctionCorrected.csv"]
    
    

    b = [0.0, 10., 15., 20., 40., 60., 80.]
    l = 180.

    fig, ax = plt.subplots(1, len(filelist), figsize = (23, 3))
    
    for i, catalog in enumerate(filelist):

        #ra_cen = np.float(((catalog.split("_"))[2])[2:])
        #dec_cen = np.float(((catalog.split("_"))[3])[3:])
        
        
        df0 = pd.read_csv(path + catalog)

        if extinctioncorr == False:
            gi = df0['gmeanpsfmag'] - df0['imeanpsfmag']
            figname = "gaia_ps1_cmatch.png"
        else:
            gi = df0['g0'] - df0['i0']
            figname = "gaia_ps1_cmatch_ExtinctionCorrected.png"
        
        filt = (gi > gilim[0]) & (gi < gilim[1]) & \
            (df0['gmeanpsfmag'] < glim[1]) & \
            (df0['gmeanpsfmag'] > glim[0])

        df = df0[filt]

        N_selected = len(df)
        print("Number of stars %i  => %i "%(len(df0), len(df)))


        print(" (bp_rp not available: %i)"%(df['bp_rp'].isnull().sum()))
        
    

        
        ax[i].plot(df['ra'], df['dec'], linestyle = "", marker = ".")
        w = 0.9

        ra_cen= np.mean(df0['ra'])
        dec_cen = np.mean(df0['dec'])

        ax[i].text(ra_cen-0.85, dec_cen-0.85, "$N(%.1f<g-i<%.1f)=%i$"%(gilim[0], gilim[1], N_selected))
        ax[i].text(ra_cen-0.85, dec_cen+0.78, "$%.0f<g<%.0f$"%(glim[0], glim[1]))
        
        ax[i].set_xlim(ra_cen - w, ra_cen + w)
        ax[i].set_ylim(dec_cen - w, dec_cen + w)
        ax[i].set_xlabel("RA")
        if i==0:
            ax[i].set_ylabel("DEC")
        ax[i].set_title(r"$(l, b)=(%.0f^{\circ}, %.0f^{\circ})$"%(l, b[i]))
        #ax[i].axis('equal')

    plt.tight_layout()
    
    plt.savefig(figname)
    plt.show()
    return


def plot_healpix_Gaia(catalog):

    import healpy as hp
    
    df = pd.read_csv(catalog)
    indx = df['healpix6'].values
    hpxmap = df['srcdens'].values
    hp.mollview(hpxmap, nest = True, coord = "E", title = "N stars per arcmin")

    plt.savefig("../figs/healpix_Gaia_test.png")

    return


def filter_data(catalog, selection = "classifier"):

	df = pd.read_csv(catalog)

	if selection == "color":
		filt = (df['g0']-df['i0'] > 0.0) & (df['g0']-df['i0']<0.4)
		#filt = (df['g0']-df['i0'] > 0.0) & (df['g0']-df['i0']<0.4) & \
		#       (df['imag']-df['iKmag']<0.05)
	elif selection == "classifier":
		filt = (df['ProbFstar']>0.5) & (df['gmag']>15.) & (df['gmag']<20.) & ((df['b']>10.) | (df['b']<-10.))
		#filt = (df['imag']-df['iKmag']<0.05) & (df['ProbFstar']>0.5) & (df['gmag']>15.) & (df['gmag']<20.) & ((df['b']>10.) | (df['b']<-10.))  
		ra = np.append(ra, df['ra'][filt])
		dec = np.append(dec, df['dec'][filt])







def plot_sky_density_healpix(catalogs, selection = "classifier", zscalelog = False):

    import healpy as hp

    ra = np.array([])
    dec = np.array([])


    for catalog in catalogs:

        
        df = pd.read_csv(catalog)


        if selection == "color":
            filt = (df['g0']-df['i0'] > 0.0) & (df['g0']-df['i0']<0.4)
            #filt = (df['g0']-df['i0'] > 0.0) & (df['g0']-df['i0']<0.4) & \
            #       (df['imag']-df['iKmag']<0.05)
        elif selection == "classifier":
            filt = df['ProbFstar']>0.5
            #filt = (df['imag']-df['iKmag']<0.05) & (df['ProbFstar']>0.5) & (df['gmag']>15.) & (df['gmag']<20.) & ((df['b']>10.) | (df['b']<-10.))  

        ra = np.append(ra, df['ra'][filt])
        dec = np.append(dec, df['dec'][filt])


    phis = ra * np.pi/180.
    thetas = -1. * dec * np.pi/180. + np.pi/2.
    nside = 64

    res = hp.nside2resol(nside, arcmin=True) / 60 
    npix = hp.nside2npix(nside) 
    print("Approximate resolution for NSIDE = %d: %.2f, Npix=%i"%(nside,res, npix))  
    indices = hp.ang2pix(nside, thetas, phis)

    npix = hp.nside2npix(nside)
    hpxmap = np.zeros(npix, dtype = np.float)
    for i in range(len(ra)):
        hpxmap[indices[i]] += 1

    
    hp.mollview(hpxmap, coord = "E", title = "N stars in HEALPix with resolution = %.1f deg"%(res), max = 300)
    # If res = 0.916 deg, the area per pixel is 0.84 deg^2. Multiplying by 1.6 with gives N per 1.3 deg^2 (PFS FoV).  


    plt.savefig("../figs/healpix_" + selection + ".png")  

    print("Total number of pixels: %i"%(len(hpxmap)))
    print("Fraction of pixels with < 5 stars: %.1f percent"%((len(hpxmap[hpxmap<5])/len(hpxmap)*100.)))

    return



def plot_sky_density_ps1(catalogs, zscalelog = False):
    
    ra = np.array([])
    dec = np.array([])
    
    
    for catalog in catalogs:
        
        df = pd.read_csv(catalog)


        coords = SkyCoord(ra = df['raMean'].values, dec = df['decMean'].values, \
unit = 'deg', frame = 'icrs')


        lgal = coords.galactic.l.value
        bgal = coords.galactic.b.value

        
        filt = ((bgal > 10.) | (bgal < -10)) & (df['g0']-df['i0']>0.0) \
& (df['g0']-df['i0']<0.4) & (df['gMeanPSFMag']<20.) #& (1./df['parallax'] > 0.2)
        ra = np.append(ra, df['raMean'][filt])
        dec = np.append(dec, df['decMean'][filt])
        #print(ra)
        #sys.exit()

        #    lab = "18 < G < 20, 0.2 < (BP-RP) < 0.8"
        #elif catalogname == "GaiaEDR3xPS1":
        #    filt = (df['dec'] > -20.) & (df['phot_g_mean_mag'] > 18.) & (df['phot_g_mean_mag'] < 20.)
        #    ra = np.append(ra, df['ra'][filt])
        #    dec = np.append(dec, df['dec'][filt])
        #    lab = "18 < G < 20, r<20., 0.0 < g-i < 0.4"
            
        #elif catalogname == "GUM":
        #    filt = (df['dec'] > -20.) #& (df['distance'] > 200.)
        #    ra = np.append(ra, df['ra'][filt])
        #    dec = np.append(dec, df['dec'][filt])
        #    lab = "18 < G < 20, 0.2 < (BP-RP) < 0.8"
    
    #ramin, ramax, decmin, decmax = get_radec_range_from_file(file)
    ramin = 0.0
    ramax = 300. 
    decmin = -40. 
    decmax = 89.
    
    wdeg = 1.
    xedges = np.arange(ramin, ramax, wdeg)
    yedges = np.arange(decmin, decmax, wdeg)
    x = ra * np.cos(dec*np.pi/180.)
    y = dec
    H, xedges, yedges = np.histogram2d(x, y, bins=(xedges, yedges))
    
   


    # Calculate the number of pixels that fall in -10<b<10

    ct = 0
    for xx in xedges:
        xbincen = xx + 0.5 * wdeg
        for yy in yedges:
            ybincen = yy + 0.5 * wdeg 
        
            c = SkyCoord(ra = xbincen, dec = ybincen, unit = 'deg', frame = 'icrs')
            bb = c.galactic.b.value
            if bb >= -10. and bb <= 10.:
                ct = ct + 1
    
   
    print("Total number of pixels:",np.size(H))
    print("Total number of pixels falling in -10<=b<=10:", ct)
    print("Number of xedges:", np.size(xedges))
    print("Number of yedges:", np.size(yedges))
    H = H.T 

    histH = np.reshape(H, np.size(H))
    
    
    frac_zerostars = ((np.size(H[H == 0])-ct) / (np.size(H)-ct)) * 100.
    
    print("Npix(0stars)/Npix(all)=%.1f"%(frac_zerostars) + r"%")
    
    frac_zerostars = ((np.size(H[H < 5])-ct) / (np.size(H)-ct)) * 100.
    
    print("Npix(<5stars)/Npix(all)=%.1f"%(frac_zerostars) + r"%")
    
    
  

    fig = plt.figure(figsize=(15,8))
    ax = fig.add_subplot(111)

    #title=r" $18<G<20$, $0.2<(BP-RP)<0.7$ (F-star candidates), distance > 500 pc $"

    if zscalelog == True:
        lH = np.log10(H)
        H = lH
        clab = r"$\log N$"
        vmax = np.min([np.max(H), np.log10(300.)])
    else:
        clab = r"N stars per 1.3 deg$^{2}$"
        
        vmax = np.min([np.max(H), 300.])
        print(vmax)
    im = ax.imshow(H, interpolation='nearest', cmap='jet',origin='lower', \
            extent=[xedges[0], xedges[-1], yedges[0], yedges[-1]], \
            vmin=0.0,vmax=vmax)


    # Plot coordinate grid 
    
    for b in np.arange(-90., 90., 20.):
        larr = np.arange(0., 360., 1.)
        barr = np.zeros_like(larr)
        barr.fill(b)
     
        gal = SkyCoord(larr[:], barr[:], frame='galactic', unit=u.deg) 
        eq = gal.transform_to('fk5')
      
     
        if b == 0:
            ls = "-"
        else:
            ls = "--"
        filt = eq.ra.value < 360.


        ax.plot(eq.ra[filt], eq.dec[filt], linestyle = "", marker = ".", ms = 2, color = "white", alpha = 0.6)
        
        # label 
        l = 180.
        gallab = SkyCoord(l, b, frame = 'galactic', unit = u.deg)
        eqlab = gallab.fk5
        ax.text(eqlab.ra.value, eqlab.dec.value, "b=%.0f"%(b), color = "white")
    

    #ax.set_title(title)
    ax.set_xlabel("RA (deg)")
    ax.set_ylabel("DEC (deg)")
    ax.set_xlim(ramin, ramax)
    ax.set_ylim(-30., 89)
    fig.colorbar(im,orientation='horizontal',label= clab)

    plt.savefig("../figs/PS1_color-selected_skydensity_equalArea.png")
   
    
    
    #fig, bx = plt.subplots(1, 1)
    
    #xmin = 0
    #xmax = 500
    #bx.hist(histH, bins = 50, range = (xmin, xmax), label = lab)
    #ymax = np.max(histH)
    #ymin = 0.
    #bx.text(xmax - (xmax - xmin) * 0.05, ymax - (ymax - ymin)*0.05, \
    #        "Npix(<5stars)/Npix(all)=%.1f"%(frac_zerostars), ha = 'right')
    
    #bx.set_xlabel(r"N stars per 1.3 deg$^{2}$")
    #bx.set_ylabel(r"N RA-DEC pixels")
    #bx.legend()
    #plt.savefig(catalogname + "_1Dhist.png")
 
   
    return()




def plot_sky_density_ps1_fstars(catalogs, zscalelog = False):
    
    ra = np.array([])
    dec = np.array([])
    
    
    for catalog in catalogs:
        
        df = pd.read_csv(catalog)

        filt = ((df['b'] > 10.) | (df['b'] < -10)) & (df['ProbFstar']>0.5) \
	& (df['gmag']>15.) & (df['gmag']<20.) #& (1./df['parallax'] > 0.2)
        ra = np.append(ra, df['ra'][filt])
        dec = np.append(dec, df['dec'][filt])
        #print(ra)
        #sys.exit()

        #    lab = "18 < G < 20, 0.2 < (BP-RP) < 0.8"
        #elif catalogname == "GaiaEDR3xPS1":
        #    filt = (df['dec'] > -20.) & (df['phot_g_mean_mag'] > 18.) & (df['phot_g_mean_mag'] < 20.)
        #    ra = np.append(ra, df['ra'][filt])
        #    dec = np.append(dec, df['dec'][filt])
        #    lab = "18 < G < 20, r<20., 0.0 < g-i < 0.4"
            
        #elif catalogname == "GUM":
        #    filt = (df['dec'] > -20.) #& (df['distance'] > 200.)
        #    ra = np.append(ra, df['ra'][filt])
        #    dec = np.append(dec, df['dec'][filt])
        #    lab = "18 < G < 20, 0.2 < (BP-RP) < 0.8"
    
    #ramin, ramax, decmin, decmax = get_radec_range_from_file(file)
    ramin = 0.0
    ramax = 360. 
    decmin = -40. 
    decmax = 89.
    
    wdeg = 1.14
    xedges = np.arange(ramin, ramax, wdeg)
    yedges = np.arange(decmin, decmax, wdeg)
    x = ra
    y = dec
    H, xedges, yedges = np.histogram2d(x, y, bins=(xedges, yedges))

    ct = 0
    for xx in xedges:
        xbincen = xx + 0.5 * wdeg
        for yy in yedges:
            ybincen = yy + 0.5 * wdeg

            c = SkyCoord(ra = xbincen, dec = ybincen, unit = 'deg', frame = 'icrs')
            bb = c.galactic.b.value
            if bb >= -10. and bb <= 10.:
                ct = ct + 1



    #print(H)
    H = H.T 





    histH = np.reshape(H, np.size(H))
    
    
    frac_zerostars = ((np.size(H[H == 0])-ct) / (np.size(H)-ct)) * 100.

    print("Npix(0stars)/Npix(all)=%.1f"%(frac_zerostars) + r"%")

    frac_zerostars = ((np.size(H[H < 5])-ct) / (np.size(H)-ct)) * 100.

    print("Npix(<5stars)/Npix(all)=%.1f"%(frac_zerostars) + r"%")



    #frac_zerostars = (np.size(H[H == 0]) / np.size(H)) * 100.
    
    #print("Npix(0stars)/Npix(all)=%.1f"%(frac_zerostars) + r"%")
    
    #frac_zerostars = (np.size(H[H < 5]) / np.size(H)) * 100.
    
    #print("Npix(<5stars)/Npix(all)=%.1f"%(frac_zerostars) + r"%")
    
    
  

    fig = plt.figure(figsize=(15,8))
    ax = fig.add_subplot(111)

    #title=r" $18<G<20$, $0.2<(BP-RP)<0.7$ (F-star candidates), distance > 500 pc $"

    if zscalelog == True:
        lH = np.log10(H)
        H = lH
        clab = r"$\log N$"
        vmax = np.min([np.max(H), np.log10(300.)])
    else:
        clab = r"N stars per 1.3 deg$^{2}$"
        
        vmax = np.min([np.max(H), 300.])
        print(vmax)
    im = ax.imshow(H, interpolation='nearest', cmap='jet',origin='lower', \
            extent=[xedges[0], xedges[-1], yedges[0], yedges[-1]], \
            vmin=0.0,vmax=vmax)


    # Plot coordinate grid 
    
    for b in np.arange(-90., 90., 20.):
        larr = np.arange(0., 360., 1.)
        barr = np.zeros_like(larr)
        barr.fill(b)
     
        gal = SkyCoord(larr[:], barr[:], frame='galactic', unit=u.deg) 
        eq = gal.transform_to('fk5')
      
     
        if b == 0:
            ls = "-"
        else:
            ls = "--"
        filt = eq.ra.value < 360.
        ax.plot(eq.ra[filt], eq.dec[filt], linestyle = "", marker = ".", ms = 2, color = "white", alpha = 0.6)
        
        # label 
        l = 180.
        gallab = SkyCoord(l, b, frame = 'galactic', unit = u.deg)
        eqlab = gallab.fk5
        ax.text(eqlab.ra.value, eqlab.dec.value, "b=%.0f"%(b), color = "white")
    

    #ax.set_title(title)
    ax.set_xlabel("RA (deg)")
    ax.set_ylabel("DEC (deg)")
    ax.set_xlim(ramin, ramax)
    ax.set_ylim(-30., 89)
    fig.colorbar(im,orientation='horizontal',label= clab)

    plt.savefig("../figs/PS1_classified_skydensity.png")
   
    
    
    #fig, bx = plt.subplots(1, 1)
    
    #xmin = 0
    #xmax = 500
    #bx.hist(histH, bins = 50, range = (xmin, xmax), label = lab)
    #ymax = np.max(histH)
    #ymin = 0.
    #bx.text(xmax - (xmax - xmin) * 0.05, ymax - (ymax - ymin)*0.05, \
    #        "Npix(<5stars)/Npix(all)=%.1f"%(frac_zerostars), ha = 'right')
    
    #bx.set_xlabel(r"N stars per 1.3 deg$^{2}$")
    #bx.set_ylabel(r"N RA-DEC pixels")
    #bx.legend()
    #plt.savefig(catalogname + "_1Dhist.png")
 
   
    return()


def plot_fstar_spec(catalog, synspecpath = "../../pfs_calibstars_data/database/Synspec"):

	df = pd.read_csv(catalog)

	#ource_extended_id,source_id,solution_id,ra,dec,barycentric_distance,pmra,pmdec,radial_velocity,mag_g,mag_bp,mag_rp,mag_rvs,v_i,mean_absolute_v,ag,av,teff,spectral_type,logg,feh,alphafe,mbol,age,mass,radius,vsini,population,has_photocenter_motion,nc,nt,semimajor_axis,eccentricity,inclination,longitude_ascending_node,orbit_period,periastron_date,periastron_argument,variability_type,variability_amplitude,variability_period,variability_phase,r_env_r_star
	teff_refs = np.array([ 5600, 6000., 6400., 6800., 7200., 7600. ], dtype=float)
	dt = 200
	founds = np.array([False, False, False, False, False, False], dtype=bool)

	plt.rcParams["font.size"] = 18
	fig, ax = plt.subplots(1, 1, figsize=(12, 8))


	for source_id, spectype, teff, logg, feh in zip(df['source_id'], df['spectral_type'], df['teff'], df['logg'], df['feh']):
		for i, teff_ref in enumerate(teff_refs):
			if founds[i]==True:
				continue
			if teff<(teff_ref + dt) and teff>(teff_ref - dt) and logg>3.5 and feh<-0.5:
				specfile = synspecpath + "/GaiaEDR3Sim%i"%(source_id) + "/Synspec_moog/GaiaEDR3Sim%i"%(source_id) + ".txt"
				if os.path.isfile(specfile):
					print(teff_ref, source_id, teff, logg, feh, spectype, specfile)
					wave, flux = np.loadtxt(specfile, usecols=(0, 1), unpack = True, skiprows=1)
					ax.plot(wave, flux, marker = "", linestyle = "-", color = cm.magma(i/len(teff_refs)), lw = 2, label=spectype[2:-1], alpha = 0.7)
					founds[i]=True
		if(len(founds[founds==False])==0):
			break

	
	ax.set_xlim(350., 1200.)
	ax.set_xlabel("Wavelength [nm]")
	ax.set_ylabel(r"Flux [erg/cm$^2$/s/$\AA$]")
	plt.legend()
	plt.savefig("../figs/Fstarspec.png")
	return


def plot_trasmission_curve():

	plt.rcParams["font.size"] = 18
	fig, ax = plt.subplots(1, 1, figsize = (12, 5))

	filelist = glob.glob("../PhotmetricFilters/*.dat")

	for filename in filelist:
		filtername = ((((filename.split("/"))[-1]).split("_"))[-1])[:-4]
		w, c = np.loadtxt(filename, usecols=(0, 1), unpack = True)
		ax.plot(w/10., c, lw = 2, label = filtername, alpha = 0.7)

	ax.set_xlim(350., 1200.)
	ax.set_xlabel("Wavelength [nm]")
	ax.set_ylabel("Transmission")
	plt.legend()
	plt.savefig("../figs/filters.png")


	return



def plot_gaia_simulation():

	path = "../Gaia/simulation/"
	filelist = ["GEDR3_l90_b20.csv", "GEDR3_l90_b40.csv","GEDR3_l90_b60.csv","GEDR3_l90_b80.csv"]
	labs = ["b=20", "b=40", "b=60", "b=80"]



	plt.rcParams["font.size"] = 18
	fig=plt.figure(figsize = (12, 12))

	for i, filename in enumerate(filelist):
		df0 = pd.read_csv(path + filename)
		filt = (df0["teff"] < 6000) & (df0["teff"]<7500) & (df0["logg"] > 3.5) & (df0["logg"]<7.)
		df = df0[filt]
		gmag = df["mag_g"]
		ax = fig.add_subplot(2, 2, (i+1))
		
		ax.hist(gmag, bins = np.arange(15., 21., 0.5)-0.25, label = labs[i])
		
		ax.legend()
		ax.set_xlabel("Gmag")
		ax.set_ylabel(r"Nstar/1deg$^2$") 
	plt.tight_layout()
	plt.savefig("../figs/Gaia_simulation_Gmaghist.png")
	return


def plot_stellarlocus():

	
	isopath = "../isochrones/PARSEC/"

	fehs = [-2.0, -1.0, 0.0 ]
	isofiles = [ "PS1_MH-2.0_-1.0_logAge6.6_10.20.dat" , "PS1_MH-1.0_-0.0_logAge6.6_10.20.dat", "PS1_MH0.0_1.0_logAge6.6_10.20.dat" ]
	lages = [10.10, 10.00, 9.70]
	lss = [":", "--", "-"]

	plt.rcParams["font.size"] = 18
	cmap = plt.get_cmap("winter")

	fig, ax = plt.subplots(3, 3, sharex = "col", sharey = "row", figsize=(15,15))

	ax[0, 1].axis("off")
	ax[0, 2].axis("off")
	ax[1, 2].axis("off")


	# plot PS1 data

	catalog = "../catalogs/PS1_highgallat/b90_grizy.csv"
	dfobs0 = pd.read_csv(catalog)

	dfobs = dfobs0.query('gMeanPSFMag!=-999. and rMeanPSFMag!=-999. and iMeanPSFMag!=-999. and zMeanPSFMag!=-999. and yMeanPSFMag!=-999.')


	gr = dfobs["gMeanPSFMag"] - dfobs["rMeanPSFMag"]
	ri = dfobs["rMeanPSFMag"] - dfobs["iMeanPSFMag"]
	iz = dfobs["iMeanPSFMag"] - dfobs["zMeanPSFMag"]
	zy = dfobs["zMeanPSFMag"] - dfobs["yMeanPSFMag"]
	col = "gray"
	ax[0, 0].plot(gr, ri, linestyle = "", marker = ".", ms = 1, color = col)
	ax[1, 0].plot(gr, iz, linestyle = "", marker = ".", ms = 1, color = col)
	ax[2, 0].plot(gr, zy, linestyle = "", marker = ".", ms = 1, color = col)
	ax[2, 1].plot(ri, zy, linestyle = "", marker = ".", ms = 1, color = col)
	ax[2, 2].plot(iz, zy, linestyle = "", marker = ".", ms = 1, color = col)
	ax[1, 1].plot(ri, iz, linestyle = "", marker = ".", ms = 1, color = col)



	fig2, bx = plt.subplots(1, 1)

	for i, feh in enumerate(fehs):
    
		df0 = pd.read_table(isopath + isofiles[i], comment="#", delimiter="\s+", usecols = ["MH", "logAge", "logTe", "logg", "label", "gP1mag", "rP1mag", "iP1mag", "zP1mag", "yP1mag"])
		filt = (df0["MH"]==np.round(feh, 1)) & ((df0["logAge"]).round(2) ==np.round(lages[i], 2)) & (df0["label"]<4)
		df = df0[filt]

	

		gr = df["gP1mag"]-df["rP1mag"]
		ri = df["rP1mag"]-df["iP1mag"]
		iz = df["iP1mag"]-df["zP1mag"]
		zy = df["zP1mag"]-df["yP1mag"]

		ffilt = (df["logTe"]>np.log10(6000.)) & (df["logTe"]<np.log10(7800.)) & (df["logg"]<5.5) & (df["logg"]>3.5)


		bx.plot(gr, df["gP1mag"], ls = lss[i], lw = 3, color = cmap(1), marker = "", label = "[Fe/H]=%.1f, log(Age)=%.2f"%(feh, lages[i]))
		bx.plot(gr[ffilt], df["gP1mag"][ffilt], ls = lss[i], lw = 5, color = cmap(250), marker = "")

		
		
		ax[0, 0].plot(gr, ri, linestyle=lss[i], marker="", linewidth=3, color = cmap(1), label="[Fe/H]=%.1f"%(feh))
		ax[0, 0].plot(gr[ffilt], ri[ffilt], linestyle=lss[i], marker="", linewidth=5, color = cmap(250))


		ax[1, 0].plot(gr, iz, linestyle=lss[i], marker="", linewidth=3, color = cmap(1), label="[Fe/H]=%.1f"%(feh))
		ax[1, 0].plot(gr[ffilt], iz[ffilt], linestyle=lss[i], marker="", linewidth=5, color = cmap(250))

		ax[2, 0].plot(gr, zy, linestyle=lss[i], marker="", linewidth=3, color = cmap(1), label="[Fe/H]=%.1f"%(feh))
		ax[2, 0].plot(gr[ffilt], zy[ffilt], linestyle=lss[i], marker="", linewidth=5, color = cmap(250))

		ax[2, 1].plot(ri, zy, linestyle=lss[i], marker="", linewidth=3, color = cmap(1), label="[Fe/H]=%.1f"%(feh))
		ax[2, 1].plot(ri[ffilt], zy[ffilt], linestyle=lss[i], marker="", linewidth=5, color = cmap(250))

		ax[2, 2].plot(iz, zy, linestyle=lss[i], marker="", linewidth=3, color = cmap(1), label="[Fe/H]=%.1f"%(feh))
		ax[2, 2].plot(iz[ffilt], zy[ffilt], linestyle=lss[i], marker="", linewidth=5, color = cmap(250))
	
		ax[1, 1].plot(ri, iz, linestyle=lss[i], marker="", linewidth=3, color = cmap(1), label="[Fe/H]=%.1f"%(feh))

		ax[1, 1].plot(ri[ffilt], iz[ffilt], linestyle=lss[i], marker="", linewidth=5, color = cmap(250))
	
		if i==2:
			ax[0, 0].legend()
	
	ax[2, 0].set_xlabel(r"$g-r$")
	ax[2, 0].set_xlim(-0.1, 1.5)


	ax[2, 1].set_xlabel(r"$r-i$")
	ax[2, 1].set_xlim(-0.1, 1.5)

	ax[2, 2].set_xlabel(r"$i-z$")
	ax[2, 2].set_xlim(-0.1, 0.9)

	ax[2, 0].set_ylabel(r"$z-y$")
	ax[2, 0].set_ylim(-0.15, 0.40)

	ax[1, 0].set_ylabel(r"$i-z$")
	ax[1, 0].set_ylim(-0.1, 0.9)

	ax[0, 0].set_ylabel(r"$r-i$")
	ax[0, 0].set_ylim(-0.1, 1.5)

	fig.subplots_adjust( 
                    wspace=0.05, 
                    hspace=0.05)

	fig.savefig("../figs/PS1_color-color.eps")

	bx.legend(prop = {'size': 10})
	bx.set_xlim(-0.1, 1.5)
	bx.set_ylim(17.5, -5)
	bx.set_xlabel(r"$g-r$")
	bx.set_ylabel(r"$M_{g}$")
	fig2.tight_layout()
	fig2.savefig("../figs/PS1_gr_g_iso.eps")
	return


def plot_hist_fprob_flags():

	catalog = "../Classified/PS1DR1_GaiaEDR3_ra0.0_1.0_dec-40.0_Classified.csv"

	df0 = pd.read_csv(catalog)

	fig, ax = plt.subplots(1, 1)
	
	bins = np.arange(14., 21., 0.5)

	filt = df0["ProbFstar"]>0.5
	df = df0[filt]
	ax.hist(df["phot_g_mean_mag"], bins = bins , label = "All", alpha=0.7)

	filt = (df["flags_dist"]==False) & (df["flags_ebv"]==0.0) 

	df_noflag = df[filt]
	ax.hist(df_noflag["phot_g_mean_mag"][filt], bins = bins, label = "No flags", alpha=0.7)

	frac_noflag = len(df_noflag)/len(df)

	ax.text(14, 500., "Unflaged Fstar fraction = %.3f"%(frac_noflag))
	ax.set_xlabel("Gmag")
	ax.set_ylabel("Nstars")

	plt.legend()
	plt.savefig("../figs/hist_fprob_flags.png")


if __name__ == "__main__":

	catalogs = glob.glob("../Classified/*.csv")
	plot_sky_density_healpix(catalogs)
	#catalog="../Gaia/simulation/GEDR3_l90_b60.csv"
	#plot_fstar_spec(catalog)
	#plot_trasmission_curve()
	#plot_segue()
	#plot_gaia_simulation()
	#plot_stellarlocus()
	#plot_hist_fprob_flags()


#plot_sky_density_healpix(catlist)


#healpixcat = "../Gaia/Simulation/gaia_source_simulation_G20_teff6000-7800_healpix.csv"
#plot_healpix_Gaia(healpixcat)



#catalogs = glob.glob("../PS1_ExtinctionCorrected/*.csv")
#plot_sky_density_ps1(catalogs, zscalelog = False)

#catalogs = glob.glob("../PS1_Classified/*.csv")
#plot_sky_density_ps1_fstars(catalogs, zscalelog = False)

#glim = [18., 20.]
#gilim = [0.0, 0.4]





#plot_gum(glim, gilim)
#plot_gaia_ps1_cmatch(glim, gilim, extinctioncorr = True)



