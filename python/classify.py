import numpy as np
import seaborn as sns
import time, datetime

from sklearn.model_selection import train_test_split
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import confusion_matrix, accuracy_score, precision_score, recall_score, f1_score

import matplotlib.pyplot as plt
import matplotlib.cm as cm

import pandas as pd

import sys, glob, os

import multiprocessing as mp



def example():

    iris_df = sns.load_dataset('iris') 
    iris_df = iris_df[(iris_df['species']=='versicolor') | \
                      (iris_df['species']=='virginica')] 


    sns.pairplot(iris_df, hue='species')
    plt.show()



    X = iris_df[['petal_length', 'petal_width']]
    Y = iris_df['species'].map({'versicolor': 0, 'virginica': 1})
    X_train, X_test, Y_train, Y_test = \
        train_test_split(X, Y, test_size = 0.2, random_state = 0)


    lr = LogisticRegression()
    lr.fit(X_train, Y_train)


    #print("coefficient = ", lr.coef_)
    #print("intercept = ", lr.intercept_)

    Y_pred = lr.predict(X_test)

    print('accuracy = ', accuracy_score(y_true=Y_test, y_pred=Y_pred))
    print('precision = ', precision_score(y_true=Y_test, y_pred=Y_pred))

    return()

def segue_ps1(extinction = True, plot = True):

    path = "../catalogs/MAST_PS1_crossmatch/"

    table = "sppParams_PS1DR1_crossmatch_GaiaEDR3_ExtinctionCorrected.csv"


    
    df00 = pd.read_csv(path + table)


    if extinction == True:
    
        filt = (df00['rMeanPSFMag']>0.) & \
            (df00['iMeanPSFMag']>0.) & (df00['zMeanPSFMag']>0.) &\
            (df00['teffadopunc']<100.0) & \
            (df00['g0']-df00['r0']>0.0) & \
            (df00['g0']-df00['r0']<1.0) & \
            (np.abs(df00['b'])>10.)
        outfig = "../figs/Segue_ps1_extinction.png"
        
    else:
        filt = (df00['gMeanPSFMag']>13.) & (df00['gMeanPSFMag']<18.) & \
            (df00['rMeanPSFMag']>0.) & \
            (df00['iMeanPSFMag']>0.) & (df00['zMeanPSFMag']>0.) & \
            (df00['teffadopunc']<100.0) & \
            (df00['gMeanPSFMag']-df00['iMeanPSFMag']>0.0) & \
            (df00['gMeanPSFMag']-df00['iMeanPSFMag']<1.5) & \
            (np.abs(df00['b'])>10.)
        
        outfig = "../figs/Segue_ps1.png"
    
    df0 = df00[filt]

    
    filt = (df0['teffadop'] > 6000.) & (df0['teffadop'] < 7800.) & \
        (df0['loggadop']>3.5) & (df0['loggadop']<5.5)  

    df0['fstars'] = filt


    df0.to_csv(path + "sspParams_PS1DR1_crossmatch_GaiaEDR3_ExtinctionCorrected_selected.csv")



    if extinction == True:
        
        df0['gr'] = df0['g0'] - df0['r0']
        df0['gi'] = df0['g0'] - df0['i0']
        df0['ri'] = df0['r0'] - df0['i0']
        df0['rz'] = df0['r0'] - df0['z0']
        df0['iy'] = df0['i0'] - df0['y0']
        df0['iz'] = df0['i0'] - df0['z0']

    else:
        
        df0['gr'] = df0['gMeanPSFMag'] - df0['rMeanPSFMag']
        df0['gi'] = df0['gMeanPSFMag'] - df0['iMeanPSFMag']
        df0['ri'] = df0['rMeanPSFMag'] - df0['iMeanPSFMag']
        df0['rz'] = df0['rMeanPSFMag'] - df0['zMeanPSFMag']
        df0['iy'] = df0['iMeanPSFMag'] - df0['yMeanPSFMag']
        df0['iz'] = df0['iMeanPSFMag'] - df0['zMeanPSFMag']        


    df = df0[['gr','gi','ri','rz','iy','iz', 'fstars']]
    #sns.pairplot(df, hue = 'fstars', \
    #             hue_order=[False, True],\
    #             markers='+',\
    #             plot_kws={'alpha': 0.8}).savefig(outfig)
    #plt.show()


    df_fstartrue = df0[df0['fstars']==True]

    df_fstartrue.to_csv(path + "fstartrue.csv")
    
    print('The number of training + test data:',len(df))

    
    X = df[['gr', 'iz']]
    Y = df['fstars'].map({False: 0, True: 1})
    X_train, X_test, Y_train, Y_test = \
        train_test_split(X, Y, test_size = 0.2, random_state = 0)
    lr = LogisticRegression()
    lr.fit(X_train, Y_train)

  
    w0 = lr.intercept_
    w1, w2 = (lr.coef_)[0]

    

    XX = np.arange(0.0, 1.0, 0.01)
    YY = -1. * (w1 / w2) * XX - w0 / w2


    if plot == True:
        fig, ax = plt.subplots(1, 1)
        filt = df['fstars'] == False
        plt.scatter(df['gr'][filt], df['iz'][filt], marker = "+", label = "non-Fstars")
        filt = df['fstars'] == True

        plt.scatter(df['gr'][filt], df['iz'][filt], marker = "+", label = "F-stars")
        plt.plot(XX, YY, linestyle = "-", color = "gray", linewidth = 2, \
                 marker = "")

        plt.xlabel("g-r")
        plt.ylabel("i-z")
        plt.legend()
        plt.savefig("../figs/gr_iz_boundary.png")

    # Plot contour
    #gi_grid = np.arange(-0.5, 1.0, 0.05)
    #rz_grid = np.arange(-0.5, 1.0, 0.05)
    #XX, YY = np.meshgrid(gi_grid, rz_grid)
    #ZZ = 1./(1. + np.exp(w0 + w1 * XX + w2 * YY))

    #plt.pcolormesh(XX, YY, ZZ, cmap = 'Oranges', shading = 'auto')
    #pp=plt.colorbar(orientation="vertical")
    #cont=plt.contour(XX, YY, ZZ,  colors=['black'])
    #plt.show()


    
    #print("coefficient = ", lr.coef_)
    #print("intercept = ", lr.intercept_)

    Y_pred = lr.predict(X_test)
    probs = lr.predict_proba(X_test)

    # Plot probability

    if plot == True:
        fig, ax = plt.subplots(1, 1)
        Z = probs[:, 1]
        sc = ax.scatter(X_test['gr'], X_test['iz'], vmin = 0.0, vmax = 1.0, c = Z, \
                        cmap = cm.seismic)

        ax.set_xlim(-0.1, 1.0)
        ax.set_ylim(-0.55, 0.9)
        ax.set_xlabel(r"$(g-r)_{0}$")
        ax.set_ylabel(r"$(i-z)_{0}$")
        cbar = plt.colorbar(sc)
        cbar.set_label("Probability of being F-star")
        plt.savefig("../figs/Segue_ps1_F-star_proba.png")
    

    
    print(probs[0, :])

    print('accuracy = ', accuracy_score(y_true=Y_test, y_pred=Y_pred))
    print('precision = ', precision_score(y_true=Y_test, y_pred=Y_pred))

    return(lr)




def Classify_GaiaPS1(lr, gmaglim = [18.,20.]):

    path = "../../../catalogs/"

    catidlist = ["RA86_4_DEC28_9", "RA96_7_DEC33_8","RA102_2_DEC35_8", "RA108_0_DEC37_6", \
                 "RA133_6_DEC41_6", "RA159_9_DEC39_6", "RA182_9_DEC32_2"]

    b = [0.0, 10., 15., 20., 40., 60., 80.]

    problims = [ 0.0, 0.5, 0.7, 0.9, 0.95 ]   # problims = 0.0 corresponds to the simple color-cut.
    
    nfstars = [[0.0 for x in range(len(b))] for i in range(len(problims))]
 
    for j, catid in enumerate(catidlist):

        df0 = pd.read_csv(path + "GEDR3_PS1DR1_" + catid + \
                      "_r0_65_ExtinctionCorrected.csv")

        
        filt = (df0["gmeanpsfmag"]> gmaglim[0]) & \
            (df0["gmeanpsfmag"]<gmaglim[1]) 
        df = df0[filt].\
            dropna(subset = ["g0", "i0", "r0", "z0"])
        print(len(df))
        
        df["gi"] = df["g0"] - df["i0"]
        df["rz"] = df["r0"] - df["z0"]

        for i, problim in enumerate(problims):

            if i == 0:
                filt = (df["gi"] > 0.0) & (df["gi"]<0.4)
                nfstars[i][j] = len(df[filt])

            else:
            
    
                X = df[['gi', 'rz']]

                Y_pred = lr.predict(X)

                probs = lr.predict_proba(X)

                probs_fstars = probs[:,1]

                filt = (probs_fstars > problim)


                nfstars[i][j] = len(Y_pred[filt])


    figname = "../../../figs/Nfstars_b_gmag%.1f-%.1f.png"%(gmaglim[0], gmaglim[1])

    fig, ax = plt.subplots(1, 1)
    cmap = plt.get_cmap("tab10")
    mks = ["x", "o", "^", "s", "*"]
    cls = [cmap(7), cmap(0), cmap(1), cmap(2), cmap(3)]
    lss = [":", "-", "-", "-", "-"]
    
    
    for i in range(len(problims)):
        if i==0:
            lab = r"$0.0<(g-i)_{0}<0.4$"
        else:
            lab = r"$P$(F-stars)$>%.2f$"%(problims[i])

        ax.plot(b, nfstars[i], linestyle = lss[i], marker = mks[i], color = cls[i], label = lab)

    ax.set_xlabel(r"Galactic latitude ($b$ [deg]) @ $l=180^{\circ}$")
    ax.set_ylabel("N stars")
    ax.set_title(r"The number of F-stars per FoV from GaiaEDR3xPS1DR1, $%.1f<g<%.1f$"%(gmaglim[0], gmaglim[1]))
    plt.legend()
    plt.savefig(figname)
    plt.show()
    print(nfstars)

    return


def Classify_PS1(lr):

    #path = "../../../catalogs/"

    path = "../PS1Gaia_ExtinctionCorrected/"

    catlist = glob.glob(path + "*_ra*.csv")


    #catidlist = ["RA86_4_DEC28_9", "RA96_7_DEC33_8","RA102_2_DEC35_8", "RA108_0_DEC37_6", \
    #             "RA133_6_DEC41_6", "RA159_9_DEC39_6", "RA182_9_DEC32_2"]

    #b = [0.0, 10., 15., 20., 40., 60., 80.]

    #problims = [ 0.0, 0.5, 0.7, 0.9, 0.95 ]   # problims = 0.0 corresponds to the simple color-cut.
    
    #nfstars = [[0.0 for x in range(len(b))] for i in range(len(problims))]
    
    start_time = time.process_time()

    ncpus = 1
    print("Number of CPUs used: ", ncpus)
    pool = mp.Pool(ncpus)

    params = [(catalog, lr) for catalog in catlist ]


    results = pool.map(Classify_PS1_sub, [param for param in params])

    pool.close()

    end_time = time.process_time()

    elapse_time = end_time - start_time
    print("CPU time = " ,elapse_time)


    #p = Pool(8)
    #for catalog in catlist:
    #    #p.apply_async(Classify_PS1_sub, args=(catalog, lr,))
    #    Classify_PS1_sub(catalog, lr)    



    #print('Waiting for all subprocess done...')
    #p.close()
    #p.join()
    #print('All subprocess done.')

    return


def Classify_PS1_sub(param):

    catalog, lr = param

    gmaglim = [0., 21.]

    catid = "_".join((((catalog.split("/"))[-1]).split("_"))[:-1])

    outfile = "../Classified/" + catid + "_Classified.csv"
    print(outfile)

    if os.path.exists(outfile):
        print(outfile, " exist!")
        return

    df0 = pd.read_csv(catalog)
    
    df0["gi"] = df0["g0"] - df0["i0"]
    df0["rz"] = df0["r0"] - df0["z0"]
    df0["gr"] = df0["g0"] - df0["r0"]
    df0["iz"] = df0["i0"] - df0["z0"]


    # Select stars brighter than the limitting magnitudes

    filt = (df0["gmag"]> gmaglim[0]) & \
         (df0["gmag"]<gmaglim[1]) & \
         (df0["rmag"]>0.) & \
         (df0["imag"]>0.) & \
         (df0["zmag"]>0.) & \
         (df0["e_gmag"]<0.01) & \
         (df0["e_rmag"]<0.01) & \
         (df0["e_imag"]<0.01) & \
         (df0["e_zmag"]<0.01) 

    df = df0[filt].dropna(subset = ["g0", "r0", "i0", "z0"])
	
    print("Number of stars original / after magnitude cuts = %i / %i "%(len(df0), len(df)))

    df["ProbFstar"] = np.array([0.0] * len(df))


    #filt = (df0["gmag"]> gmaglim[0]) & \
    #     (df0["gmag"]<gmaglim[1]) & \
         #(df0["rmag"]>0) & \
         #(df0["imag"]>0) & \
         #(df0["zmag"]>0) & \
         #(df0["e_gmag"]<0.01) & \
         #(df0["e_rmag"]<0.01) & \
         #(df0["e_imag"]<0.01) & \
         #(df0["e_zmag"]<0.01) & \
    #     (df0["gr"]<1.) & \
    #     (df0["gr"]>0.) & \

    # Adopt color selection
    filt_color = (df["gr"]<1.) & (df["gr"]>0.)

    df_c1 = df[filt_color]
    df_c2 = df[~filt_color]

    X = df_c1[['gr', 'iz']]
    Y_pred = lr.predict(X)
    probs = lr.predict_proba(X)
    probs_fstars = probs[:,1]

    df_c1['ProbFstar'] = probs_fstars

    dfout = pd.concat([df_c1, df_c2])

    dfout.to_csv(outfile)

    return



def highprobFstars(filelist, threshold):


    for catfile in filelist:

        PS1_ids = ()
        Gaia_ids =()

        # Coordinate 
        ras = ()
        decs = ()
        epochs = ()

        # Astrometry 
        plxs = ()
        plx_errors = ()
        pmras = ()
        pmra_errors = ()
        pmdecs = ()
        pmdec_errors = ()
  


        # Magnitudes
        gmag = ()
        rmag = ()
        imag = ()
        zmag = ()
        ymag = ()

        # Flux
        gFluxJy =()
        rFluxJy =()
        iFluxJy =()
        zFluxJy =()
        yFluxJy =()

    



        flags_dists = ()
        flags_ebvs = ()
        probfstars = ()

        
        df0 = pd.read_csv(catfile)


        filt = (df0["ProbFstar"] > threshold)
        objid = (df0["objID"][filt]).astype(str)
        source_id = (df0["source_id"][filt]).astype(str)
        ra = df0["ra"][filt]
        dec = df0["dec"][filt]
        epoch =df0["ref_epoch"][filt] 
        
        flags_dist = df0["flags_dist"][filt]
        flags_ebv = df0["flags_ebv"][filt] 
        probfstar = df0["ProbFstar"][filt]
        g = df0["gmag"][filt]
        r = df0["rmag"][filt]
        i = df0["imag"][filt]       
        z = df0["zmag"][filt]
        y = df0["ymag"][filt]        



        PS1_ids = np.hstack((PS1_ids, objid))
        ras = np.hstack((ras, ra))
        decs = np.hstack((decs, dec))
        epochs = np.hstack((epochs, epoch))
        Gaia_ids = np.hstack((Gaia_ids, source_id))   


        flags_dists = np.hstack((flags_dists, flags_dist))
        flags_ebvs = np.hstack((flags_ebvs, flags_ebv))
        probfstars = np.hstack((probfstars, probfstar))

        gmag = np.hstack((gmag, g))
        rmag = np.hstack((rmag, r))
        imag = np.hstack((imag, i))
        zmag = np.hstack((zmag, z))
        ymag = np.hstack((ymag, y))


        
        f0 = 3631. #[Jy]
        gFluxJy = np.hstack((gFluxJy, f0*10**(-0.4*g)))
        rFluxJy = np.hstack((rFluxJy, f0*10**(-0.4*r)))
        iFluxJy = np.hstack((iFluxJy, f0*10**(-0.4*i)))
        zFluxJy = np.hstack((zFluxJy, f0*10**(-0.4*z)))
        yFluxJy = np.hstack((yFluxJy, f0*10**(-0.4*y)))

        plx = df0["parallax"][filt]
        plx_e = df0["parallax_error"][filt]
        pmra = df0["pmra"][filt]
        pmra_e = df0["pmra_error"][filt]
        pmdec = df0["pmdec"][filt]
        pmdec_e = df0["pmdec_error"][filt]
        
        plxs = np.hstack((plxs, plx))
        plx_errors = np.hstack((plx_errors, plx_e))
        pmras = np.hstack((pmras, pmra))
        pmra_errors = np.hstack((pmra_errors, pmra_e))
        pmdecs = np.hstack((pmdecs, pmdec))
        pmdec_errors = np.hstack((pmdec_errors, pmdec_e))
             
     
        data = {"PS1_objid":PS1_ids, "GaiaEDR3_sourceid":Gaia_ids, "ra":ras, "dec":decs, "epoch":epochs, "parallax":plxs, "parallax_error":plx_errors, "pmra":pmras, "pmra_error":pmra_errors, "pmdec":pmdecs, "pmdec_error":pmdec_errors, "gPS1":gmag, "rPS1":rmag, "iPS1":imag, "zPS1":zmag, "yPS1":ymag, "gFluxJy":gFluxJy, "rFluxJy":rFluxJy, "iFluxJy":iFluxJy, "zFluxJy":zFluxJy, "yFluxJy":yFluxJy, "flags_dist":flags_dists.astype(bool), "flags_ebv":flags_ebvs.astype(bool), "probfstar":probfstars}
        df = pd.DataFrame(data)
        outcatalog = "../TargetList/Fstar_v1.0_probfstar%.1f"%(threshold) + "_" + (((catfile.split("/"))[-1]).strip("PS1DR1_GaiaEDR3_")).strip("_Classified.csv") + ".csv"
        df.to_csv(outcatalog, index = False)

def produce_readme():

	samplefile = "../TargetList/Fstar_v1.0_part.csv"
	df = pd.read_csv(samplefile)
	colnames = df.columns
	
	line = "="*100
	empty = " "*100 
	empty5 = "\n"*5
	empty10 = "\n"*10
	empty3 = "\n"*3

	line2 = " --------------------------------------"
	title = "        Candidate F-type stars"
	virsion = "            Virsion: 1.0"
	dt_now = datetime.datetime.now()
	datestring = "     Table created on %i/%i/%i"%(dt_now.day, dt_now.month, dt_now.year)	
	author = "      Authors: PFS obsproc team"
	section = " File description"


	f=open("../TargetList/Readme", "w")
	for string in empty, line2, empty, title, empty, virsion, empty, datestring, empty, author, empty, line2, empty3, section, empty:
		f.write(string + "\n")

	for i, column in enumerate(colnames):
		if column == "PS1_objid":
			unit = " - "
			description = "PanStarrs1 DR1 objID"
		elif column == "GaiaEDR3_sourceid":
			unit = " - "
			description = "Gaia EDR3 source_id"
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
		elif column == "gPS1":
			unit = "mag"
			description = "PS1 magnitude in g-band"
		elif column == "rPS1": 
			unit = "mag"
			description = "PS1 magnitude in r-band"
		elif column == "iPS1":
			unit = "mag"
			description = "PS1 magnitude in i-band"
		elif column == "zPS1":
			unit = "mag"
			description = "PS1 magnitude in z-band"
		elif column == "yPS1":
			unit = "mag"
			description = "PS1 magnitude in y-band"
		elif column == "gFluxJy":
			unit = "Jy"
			description = "g-band flux"
		elif column == "rFluxJy":
			unit = "Jy"
			description = "r-band flux"
		elif column == "iFluxJy":
			unit = "Jy"
			description = "i-band flux"
		elif column == "zFluxJy":
			unit = "Jy"
			description = "z-band flux"
		elif column == "yFluxJy":
			unit = "Jy"
			description = "y-band flux"		
		elif column == "flags_dist":
			unit = " - "
			description = "Distance uncertanty flag, True if parallax_error/parallax > 0.2"
		elif column == "flags_ebv":
			unit = " - "
			description = "E(B-V) uncertainty flag, True if E(B-V) uncertainty is greater than 20%"
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
	note1 = "  - Recommended F-star candidates satisfy flags_dist=False, flags_ebv=False and probfstar>0.5."
	note2 = "  - For ra, dec, epoch, parallax, pmra, pmdec and their associated uncertainties." 
	note3 = "    see the Gaia documentation at:"
	note4 = "    https://gea.esac.esa.int/archive/documentation/GEDR3/Gaia_archive/chap_datamodel/sec_dm_main_tables/ssec_dm_gaia_source.html#gaia_source-ra"

	for string in empty3, section, empty, note1, note2, note3, note4, empty3:
		f.write(string + "\n")


	f.close()


if __name__=="__main__":


	#lr = segue_ps1(extinction = True, plot = True)
	#Classify_PS1(lr)


	filelist = glob.glob("../Classified/*.csv")
	threshold = 0.5
	highprobFstars(filelist, threshold)

	#produce_readme()

