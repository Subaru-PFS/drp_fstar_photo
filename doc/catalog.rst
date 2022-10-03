# The catalog of flux standard stars 


This catalog is created to ensure the quality of flux calibration in the PFS 2D data reduction pipeline. 

The requrement of the catalog can be summarized as below: 


* To ensure the accurate evaluation of continuum level, the selected stars should ideally be F-type main-sequense stars, 
which typically show a smaller number of absorption features than other types.
* The stars should be brighter than g~20 (TBC) so that a sufficient signal-to-noise is reached within a single exposure (15 miniutes) with the low-resolution mode.  
* The stars should be sufficiently numerous and homogeneously distributed 
to evaluate the throughput variation accros a PFS field-of-view. 





## Base catalogs


* Gaia DR3

** The entire catalog has been downloaded from Gaia archive website and installed to a local database.


* PanStarr1 DR2 

** Downloaded through MAST casjobs query for a slice of one degree in light accention (ra = 359-360 in the following example). 


```

select m.objID, objName, raMean as RAJ2000, decMean as DEJ2000, l, b, gMeanPSFMag as gmag, gMeanPSFMagErr as e_gmag, gFlags, rMeanPSFMag as rmag, rMeanPSFMagErr as e_rmag, rFlags, iMeanPSFMag as imag, iMeanPSFMagErr as e_imag, iFlags, zMeanPSFMag as zmag, zMeanPSFMagErr as e_zmag, zFlags, yMeanPSFMag as ymag, yMeanPSFMagErr as e_ymag, yFlags, ng, nr, ni, nz, ny, objInfoFlag, qualityFlag, psc.ps_score into mydb.RA359_360_gmag135_210 from MeanObjectView m join HLSP_PS1_PSC.pointsource_scores psc on psc.objid=m.objID where raMean between 359 and 360 and (b<-10 or b>10) and gMeanPSFMag<21. and gMeanPSFMag>13.5 and ng>0 and nr>0 and ni>0 and nz>0

```

* Gaia x PanStarrs1 cross match

** Cross-matched table (list of Gaia ID and corresponding PanStarrs1 ID)


## Create a combined catalog

Total number of entries:  276780591






## Selection methods

### Logistic regression 

	Input: g, r, i, z, extinction map
    Output: Probability of being an F-type star 


### Stellar parameter estimates based on MCMC 

	Input: g, r, i, z, parallax, extinction map, stellar isochrone models
	Output: Posterior probability distribution of Teff and logg 


### Stellar parameter estimates based on the brutus code 

* Input: g, r, i, z, parallax, stellar isochrone models
* Output: Posterior probability distribution of Teff and logg 


* The algorithm in brutus: 
  
  ** Interporate a grid of bolometric correction for a given set of labels: 
      temperature, logg, feh, afe, av, rv
  ** The weights and bias for the interporation are pre-computed. The interporation 
     is therefore very fast 
  ** 

  



## Properties of the selected stars

### Sky distribution 






### Number density 

### Magnitude distribution 


## Preparation of specific run

### Sep.21-27, 2022

RA: 18-24h, 0-6h

276780591










