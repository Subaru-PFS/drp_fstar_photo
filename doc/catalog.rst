The catalog of flux standard stars for PFS
==================================================








Catalog description
---------------------------------------------------

The catalog of flux standard stars for PFS contains stars classified into spectral types 
based on photometric data from PanSTARRS1 DR2 and astrometric data from Gaia DR3. 





Flux standard stars for PFS
--------------------------------------------------
The photometric selection is made to satisfy following requirements.
* To ensure the accurate evaluation of continuum level, the candidate standard stars should ideally be F-type main-sequense stars, 
which typically show a smaller number of strong absorption features than other types.
* The stars should be brighter than g~20 (TBC) so that a sufficient signal-to-noise is reached within a single exposure (15 miniutes) 
taken with the low-resolution mode.  
* The stars should be sufficiently numerous and homogeneously distributed on the sky 
to evaluate the throughput variation accros the PFS field-of-view. 



Base catalogs
---------------

1. Gaia DR3

The entire catalog has been downloaded from Gaia archive website and installed to a local database.


2. PanSTARRS1 DR2 

Downloaded through MAST casjobs query for a slice of one degree in light accention (ra = 359-360 in the following example). 


```
select m.objID, objName, raMean as RAJ2000, decMean as DEJ2000, l, b, gMeanPSFMag as gmag, gMeanPSFMagErr as e_gmag, gFlags, rMeanPSFMag as rmag, rMeanPSFMagErr as e_rmag, rFlags, iMeanPSFMag as imag, iMeanPSFMagErr as e_imag, iFlags, zMeanPSFMag as zmag, zMeanPSFMagErr as e_zmag, zFlags, yMeanPSFMag as ymag, yMeanPSFMagErr as e_ymag, yFlags, ng, nr, ni, nz, ny, objInfoFlag, qualityFlag, psc.ps_score into mydb.RA359_360_gmag135_210 from MeanObjectView m join HLSP_PS1_PSC.pointsource_scores psc on psc.objid=m.objID where raMean between 359 and 360 and (b<-10 or b>10) and gMeanPSFMag<21. and gMeanPSFMag>13.5 and ng>0 and nr>0 and ni>0 and nz>0
```


Median 50% completeness of PFS photometry:  g=23.2, r=23.2, i=23.1, z=22.3 and y=21.2 with significant variation accross the sky.
(https://outerspace.stsci.edu/display/PANSTARRS/PS1+Photometric+Depth) 




3. Gaia x PanSTARRS1 cross match

Cross-matched table (list of Gaia ID and corresponding PanSTARRS1 ID) from Gaia archive


4. Create a combined catalog


```
   count
-----------
 201818351
(1 row)
```






Quality cuts
-------------------------------------

Inorder to ensure the quality of the photometric data, we adopt the following quality cuts. 

- g>14 to avoid saturated objects.
- Non-negative values of uncertainties in PS1 magnitudes.  
- ps_score > 0.8 : To ensure an object is a point source according to the classification scheme of Tachibana+18 (https://iopscience.iop.org/article/10.1088/1538-3873/aae3d9). 
 


```
   count   
-----------
 185990390
(1 row)

```
(92% of the original cross-matched catalog.)















Selection methods
-------------------------------------



1. Logistic regression 

	Input: g, r, i, z, extinction map
    Output: Probability of being an F-type star 



2. Stellar parameter estimates based on the brutus code

 
* Input: g, r, i, z, parallax, stellar isochrone models
* Output: Posterior probability distribution of Teff and logg 


* The algorithm in brutus: 
  
  ** Interporate a grid of bolometric correction for a given set of labels: 
      temperature, logg, feh, afe, av, rv
  ** The weights and bias for the interporation are pre-computed. The interporation 
     is therefore very fast 
  ** 

 

3. Stellar parameter estimates based on MCMC

    Input: g, r, i, z, parallax, extinction map, stellar isochrone models
    Output: Posterior probability distribution of Teff and logg





Properties of the selected stars
-----------------------------------------




### Sky distribution 






### Number density 

### Magnitude distribution 


## Preparation of specific run


Reference
------------------------------

https://docs.github.com/ja/get-started/writing-on-github/getting-started-with-writing-and-formatting-on-github/basic-writing-and-formatting-syntax










