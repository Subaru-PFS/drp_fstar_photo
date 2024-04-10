                                                                                                    
# The catalog of flux calibration stars for PFS observations
                                                                                                    
- Virsion: 3.3
                                                                                                    
- Table created on 10/4/2024
                                                                                                    
- Authors: PFS obsproc team
                                                                                                    
- For more details, please see the [documentation](doc/catalog.md).
                                                                                                    
 --------------------------------------




## File description
                                                                                                    
 1. obj_id               -        Gaia source ID
 2. catalog              -        Catalog name
 3. ra                  deg       Right ascension in ICRS at the reference epoch
 4. dec                 deg       Declination in ICRS at the reference epoch
 5. epoch               year      The reference epoch for ra, dec, parallax, pmra, and pmdec
 6. parallax            mas       Parallax
 7. parallax_error      mas       Standard error of parallax
 8. pmra                mas/yr    Proper motion in right ascension direction
 9. pmra_error          mas/yr    Standard error of pmra
 10. pmdec              mas/yr    Proper motion in declination direction
 11. pmdec_error        mas/yr    Standard error of pmdec
 12. teff_gspphot       K         Effective temperature from GSP-Phot Aeneas best library using BP/RP spectra
 13. teff_gspphot_lower K         Lower confidence level (16%) of effective temperature from GSP-Phot Aeneas best library using BP/RP spectra
 14. teff_gspphot_upper K         Upper confidence level (84%) of effective temperature from GSP-Phot Aeneas best library using BP/RP spectra
 15. fstar_gaia          -        True if teff_gspphot is between 6000K and 7500K
 16. psf_flux_g         nJy       g-band flux
 17. psf_flux_r         nJy       r-band flux
 18. psf_flux_i         nJy       i-band flux
 19. psf_flux_z         nJy       z-band flux
 20. psf_flux_y         nJy       y-band flux
 21. psf_flux_error_g   nJy       g-band flux error
 22. psf_flux_error_r   nJy       r-band flux error
 23. psf_flux_error_i   nJy       i-band flux error
 24. psf_flux_error_z   nJy       z-band flux error
 25. psf_flux_error_y   nJy       y-band flux error
 26. psf_mag_g          mag       g-band magnitude
 27. psf_mag_r          mag       r-band magnitude
 28. psf_mag_i          mag       i-band magnitude
 29. psf_mag_z          mag       z-band magnitude
 30. psf_mag_y          mag       y-band magnitude
 31. psf_mag_error_g    mag       g-band magnitude error
 32. psf_mag_error_r    mag       r-band magnitude error
 33. psf_mag_error_i    mag       i-band magnitude error
 34. psf_mag_error_z    mag       z-band magnitude error
 35. psf_mag_error_y    mag       y-band magnitude error
 36. filter_g            -        g-band filter
 37. filter_r            -        r-band filter
 38. filter_i            -        i-band filter
 39. filter_z            -        z-band filter
 40. filter_y            -        y-band filter
 41. prob_f_star         -        Probability of being an F-type star
 42. teff_brutus        K         Effective temperature
 43. teff_brutus_low    K         16% quantile of the teff posterior distribution
 44. teff_brutus_high   K         84% quantile of the teff posterior distribution
 45. logg_brutus        dex       Surface gravity
 46. logg_brutus_low    dex       16% quantile of the logg posterior distribution
 47. logg_brutus_high   dex       84% quantile of the logg posterior distribution




 Notes
                                                                                                    
  - For ra, dec, epoch, parallax, pmra, pmdec, teff_gspphoto and their associated uncertainties.
    see the Gaia archive website and the documentation.




