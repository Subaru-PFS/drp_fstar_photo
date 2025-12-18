                                                                                                    
# The catalog of flux calibration stars for PFS observations
- Virsion: 3.4
- The version used in the latest run (September.2025 run): 3.4.
- Authors: PFS obsproc team
- For more details, please see the [documentation](doc/catalog.md).
                                                                                                    
## Sky coverage
- 0<RA<360 deg, DEC>-30 deg and |b|>10 deg based on the PanStarrs1 DR2 catalog.
- Other fields are supplemented by stars from Gaia DR3, based on `teff_gspphot`.
 --------------------------------------




## File description
                                                                                                    
| name					 | unit		 | description															  |
|-----------------------|-----------|------------------------------------------------------------------------|
| `obj_id`              |  -        | Gaia source ID                                                         |
| `catalog`             |  -        | Catalog name                                                           |
| `ra`                  | deg       | Right ascension in ICRS at the reference epoch                         |
| `dec`                 | deg       | Declination in ICRS at the reference epoch                             |
| `epoch`               | year      | The reference epoch for ra, dec, parallax, pmra, and pmdec             |
| `parallax`            | mas       | Parallax                                                               |
| `parallax_error`      | mas       | Standard error of parallax                                             |
| `pmra`                | mas/yr    | Proper motion in right ascension direction                             |
| `pmra_error`          | mas/yr    | Standard error of pmra                                                 |
| `pmdec`               | mas/yr    | Proper motion in declination direction                                 |
| `pmdec_error`         | mas/yr    | Standard error of pmdec                                                |
| `teff_gspphot`        | K         | Effective temperature from GSP-Phot Aeneas best library using BP/RP spectra |
| `teff_gspphot_lower`  | K         | Lower confidence level (16%) of effective temperature from GSP-Phot Aeneas best library using BP/RP spectra |
| `teff_gspphot_upper`  | K         | Upper confidence level (84%) of effective temperature from GSP-Phot Aeneas best library using BP/RP spectra |
| `fstar_gaia`          |  -        | True if teff_gspphot is between 6000K and 7500K                        |
| `psf_flux_g`          | nJy       | g-band flux                                                            |
| `psf_flux_r`          | nJy       | r-band flux                                                            |
| `psf_flux_i`          | nJy       | i-band flux                                                            |
| `psf_flux_z`          | nJy       | z-band flux                                                            |
| `psf_flux_y`          | nJy       | y-band flux                                                            |
| `psf_flux_error_g`    | nJy       | g-band flux error                                                      |
| `psf_flux_error_r`    | nJy       | r-band flux error                                                      |
| `psf_flux_error_i`    | nJy       | i-band flux error                                                      |
| `psf_flux_error_z`    | nJy       | z-band flux error                                                      |
| `psf_flux_error_y`    | nJy       | y-band flux error                                                      |
| `psf_mag_g`           | mag       | g-band magnitude                                                       |
| `psf_mag_r`           | mag       | r-band magnitude                                                       |
| `psf_mag_i`           | mag       | i-band magnitude                                                       |
| `psf_mag_z`           | mag       | z-band magnitude                                                       |
| `psf_mag_y`           | mag       | y-band magnitude                                                       |
| `psf_mag_error_g`     | mag       | g-band magnitude error                                                 |
| `psf_mag_error_r`     | mag       | r-band magnitude error                                                 |
| `psf_mag_error_i`     | mag       | i-band magnitude error                                                 |
| `psf_mag_error_z`     | mag       | z-band magnitude error                                                 |
| `psf_mag_error_y`     | mag       | y-band magnitude error                                                 |
| `filter_g`            |  -        | g-band filter                                                          |
| `filter_r`            |  -        | r-band filter                                                          |
| `filter_i`            |  -        | i-band filter                                                          |
| `filter_z`            |  -        | z-band filter                                                          |
| `filter_y`            |  -        | y-band filter                                                          |
| `prob_f_star`         |  -        | Probability of being an F-type star                                    |
| `teff_brutus`         | K         | Effective temperature                                                  |
| `teff_brutus_low`     | K         | 16% quantile of the teff posterior distribution                        |
| `teff_brutus_high`    | K         | 84% quantile of the teff posterior distribution                        |
| `logg_brutus`         | dex       | Surface gravity                                                        |
| `logg_brutus_low`     | dex       | 16% quantile of the logg posterior distribution                        |
| `logg_brutus_high`    | dex       | 84% quantile of the logg posterior distribution                        |

## Usage

 Please select flux standards based on  the following criteria. Adgust the `prob_f_star` threshold to 
achieve the desired distribution of flux standards within the target field.

```sql
(prob_f_star > 0.5) OR (fstar_gaia = True)
```


#### Notes
                                                                                                    
  For ra, dec, epoch, parallax, pmra, pmdec, teff_gspphoto and their associated uncertainties, see [Gaia archive website](https://gea.esac.esa.int/archive/) and the [Gaia DR3 documentation](https://gea.esac.esa.int/archive/documentation/GDR3/).




