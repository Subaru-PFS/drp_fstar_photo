                                                                                                    
# PFS Flux Calibration Star Catalog

**Latest Version:** 3.5  
**Previous Version (Jan 2025 run):** 3.4  
**Authors:** PFS ObsProc Team

For detailed methodology and selection criteria, see the [documentation](doc/catalog.md).
                                                                                                    

## Sky Coverage

- $0 < \mathrm{RA} < 360$ deg, $\mathrm{DEC} > -30$ deg, and $|b| > 10$ deg (from Pan-STARRS1 DR2)
- Other regions are supplemented with stars from Gaia DR3, using `teff_gspphot`.




## Catalog File Description

| Name                | Unit      | Description |
|---------------------|-----------|-------------|
| `obj_id`            | -         | Gaia source ID |
| `catalog`           | -         | Catalog name |
| `ra`                | deg       | Right ascension (ICRS, reference epoch) |
| `dec`               | deg       | Declination (ICRS, reference epoch) |
| `epoch`             | year      | Reference epoch for astrometry |
| `parallax`          | mas       | Parallax |
| `parallax_error`    | mas       | Parallax error |
| `pmra`              | mas/yr    | Proper motion in RA |
| `pmra_error`        | mas/yr    | Proper motion error (RA) |
| `pmdec`             | mas/yr    | Proper motion in Dec |
| `pmdec_error`       | mas/yr    | Proper motion error (Dec) |
| `teff_gspphot`      | K         | Effective temperature (Gaia GSP-Phot, BP/RP spectra) |
| `teff_gspphot_lower`| K         | 16% confidence lower bound (GSP-Phot) |
| `teff_gspphot_upper`| K         | 84% confidence upper bound (GSP-Phot) |
| `fstar_gaia`        | -         | True if $6000 < T_{\rm eff} < 7500$ K (Gaia) |
| `psf_flux_g`        | nJy       | $g$-band flux |
| `psf_flux_r`        | nJy       | $r$-band flux |
| `psf_flux_i`        | nJy       | $i$-band flux |
| `psf_flux_z`        | nJy       | $z$-band flux |
| `psf_flux_y`        | nJy       | $y$-band flux |
| `psf_flux_error_g`  | nJy       | $g$-band flux error |
| `psf_flux_error_r`  | nJy       | $r$-band flux error |
| `psf_flux_error_i`  | nJy       | $i$-band flux error |
| `psf_flux_error_z`  | nJy       | $z$-band flux error |
| `psf_flux_error_y`  | nJy       | $y$-band flux error |
| `psf_mag_g`         | mag       | $g$-band magnitude |
| `psf_mag_r`         | mag       | $r$-band magnitude |
| `psf_mag_i`         | mag       | $i$-band magnitude |
| `psf_mag_z`         | mag       | $z$-band magnitude |
| `psf_mag_y`         | mag       | $y$-band magnitude |
| `psf_mag_error_g`   | mag       | $g$-band magnitude error |
| `psf_mag_error_r`   | mag       | $r$-band magnitude error |
| `psf_mag_error_i`   | mag       | $i$-band magnitude error |
| `psf_mag_error_z`   | mag       | $z$-band magnitude error |
| `psf_mag_error_y`   | mag       | $y$-band magnitude error |
| `filter_g`          | -         | $g$-band filter |
| `filter_r`          | -         | $r$-band filter |
| `filter_i`          | -         | $i$-band filter |
| `filter_z`          | -         | $z$-band filter |
| `filter_y`          | -         | $y$-band filter |
| `prob_f_star`       | -         | Probability of being an F-type star |
| `teff_brutus`       | K         | Effective temperature (brutus SED fit) |
| `teff_brutus_low`   | K         | 16% quantile of $T_{\rm eff}$ posterior |
| `teff_brutus_high`  | K         | 84% quantile of $T_{\rm eff}$ posterior |
| `logg_brutus`       | dex       | Surface gravity (brutus SED fit) |
| `logg_brutus_low`   | dex       | 16% quantile of $\log g$ posterior |
| `logg_brutus_high`  | dex       | 84% quantile of $\log g$ posterior |
| `is_gc_neighbor`    | -         | True if within 10× half-light radius of a globular cluster |
| `is_dense_region`   | -         | True if in a crowded field |


## Usage

The catalogs are stored in a database located at a Subaru Telescope server. 
To select flux standards from the database, use the following criteria. Adjust the `prob_f_star` threshold as needed to achieve the desired spatial distribution in your target field.

**Basic selection:**

```sql
(prob_f_star > 0.5) OR (is_fstar_gaia = True)
```

**Exclude crowded regions and globular cluster neighbors:**

```sql
((prob_f_star > 0.5) OR (is_fstar_gaia = True)) AND (is_gc_neighbor = False) AND (is_dense_region = False)
```

---

### Notes

- For details on astrometric and photometric parameters (e.g., `ra`, `dec`, `epoch`, `parallax`, `pmra`, `pmdec`, `teff_gspphot`), see the [Gaia archive](https://gea.esac.esa.int/archive/) and [Gaia DR3 documentation](https://gea.esac.esa.int/archive/documentation/GDR3/).
- The catalog of Galactic globular clusters are taken from [Globular cluster database](https://people.smp.uq.edu.au/HolgerBaumgardt/globular/). See [Baumgardt & Vasiliev (2021)](https://ui.adsabs.harvard.edu/abs/2021MNRAS.505.5957B/abstract) for details.

---




