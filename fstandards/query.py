import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib
from pylab import meshgrid 
import astropy.coordinates as coord
import astropy.units as u
from astropy.io import ascii
from astropy.coordinates import SkyCoord, FK5
from astropy.table import Table

import sys
import re
import pylab
import json
import requests

try: # Python 3.x
    from urllib.parse import quote as urlencode
    from urllib.request import urlretrieve
except ImportError:  # Python 2.x
    from urllib import pathname2url as urlencode
    from urllib import urlretrieve

try: # Python 3.x
    import http.client as httplib 
except ImportError:  # Python 2.x
    import httplib   


from astroquery.gaia import Gaia

from astroquery.mast import Observations
from astroquery.mast import Catalogs




#Gaia.login()


def query_gaiaedr3_pointings(coords, frame, radius, output_dir):


    
    Gaia.ROW_LIMIT = -1
    Gaia.MAIN_GAIA_TABLE = "gaiaedr3.gaia_source"
    

    for coord in coords:
        
        c = SkyCoord(coord[0], coord[1], frame = frame, unit = "deg")


        radius = u.Quantity(radius, u.deg)
        j = Gaia.cone_search_async(c, radius)
    

        r = j.get_results()
   
        outpath = output_dir + \
            '/GEDR3_RA%.0f_DEC%.0f_r%.0f.csv'%\
            (c.icrs.ra.degree, c.icrs.dec.degree, (radius.value)*60)
        ascii.write(r, outpath, format = 'csv', overwrite = True)

    return()

def query_gaiaedr3simulation_pointing(coord, frame):

	#Gaia.MAIN_GAIA_TABLE = "gaiaedr3.gaia_universe_model"

	c = SkyCoord(coord[0], coord[1], frame = frame, unit = "deg")

	width = u.Quantity(1., u.deg)
	height = u.Quantity(1., u.deg)
	
	#radius = u.Quantity(radius, u.deg)

	print("Query for RA=%.5f DEC=%.5f\n"%(c.icrs.ra.degree, c.icrs.dec.degree))
	Gaia.ROW_LIMIT = -1
	Gaia.MAIN_GAIA_TABLE = "gaiaedr3.gaia_universe_model"

	r = Gaia.query_object_async(coordinate=c, width = width, height = height)
	#query = 'SELECT * FROM gaiaedr3.gaia_universe_model '\
	#	+ 'WHERE 1=CONTAINS('\
	#	+ 'POINT("ICRS", %.5f, %.5f), '%(c.icrs.ra.degree, c.icrs.dec.degree) \
	#	+ 'CIRCLE("ICRS", ra, dec, %.5f))'%(radius)
  
	#query = 'SELECT TOP 10 * FROM gaiaedr3.gaia_universe_model'
    #j = Gaia.launch_job_async(query = query, output_file = outpath, \
    #                          dump_to_file = True, output_format = 'csv')

	#j = Gaia.launch_job_async(query = query)
	#print(j)
	#r = j.get_results()
	
	return(r)


def get_gaiaedr3sim_galactic():

	# Galactic coordinates:
	coords = [(90., 20.), (90., 40.), (90., 60), (90., 80)]

	frame = "galactic"

	for coord in coords:

		r = query_gaiaedr3simulation_pointing(coord, frame) 


		output_dir = "../Gaia/simulation"
		outpath = output_dir + \
            '/GEDR3_l%.0f_b%.0f.csv'%\
            (coord[0], coord[1])
		ascii.write(r, outpath, format = 'csv', overwrite = True)

	return()





def query_gaiaedr3_area(crange, Glim, output_dir):


    ra_min, ra_max = crange[0]
    dec_min, dec_max = crange[1]


    query = 'SELECT * FROM gaiaedr3.gaia_source as g '\
        + 'WHERE (ra > %.2f AND ra < %.2f AND '%(ra_min, ra_max) \
        + 'dec > %.2f AND dec < %.2f AND '%(dec_min, dec_max) \
        + 'phot_g_mean_mag > %.1f AND phot_g_mean_mag < %.1f)'\
        %(Glim[0],Glim[1])

    
    catalog_name = "gedr3_ra%.2f-%.2f_dec%.2f-%.2f_glim%.0f-%.0f.csv"\
        %(ra_min, ra_max, dec_min, dec_max, Glim[0]*10, Glim[1]*10)
    
    Gaia.ROW_LIMIT = -1
    j = Gaia.launch_job_async(query = query, output_file = output_dir + \
                              "/" + catalog_name, \
                              dump_to_file = True, output_format = 'csv')
    #r = j.get_results()

    #catalog_name = "gedr3_ra%.2f-%.2f_dec%.2f-%.2f_glim%.0f-%.0f.csv"\
    #    %(ra_min, ra_max, dec_min, dec_max, Glim[0]*10, Glim[1]*10)
    
    #outpath = output_dir + "/" + catalog_name

    #ascii.write(r, outpath, format = 'csv', overwrite = True)
    
    return(catalog_name)


def query_ps1_area(crange, glim, output_dir):


    # Check the latest information here
    #       https://catalogs.mast.stsci.edu/docs/panstarrs.html

    ra_min, ra_max = crange[0]
    dec_min, dec_max = crange[1]
    
    #catalog_data = Catalogs.query_criteria(raMean = [("gte", ra_min), ("lte", ra_max)],
    #                                       decMean = [("gte", dec_min), ("lte", dec_max)], \
    #                                       catalog="Panstarrs", \
    #                                       data_release = "dr2", table = "mean", \
    #                                       ng = [("gte", 2)], \
    #                                       nr = [("gte", 2)], \
    #                                       ni = [("gte", 2)], \
    #                                       nz = [("gte", 2)], \
    #                                       ny = [("gte", 2)], \
    #                                         gmeanpsfmag = [("lte", glim[1]), ("gte", glim[0])],
    #                                       #columns = ["objName", "objID"],\
    #                                       pagesize = 500000)
    catalog_data = Catalogs.query_criteria(raMean = [("gte", ra_min), ("lt", ra_max)],
                                           decMean = [("gte", dec_min), ("lt", dec_max)], \
                                           catalog="Panstarrs", \
                                           data_release = "dr1", table = "mean", \
                                           #ng = [("gte", 2)], \
                                           #nr = [("gte", 2)], \
                                           #ni = [("gte", 2)], \
                                           #nz = [("gte", 2)], \
                                           #ny = [("gte", 2)], \
                                           gMeanPSFMag = [("lte", glim[1]), ("gte", glim[0])], \
                                           overwrite = True)
                                           #columns = ["objName", "objID"],\
                                           #pagesize = 500000)

    catalog_name = 'ps1dr1_ra%.2f-%.2f_dec%.2f-%.2f_glim%.0f-%.0f.csv'%\
    (ra_min, ra_max, dec_min, dec_max, glim[0]*10, glim[1]*10)
    
    outpath = output_dir + '/' + catalog_name
    ascii.write(catalog_data, outpath, format = 'csv', overwrite = True)


    return(catalog_name)


def cmatch_gaia_ps1(table_gaia, table_ps1, output_dir):

    #Gaia.login()

    table_name_gaia = (((table_gaia.strip(".csv").split("/"))[-1]).lower()).replace('-','_')
    table_name_gaia = table_name_gaia.replace('.', '_')
    
    job = Gaia.upload_table(upload_resource = table_gaia, 
                            table_name=table_name_gaia, format="csv")

    table_name_ps1 = (((table_ps1.strip(".csv").split("/"))[-1]).lower()).replace('-','_')
    table_name_ps1 = table_name_ps1.replace('.', '_')

    
    job = Gaia.upload_table(upload_resource = table_ps1, 
                            table_name=table_name_ps1, format="csv")


    
    
    query = ("select * from user_mishig01." + table_name_gaia + " as ga "
             "inner join gaiaedr3.panstarrs1_best_neighbour as gp "
             "on ga.source_id = gp.source_id "
             "inner join user_mishig01." + table_name_ps1 + " as ps "
             "on gp.original_ext_source_id = ps.objID ")
 

    Gaia.ROW_LIMIT = -1

    outpath = output_dir + "/" + table_name_gaia.replace("gedr3","gedr3_ps1dr1") + ".csv"

    job_cmatch = Gaia.launch_job_async(query=query, output_file = outpath, dump_to_file = True, output_format = 'csv')
    #results_cmatch = job_cmatch.get_results()
    
    job_del = Gaia.delete_user_table(table_name_gaia)
    job_del = Gaia.delete_user_table(table_name_ps1) 
    return()


def ps1cone(ra,dec,radius,table="mean",release="dr1",format="csv",columns=None,
           baseurl="https://catalogs.mast.stsci.edu/api/v0.1/panstarrs", verbose=False,
           **kw):
    """Do a cone search of the PS1 catalog
    
    Parameters
    ----------
    ra (float): (degrees) J2000 Right Ascension
    dec (float): (degrees) J2000 Declination
    radius (float): (degrees) Search radius (<= 0.5 degrees)
    table (string): mean, stack, or detection
    release (string): dr1 or dr2
    format: csv, votable, json
    columns: list of column names to include (None means use defaults)
    baseurl: base URL for the request
    verbose: print info about request
    **kw: other parameters (e.g., 'nDetections.min':2)
    """
    
    data = kw.copy()
    data['ra'] = ra
    data['dec'] = dec
    data['radius'] = radius
    return ps1search(table=table,release=release,format=format,columns=columns,
                    baseurl=baseurl, verbose=verbose, **data)


def ps1area(ramin, ramax, decmin, decmax, \
            table="mean",release="dr1",format="votable",columns=None,
           baseurl="https://catalogs.mast.stsci.edu/api/v0.1/panstarrs", verbose=False,
           **kw):
    """Do a cone search of the PS1 catalog
    
    Parameters
    ----------
    ra (float): (degrees) J2000 Right Ascension
    dec (float): (degrees) J2000 Declination
    radius (float): (degrees) Search radius (<= 0.5 degrees)
    table (string): mean, stack, or detection
    release (string): dr1 or dr2
    format: csv, votable, json
    columns: list of column names to include (None means use defaults)
    baseurl: base URL for the request
    verbose: print info about request
    **kw: other parameters (e.g., 'nDetections.min':2)
    """
    
    data = kw.copy()
    data['raMean.gt'] = ramin
    data['raMean.lt'] = ramax
    data['decMean.gt'] = decmin
    data['decMean.lt'] = decmax

    data['pagesize.eq'] = 500000
    return ps1search(table=table,release=release,format=format,columns=columns,
                    baseurl=baseurl, verbose=verbose, **data)




def ps1search(table="mean",release="dr1",format="votable",columns=None,
           baseurl="https://catalogs.mast.stsci.edu/api/v0.1/panstarrs", verbose=False,
           **kw):
    """Do a general search of the PS1 catalog (possibly without ra/dec/radius)
    
    Parameters
    ----------
    table (string): mean, stack, or detection
    release (string): dr1 or dr2
    format: csv, votable, json
    columns: list of column names to include (None means use defaults)
    baseurl: base URL for the request
    verbose: print info about request
    **kw: other parameters (e.g., 'nDetections.min':2).  Note this is required!
    """
    
    data = kw.copy()
    if not data:
        raise ValueError("You must specify some parameters for search")
    checklegal(table,release)
    if format not in ("csv","votable","json"):
        raise ValueError("Bad value for format")
    url = "{baseurl}/{release}/{table}.{format}".format(**locals())
    if columns:
        # check that column values are legal
        # create a dictionary to speed this up
        dcols = {}
        for col in ps1metadata(table,release)['name']:
            dcols[col.lower()] = 1
        badcols = []
        for col in columns:
            if col.lower().strip() not in dcols:
                badcols.append(col)
        if badcols:
            raise ValueError('Some columns not found in table: {}'.format(', '.join(badcols)))
        # two different ways to specify a list of column values in the API
        # data['columns'] = columns
        data['columns'] = '[{}]'.format(','.join(columns))

    # either get or post works
    #    r = requests.post(url, data=data)
    r = requests.get(url, params=data)

    
    if verbose:
        print(r.url)
    r.raise_for_status()
    if format == "json":
        return r.json()
    else:
        return r.text


def checklegal(table,release):
    """Checks if this combination of table and release is acceptable
    
    Raises a VelueError exception if there is problem
    """
    
    releaselist = ("dr1", "dr2")
    if release not in ("dr1","dr2"):
        raise ValueError("Bad value for release (must be one of {})".format(', '.join(releaselist)))
    if release=="dr1":
        tablelist = ("mean", "stack")
    else:
        tablelist = ("mean", "stack", "detection")
    if table not in tablelist:
        raise ValueError("Bad value for table (for {} must be one of {})".format(release, ", ".join(tablelist)))


def ps1metadata(table="mean",release="dr1",
           baseurl="https://catalogs.mast.stsci.edu/api/v0.1/panstarrs"):
    """Return metadata for the specified catalog and table
    
    Parameters
    ----------
    table (string): mean, stack, or detection
    release (string): dr1 or dr2
    baseurl: base URL for the request
    
    Returns an astropy table with columns name, type, description
    """
    
    checklegal(table,release)
    url = "{baseurl}/{release}/{table}/metadata".format(**locals())
    r = requests.get(url)
    r.raise_for_status()
    v = r.json()
    # convert to astropy table
    tab = Table(rows=[(x['name'],x['type'],x['description']) for x in v],
               names=('name','type','description'))
    return tab


def uploadPS1_crossmatchGaiaEDR3(PS1catalog):

    Gaia.login()

    #table_name_gaia = table_gaia.strip(".csv")
    #job = Gaia.upload_table(upload_resource = table_gaia,
    #                        table_name=table_name_gaia, format="csv")

    table_name_ps1 = PS1catalog.strip(".csv")
    job = Gaia.upload_table(upload_resource = PS1catalog,
                            table_name=table_name_ps1, format="csv")




    query = ("select * from user_mishigak." + table_name_ps1 + " as ps "
             "inner join gaiaedr3.gaia_source as ga "
             "on 1 = coontains (point('', raMean, decMean), circle('', ra, dec, 1./3600.)) " )
           


    Gaia.ROW_LIMIT = -1
    job_cmatch = Gaia.launch_job_async(query=query, dump_to_file=True, output_format='csv')

    print(job_cmatch)
    results_cmatch = job_cmatch.get_results()



    #outpath = "../PS1_ExtinctionCorrected_Gaia/" + table_name_ps1.replace("miho","PS1") + ".csv"

    #ascii.write(results_cmatch, outpath, format = 'csv')



    return


def query_gaia_ps1_crossmatch():

    #Gaia.login(user='', password='')
    
    #Gaia.login()

    catalog_dir = \
        "../PS1DR1_Gaia"
    output_dir = catalog_dir 


    glim = [0., 28.]
    Glim = [0., 28.]



    ramin = 0.
    ramax = 160.
    decmin = -20.
    decmax = 60.


    imax = 320
    for i in range(0, imax):

        rastep = (ramax - ramin)/imax

        ra1 = ramin + i * rastep
        ra2 = ra1 + rastep
    
        jmax = 80
        for j in range(0, jmax):
            decstep = (decmax - decmin)/jmax
        
            dec1 = decmin + j * decstep
            dec2 = dec1 + decstep
    
            crange = [[ra1, ra2], [dec1, dec2]]

        
    
            catalog_name_ps1 = query_ps1_area(crange, glim, catalog_dir)
        
    
            catalog_name_gaia = query_gaiaedr3_area(crange, Glim, catalog_dir)
        
            table_gaia = catalog_dir + "/" + catalog_name_gaia 


            table_ps1 = catalog_dir + "/" + catalog_name_ps1

            print("Crossmatching " + catalog_name_gaia + " x " \
                  + catalog_name_ps1 )
        
        
            r= cmatch_gaia_ps1(table_gaia, table_ps1, output_dir)

          
    
    return

if __name__ == "__main__":

	get_gaiaedr3sim_galactic()




