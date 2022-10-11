import pandas as pd
import numpy as np
import uuid
import sqlalchemy 
import glob, os
import time
import d6tstack.utils
import d6tstack.combine_csv
import psycopg2


def apply(df0):
    
    df = df0.filter(['objID', 'RAJ2000', 'DEJ2000', 'Ng', 'Nr', 'Ni', 'Nz', 'Ny', \
            'gmag', 'e_gmag', 'gKmag', 'e_gKmag', 'rmag', 'e_rmag', \
            'rKmag', 'e_rKmag', 'imag', 'e_imag', 'iKmag', 'e_iKmag', \
            'zmag', 'e_zmag', 'zKmag', 'e_zKmag', 'ymag', 'e_ymag', \
            'yKmag', 'e_yKmag'], axis = 1)

    df = df.fillna(-999.)  
    return(df)


def process_gaia_source(df0):
 
    df = df0.filter(['source_id', 'ref_epoch', 'ra', 'ra_error', 'dec', 'dec_error', \
     'parallax', 'parallax_error', 'pmra', 'pmra_error', 'pmdec','pmdec_error', 'ruwe', 'duplicated_source', \
      'phot_g_mean_flux', 'phot_g_mean_flux_error', 'phot_g_mean_mag', 'phot_bp_mean_flux', \
       'phot_bp_mean_flux_error', 'phot_bp_mean_mag', 'phot_rp_mean_flux', 'phot_rp_mean_flux_error', \
        'phot_rp_mean_mag', 'bp_rp', 'dr2_radial_velocity', 'dr2_radial_velocity_error', 'dr2_rv_template_teff', \
         'dr2_rv_template_logg', 'dr2_rv_template_fe_h', 'l', 'b'])
  

    return(df)


def register_ps1(ps1dir):


    cfg_uri_psql = 'postgresql+psycopg2://pfs:fstandard@localhost:5432/gaia_ps1'

    sqlengine = sqlalchemy.create_engine(cfg_uri_psql)
   
    cfg_fnames = list(glob.glob(ps1dir + "/PS1DR1_ra*.csv"))
    
    d6tstack.combine_csv.CombinerCSV(cfg_fnames, apply_after_read = apply, add_filename = False).\
            to_psql_combine(cfg_uri_psql, 'ps1dr1', if_exists = 'replace')
    
    print(pd.read_sql_table('ps1dr1',sqlengine).head())



def register_gaia_ps1_bestneighbour(gaiaps1dir):


    cfg_uri_psql = 'postgresql+psycopg2://pfs:fstandard@localhost:5432/gaia_ps1'

    sqlengine = sqlalchemy.create_engine(cfg_uri_psql)

    cfg_fnames = list(glob.glob(gaiaps1dir + "/*.csv"))

    d6tstack.combine_csv.CombinerCSV(cfg_fnames, add_filename = False).\
            to_psql_combine(cfg_uri_psql, 'gaia_ps1_bestneighbour', if_exists = 'replace')

    print(pd.read_sql_table('gaia_ps1_bestneighbour',sqlengine).head())

def register_gaia(gaiadir):

    cfg_uri_psql = 'postgresql+psycopg2://pfs:fstandard@localhost:5432/gaia_ps1'

    sqlengine = sqlalchemy.create_engine(cfg_uri_psql)

    cfg_fnames = list(glob.glob(gaiadir + "/GaiaSource_*.csv"))

 
    d6tstack.combine_csv.CombinerCSV(cfg_fnames, apply_after_read = process_gaia_source, add_filename = False).\
            to_psql_combine(cfg_uri_psql, 'gaia_source', if_exists = 'replace')

    print(pd.read_sql_table('gaia_source',sqlengine).head())


def get_connection():
    
	dsn = "dbname=gaia_ps1 host=133.40.210.124 user=pfs password=fstandard"
    #dsn = os.environ.get('postgresql://pfs:fstandard@localhost:5432/gaia_ps1')
	return psycopg2.connect(dsn)


def query_ps1():
	
	with get_connection() as conn:
		with conn.cursor() as cur: 
			#cur.execute('CREATE TABLE test_match AS SELECT gp.*, ps.* FROM gaia_ps1_bestneighbour gp JOIN ps1dr1 ps ON gp.original_ext_source_id=ps."objID"')
			#cur.execute('SELECT gp.* FROM gaia_ps1_bestneighbour gp')
			#cur.execute('SET max_parallel_workers_per_gather=16')
			cur.execute('CREATE TABLE gaia_ps1_crossmatch AS SELECT gp.*, ps.*, ref_epoch, ra, dec, parallax, parallax_error, pmra, pmra_error, pmdec, pmdec_error, ruwe, phot_g_mean_mag, phot_bp_mean_mag, phot_rp_mean_mag, bp_rp, dr2_radial_velocity, dr2_radial_velocity_error, l, b FROM gaia_ps1_bestneighbour gp JOIN ps1dr1 ps ON gp.original_ext_source_id=ps."objID" JOIN gaia_source ga ON ga.source_id=gp.source_id')
			#colnames = [col.name for col in cur.description]
			print("FINISHED.")

#gaiaps1dir = "../Gaia/gaia_PS1_crossmatch"
#register_gaia_ps1_bestneighbour(gaiaps1dir)


#gaiadir = "../Gaia/gaia_source"
#register_gaia(gaiadir)

#ps1dir = "../PS1_Gaia_crossmatch_stilts"
#register_ps1(ps1dir)

if __name__ == "__main__":
	#ps1dir = "../PS1_Gaia_crossmatch_stilts"
	#register_ps1(ps1dir)

	query_ps1()



#sqlengine = sqlalchemy.create_engine(cfg_uri_psql)

#start_time = time.time()

#df.to_sql('benchmark',sqlengine,if_exists='replace')
#print("--- %s seconds ---" % (time.time() - start_time))

#start_time = time.time()
#d6tstack.utils.pd_to_psql(df, cfg_uri_psql, 'benchmark', if_exists='replace')
#print("--- %s seconds ---" % (time.time() - start_time))



#data = {'col1':np.arange(0, 10, 1), 'col2':np.arange(10, 20, 1)}

#df = pd.DataFrame(data)

#connection_config = {
#    'user': 'miho',
#    'password': '',
#    'host': 'localhost',
#    'port': '5432', 
#    'database': 'gaia_ps1'
#}

#engine = create_engine('postgresql://{user}:{password}@{host}:{port}/{database}'.format(**connection_config))

#df.to_sql('test_new', con=engine, if_exists='replace', index=True)
