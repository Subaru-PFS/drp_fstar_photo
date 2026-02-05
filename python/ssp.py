import numpy as np
import os, sys


code_dir = os.path.abspath(os.getcwd()) + "/../"
sys.path.append(code_dir)

import fstandards as fs

from fstandards.database import *







def select_msto_sspfield():


	start_time = time.process_time()

	#ncpus = 30
	#print("Number of CPUs used: ", ncpus)
	#pool = mp.Pool(ncpus)

	rastep = 10.

	ramins = np.array([330., 340., 350., 0., 10., 20., 30., 40.])

	regions = [(ramin, ramin+rastep, -1., 7) for ramin in ramins ]
	

	for region in regions:
		df = query_ps1_region(region[0], region[1], region[2], region[3])	
		filt = (df['parallax'] < 0.2) & (df['gmag']< 20.) & (df['gmag']-df['rmag']>0.1) & (df['gmag']-df['rmag']<0.5) 

		print(region, "Target density (per square degree) = %f"%(len(df[filt])/(8.*rastep)))

	end_time = time.process_time()

	elapse_time = end_time - start_time
	print("CPU time = " ,elapse_time)

	return()



def query_ps1_region(ramin, ramax, decmin, decmax):

    with get_connection() as conn:
        with conn.cursor() as cur:
            #cur.execute('CREATE TABLE test_match AS SELECT gp.*, ps.* FROM gaia_ps1_bestneighbour gp JOIN ps1dr1 ps ON gp.original_ext_source_id=ps."objID"')
            #cur.execute('SELECT gp.* FROM gaia_ps1_bestneighbour gp')
            cur.execute('SET max_parallel_workers_per_gather=8')
            df = pd.read_sql(sql='SELECT gp."objID", source_id, ref_epoch, ra, dec, l, b, parallax, phot_g_mean_mag, phot_bp_mean_mag, phot_rp_mean_mag, gmag, e_gmag, rmag, e_rmag, imag, e_imag, zmag,e_zmag, ymag, e_ymag FROM gaia_ps1_crossmatch gp WHERE ra>=%f AND ra<%f AND dec>=%f AND dec<%f;'%(ramin, ramax, decmin, decmax), con=conn)

            df.to_csv("test.csv")
            #colnames = [col.name for col in cur.description]
            print("FINISHED.")
    return(df)


if __name__ == "__main__":
	select_msto_sspfield()


