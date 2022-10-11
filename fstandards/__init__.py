import os, sys

code_dir = os.path.abspath(os.getcwd()) + "/../python"
#code_dir = "/Users/ishigakimiho/PFS/OpenUse/pfs_fstandard_preselection/python"
sys.path.append(code_dir)




import fstandards as fs

from fstandards.query import * 
from fstandards.extinction import *
from fstandards.plot import *

