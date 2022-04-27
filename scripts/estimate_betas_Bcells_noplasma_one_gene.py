import os
import sys
import scanpy as sc
import pandas as pd
import xarray as xr
from numpy import ones
from pandas_plink import read_plink1_bin
from numpy.linalg import cholesky
import time
from limix.qc import quantile_gaussianize

from cellregmap import estimate_betas


arg = {}

# chrom
arg["chrom"] = str(sys.argv[1])

# gene index
arg["i"] = int(sys.argv[2])

mydir = "/share/ScratchGeneral/anncuo/OneK1K/"
input_files_dir = mydir+"input_files_CellRegMap/"


######################################
###### sample mapping file (SMF) #####
######################################

## this file will map cells to donors
## in this case, it is limited to B cells only (no plasma)
sample_mapping_file = input_files_dir+"smf_Bcells_noplasma.csv"
sample_mapping = pd.read_csv(sample_mapping_file, dtype={"individual_long": str, "genotype_individual_id": str, "phenotype_sample_id": str}, index_col=0)
