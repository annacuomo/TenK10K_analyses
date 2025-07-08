import os
import sys
import numpy as np
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
## here, Monocytes only
sample_mapping_file = input_files_dir+"smf_monocytes.csv"
sample_mapping = pd.read_csv(sample_mapping_file, dtype={"individual_long": str, "genotype_individual_id": str, "phenotype_sample_id": str}, index_col=0)

## extract unique individuals
donors0 = sample_mapping["genotype_individual_id"].unique()
donors0.sort()
print("Number of unique donors: {}".format(len(donors0)))

# Filter on specific gene-SNP pairs
# eQTL from Monocytes (Mono NC + Mono C) that were significantly sex-biased
mono_eqtl_file = input_files_dir+"sex_beta_fvf_Monocyte.csv"
mono_eqtl = pd.read_csv(mono_eqtl_file, index_col = 0)

genes = mono_eqtl[mono_eqtl['chrom']==int(arg["chrom"])]['feature'].unique()

##########################################
###### check if file already exists ######
##########################################

gene_name = genes[arg["i"]]

folder = mydir + "CRM_interaction/Monocytes_Mono_eQTLs/sex_interactions/betas/"
outfilename = f"{folder}{gene_name}"
print(outfilename)

outfilename_betaGxC = outfilename+"_betaGxC.csv"

if os.path.exists(outfilename_betaGxC):
    print("File already exists, exiting")
    sys.exit()

######################################
############ kinship file ############
######################################

## read in GRM (genotype relationship matrix; kinship matrix)
kinship_file = input_files_dir+"grm_wide.csv"
K = pd.read_csv(kinship_file, index_col=0)
K.index = K.index.astype('str')
assert all(K.columns == K.index) #symmetric matrix, donors x donors

K = xr.DataArray(K.values, dims=["sample_0", "sample_1"], coords={"sample_0": K.columns, "sample_1": K.index})
K = K.sortby("sample_0").sortby("sample_1")
donors = sorted(set(list(K.sample_0.values)).intersection(donors0))
print("Number of donors after kinship intersection: {}".format(len(donors)))

## subset to relevant donors
K = K.sel(sample_0=donors, sample_1=donors)
assert all(K.sample_0 == donors)
assert all(K.sample_1 == donors)

## and decompose such as K = hK @ hK.T (using Cholesky decomposition)
hK = cholesky(K.values)
hK = xr.DataArray(hK, dims=["sample", "col"], coords={"sample": K.sample_0.values})
assert all(hK.sample.values == K.sample_0.values)

del K
print("Sample mapping number of rows BEFORE intersection: {}".format(sample_mapping.shape[0]))
## subsample sample mapping file to donors in the kinship matrix
sample_mapping = sample_mapping[sample_mapping["genotype_individual_id"].isin(donors)]
print("Sample mapping number of rows AFTER intersection: {}".format(sample_mapping.shape[0]))

## use sel from xarray to expand hK (using the sample mapping file)
hK_expanded = hK.sel(sample=sample_mapping["genotype_individual_id"].values)
assert all(hK_expanded.sample.values == sample_mapping["genotype_individual_id"].values)

######################################
############ genotype file ###########
######################################

## read in genotype file (plink format)
plink_folder = mydir + "plink_files/"
plink_file = plink_folder+"plink_chr"+str(arg["chrom"])+".bed"
G = read_plink1_bin(plink_file)

leads = mono_eqtl[mono_eqtl['feature']==gene_name]['snp_id'].unique()
G_sel = G[:,G['snp'].isin(leads)]

# expand out genotypes from cells to donors (and select relevant donors in the same step)
G_expanded = G_sel.sel(sample=sample_mapping["individual_long"].values)
#assert all(hK_expanded.sample.values == G_expanded.sample.values)

del G

######################################
########### phenotype file ###########
######################################

# open anndata 
my_file = mydir + "expression_objects/sce"+str(arg["chrom"])+".h5ad"
adata = sc.read(my_file)
# sparse to dense
mat = adata.raw.X.todense()
# make pandas dataframe
mat_df = pd.DataFrame(data=mat.T, index=adata.raw.var.index, columns=adata.obs.index)
# turn into xr array
phenotype = xr.DataArray(mat_df.values, dims=["trait", "cell"], coords={"trait": mat_df.index.values, "cell": mat_df.columns.values})
phenotype = phenotype.sel(cell=sample_mapping["phenotype_sample_id"].values)

del mat
del mat_df

######################################
############ context file ############
######################################

# cells (Monocytes only) by PCs + sex + age
C_file = input_files_dir+"PCs_sex_age_monocytes.csv"
C = pd.read_csv(C_file, index_col = 0)
C = xr.DataArray(C.values, dims=["cell", "pc"], coords={"cell": C.index.values, "pc": C.columns.values})
C = C.sel(cell=sample_mapping["phenotype_sample_id"].values)
assert all(C.cell.values == sample_mapping["phenotype_sample_id"].values)

# C_gauss = quantile_gaussianize(C)

######################################
########### Prepare model ############
######################################

n_cells = phenotype.shape[1]
W = ones((n_cells, 1)) # just intercept as covariates

# select gene
y = phenotype.sel(trait=gene_name)

y = quantile_gaussianize(y)
y = y.values.reshape(y.shape[0],1)

#del phenotype

GG = G_expanded.values

del G_sel

# get MAF
MAF_dir = mydir + "snps_with_maf_greaterthan0.05/"
myfile = MAF_dir+"chr"+str(arg['chrom'])+".SNPs.txt"
df_maf = pd.read_csv(myfile, sep="\t")

snps = G_expanded["snp"].values
mafs = np.array([])
for snp in snps:
    mafs = np.append(mafs, df_maf[df_maf["SNP"] == snp]["MAF"].values)

print("Running for gene {}".format(gene_name))

betas = estimate_betas(y=y, W=W, E=C.values[:,0:1], E1=C.values[:,0:12], E2=C.values[:,0:12], G=GG, hK=hK_expanded, maf=mafs)

beta_G = betas[0]
beta_GxC = betas[1][0]

beta_G_df = pd.DataFrame({"chrom":G_expanded.chrom.values,
               "betaG":beta_G,
               "variant":G_expanded.snp.values})

beta_G_df.to_csv(outfilename+"_betaG.csv")

cells = phenotype["cell"].values
snps = G_expanded["variant"].values

beta_GxC_df = pd.DataFrame(data = beta_GxC, columns = snps, index = cells) 
beta_GxC_df.head()

beta_GxC_df.to_csv(outfilename_betaGxC)
