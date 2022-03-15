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

from cellregmap import run_association_fast


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
## in this case, it is limited to CD4 positive T cells only 
sample_mapping_file = input_files_dir+"smf_Tcells.csv"
sample_mapping = pd.read_csv(sample_mapping_file, dtype={"individual_long": str, "genotype_individual_id": str, "phenotype_sample_id": str}, index_col=0)

## extract unique individuals
donors0 = sample_mapping["genotype_individual_id"].unique()
donors0.sort()
print("Number of unique donors: {}".format(len(donors0)))

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

genes = phenotype.trait.values

##########################################
###### check if file already exists ######
##########################################

gene_name = genes[arg["i"]]

folder = mydir + "CRM_association/all_Tcells/"
outfilename = f"{folder}{gene_name}.tsv"
print(outfilename)

if os.path.exists(outfilename):
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
plink_file = plink_folder+"plink_chr"+str(chrom)+".bed"
G = read_plink1_bin(plink_file)

def cis_snp_selection(feature_id, annotation_df, G, window_size):
    feature = annotation_df.query("gene_name==\"{}\"".format(feature_id)).squeeze()
    chrom = str(feature['seqid'])
    start = feature['start']
    end = feature['end']
    # make robust to features self-specified back-to-front
    lowest = min([start,end])
    highest = max([start,end])
    # for cis, we sequentially add snps that fall within each region
    G = G.where((G.chrom == str(chrom)) & (G.pos > (lowest-window_size)) & (G.pos < (highest+window_size)), drop=True)
    return G


# gene annotation file linking gene to genomic position
annotation_file = mydir+"GeneLocations.tsv"
anno_df = pd.read_csv(annotation_file, sep="\t", index_col=0)

# window size (cis)
w = 100000

G_sel = cis_snp_selection(gene_name, anno_df, G, w)

# expand out genotypes from cells to donors (and select relevant donors in the same step)
G_expanded = G_sel.sel(sample=sample_mapping["individual_long"].values)
assert all(hK_expanded.sample.values == G_expanded.sample.values)

del G

######################################
############ context file ############
######################################

# cells by PCs (10)
C_file = input_files_dir+"PCs.csv.pkl"
C = pd.read_pickle(C_file)
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

del phenotype

GG = G_expanded.values

del G_sel

print("Running for gene {}".format(gene_name))

pvals = run_association_fast(y, W, C.values[:,0:10], G=GG, hK=hK_expanded)[0]

pv = pd.DataFrame({"chrom":G_expanded.chrom.values,
               "pv":pvals,
               "variant":G_expanded.snp.values})
pv.to_csv(outfilename)