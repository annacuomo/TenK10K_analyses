import os
import pandas as pd
import scanpy as sc

# phenotype filename (gene expression, all genes for a given chromosome)
chrom = 1
phenotypeFile = f'/share/ScratchGeneral/anncuo/OneK1K/Sex_interactions/Monocytes/input_files/phenotypes_chr{chrom}.tsv'

outdir = '/share/ScratchGeneral/anncuo/OneK1K/Sex_interactions/Monocytes/input_files/'

if not os.path.exists(outdir):
    os.makedirs(outdir, exist_ok=True)

# expression file (H5 format)
mydir = '/share/ScratchGeneral/anncuo/OneK1K/expression_objects/'
phenotype_input_file_h5 = mydir + f'sce{chrom}.h5ad'

# open and extract using scanpy
adata = sc.read(phenotype_input_file_h5)
# sparse to dense
mat = adata.raw.X.todense()
# make pandas dataframe
mat_df = pd.DataFrame(
    data=mat.T, index=adata.raw.var.index, columns=adata.obs.index
)

# save
mat_df.to_csv(phenotypeFile)
