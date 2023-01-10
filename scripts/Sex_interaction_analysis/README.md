# Analysis of sex-interacting eQTLs in Monocytes

Focusing on Monocytes from the [OneK1K](https://onek1k.org/) data, testing for interactions between genotypes (at common loci, freq>5%) and chromosomal sex on single-cell gene expression.

Using [Limix_QTL](https://github.com/single-cell-genetics/limix_qtl) a wrapper around [LIMIX](https://github.com/limix/glimix-core), specifically the [interaction runner](https://github.com/single-cell-genetics/limix_qtl/blob/master/Limix_QTL/run_interaction_QTL_analysis.py).

## Preprocessing

* Make chunks from gene annotation file [script](create_chunks.R) - note: this works as a script, not from a notebook.
* Take other input files from CellRegMap run
  * genotypes (plink format, one per chromosome, MAF>5%)
  * kinship (GRM, wide format)
  * covariates (PCs 1:20, sex)
  * sample mapping file (linking cells to individuals)
  * feature variant filter (only SNP-gene pairs found to be eQTL in monocytes [1])
* make phenotype file running [this script](make_phenotype.py) - note: one per chromosome

## Limix QTL wrapper installation

Follow instructions from [here](https://github.com/single-cell-genetics/limix_qtl/wiki/Installation#installation-using-conda) to install using conda

## Snakemake to run interaction eQTL mapping (single cell version)

* [This script](snakemake_sex_interaction_singlecells.smk) is the actual snakemake,
* [this one](snakemake_runner.sh) runs it using PBS (Portable Batch System, qsub), the high performance computing (HPC) cluster at Garvan.


## References

[1] Yazar et al Science 2022
