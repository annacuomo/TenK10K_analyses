# Analysis of sex-interacting eQTLs in Monocytes

Focusing on Monocytes from the [OneK1K](https://onek1k.org/) data, testing for interactions between genotypes (at common loci, freq>5%) and chromosomal sex on single-cell gene expression.

Using [Limix_QTL](https://github.com/single-cell-genetics/limix_qtl) a wrapper around [LIMIX](https://github.com/limix/glimix-core), specifically the [interaction runner](https://github.com/single-cell-genetics/limix_qtl/blob/master/Limix_QTL/run_interaction_QTL_analysis.py).

## Preprocessing

* Make chunks from gene annotation file [script](create_chunks.R) - note: this works as a script, not from a notebook.
* Create other input files


## Snakemake to run interaction eQTL mapping (single cell version)

* [This script](snakemake_sex_interaction_singlecells.smk) is the actual snakemake,
* [this one](snakemake_runner.sh) runs it using PBS (Portable Batch System, qsub), the high performance computing (HPC) cluster at Garvan.
