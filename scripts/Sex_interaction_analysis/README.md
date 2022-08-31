## Snakemake to run interaction eQTL mapping 

### Analysis of sex-eQTLs in Monocytes

Focusing on Monocytes from the [OneK1K] data, testing for interactions between genotypes (at common loci) and sex on single-cell expression

Using [Limix_QTL](https://github.com/single-cell-genetics/limix_qtl) a wrapper around [LIMIX](https://github.com/limix/glimix-core), specifically the [interaction runner](https://github.com/single-cell-genetics/limix_qtl/blob/master/Limix_QTL/run_interaction_QTL_analysis.py).

[This script](snakemake_sex_interaction.smk) is the actual snakemake, [this one](snakemake_runner.sh) runs it using PBS (Portable Batch System, qsub), the high performance computing (HPC) cluster at Garvan.