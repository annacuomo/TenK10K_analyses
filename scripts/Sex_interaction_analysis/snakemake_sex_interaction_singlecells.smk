"""
Snakefile for monocytes (OneK1K v1) across subtypes (CD14, CD16) - single-cell version
at first, eQTLs only (identified in any Monocyte sub cell type, OneK1K original paper)
15 expression PCs as covariates (fixed effects)
(expanded) kinship to account for multiple cells from the same individual

Author: Anna Cuomo
Affiliation: Garvan Institute (formerly EMBL-EBI and Wellcome Sanger Institute)
Date: Tuesday 10th January 2023
# Run as: bash snakemake_runner.sh
"""

import glob
import os
from subprocess import run
import pandas as pd
import re
from os.path import join

shell.prefix("set -euo pipefail;")

def _multi_arg_start(flag, files):
    flag += " "
    return " ".join(flag + f for f in files)

def _multi_arg_end(flag, files):
    flag = " "+flag
    return " ".join(f + flag for f in files)

def _multi_arg_both_ends(flag1, flag2, files):
    flag1 += " "
    flag2 = " "+flag2
    return " ".join(flag1 + f + flag2 for f in files)

def flattenChunk(chunk):
    return chunk.replace(":", "_").replace("-", "_")

def extendChunk(chunk):
    relChunk = chunk.pop()
    chunkSplitted = relChunk.split("_")
    return chunkSplitted[0]+":"+chunkSplitted[1]+"-"+chunkSplitted[2]


# Variables
chunkFile = '/share/ScratchGeneral/anncuo/OneK1K/Sex_interactions/chunks.txt'
genotypeFile = '/share/ScratchGeneral/anncuo/OneK1K/plink_files/plink_chr1'
annotationFile = '/share/ScratchGeneral/anncuo/OneK1K/Sex_interactions/LCL.featureCounts.features.tsv'
phenotypeFile = '/share/ScratchGeneral/anncuo/OneK1K/Sex_interactions/Monocytes/input_files/phenotypes_chr1.tsv'
kinshipFiles = '/share/ScratchGeneral/anncuo/OneK1K/input_files_CellRegMap/grm_wide.tsv'
featureVariantFile = '/share/ScratchGeneral/anncuo/OneK1K/input_files_CellRegMap/fvf_Monocyte_eqtls.tsv'
sampleMappingFile = '/share/ScratchGeneral/anncuo/OneK1K/input_files_CellRegMap/smf_monocytes.tsv'
covariateFile = '/share/ScratchGeneral/anncuo/OneK1K/input_files_CellRegMap/PCs_sex_monocytes.tsv'
interactionTerm = 'sex'
numberOfPermutations = '1000'
minorAlleleFrequency = '0.05'
hwe = '0.000001'
callRate = '1'
windowSize = '250000'
blockSize = '15000'
outputFolder = '/share/ScratchGeneral/anncuo/OneK1K/Sex_interactions/Monocytes/results_singlecells/'

print(outputFolder)

finalQtlRun = '/share/ScratchGeneral/anncuo/OneK1K/Sex_interactions/Monocytes/results_singlecells/top_qtl_results_all.txt'
finalQtlRun1 = '/share/ScratchGeneral/anncuo/OneK1K/Sex_interactions/Monocytes/results_singlecells/qtl_results_all.txt'

print(finalQtlRun)
print(finalQtlRun1)

with open(chunkFile,'r') as f:
    chunks = [x.strip() for x in f.readlines()]

print(chunks[1:5])
    
qtlOutput = []
for chunk in chunks:
    processedChunk = flattenChunk(chunk)
    processedChunk=expand(outputFolder+'{chunk}.finished',chunk=processedChunk )
    qtlOutput.append(processedChunk)

print(qtlOutput[1:5])
    
## flatten these lists
qtlOutput = [filename for elem in qtlOutput for filename in elem]

rule all:
    input:
        qtlOutput, finalQtlRun, finalQtlRun1

rule run_qtl_mapping:
    input:
        af = annotationFile,
        pf = phenotypeFile,
        cf = covariateFile,
        kf = kinshipFiles,
        smf = sampleMappingFile,
        fvf = featureVariantFile
    output:
        '/share/ScratchGeneral/anncuo/OneK1K/Sex_interactions/Monocytes/results_singlecells/{chunk}.finished'
    params:
        gen = genotypeFile,
        od = outputFolder,
        it = interactionTerm,
        np = numberOfPermutations,
        maf = minorAlleleFrequency,
        hwe = hwe,
        cr = callRate,
        w = windowSize,
        bs = blockSize
    run:
        chunkFull = extendChunk({wildcards.chunk})
        shell(
            "/share/ScratchGeneral/anncuo/jupyter/conda_notebooks/envs/limix_qtl/bin/python /share/ScratchGeneral/anncuo/github_repos/limix_qtl/Limix_QTL/run_interaction_QTL_analysis.py "
            "--plink {params.gen} "
            " -af {input.af} "
            " -pf {input.pf} "
            " -cf {input.cf} "
            " -od {params.od} "
            " -smf {input.smf} "
            " -fvf {input.fvf} "
            " -rf {input.kf} "
            " -gr {chunkFull} "
            " -i {params.it} "
            " -np {params.np} "
            " -maf {params.maf} "
            " -hwe {params.hwe} "
            " -cr {params.cr} "
            " -c -gm standardize "
            " -w {params.w} "
            " --block_size {params.bs} ")
        shell("touch {output}")


rule aggregate_qtl_results:
    input:
        IF = outputFolder,
        OF = outputFolder,
        finalFiles = qtlOutput
    output:
        '/share/ScratchGeneral/anncuo/OneK1K/Sex_interactions/Monocytes/results_singlecells/top_qtl_results_all.txt'
    run:
        shell(
            "/share/ScratchGeneral/anncuo/jupyter/conda_notebooks/envs/limix_qtl/bin/python /share/ScratchGeneral/anncuo/github_repos/limix_qtl/Limix_QTL/post-processing/minimal_postprocess.py "
            "-id {input.IF} "
            " -od {input.OF} "
            " -sfo -tfb ")

rule aggregate_qtl_results_all:
    input:
        IF = outputFolder,
        OF = outputFolder,
        finalFiles = qtlOutput
    output:
        "/share/ScratchGeneral/anncuo/OneK1K/Sex_interactions/Monocytes/results_singlecells/qtl_results_all.txt"
    run:
        shell(
            "/share/ScratchGeneral/anncuo/jupyter/conda_notebooks/envs/limix_qtl/bin/python /share/ScratchGeneral/anncuo/github_repos/limix_qtl/Limix_QTL/post-processing/minimal_postprocess.py "
            "-id {input.IF} "
            " -od {input.OF} "
            " -sfo -mrp 1 ")
