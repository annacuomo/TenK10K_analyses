## Genes to use as reference when mapping reads from scRNA-seq

### The concern

I was looking for expression information for a specific gene (IGLL5, Ensembl ID: ENSG00000254709). 
This genes is pretty highly expressed (in B cells) and indeed was found [to have eQTLs in the original study](https://onek1k.org/dashboard?type=cd4et%2Ccd4nc%2Ccd4sox4%2Ccd8nc%2Ccd8et%2Ccd8s100b%2Cnk%2Cnkr%2Cbmem%2Cbin%2Cplasma%2Cmonoc%2Cmononc%2Cdc&resultSet=esnp&search=igll5).

However, when I try to look it up in the new expression objects obtained using a new alignment (from GRCh37 to GRCh38), the gene appears to be gone.

### The solution

Our VP of comp bio got back to me. See below his feedback and let me know if you need more input. 
"Anything annotated as a readthrough transcript is removed without further consideration, so that's probably the direct reason. In Gencode 42 IGLL5 does not have the read-through annotation, so I'm guessing it was a mistake in build 32. Probably the best thing to do would be to build a new reference using Gencode 42, based on the instructions here: https://support.10xgenomics.com/single-cell-gene-expression/software/release-notes/build#GRCh38_2020A
We'll probably update to the latest Gencode next year.
IGLL5 is a special IGL transcript that doesn't require any VDJ recombination in order to be expressed, but uses the same IGLJ1 and IGLC1 that could be put into a normal IGL VDJ recombined transcript. So for a 3' read you probably won't be able to tell the difference between a recombined IGL transcript using IGLC1 and a IGLL5 transcript -- you'll likely get an ambiguous gene assignment, which wouldn't be counted into the final matrix. You might want to remove either IGLL5 or IGLJ1 and IGLC1 in order to force a count. The 5' region of IGLL5 is unique so for 5' you probably will get correct counting."

### The plan

I downloaded the [Gencode release 42](https://www.gencodegenes.org/human/) - this is reference GRCh38 (specifically GRCh38.p13).

Reference building 10X Genomics guidelines: https://support.10xgenomics.com/single-cell-gene-expression/software/release-notes/build#GRCh38_2020A
