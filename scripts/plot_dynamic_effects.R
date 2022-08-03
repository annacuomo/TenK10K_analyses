library(data.table)
library(dplyr)
library(ggplot2)
library(ggthemes)
library(Seurat)

# example SNP - gene pair
rsid <- "rs2848626" # 11:57283988
genename <- "SELL"

## x axis: pseudotime
pseudotime_file = "/share/ScratchGeneral/anncuo/OneK1K/slingshot_pseudotime.RDS"
pt = as.data.frame(readRDS(pseudotime_file))

## colour (i.e., stratified by): genotype
# Get the snpid
filename = "/share/ScratchGeneral/anncuo/OneK1K/Seyhan_scripts/hrc_ids_all.rds"
df_hrc <- readRDS(filename)
eSNP_df <- as.data.frame(df_hrc) %>% filter(ID==rsid)
eSNP <- eSNP_df$snpid
chrNumber <- gsub(":.*","", eSNP)

# Prepare genotype file
genotype_dir <- '/share/ScratchGeneral/anncuo/OneK1K/Seyhan_scripts/Genotype_Files'
genotype <- fread(sprintf("%s/genotype_chr%s.tsv", genotype_dir, chrNumber))
genotype_df <- genotype %>% select("sampleid", all_of(eSNP))
colnames(genotype_df)[2] <- "genotype"
genotype_df$individual = genotype_df$sampleid

## y axis: single-cell expression
expr_filename = "/share/ScratchGeneral/anncuo/OneK1K/Bcells_sce.rds"
expr_sce = readRDS(expr_filename)
expr_df <- data.frame(gene = FetchData(expr_sce, vars = genename)) 
expr_df$barcode = rownames(expr_df)
colnames(expr_df)[1] = "gene"

### merge dataframes
pt_geno = left_join(genotype_df, pt, by="individual")
pt_geno_expr = left_join(expr_df, pt_geno, by="barcode")

df_snps <- fread(sprintf("/share/ScratchGeneral/anncuo/OneK1K/snps_with_maf_greaterthan0.05/chr%s.SNPs.txt", chrNumber))
eSNP_letter <- gsub("_.*", "", eSNP)
df_snps <- df_snps %>% filter(SNP==eSNP_letter)
A1 <- df_snps$'REF(0)'
A2 <- df_snps$'ALT(1)'
pt_geno_expr$genotype <- factor(pt_geno_expr$genotype, labels=c(paste0(A1,A1), paste0(A1,A2), paste0(A2,A2)))

## save plot
fig_dir = "/share/ScratchGeneral/anncuo/OneK1K/CRM_interaction/Bcells_Bcell_eQTLs/Figures/pseudo_gene_by_genotype/"
pdf(paste0(fig_dir,genename,"-",rsid,".pdf"), width=10, height=6)
myplot <- ggplot(pt_geno_expr, aes(x=pseudotime, y=gene, colour=as.factor(genotype))) + geom_point() + 
     stat_smooth(se=F, linetype = 2, aes(group=as.factor(genotype), colour=as.factor(genotype))) + 
     scale_color_canva(palette = "Warm and cool") + theme_classic() +
     theme(text = element_text(size=20)) + ylab(genename) + labs(colour=rsid)
print(myplot)
dev.off()
