library(data.table)
library(dplyr)
library(ggplot2)
library(ggthemes)
library(Seurat)
library(qvalue)

# # example SNP - gene pair
# rsid <- "rs2848626" # 11:57283988
# genename <- "SELL"

## start by loading files that can be opened just once
pseudotime_file = "/share/ScratchGeneral/anncuo/OneK1K/slingshot_pseudotime.RDS"
pt = as.data.frame(readRDS(pseudotime_file))

filename = "/share/ScratchGeneral/anncuo/OneK1K/Seyhan_scripts/hrc_ids_all.rds"
df_hrc <- readRDS(filename)

expr_filename = "/share/ScratchGeneral/anncuo/OneK1K/Bcells_sce.rds"
expr_sce = readRDS(expr_filename)

# get significant SNP-gene pairs instead
all_results = read.csv("/share/ScratchGeneral/anncuo/OneK1K/CRM_interaction/Bcells_Bcell_eQTLs/summary.csv", row.names=1)
all_results$qv = qvalue(all_results$pv_raw)$qvalues
sign_results = all_results[all_results$qv<0.05,]

for (i in 1:nrow(sign_results)){
    print(round(i/nrow(sign_results)*100), digits=2)
    genename = sign_results$gene[i]
    snp = sign_results$snp_id[i]
    rsid = df_hrc[grep(snp,df_hrc$snpid),]$ID
    # check if file already exists
    fig_dir = "/share/ScratchGeneral/anncuo/OneK1K/CRM_interaction/Bcells_Bcell_eQTLs/Figures/pseudo_gene_by_genotype/"
    filename = paste0(fig_dir,genename,"-",rsid,".pdf")
    if (file.exists(filename)){next}
    # Get the snpid
    eSNP_df <- as.data.frame(df_hrc) %>% filter(ID==rsid)
    eSNP <- eSNP_df$snpid
    chrNumber <- gsub(":.*","", eSNP)
    # Prepare genotype file
    genotype_dir <- '/share/ScratchGeneral/anncuo/OneK1K/Seyhan_scripts/Genotype_Files'
    genotype <- fread(sprintf("%s/genotype_chr%s.tsv", genotype_dir, chrNumber))
    genotype_df <- genotype %>% select("sampleid", all_of(eSNP))
    colnames(genotype_df)[2] <- "genotype"
    genotype_df$individual = genotype_df$sampleid
    # Prepare gene expression
    expr_df <- data.frame(gene = FetchData(expr_sce, vars = genename)) 
    expr_df$barcode = rownames(expr_df)
    colnames(expr_df)[1] = "gene"
    ### merge dataframes
    pt_geno = left_join(genotype_df, pt, by="individual")
    pt_geno_expr = left_join(expr_df, pt_geno, by="barcode")
    # get ref and alt
    df_snps <- fread(sprintf("/share/ScratchGeneral/anncuo/OneK1K/snps_with_maf_greaterthan0.05/chr%s.SNPs.txt", chrNumber))
    eSNP_letter <- gsub("_.*", "", eSNP)
    df_snps <- df_snps %>% filter(SNP==eSNP_letter)
    A1 <- df_snps$'REF(0)'
    A2 <- df_snps$'ALT(1)'
    pt_geno_expr$genotype <- factor(pt_geno_expr$genotype, labels=c(paste0(A1,A1), paste0(A1,A2), paste0(A2,A2)))
    ## save plot
    pdf(filename, width=10, height=6)
    myplot <- ggplot(pt_geno_expr, aes(x=pseudotime, y=gene, colour=as.factor(genotype))) + geom_point() + 
        stat_smooth(se=F, linetype = 2, aes(group=as.factor(genotype), colour=as.factor(genotype))) + 
        scale_color_canva(palette = "Art history inspired") + theme_classic() +
        theme(text = element_text(size=20)) + ylab(genename) + labs(colour=rsid)
    print(myplot)
    dev.off()
}
