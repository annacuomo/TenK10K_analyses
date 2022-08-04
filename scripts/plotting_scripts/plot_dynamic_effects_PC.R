library(data.table)
library(dplyr)
library(ggplot2)
library(ggthemes)
library(Seurat)
library(qvalue)

# # example SNP - gene pair
# rsid <- "rs2848626" # 11:57283988
# genename <- "SELL"
args = commandArgs(trailingOnly=TRUE)

pc_number <- args[1]
pc = paste0("PC_",pc_number)
print(paste0("Running for PC",pc))

## start by loading files that can be opened just once
# chage this to be PCs (context file?)
PCs_filename = "/share/ScratchGeneral/anncuo/OneK1K/CellRegMap_input_files/all_B_cells/PCs_Bcells.csv"
PCs_df = read.csv(PCs_filename)
colnames(PCs_df)[1] = "barcode"
print(paste0("Number of cells:", ncol(PCs_df)-1))

# match cells to individuals
smf_filename = "/share/ScratchGeneral/anncuo/OneK1K/CellRegMap_input_files/all_B_cells/smf_Bcells.csv"
smf = read.csv(smf_filename, row.names=1)
smf$barcode = smf$phenotype_sample_id
smf$individual = smf$individual_long
PCs_donor = left_join(PCs_df[,1:11], smf[,c("barcode","individual")], by="barcode")

filename = "/share/ScratchGeneral/anncuo/OneK1K/Seyhan_scripts/hrc_ids_all.rds"
df_hrc <- readRDS(filename)

print("Load Seurat object")
expr_filename = "/share/ScratchGeneral/anncuo/OneK1K/Bcells_sce.rds"
expr_sce = readRDS(expr_filename)

print("Load CellRegMap results")
# get significant SNP-gene pairs instead
all_results = read.csv("/share/ScratchGeneral/anncuo/OneK1K/CRM_interaction/Bcells_Bcell_eQTLs/summary.csv", row.names=1)
all_results$qv = qvalue(all_results$pv_raw)$qvalues
sign_results = all_results[all_results$qv<0.05,]

print("Start plotting loop")
for (i in 1:nrow(sign_results)){
    print(c(i,round(i/nrow(sign_results)*100, digits=2)))
    genename = sign_results$gene[i]
    snp = sign_results$snp_id[i]
    rsid = df_hrc[grep(snp,df_hrc$snpid),]$ID[1]
    # check if file already exists
    fig_dir = paste0("/share/ScratchGeneral/anncuo/OneK1K/CRM_interaction/Bcells_Bcell_eQTLs/Figures/PC",pc,"_gene_by_genotype/")
    filename = paste0(fig_dir,genename,"-",rsid,".pdf")
    if (file.exists(filename)){next}
    # Get the snpid
    eSNP_df <- as.data.frame(df_hrc) %>% filter(ID==rsid)
    eSNP <- eSNP_df$snpid[1]
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
    pc_geno = left_join(genotype_df, PCs_donor, by="individual")
    pc_geno_expr = left_join(expr_df, pc_geno, by="barcode")
    # get ref and alt
    df_snps <- fread(sprintf("/share/ScratchGeneral/anncuo/OneK1K/snps_with_maf_greaterthan0.05/chr%s.SNPs.txt", chrNumber))
    eSNP_letter <- gsub("_.*", "", eSNP)
    df_snps <- df_snps %>% filter(SNP==eSNP_letter)
    A1 <- df_snps$'REF(0)'
    A2 <- df_snps$'ALT(1)'
    if (2 %in% unique(pc_geno_expr$genotype) == F){pc_geno_expr$genotype <- factor(pc_geno_expr$genotype, labels=c(paste0(A1,A1), paste0(A1,A2)))}
    else {pc_geno_expr$genotype <- factor(pc_geno_expr$genotype, labels=c(paste0(A1,A1), paste0(A1,A2), paste0(A2,A2)))}
    ## save plot
    pdf(filename, width=10, height=6)
    myplot <- ggplot(pc_geno_expr, aes(x=pc, y=gene, colour=as.factor(genotype))) + geom_point() + 
        stat_smooth(se=F, linetype = 2, aes(group=as.factor(genotype), colour=as.factor(genotype))) + 
        scale_color_canva(palette = "Art history inspired") + theme_classic() +
        theme(text = element_text(size=20)) + ylab(genename) + labs(colour=rsid)
    print(myplot)
    dev.off()
}
