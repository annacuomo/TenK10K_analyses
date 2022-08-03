library(dplyr)
library(ggplot2)
library(Seurat)

sceb_filename = "/share/ScratchGeneral/anncuo/OneK1K/b_cells_phate_slingshot.RDS"
sceb = readRDS(sceb_filename)

phate_df = as.data.frame(Embeddings(sceb, reduction = "phate"))
df_to_plot = phate_df

summary_betaGxC_file = "/share/ScratchGeneral/anncuo/OneK1K/CRM_interaction/Bcells_noplasma_Bcell_eQTLs/betas/summary_betaGxC.csv"
df1 = read.csv(summary_betaGxC_file, row.names=1)

for (gene in colnames(df1)){
    print(gene)
    if (gene == "barcode"){next}
    # check if file already exists
    fig_dir = "/share/ScratchGeneral/anncuo/OneK1K/CRM_interaction/Bcells_Bcell_eQTLs/Figures/phate_by_betaGxC/"
    filename = paste0(fig_dir,gene,".pdf")
    if (file.exists(filename)){next}
    df2 = inner_join(df_to_plot, as.data.frame(df1[,c("barcode",gene)]), by="barcode")
    colnames(df2)[ncol(df2)] = "beta"
    ## save plot
    pdf(filename, width=8, height=6)
    myplot <- ggplot(df2, aes(x = PHATE_2, y = PHATE_1, colour = beta)) + geom_point(alpha=0.2) + 
        theme_classic() + scale_colour_gradientn(colors = brewer.pal(9,"Spectral")) + ggtitle(gene))
    print(myplot)
    dev.off()
}