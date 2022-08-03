library(ggplot2)
library(dplyr)

pseudotime_file = "/share/ScratchGeneral/anncuo/OneK1K/slingshot_pseudotime.RDS"
pt = readRDS(pseudotime_file)

summary_betaGxC_file = "/share/ScratchGeneral/anncuo/OneK1K/CRM_interaction/Bcells_noplasma_Bcell_eQTLs/betas/summary_betaGxC.csv"
df1 = read.csv(summary_betaGxC_file, row.names=1)
df1$barcode = rownames(df1)

for (gene in colnames(df1)){
    # check if file already exists
    fig_dir = "/share/ScratchGeneral/anncuo/OneK1K/CRM_interaction/Bcells_Bcell_eQTLs/Figures/pseudotime_vs_betaGxC/"
    filename = paste0(fig_dir,gene,".pdf")
    if (file.exists(filename)){next}
    df2 = inner_join(pt, as.data.frame(df1[,c("barcode",gene)]), by="barcode")
    colnames(df2)[ncol(df2)] = "betaGxC"
    ## save plot
    pdf(filename, width=8, height=6)
    myplot <- ggplot(df2, aes(x = pseudotime, y = betaGxC)) + geom_point(alpha=0.2) + theme_classic() +
        ggtitle(gene)+ stat_smooth(se = F, linetype=2, col="darkgrey")
    print(myplot)
    dev.off()
}