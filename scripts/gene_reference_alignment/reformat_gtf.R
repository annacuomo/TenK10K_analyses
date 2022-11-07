# Load R packages
library(plyr)

# Get object (Garvan HPC cluster)
gene_ref_folder <- "/share/ScratchGeneral/anncuo/OneK1K/gene_ref/"
# gtf_file <- paste0(gene_ref_folder, "gencode.v42.basic.annotation.gtf")
gtf_file <- paste0(gene_ref_folder, "gencode.v42.chr_patch_hapl_scaff.annotation.gtf")

# Load object (Gencode v42)
gtf <- read.table(gtf_file, header = FALSE, sep = "\t")

# Turn into data frame
gtf_df <- as.data.frame(gtf)

# Select genes only (no transcripts, exons)
gtf_genes_df <- gtf_df[gtf_df$V3 == "gene", ]
print(paste0("Total number of genes: ", nrow(gtf_genes_df)))

# Rename columns
colnames(gtf_genes_df) <- c("seqname", "source", "feature",
    "start", "end", "score", "strand", "frame", "attribute")

# Extract attributes
attributes_df <- data.frame()
for (i in seq_len(nrow(gtf_genes_df))){
    line <- gtf_genes_df$attribute[i]
    elems <- unlist(strsplit(line, "; "))
    mat <- matrix(unlist(strsplit(elems, " ")), nrow = 2, ncol = length(elems))
    df_curr <- as.data.frame(mat)
    colnames(df_curr) <- mat[1, ]
    df_curr <- df_curr[-1, ]
    attributes_df <- rbind.fill(attributes_df, df_curr)
}

# Add reformatted attributes to data frame
full_df <- cbind(gtf_genes_df, attributes_df)
full_df$attribute <- c()

# Add columns Ensembl Gene ID (no version)
full_df$ensembl_gene_id <- gsub("\\..*", "", full_df$gene_id)

# Clean attribute columns
full_df$tag <- gsub(";", "", full_df$tag)
full_df$level <- gsub(";", "", full_df$level)
full_df$havana_gene <- gsub(";", "", full_df$havana_gene)

# Save formatted object
csv_filename <- gsub(".gtf", ".genesonly.csv", gtf_file)
write.csv(full_df, csv_filename)




