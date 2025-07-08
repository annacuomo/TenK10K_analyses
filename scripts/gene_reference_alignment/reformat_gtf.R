# Load R packages
library(plyr)

# Get object (Garvan HPC cluster)
# downloaded from https://www.gencodegenes.org/human/release_{}.html
gene_ref_folder <- "/share/ScratchGeneral/anncuo/OneK1K/gene_ref/"
gtf_file <- paste0(gene_ref_folder, "gencode.v42.basic.annotation.gtf")  # v42
gtf_file <- paste0(gene_ref_folder, "gencode.v38.annotation.gtf")        # v38

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


# Select transcripts only
gtf_t_df <- gtf_df[gtf_df$V3 == "transcript", ]
print(paste0("Total number of transcripts: ", nrow(gtf_t_df)))

# Rename columns
colnames(gtf_t_df) <- c("seqname", "source", "feature",
    "start", "end", "score", "strand", "frame", "attribute")

# Extract attributes
attributes_df <- data.frame()
for (i in seq_len(nrow(gtf_t_df))){
    line <- gtf_t_df$attribute[i]
    elems <- unlist(strsplit(line, "; "))
    mat <- matrix(unlist(strsplit(elems, " ")), nrow = 2, ncol = length(elems))
    df_curr <- as.data.frame(mat)
    colnames(df_curr) <- mat[1, ]
    n <- length(colnames(df_curr)[grep("tag", colnames(df_curr))]) # split tags
    colnames(df_curr)[grep("tag", colnames(df_curr))] <- paste0("tag", 1:n)
    df_curr <- df_curr[-1, ]
    attributes_df <- rbind.fill(attributes_df, df_curr)
}

# Add reformatted attributes to data frame
full_df <- cbind(gtf_t_df, attributes_df)
full_df$attribute <- c()

# Add columns Ensembl Gene ID (no version)
full_df$ensembl_gene_id <- gsub("\\..*", "", full_df$gene_id)

# Clean attribute columns
full_df$tag1 <- gsub(";", "", full_df$tag1)
full_df$tag2 <- gsub(";", "", full_df$tag2)
full_df$level <- gsub(";", "", full_df$level)
full_df$havana_gene <- gsub(";", "", full_df$havana_gene)

# Save formatted object
csv_filename <- gsub(".gtf", ".transcriptsonly.csv", gtf_file)
write.csv(full_df, csv_filename)



