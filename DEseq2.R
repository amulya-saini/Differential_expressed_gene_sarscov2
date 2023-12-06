# Specifying the path to the miRNA count file
miRNA_path <- "C:\\Users\\saini\\Downloads\\all_samples_miRNA_count.txt"

# reading the file into a dataframe
miRNA_df <- read.table(miRNA_path, header = TRUE, sep = "\t")

# columns to be dropped
columns_to_drop <- c("Chr", "Start", "End", "Strand", "Length")

# Using square brackets to drop the specified columns
miRNA_df <- miRNA_df[, !colnames(miRNA_df) %in% columns_to_drop]

# trimming the column name to only have ascension number as sample name
for (i in 2:ncol(miRNA_df)) {  # Start from the second column, skipping the first one
  colnames(miRNA_df)[i] <- sub(".*\\.(SRR[0-9]+)_sorted\\.bam", "\\1", colnames(miRNA_df)[i])
}

# Specifying the path to the gene count file
Gene_path <- "C:\\Users\\saini\\Downloads\\all_samples_gene_count.txt"

# reading the file into a dataframe
gene_df <- read.table(Gene_path, header = TRUE, sep = "\t")

# columns to be dropped
columns_to_drop <- c("Chr", "Start", "End", "Strand", "Length")

# Using square brackets to drop the specified columns
gene_df <- gene_df[, !colnames(gene_df) %in% columns_to_drop]

# trimming the column name to only have ascension number as sample name
for (i in 2:ncol(gene_df)) {  # Start from the second column, skipping the first one
  colnames(gene_df)[i] <- sub(".*\\.(SRR[0-9]+)_sorted\\.bam", "\\1", colnames(gene_df)[i])
}

# combining the dataframes
deg_df <- rbind(miRNA_df, gene_df)

# converting the first column to row names
deg_df <- deg_df %>% remove_rownames %>% column_to_rownames(var="Geneid")

View(deg_df)

# Store row names
rownames_deg_df <- rownames(deg_df)

# Convert the columns to numeric
deg_df <- sapply(deg_df, as.numeric)

# Restore row names
rownames(deg_df) <- rownames_deg_df

# sum of all the counts in each row
row_sums <- rowSums(deg_df)

# dropping gene for which the sum of counts is 0
deg_df <- deg_df[row_sums > 0, ]

# view
View(deg_df)

# Specifying the path to the miRNA count file
meta_path <- "C:\\Users\\saini\\Downloads\\SraRunTable.txt"

# reading the file into a dataframe
metadata <- read.table(meta_path, header = TRUE, sep = ",")

# converting the first column to row names
metadata <- metadata %>% remove_rownames %>% column_to_rownames(var="Run")
View(metadata)

# loading the required libraries
library( "DESeq2" )
library(ggplot2)
library(tidyverse)

# making sure the row names in deg_df matches to column names in metadata
all(colnames(deg_df) %in% rownames(metadata))

# are they in the same order?
all(colnames(deg_df) == rownames(metadata))

# Creating DESeqDataSetFromMatrix
dds <- DESeqDataSetFromMatrix(countData = deg_df,
                              colData = metadata,
                              design = ~ infection + Time)

# Running DESeq analysis
dds <- DESeq(dds)

# Check the levels of the "Time" factor
levels(dds$Time)

########### DEGs for Sarcov2 at different time points ###########

# Differential expression analysis
results_timesarcov2 <- results(dds, contrast = c("Time", "24h", "72h"), name = "SARS-CoV-2 infected")

# Order the results by adjusted p-value
results_timesarcov2 <- results_timesarcov2[order(results_timesarcov2$padj),]

# Create a data frame with the ordered results
results_timesarcov2 <- data.frame(results_timesarcov2)

# Print significant genes (padj < 0.05) with row names
significant_genes_ts <- results_timesarcov2[which(results_timesarcov2$padj < 0.05), ]

# getting DEGs into a list
significant_genes_ts_list <- paste(as.list(rownames(significant_genes_ts)), collapse = ",")
print(significant_genes_ts_list)

########### DEGs between control vs SARS-CoV-2 ###########

# Differential expression analysis 
results_control_vs_SARS <- results(dds, contrast = c("infection", "untreated (control)", "SARS-CoV-2 infected"))

# Order the results by adjusted p-value
results_control_vs_SARS <- results_control_vs_SARS[order(results_control_vs_SARS$padj),]

# Create a data frame with the ordered results
results_control_vs_SARS <- data.frame(results_control_vs_SARS)

# Print significant genes (padj < 0.05) with row names
significant_genes <- results_control_vs_SARS[which(results_control_vs_SARS$padj < 0.05), ]
View(significant_genes)

# getting DEGs into a list
significant_genes_list <- paste(as.list(rownames(significant_genes)), collapse = ",")
print(significant_genes_list)
