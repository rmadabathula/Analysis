
library(BiocManager)
BiocManager::install("GenomicFeatures")
library(GenomicFeatures)
install.packages("pheatmap")

calcTPM <- function(counts, lengths) {
  lengths_kb <- lengths / 1000
  rate <- counts / lengths_kb
  tpm <- rate / sum(rate) * 1e6
    return(tpm)                       ### function to calculate Tpm ###
}


count_file <- read.csv(file = "Countswatermelon.csv",
                         row.names = 1, header = TRUE
                         )
gtf_file <- "97103_v2.5.gff3.gz"

txdb <-makeTxDbFromGFF(gtf_file, format =  "gff3")
txdb


exons_by_gene <- exonsBy(txdb, by = "gene")
gene_lengths <- vapply(exons_by_gene, function(gr) sum(width(reduce(gr))), numeric(1))    ### summing up the exons in a gene & determining the length of each gene ####


common_genes <- intersect(rownames(count_file), names(gene_lengths))

if(length(common_genes) == 0) {
  stop("No Overlapping gene IDs between the count table and the GTF-derived gene lengths.")
}

counts_data <- count_file[common_genes, ]
gene_lengths <- gene_lengths[common_genes]

attr(count_file, "gene_lengths") <- gene_lengths

counts_data[] <- lapply(counts_data, function(x) as.numeric(as.character(x)))

tpm_matrix <- apply(counts_data, 2, calcTPM, lengths = gene_lengths)

tpm_matrix <- as.data.frame(tpm_matrix)

write.table(tpm_matrix, file = "TPM_values.txt", sep = "\t", quote = FALSE, col.names = NA)

View(tpm_matrix)


genes_of_interest <- c("AT1G53430", "AT1G53440","AT5G01950")
genes_of_interest

heatmap_data <- tpm_matrix[rownames(tpm_matrix) %in% genes_of_interest, ]

if(nrow(heatmap_data) == 0) {
  stop("None of the specified gene of interst were found in the TPM matrix")
}

log_tpm <- log2(heatmap_data + 1)

metadata <- read.csv("Sample_Metadata.csv", row.names = 1)
library(pheatmap)
pheatmap(log_tpm,
         cluster_rows = TRUE, cluster_cols = FALSE, annotation_col = metadata,
         main = "Heatmap of TPM values for Selected Arabidopsis genes")


                        





# ============================
# PCA for TPM (UG/SG/LG × Control/WF) with ellipses + labels
# Input: TPM_filtered.csv
# Rows = genes, Columns = samples (UG_C1... LGWF_T3)
# ============================

# 0) Packages
# install.packages(c("ggplot2","ggrepel"))  # run once if needed
library(ggplot2)
library(ggrepel)

# 1) Load TPM (genes as rows, samples as columns)
tpm <- read.csv("TPM_filtered.csv", row.names = 1, check.names = FALSE)

# 2) (Optional but recommended) filter low-expression genes
# Keep genes with average TPM > 1 across all samples
tpm <- tpm[rowMeans(tpm) > 1, , drop = FALSE]

# 3) Log-transform
tpm_log <- log2(tpm + 1)

# 4) Transpose for PCA (rows = samples, cols = genes)
tpm_t <- t(tpm_log)

# 5) PCA
pca <- prcomp(tpm_t, center = TRUE, scale. = TRUE)

# 6) Variance explained (for axis labels)
pvar <- (pca$sdev^2) / sum(pca$sdev^2)
pc1_lab <- paste0("PC1 (", round(pvar[1] * 100, 2), "%)")
pc2_lab <- paste0("PC2 (", round(pvar[2] * 100, 2), "%)")

# 7) Build PCA dataframe (PC1/PC2)
pca_df <- as.data.frame(pca$x[, 1:2])
pca_df$Sample <- rownames(pca_df)

# 8) Parse group info from your sample names
# Your columns look like:
# UG_C1, UG_C2, UG_C3
# SG_C1, SG_C2, SG_C3
# LG_C1, LG_C2, LG_C3
# UGWF_T1, UGWF_T2, UGWF_T3
# SGWF_T1, SGWF_T2, SGWF_T3
# LGWF_T1, LGWF_T2, LGWF_T3

# Genotype: UG / SG / LG
pca_df$Genotype <- ifelse(grepl("^UG", pca_df$Sample, ignore.case = TRUE), "UG",
                          ifelse(grepl("^SG", pca_df$Sample, ignore.case = TRUE), "SG",
                                 ifelse(grepl("^LG", pca_df$Sample, ignore.case = TRUE), "LG", NA)))

# Condition: Control (C) vs WF
pca_df$Condition <- ifelse(grepl("WF", pca_df$Sample, ignore.case = TRUE), "WF", "Control")

# Combined group
pca_df$Group <- paste0(pca_df$Genotype, "-", pca_df$Condition)

# Set legend order (optional)
pca_df$Group <- factor(
  pca_df$Group,
  levels = c("UG-Control","UG-WF","SG-Control","SG-WF","LG-Control","LG-WF")
)

# 9) Plot (paper-style: colors + ellipses + repel labels)
p <- ggplot(pca_df, aes(x = PC1, y = PC2, color = Group)) +
  geom_point(size = 3) +
  stat_ellipse(aes(fill = Group), geom = "polygon", alpha = 0.20, color = NA) +
  geom_text_repel(aes(label = Sample), size = 3, max.overlaps = 50) +
  labs(title = "2D PCA Plot", x = pc1_lab, y = pc2_lab) +
  theme_classic() +
  theme(
    plot.title = element_text(hjust = 0.5),
    legend.title = element_blank()
  )

print(p)

# 10) (Optional) Save figure
ggsave("PCA_TPM_UG_SG_LG_Control_WF.png", plot = p, width = 8, height = 6, dpi = 300)

# 11) (Optional) See PCA summary in console
summary(pca)

