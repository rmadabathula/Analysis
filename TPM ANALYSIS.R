
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





