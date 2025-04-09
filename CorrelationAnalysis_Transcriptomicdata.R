library(tidyverse)
library(DESeq2)

##1.Load Count Matrix & Sample Metadata
count_data <- read.csv("Countswatermelon.csv", row.names = 1)
head(count_data)
sample_info <- read.csv("Sampleinfo.csv", row.names = 1)

###2.Reorder count data to match sample info
count_data <- count_data[, rownames(sample_info)]

###3.Create DESeq2 Object & Apply VST Normalization
dds <- DESeqDataSetFromMatrix(countData = count_data,
                              colData = sample_info,
                              design = ~Condition)
vsd <- vst(dds, blind = TRUE)
vsd_mat <- assay(vsd)

###4.Calculate Pearson Correlation Matrix
cor_matrix <- cor(vsd_mat, method = "pearson")

install.packages("reshape2")
library(reshape2) 

###5.Prepare Lower Triangle Correlation Matrix
cor_melt <- melt(cor_matrix)

###6.Remove diagonal and upper triangle
cor_melt <- cor_melt[as.character(cor_melt$Var1) != as.character(cor_melt$Var2), ]
cor_melt <- cor_melt[as.numeric(factor(cor_melt$Var1, levels = rownames(cor_matrix))) >
                       as.numeric(factor(cor_melt$Var2, levels = rownames(cor_matrix))), ]

install.packages("RColorBrewer")
library(RColorBrewer)  

###7.Plot Lower Triangle Heatmap
heatmap_plot <- ggplot(cor_melt, aes(x = Var2, y = Var1, fill = value)) +
  geom_tile(color = "white") +
  scale_fill_gradientn(colors = rev(brewer.pal(9, "RdYlGn")),
                       limits = c(0, 1),
                       name = "Pearson\nCorrelation") +
  scale_y_discrete(limits = rev(levels(factor(cor_melt$Var1)))) +  # Flip Y-axis
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 9),
        axis.text.y = element_text(size = 9)) +
  labs(x = NULL, y = NULL)

print(heatmap_plot)

ggsave("transcriptomic_correlation_heatmap_fixed_highres.png", 
       plot = heatmap_plot, 
       width = 12, 
       height = 10, 
       dpi = 600)
