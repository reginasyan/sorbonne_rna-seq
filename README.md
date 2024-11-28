TRANSCRIPTOMIC DSeq

# Install necessary packages if you haven't already
install.packages("cowplot")
install.packages("gplots")
devtools::install_github("nicolash2/ggvenn")
BiocManager::install("DESeq2")  # if DESeq2 is not already installed
install.packages("dplyr")
install.packages("ggplot2")

# Load libraries
library("matrixTests")
library("ggrepel")
library("enrichR")
library("cowplot")
library("ggplot2")
library("MASS")
library("ggvenn")
library("gplots")
library("DESeq2")
library("dplyr")

# Set working direction
setwd("C:/Users/theop/Desktop/Master I2S/Transcriptomic/Partiels/")

#Creation Dataframe
data = read.delim("./02_data/project-SARSHuman.txt" ,sep = "\t")


Defining variables and prepare the data

# Define the timepoints for each condition
sars_timepoints <- c(0, 12, 24, 36, 48, 60, 72, 84, 94) 
h1n1_timepoints <- c(0, 6, 12, 18, 24, 36, 48)           

# List of all the conditions you have
conditions <- c("SARS.Cov", "SARS.dORF6", "SARS.BatSRBD", "H1N1", "mock")

# Sample names
sample_names <- colnames(data)  # This should be your sample names like "H1N1_00h_1"

# Extract condition and timepoint information from sample names
# Create a vector to store the conditions
condition <- sapply(sample_names, function(x) strsplit(x, "_")[[1]][1])  # First part is condition (e.g., H1N1)
# Create a vector to store the timepoints
timepoint <- sapply(sample_names, function(x) strsplit(x, "_")[[1]][2])  # Second part is timepoint (e.g., 00h)

# Remove the 'h' from timepoints and convert to numeric values
timepoint_numeric <- as.numeric(gsub("h", "", timepoint))

# Create the metadata dataframe with replicates grouped by condition and timepoint
metadata <- data.frame(
  condition = factor(condition, levels = conditions),  # Condition as a factor
  timepoint = timepoint_numeric  # Timepoint as numeric
)

# Ensure that metadata has the same number of rows as the columns in the data
rownames(metadata) <- sample_names

Preparation of the DSeq Analysis + DSeq analysis

# Liste des conditions et points de temps disponibles pour chaque condition
conditions_to_compare <- c("SARS.Cov", "SARS.dORF6", "SARS.BatSRBD", "H1N1")
mock_timepoints <- c(0, 12, 24, 36, 48, 60, 72, 84, 94)  # Points de temps pour le groupe contrôle (mock)

# Créer une liste vide pour stocker les résultats
results_list <- list()

# Creating DESeq2 dataset from count data and metadata
dds <- DESeqDataSetFromMatrix(countData = round(data), colData = metadata, design = ~ condition + timepoint)

# Check the DESeq2 object
dds

# Run the DESeq analysis
dds <- DESeq(dds)

#Comparing the data to control####
contrast = c("variable", "level1", "level2")

# Create a vector of conditions to compare with Control
conditions_to_compare <- c("SARS.Cov", "SARS.dORF6", "SARS.BatSRBD", "H1N1")


Defining the limits of the axis for the plots

# Step 1 : calculate the limitations of the axis
all_log2fc <- c()
all_neglog10pval <- c()

for (condition in conditions_to_compare) {
  for (time in unique(metadata$timepoint)) {
    condition_samples <- grep(paste0(condition, "_", time, "h"), sample_names)
    mock_samples <- grep(paste0("mock_", time, "h"), sample_names)
    
    if (length(condition_samples) > 0 && length(mock_samples) > 0) {
      res <- results(dds, contrast = c("condition", condition, "mock"))
      if (!is.null(res)) {
        res_timepoint <- res[metadata$timepoint == time, ]
        if (!is.null(res_timepoint)) {
          volcano_data <- as.data.frame(res_timepoint)
          all_log2fc <- c(all_log2fc, volcano_data$log2FoldChange)
          all_neglog10pval <- c(all_neglog10pval, -log10(volcano_data$pvalue))
        }
      }
    }
  }
}

# Definind the limits of the axis
x_range <- range(all_log2fc, na.rm = TRUE)
y_range <- range(all_neglog10pval, na.rm = TRUE)


Generating the plots with a loop :

# Step 2 : Loop
for (condition in conditions_to_compare) {
  for (time in unique(metadata$timepoint)) {
    condition_samples <- grep(paste0(condition, "_", time, "h"), sample_names)
    mock_samples <- grep(paste0("mock_", time, "h"), sample_names)
    
    if (length(condition_samples) > 0 && length(mock_samples) > 0) {
      res <- results(dds, contrast = c("condition", condition, "mock"))
      if (!is.null(res)) {
        res_timepoint <- res[metadata$timepoint == time, ]
        if (!is.null(res_timepoint)) {
          volcano_data <- as.data.frame(res_timepoint)
          volcano_data$gene <- rownames(volcano_data)
          volcano_data$category <- with(volcano_data, 
                                        ifelse(log2FoldChange < 0 & padj < 0.05, "Down-regulated + significant",
                                               ifelse(log2FoldChange < 0 & padj >= 0.05, "Down-regulated + non-significant",
                                                      ifelse(log2FoldChange > 0 & padj < 0.05, "Up-regulated + significant",
                                                             "Up-regulated + non-significant"))))
          
          # Identifier les gènes à annoter
          top_up_genes <- volcano_data %>%
            filter(log2FoldChange > 0 & padj < 0.05) %>% # Significant or no
            arrange(desc(log2FoldChange), padj) %>% # Sort the genes from the lower to the higher in terms of expression
            head(5)
          
          top_down_genes <- volcano_data %>%
            filter(log2FoldChange < 0 & padj < 0.05) %>%
            arrange(log2FoldChange, padj) %>%
            head(5)
          
          plot_title <- paste(condition, time, "h vs Control", sep = " ")
          p <- ggplot(volcano_data, aes(x = log2FoldChange, y = -log10(pvalue), color = category)) +
            geom_point(alpha = 0.6) +
            scale_color_manual(values = c(
              "Down-regulated + significant" = "blue",
              "Down-regulated + non-significant" = "lightblue",
              "Up-regulated + significant" = "red",
              "Up-regulated + non-significant" = "pink"
            )) +
            xlim(x_range) + ylim(y_range) +  # Appliquer les limites fixes
            theme_bw() +
            labs(
              title = plot_title,
              x = "Log2 Fold Change",
              y = "-Log10 P-value",
              color = "Category"
            ) +
            theme(
              plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
              legend.title = element_text(size = 10),
              legend.text = element_text(size = 8)
            ) +
            # Annotation des gènes up-regulated (en vert foncé)
            geom_text_repel(
              data = top_up_genes,
              aes(label = gene),
              color = "darkgreen",  # Couleur spécifique pour les up-regulated
              size = 3,
              max.overlaps = 10
            ) +
            # Annotation des gènes down-regulated (en violet)
            geom_text_repel(
              data = top_down_genes,
              aes(label = gene),
              color = "darkorange1",  # Couleur spécifique pour les down-regulated
              size = 3,
              max.overlaps = 10
            )
          
          save_dir <- "C:/Users/theop/Desktop/Master I2S/Transcriptomic/Partiels/03_figures"
          condition_dir <- file.path(save_dir, condition)
          if (!dir.exists(condition_dir)) {
            dir.create(condition_dir)
          }
          plot_filename <- file.path(condition_dir, paste0("volcano_", condition, "_", time, "h_vs_Control.pdf"))
          ggsave(plot_filename, plot = p, width = 8, height = 6)
          cat("Volcano plot saved for", condition, "at", time, "h as", plot_filename, "\n")
        }
      }
    }
  }
}



