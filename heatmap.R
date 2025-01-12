#libraries
library("matrixTests") 
library("ggrepel") 
library("DESeq2") 
library("dplyr")
library('ComplexHeatmap')
library("circlize")


setwd("F:/Documents/Rstudio/Project-SARSHuman/")

#Creation Dataframe
data = read.delim("./01_data/project-SARSHuman.txt" ,sep = "\t")

idx_h1n1 = grepl("H1N1",colnames(data))
idx_sars_bat = grepl("SARS.BatSRBD",colnames(data))
idx_sars_cov = grepl("SARS.Cov",colnames(data))
idx_sars_dorf = grepl("SARS.dORF6",colnames(data))
idx_mock = grepl("mock",colnames(data))


data_h1n1 = data[, c(which(idx_h1n1), which(idx_mock))]
data_sars_bat = data[, c(which(idx_sars_bat), which(idx_mock))]
data_sars_cov = data[, c(which(idx_sars_cov), which(idx_mock))]
data_sars_dorf = data[, c(which(idx_sars_dorf), which(idx_mock))]

#Gene differential expression analysis
DE = function(data, dataset_name){
  sample_names <- colnames(data)
  timepoint <- sapply(sample_names, function(x) strsplit(x, "_")[[1]][2]) # Second part is timepoint (e.g., 00h)
  
  metadata <- data.frame(timepoint)
  metadata$type<- as.vector(sapply(sample_names,function(x) paste(strsplit(x,"\\_")[[1]][1],collapse="_")))
  metadata$type <- factor(metadata$type)
  metadata$timepoint <- factor(metadata$timepoint)
  metadata$type <- relevel(metadata$type, ref = "mock")
  rownames(metadata) <- sample_names  
  dds <- DESeqDataSetFromMatrix(countData = round(data), colData = metadata, design = ~type)
  dds <- DESeq(dds)
  colnames(dds)
  result_table <- data.frame()
  
  # Loop Through Timepoints
  for (tp in unique(metadata$timepoint[metadata$type == dataset_name])) {
    dds_subset <- dds[, metadata$timepoint == tp]
    dds_subset <- DESeq(dds_subset)
    # Perform the comparison
    res <- results(dds_subset, contrast = c("type", dataset_name, "mock"))
    # Extract significant genes (adjusted p-value < 0.05 and |log2FoldChange| > 0 for heatmap)
    res <- res[!is.na(res$padj) & !is.na(res$log2FoldChange), ]
    sig_genes <- res[res$padj < 0.05 & abs(res$log2FoldChange) > 0, ]
    sig_genes$timepoint <- tp
    sig_genes$dataset <- dataset_name
    result_table <- rbind(result_table, data.frame(
      gene = rownames(sig_genes),
      fold_change = sig_genes$log2FoldChange,
      p_value_adj = sig_genes$padj,
      timepoint = sig_genes$timepoint,
      dataset = sig_genes$dataset
    ))
  }
  return(result_table)
}

genes_h1n1 = DE(data_h1n1,"H1N1" )
genes_sars_bat = DE(data_sars_bat,"SARS.BatSRBD" )
genes_sars_cov = DE(data_sars_cov,"SARS.Cov" )
genes_sars_dorf = DE(data_sars_dorf,"SARS.dORF6" )

all_results = list(genes_h1n1, genes_sars_bat, genes_sars_cov, genes_sars_dorf)


# Combine all tables into one
combined_table <- do.call(rbind, all_results)
combined_table$direction = ifelse(combined_table$fold_change > 0, "Upregulated", "Downregulated")


#Create amtrix data for heatmaps
heatmap_data <- dcast(combined_table, gene ~ timepoint + dataset, value.var = "fold_change")
# Convert to matrix format (row names = genes)
rownames(heatmap_data) <- heatmap_data$gene
heatmap_data <- as.matrix(heatmap_data[, -1])  # Remove the 'gene' column

heatmap_data[is.na(heatmap_data)] <- 0  # Replace NAs with 0 -> neutral fold change for non stastistically enriched

heatmap_data_filtered <- heatmap_data[apply(heatmap_data, 1, sd) > 0.5, ]


#Groups of viruses for heatmap generation
virus_groups <- list(
  H1N1 = grep("H1N1", colnames(heatmap_data), value = TRUE),
  SARS_Cov = grep("SARS.Cov", colnames(heatmap_data), value = TRUE),
  SARS_dORF6 = grep("SARS.dORF6", colnames(heatmap_data), value = TRUE),
  SARS_BatSRBD = grep("SARS.BatSRBD", colnames(heatmap_data), value = TRUE)
)

# Create a consistent row order using hierarchical clustering
row_order <- hclust(dist(heatmap_data))$order

col_fun <- colorRamp2(c(-10, 0, 10), c("darkblue", "#f7f7f7", "#d73027"))

# Create individual heatmaps for each virus group
heatmaps <- lapply(names(virus_groups), function(virus) {
  
  #Subset specific virus data
  group_data <- heatmap_data[, virus_groups[[virus]], drop = FALSE] # Filter columns

  range <- range(group_data, na.rm = TRUE)
  
  # Update column names to show only the timepoints
  colnames(group_data) <- sapply(colnames(group_data), function(x) strsplit(x, "_")[[1]][1])
  
  # Create a bottom annotation with the virus name
  bottom_anno <- HeatmapAnnotation(
    Virus = rep(virus, ncol(group_data)),  # Repeat the virus name for each column
    col = list(Virus = c(H1N1 = "gray", SARS_Cov = "orange", SARS_dORF6 = "green", SARS_BatSRBD = "purple")),
    annotation_name_side = 'right',
    show_annotation_name = FALSE# Show the annotation title on the left
  )
  
  #Annotation of overall regulation direction
  row_anno <- rowAnnotation(
    Regulation = ifelse(rowMeans(group_data) > 0, "Up", "Down"),  # Pass the vector, not a data frame
    col = list(Regulation = c("Up" = "#d73027", "Down" = 'darkblue')),  #Define colors for the annotation
    show_legend = FALSE
  )
  
  Heatmap(
    group_data[row_order, ], # Apply consistent row order
    name = virus,
    col = col_fun,
    cluster_rows = FALSE,          # Disable row clustering to keep order consistent  
    cluster_columns = TRUE,  # Cluster only columns for each virus
    show_row_names = FALSE,  # Hide row names (optional)
    show_heatmap_legend = FALSE,
    bottom_annotation = bottom_anno,
    right_annotation = row_anno,
    )
})

# Create a shared legend
fold_change_legend <- Legend(
  title = "Fold Change", 
  col_fun = col_fun, 
  at = c(-10, 0, 10), 
  labels = c("-10", "0", "10")
)

# Create a shared legend for Direction
direction_legend <- Legend(
  title = "Overall regulation",
  labels = c("Up", "Down"),
  legend_gp = gpar(fill = c("#d73027", "darkblue"))
)

# Combine the two legends vertically
combined_legends <- packLegend(
  fold_change_legend,
  direction_legend,
  direction = "vertical"
)


# Combine all heatmaps into a single figure
pdf("./03_plots/combined_heatmaps.pdf", width = 12, height = 8)
draw(Reduce(`+`, heatmaps),annotation_legend_list = list(combined_legends))
dev.off()

