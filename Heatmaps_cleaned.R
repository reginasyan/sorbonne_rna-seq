#loading Libraries
library("ggrepel") 
library("DESeq2") 
library("dplyr")
library('ComplexHeatmap')
library("circlize")
library("reshape2")


setwd("F:/Documents/Rstudio/Project-SARSHuman/")

#Creation Dataframe from raw data
data = read.delim("./01_data/project-SARSHuman.txt" ,sep = "\t")


###############################
###DESeq Enrichment analysis###
###############################

# Define the timepoints for each condition
sars_timepoints <- c(0, 12, 24, 36, 48, 60, 72, 84, 96) 
h1n1_timepoints <- c(0, 6, 12, 18, 24, 36, 48)           

# List of all the conditions you have and a vector without control
conditions <- c("SARS.Cov", "SARS.dORF6", "SARS.BatSRBD", "H1N1", "mock")
conditions_to_compare <- c("SARS.Cov", "SARS.dORF6", "SARS.BatSRBD", "H1N1")

# Sample names
sample_names <- colnames(data)  # This should be your sample names like "H1N1_00h_1"

# Create a vector to store for every conditions as name, the virus tested
condition <- sapply(sample_names, function(x) strsplit(x, "_")[[1]][1])  # First part is condition (e.g., H1N1)
# Create a vector to store fore every condition as name, the timepoints
timepoint <- sapply(sample_names, function(x) strsplit(x, "_")[[1]][2])  # Second part is timepoint (e.g., 00h)
timepoint_numeric <- as.numeric(gsub("h", "", timepoint)) # Remove the 'h' from timepoints and convert to numeric values

# Create the metadata dataframe with replicates grouped by condition and timepoint
metadata <- data.frame(
  condition = factor(condition, levels = conditions),  #The vector Conditions is used as a factor
  timepoint = timepoint_numeric # Timepoint as numeric
)

# Column condition.timepoint allow to factor the 3 replicats of a condition for each timepoint
metadata$condition_timepoint = interaction(metadata$condition, metadata$timepoint, sep = ".")

# Creating DESeq2 dataset from count data and metadata
dds <- DESeqDataSetFromMatrix(countData = round(data), colData = metadata, design = ~ condition_timepoint) #Allows to compare between different condition and adress each timepoint

# Run the DESeq analysis
dds <- DESeq(dds)

######################################
###Creation of enrichment Dataframe###
######################################

# We create empty dataframe with genes as row names
result_table <- data.frame(row.names = rownames(dds))

# we loop so that every enriched genes appears as 1 or -1 and others as 0 in every condition
for(conds in conditions_to_compare) {
  time_vect <- c() # We ensure time points exists for each condition with following if()
  
  if(conds %in% c("SARS.Cov", "SARS.dORF6", "SARS.BatSRBD")) {
    time_vect = sars_timepoints
  }else{
    time_vect = h1n1_timepoints
  }
  
  for(time in time_vect) { 
    temp_data <- data.frame()
    
    #The comparaison is made bewteen a condition i.e SARS.Cov_24h and Mock_24h for example###
    comp <- c("condition_timepoint", paste0(conds, ".", time), paste0("mock.", time)) #Creates the comparison vector for the given condition and time point to give to DESeq2 function
    
    res <- results(dds, contrast = comp) #Extract DESeq results from dds for this comparison
    temp_data <- as.data.frame(res) # Convert DESeq result to dataframe to extract the relevant info
    
    print(paste0(conds,'_',time,'h')) # Keeps track of function
    
    genes <- res[!is.na(res$padj) & !is.na(res$log2FoldChange), ]
    genes$timepoint <- time
    genes$dataset <- conds
    result_table <- rbind(result_table, data.frame(
      gene = rownames(genes),
      fold_change = genes$log2FoldChange,
      p_value_adj = genes$padj,
      timepoint = genes$timepoint,
      dataset = genes$dataset)
    )
  }
}


#write.table(result_table,file = "./01_data/DESEQ_enrichment_heatmap.txt",sep = "\t",quote = FALSE)

DFE <- result_table

#DFE <- read.delim("./01_data/DESEQ_enrichment_heatmap.txt" ,sep = "\t")

DFE$timepoint[DFE$timepoin == '0h'] = '00h'
DFE$timepoint[DFE$timepoin == '6h'] = '06h'

DFE$direction = ifelse((DFE$fold_change > 1) & (DFE$p_value_adj < 0.05), "Upregulated", 
                       ifelse((DFE$fold_change < -1) & (DFE$p_value_adj < 0.05),"Downregulated",
                              ifelse(DFE$p_value_adj < 0.05,'psignificant',"NonSignificant")))


#Subsets significant genes from DFE
Sig_genes = subset(DFE,DFE$direction != 'NonSignificant')
write.table(result_table,file = "./01_data/sig_genes_est.txt",sep = "\t",quote = FALSE)

#Create matrix data for heatmaps
heatmap_data_3 <- dcast(DFE, gene ~ timepoint + dataset, value.var = "fold_change")
# Convert to matrix format (row names = genes)
rownames(heatmap_data_3) <- heatmap_data_3$gene
heatmap_data_3 <- as.matrix(heatmap_data_3[, -1])  # Remove the 'gene' column

heatmap_data_3[is.na(heatmap_data_3)] <- 0 


##############################
#####Multiple Virus Heatmap###
##############################

#Groups of viruses for heatmap generation

virus_groups_data3 <- list(
  H1N1 = grep("H1N1", colnames(heatmap_data_3), value = TRUE),
  SARS_Cov = grep("SARS.Cov", colnames(heatmap_data_3), value = TRUE),
  SARS_dORF6 = grep("SARS.dORF6", colnames(heatmap_data_3), value = TRUE),
  SARS_BatSRBD = grep("SARS.BatSRBD", colnames(heatmap_data_3), value = TRUE)
)


# Create a consistent row order using hierarchical clustering
row_order <- hclust(dist(heatmap_data_3))$order

#Shared color function
col_fun <- colorRamp2(c(-2, 0, 2), c("darkblue", "#f7f7f7", "#d73027"))

# Create individual heatmaps for each virus group
heatmaps <- lapply(names(virus_groups_data3), function(virus) {
  print(virus)
  #Subset specific virus data
  group_data <- heatmap_data_3[, virus_groups_data3[[virus]], drop = FALSE] # Filter columns
  
  #Calculate range of data for auto scalling of log fold colors
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
  
  #Calculate the mean value of only p significant genes on each row which are logfold significant (Avoids Zero values in log fold mean)
  bool_data_mean <- replace(group_data,group_data > 1, 1)
  bool_data_mean <- replace(group_data,group_data < -1, -1)
  bool_data_mean <- replace(group_data,group_data > -1 & group_data < 1, 0)
  mean_data <- rowMeans(bool_data_mean,na.rm=TRUE)
  
  #mean_data <- rowMeans(replace(bool_data_mean, bool_data_mean == 0, NA), na.rm = TRUE)
  mean_data[is.na(mean_data)] <- 0 
  
  #Compute regulation vector based on logFold gatting if a gene is in average upregulated or downregulated for all timestamp (At least a single gene is upregulated in a single timepoint for a single virus)
  #This code is structured to put a threshold for two different series of timepoints
  regulation_vector <- if (virus == 'H1N1') {
    ifelse(mean_data > (1/7), "Up",
           ifelse(mean_data < -(1/7), "Down", "Ns"))
  } else {
    ifelse(mean_data > (1/9), "Up",
           ifelse(mean_data < -(1/9), "Down", "Ns"))
  }
  print(mean_data[mean_data>1])
  
  #Annotation of overall regulation direction
  row_anno <- rowAnnotation(
    Regulation = regulation_vector,
    col = list(Regulation = c("Up" = "#d73027", "Down" = 'darkblue',"Ns"='grey80')),  #Define colors for the annotation
    show_legend = FALSE
  )
  
  Heatmap(
    group_data[row_order, ], # Apply consistent row order
    name = virus,
    col = col_fun,
    cluster_rows = TRUE,
    show_row_dend = FALSE, # Disable row clustering to keep order consistent  
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
  at = c(-2, 0, 2), 
  labels = c("-2", "0", "2")
)

# Create a shared legend for Direction
direction_legend <- Legend(
  title = "Overall regulation",
  labels = c("Up","Non significant", "Down"),
  legend_gp = gpar(fill = c("#d73027","grey80", "darkblue"))
)

# Combine the two legends vertically
combined_legends <- packLegend(
  fold_change_legend,
  direction_legend,
  direction = "vertical"
)


# Combine all heatmaps into a single figure
pdf("./03_plots/Corrected_heatmaps/sig_genes/combined_heatmaps_4viruses_scale2_pathway_dev2.pdf", width = 12, height = 8)
draw(Reduce(`+`, heatmaps),annotation_legend_list = list(combined_legends))
dev.off()


#######################################
##HEATMAP of dorf and BAT  clustering##
#######################################

# Subset the columns for dORF6 and BatSRBD
dORF6_cols  <- virus_groups_data3[["SARS_dORF6"]]
batSRBD_cols <- virus_groups_data3[["SARS_BatSRBD"]]

subset_data <- heatmap_data_3[, c(dORF6_cols, batSRBD_cols), drop = FALSE]

col_virus <- data.frame(
  Virus = c(
    rep("SARS.dORF6", length(dORF6_cols)),
    rep("SARS.BatSRBD", length(batSRBD_cols))
  )
)

rownames(col_virus) <- colnames(subset_data)
colnames(subset_data) <- sapply(colnames(subset_data), function(x) strsplit(x, "_")[[1]][1])

# Use the same color function you've been using
col_fun <- colorRamp2(c(-2, 0, 2), c("darkblue", "#f7f7f7", "#d73027"))

# Define color mapping for the column annotation
col_anno <- HeatmapAnnotation(
  df = col_virus,
  col = list(
    Virus = c("SARS.dORF6"="green", "SARS.BatSRBD"="purple")
  ),
  show_legend = FALSE
)

# Create a shared legend
fold_change_legend <- Legend(
  title = "Fold Change", 
  col_fun = col_fun, 
  at = c(-2, 0, 2), 
  labels = c("-2", "0", "2")
)

# Create a separate legend for the Virus annotation
virus_legend <- Legend(
  title = "Virus",
  at = c("SARS.dORF6", "SARS.BatSRBD"),
  legend_gp = gpar(fill = c("green", "purple"))
)

# Pack the two legends vertically (stacked)
combined_legends <- packLegend(
  fold_change_legend,
  virus_legend,
  direction = "vertical"
)


heatmap_dorf_bat <- Heatmap(
  subset_data,
  name = "FoldChange",       # or any label
  col  = col_fun,
  cluster_rows = TRUE,
  show_row_dend = FALSE,    
  cluster_columns = TRUE,    # Let Heatmap cluster the columns
  top_annotation  = col_anno, 
  show_row_names  = FALSE,   # optional
  show_column_names = TRUE,  
  show_heatmap_legend = FALSE,
  column_names_side = "top"
)

# Draw the heatmap
pdf("./03_plots/Corrected_heatmaps/sig_genes/dORF6_vs_BatSRBD_ClusteredRow_data3_Scale2.pdf", width = 8, height = 10)
draw(
  heatmap_dorf_bat,
  annotation_legend_side = "right",
  annotation_legend_list = list(combined_legends)
)
dev.off()


############################
#Heatmap Based on timepoint#
############################

# Extract timepoints from column names
all_colnames <- colnames(heatmap_data_3)
timepoints <- sapply(all_colnames, function(x) strsplit(x, "_")[[1]][1])

# Get unique timepoints in ascending order
unique_timepoints <- unique(timepoints)

#Filter out unique H1N1 timepoints
unique_timepoints <- unique_timepoints[unique_timepoints != "06h"]
unique_timepoints <- unique_timepoints[unique_timepoints != "18h"]

col_fun <- colorRamp2(c(-5, 0, 5), c("darkblue", "#f7f7f7", "#d73027"))


#Heatmap function
make_timepoint_heatmap <- function(tp, full_data, col_fun) {
  # Identify columns that match this timepoint
  cols_for_tp <- grep(paste0("^", tp, "_"), colnames(full_data), value = TRUE)
  
  # Subset the heatmap_data for these columns
  tp_data <- full_data[, cols_for_tp, drop = FALSE]
  
  
  # Create a column annotation, labeling the virus for each column
  virus_labels <- sapply(colnames(tp_data), function(x) strsplit(x, "_")[[1]][2])
  col_anno_df <- data.frame(Virus = virus_labels)
  rownames(col_anno_df) <- colnames(tp_data)
  
  # A color map for each virus (adjust or expand as needed)
  virus_colors <- c(
    "H1N1"       = "gray", 
    "SARS.Cov"   = "orange", 
    "SARS.dORF6" = "green",
    "SARS.BatSRBD" = "purple"
  )
  
  # Build the annotation
  col_anno <- HeatmapAnnotation(
    df = col_anno_df,
    col = list(Virus = virus_colors),
    show_legend = FALSE,
    show_annotation_name = FALSE
  )
  
  # Create the Heatmap object
  Heatmap(
    tp_data,
    column_title = tp,           # The label (timepoint) for this heatmap
    column_title_side = "bottom",
    col = col_fun,           # color scale from your existing col_fun
    cluster_rows = TRUE,
    show_row_dend = FALSE,# cluster rows
    cluster_columns = TRUE,  # cluster columns
    show_row_names = FALSE,
    show_column_names = FALSE,
    top_annotation = col_anno,
    show_heatmap_legend = FALSE,  # we'll add a shared legend later
    column_names_side = "top"     # put column names on top
  )
}

# Make sure you have your color function defined for the fold-change
# For example:
# col_fun <- colorRamp2(c(-2, 0, 2), c("darkblue", "#f7f7f7", "#d73027"))

# Build one heatmap per timepoint
timepoint_heatmaps <- lapply(unique_timepoints, function(tp) {
  make_timepoint_heatmap(tp, heatmap_data_3, col_fun)
})

# Fold-change legend
fold_change_legend <- Legend(
  title = "Fold Change", 
  col_fun = col_fun, 
  at = c(-5, 0, 5),  # or c(-5, 0, 5) if using a wider range
  labels = c("-5", "0", "5")
)

# Virus legend
virus_legend <- Legend(
  title = "Virus",
  at = c("H1N1", "SARS.Cov", "SARS.dORF6", "SARS.BatSRBD"),
  legend_gp = gpar(fill = c("gray", "orange", "green", "purple"))
)
# Pack them vertically
combined_legends <- packLegend(
  fold_change_legend,
  virus_legend,
  direction = "vertical"
)

pdf("./03_plots/Corrected_heatmaps/sig_genes/AllViruses_SideBySide_Timepoints_Clustered_data3_Scale2.pdf", width = 14, height = 8)

# Combine the heatmaps side by side using Reduce(`+`, ...)
draw(
  Reduce(`+`, timepoint_heatmaps),
  annotation_legend_side = "right",
  annotation_legend_list = list(combined_legends)
)

dev.off()
