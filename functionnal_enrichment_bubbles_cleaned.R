#libraries
library("ggrepel") 
library("enrichR") 
library("ggplot2") 
library("DESeq2") 
library("dplyr")
library("reshape2")
library("stringr")


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

###############################
###DESeq Enrichment analysis###
###############################

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
    # Extract significant genes (adjusted p-value < 0.05 and |log2FoldChange| > 1)
    res <- res[!is.na(res$padj) & !is.na(res$log2FoldChange), ]
    sig_genes <- res[res$padj < 0.05 & abs(res$log2FoldChange) > 1, ]
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

# Save the combined table to CSV
write.csv(combined_table, "./01_data/combined_results.csv", row.names = FALSE)

#combined_table = read.delim("./01_data/combined_results.csv" ,sep = ",")

##########################
###Function DEF section###
##########################


buble_plot = function(data_enriched,subtitle="",database_slc='') {
  #Needs to be feed combined_table in data_enr argument
  
  #Placeholder DataFrame for empty data
  placeholder = data.frame(Term = character(0),
                           regulation = character(0),
                           gene_count = numeric(0),
                           gene_ratio = numeric(0),
                           Adjusted.P.value = numeric(0)
                           )
  
  virus = strsplit(subtitle,"_")[[1]][1]
  time = strsplit(subtitle,"_")[[1]][2]
    
  enr_data = subset(combined_table,timepoint==time & dataset == virus, select = c(gene,direction))
  
  genenames_up = enr_data$gene[enr_data$direction == 'Upregulated']
  genenames_down = enr_data$gene[enr_data$direction == 'Downregulated']

  if(length(genenames_up) > 0) {
    enriched_up = enrichr(genenames_up, database_slc)
    up = enriched_up[[database_slc]]
    up = up[up$P.value < 0.05,] # Subset the vector for only the significant patwhays
    
  }else{
    up = placeholder
  }
  
  if(length(genenames_down) > 0) {
    enriched_down = enrichr(genenames_down, database_slc)
    down = enriched_down[[database_slc]]
    down = down[down$P.value < 0.05,] # Subset the vector for only the significant patwhays
  }else{
    down = placeholder
  }

  
  #Checks if the upregulated subset exists
  if(nrow(up) > 0) {
    up$regulation = 'Up-Regulated'
    up$gene_count = str_count(up$Genes, ";") + 1
    up$overlap_int = as.numeric(sapply(strsplit(up$Overlap,split='/'), `[`, 1))
    up$totalpath_int = as.numeric(sapply(strsplit(up$Overlap,split='/'), `[`, 2))
    up$gene_ratio = up$overlap_int / up$totalpath_int  #Gene ratio is calculated by number of overlapping gene in a pathway divided by overall imput number
    up = up[order(up$Adjusted.P.value,decreasing = FALSE),] #Sort By #Adjusted.P.value
    
    ifelse(nrow(down)>0, up <- up[!up$Term %in% down$Term, ],NA)   #Excluding overlapping pathways
    ifelse(nrow(down)>0, up <- up[1:10,], up <- up[1:20,]) #Select enriched Pathways
    
  } else { 
    #Placeholder Data
    up <- placeholder
  }
  
  #Checks if donwregulated subset exists
  if(nrow(down) > 0) {
    down$regulation = 'Down-Regulated'
    down$gene_count = str_count(down$Genes, ";") + 1
    down$overlap_int = as.numeric(sapply(strsplit(down$Overlap,split='/'), `[`, 1))
    down$totalpath_int = as.numeric(sapply(strsplit(down$Overlap,split='/'), `[`, 2))
    down$gene_ratio = down$overlap_int / down$totalpath_int  #Gene ratio is calculated by number of overlapping gene in a pathway divided by overall imput number
    down = down[order(down$Adjusted.P.value,decreasing = FALSE),] #Sort By #Adjusted.P.value
    
    ifelse(nrow(up)>0,down <- down[!down$Term %in% up$Term, ],NA)   #Excluding overlapping pathways
    ifelse(nrow(up)>0, down <- down[1:10,], down <- down[1:20,]) #Select enriched Pathways
    
  } else {
    #Placeholder Data
    down <- placeholder
  }
  
  buble <- rbind(up,down)
  
  #Omit NA rows if less than 10 pathways enriched for each
  buble <- na.omit(buble)
  
  # In case Both Up and Down are empty and not significant
  if(nrow(buble) == 0) {
    message("No significant pathways found for either up- or down-regulated genes.")
    return(ggplot() + ggtitle(paste(subtitle,'enrichment vs mock',"No significant pathways to display")))
  }
  
  
  buble <- buble[order(buble$gene_ratio,decreasing = FALSE),] #Ensures Buble is also sorted by gene_ratio
  
  # Order by gene_ratio
  buble <- buble %>% mutate(Term_wrapped = str_wrap(Term, width = 40),
                            Term_wrapped = factor(Term_wrapped, levels = unique(Term_wrapped))
  )
  
  plot = ggplot(buble, aes(x = gene_ratio, y = Term_wrapped, size = gene_count, color = Adjusted.P.value)) +
    geom_point() +
    scale_color_gradient(low = "red", high = "blue") +
    facet_wrap(~regulation) +
    theme_bw() +
    labs(
      x = "Gene Ratio",
      y = "",
      size = "Gene Count",
      color = "Adjusted P-value",
      title = paste(subtitle,'enrichment vs mock')
    ) + 
    theme(
      legend.position = "right",  
      legend.title = element_text(size = 10),
      legend.text = element_text(size = 8),
      axis.text.y = element_text(size = 8, lineheight = 0.8),  # Adjust y-axis text
      plot.subtitle = element_text(face = "italic"),
      plot.margin = margin(t = 10,
                           r = 10,
                           b = 10,
                           l = 10))
  return(plot)
}


buble_plot_simple = function(data_enriched,subtitle="",database_slc='') {
  #Needs to be feed combined_table in data_enr argument
  #A simple version of the buble plot version without Up/Down regulation differenciation
  
  #Placeholder DataFrame for empty data
  placeholder = data.frame(Term = character(0),
                           regulation = character(0),
                           gene_count = numeric(0),
                           gene_ratio = numeric(0),
                           Adjusted.P.value = numeric(0)
  )
  
  virus = strsplit(subtitle,"_")[[1]][1]
  time = strsplit(subtitle,"_")[[1]][2]
  
  enr_data = subset(combined_table,timepoint==time & dataset == virus, select = c(gene,direction))
  
  genes = enr_data$gene
  
  if(length(genes) > 0) {
    enr = enrichr(genes, database_slc)
    enriched = enr[[database_slc]]
    enriched = enriched[enriched$P.value < 0.05,] # Subset the vector for only the significant patwhays
    
  }else{
    enriched = placeholder
  }
  
  #Checks if the upregulated subset exists
  if(nrow(enriched) > 0) {
    enriched$regulation = 'Enriched_Pathways'
    enriched$gene_count = str_count(enriched$Genes, ";") + 1
    enriched$overlap_int = as.numeric(sapply(strsplit(enriched$Overlap,split='/'), `[`, 1))
    enriched$totalpath_int = as.numeric(sapply(strsplit(enriched$Overlap,split='/'), `[`, 2))
    enriched$gene_ratio = enriched$overlap_int / enriched$totalpath_int #16599 #Gene ratio is calculated by number of overlapping gene in a pathway divided by overall imput number
    enriched = enriched[order(enriched$P.value,decreasing = FALSE),] #Sort By gene ratio
    
    enriched <- enriched[1:20,] #Select most enriched Pathways
   
  } else { 
    #Placeholder Data
    enriched <- placeholder
  }
  
  buble <- enriched
  
  #Omit NA rows if less than 10 pathways enriched for each
  buble <- na.omit(buble)
  
  # In case Both Up and Down are empty and not significant
  if(nrow(buble) == 0) {
    message("No significant pathways found for either up- or down-regulated genes.")
    return(ggplot() + ggtitle(paste(subtitle,'enrichment vs mock',"No significant pathways to display")))
  }
  
  
  buble <- buble[order(buble$gene_ratio,decreasing = FALSE),] #Ensures Buble is also sorted by gene_ratio
  
  # Order by gene_ratio
  buble <- buble %>% mutate(Term_wrapped = str_wrap(Term, width = 40),
                            Term_wrapped = factor(Term_wrapped, levels = unique(Term_wrapped))
  )
  
  plot = ggplot(buble, aes(x = gene_ratio, y = Term_wrapped, size = gene_count, color = Adjusted.P.value)) +
    geom_point() +
    scale_color_gradient(low = "red", high = "blue") +
    facet_wrap(~regulation) +
    theme_bw() +
    labs(
      x = "Gene Ratio",
      y = "",
      size = "Gene Count",
      color = "Adjusted P-value",
      title = paste(subtitle,'enrichment vs mock')
    ) + 
    theme(
      legend.position = "right",  # Place legends to the right
      legend.title = element_text(size = 10),
      legend.text = element_text(size = 8),
      axis.text.y = element_text(size = 8, lineheight = 0.8),  # Adjust y-axis text
      plot.subtitle = element_text(face = "italic"),
      plot.margin = margin(t = 10,
                           r = 10,
                           b = 10,
                           l = 10))
  return(plot)
}


##########################
###BUBLE PLOT CALL LOOP###
##########################

#Createas a list of condition_timpeoint from colnames(data)
sample_names <- colnames(data)

char_array = sample_names
a = data.frame("data"=char_array,"data2"= 1:length(sample_names))
a$data = substr(a$data,1,nchar(a$data)-2)

condition_a = unique(a$data)

#Lists of conds ie H1N1_12h
condition_list <- condition_a[!grepl("mock_", condition_a)]

condition_list <- c("H1N1_12h","H1N1_48h","SARS.dORF6_12h","SARS.dORF6_72h")

#Vector of selectable databases for enrichR
dbs <- c("GO_Biological_Process_2021",
         "BioPlanet_2019",
         "Reactome_2022",
         "KEGG_2021_Human",
         "ENCODE_and_ChEA_Consensus_TFs_from_ChIP-X")

dbs <- c("BioPlanet_2019")

for(database in dbs){
  print(paste('Parsing results for :',database))
  
  plt_buble_list = list() #List of buble plots
  #Open pdf file to write to
  pdf(file = paste("./03_plots/bubble2/enrichment_",database,"_VSmock_bubblesPlot_Ajustedpval_nopgate_foldchangegen11.pdf",sep="" ), #Change Filename based on function used
      width = 7.5,
      height = 6)
  
  #Looping
  for(condition in condition_list) {
    print(condition)
    
    plot = buble_plot(combined_table,condition,database)
    plot(plot)
    plt_buble_list[[condition]] = plot #Saves plot to the environement
    
  }
  dev.off()
}

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
    genes <- genes[genes$padj < 0.05 & abs(genes$log2FoldChange) > 0, ]
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

# Combine all tables into one
combined_table <- result_table
combined_table$direction = ifelse(combined_table$fold_change > 0, "Upregulated", "Downregulated")

