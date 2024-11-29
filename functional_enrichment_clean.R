# Install necessary packages if you haven't already
install.packages("cowplot")
install.packages("gplots")
devtools::install_github("nicolash2/ggvenn")
BiocManager::install("DESeq2")  # if DESeq2 is not already installed
install.packages("dplyr")
install.packages("ggplot2")
install.packages("stringr")

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
library("stringr")

# Set working direction
setwd("D:/Documents/Rstudio/Project-SARSHuman/")

#Creation Dataframe
data = read.delim("./01_data/project-SARSHuman.txt" ,sep = "\t")

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

#Create a vector of conditions to compare with Control
conditions_to_compare <- c("SARS.Cov", "SARS.dORF6", "SARS.BatSRBD", "H1N1")

#Functional_enr database
func_enr_data <- data.frame(row.names = rownames(dds))

#We create a dataframe containning every enriched genes as 1 and others as 0 in every condition for the enrichedR package

for (condition in conditions_to_compare) {
  for (time in unique(metadata$timepoint)) {
    condition_samples <- grep(paste0(condition, "_", time, "h"), sample_names)
    mock_samples <- grep(paste0("mock_", time, "h"), sample_names)
    if (length(condition_samples) > 0 && length(mock_samples) > 0) {
      res <- results(dds, contrast = c("condition", condition, "mock"))
      if (!is.null(res)) {
        res_timepoint <- res[metadata$timepoint == time, ]
        if (!is.null(res_timepoint)) {
          temp_data <- as.data.frame(res_timepoint)
          print(paste(condition,time))
          temp_data$enriched = 0
          temp_data$enriched[temp_data$pvalue<0.05] = 1
          
          df = temp_data['enriched']
          names(df)[names(df) == 'enriched'] <- paste(condition,'_',time,'h',sep='')
          
          func_enr_data <- merge(func_enr_data, df, by = "row.names", all = TRUE)
          
          rownames(func_enr_data) <- func_enr_data$Row.names #Magic trick
          func_enr_data$Row.names <- NULL
        }
      }
    }
  }
}



#Gene count ggplot TEST###
fe_plot = function(genenames,subtitle="",database_slc=''){
  
  #Vector of selectable databases
  dbs <- c("GO_Molecular_Function_2023",
           "GO_Cellular_Component_2023",
           "GO_Biological_Process_2023",
           "Human_Gene_Atlas",
           "BioPlanet_2019")
  
  enriched <- enrichr(genenames, dbs)
  
  fe_df = enriched[[database_slc]]
  fe_df = fe_df[order(fe_df$Adjusted.P.value,decreasing = FALSE),]
  fe_df = fe_df[1:15,] #select 15 most significant pathway enriched
  
  fe_df$Term = factor(fe_df$Term,levels = rev(fe_df$Term))
  fe_df$label <- str_wrap(fe_df$Term, width = 40, exdent = 5,whitespace_only=FALSE) #str_trunc(fe_df$Term, width = 40)  #
  
  
  fe_df$gene_count = ifelse(fe_df$P.value < 0.05, str_count(fe_df$Genes, ";"), 0)
  fe_df$Adjusted.P.value = fe_df$Adjusted.P.value * 0.01
  
  plot = ggplot() +
    labs(title="",
         subtitle = paste(subtitle,'enrichment vs mock'), #Change paste statement for regina's data
         tag="") +
    geom_col(data=fe_df,aes(y=Term,
                            x= pmin(pmax(gene_count, 0), 30), #Squishes out of bounds Values
                            fill=Adjusted.P.value)) +
    geom_point(
      data = fe_df %>% filter(gene_count < 0 | gene_count > 30),  # Filter out-of-bound values
      aes(x = pmin(pmax(gene_count, 0), 30), y = Term),  # Plot at squished positions
      shape = 4,  # Cross marker for special ticks
      size = 3,
      color = "gray40"
    )+
    scale_fill_gradient(low = "red", high="blue",
                        oob = scales::squish,  # Include out-of-bounds values
                        limits = c(0.001, 0.01),
                        name = "Adj. P-value\n") +
    scale_y_discrete(labels = fe_df$label,
                     expand = expansion(mult = c(0.05, 0.05))) +
    scale_x_continuous(limits = c(0,30)) +
    xlab("Gene Count")+
    ylab("") +
    theme_bw() +
    theme(plot.subtitle = element_text(face = "italic"),
          legend.position = "right",
          legend.title = element_text(size = 10),
          legend.text = element_text(size = 8),
          axis.text.y = element_text(size = 6,
                                     lineheight = 0.8),
          plot.margin = margin(t = 5,
                               r = 30,
                               b = 10,
                               l = 10),
          axis.ticks.length = unit(0.2, "cm"))
}

fe_plot_countsorted = function(genenames,subtitle="",database_slc=''){
  
  #Vector of selectable databases
  dbs <- c("GO_Molecular_Function_2023",
           "GO_Cellular_Component_2023",
           "GO_Biological_Process_2023",
           "Human_Gene_Atlas",
           "BioPlanet_2019")
  
  enriched <- enrichr(genenames, dbs)
  
  fe_df = enriched[[database_slc]]
  
  fe_df$gene_count = ifelse(fe_df$P.value < 0.05, str_count(fe_df$Genes, ";"), 0)
  fe_df$Adjusted.P.value = fe_df$Adjusted.P.value * 0.01
  
  fe_df = fe_df[order(fe_df$gene_count,decreasing = TRUE),]
  fe_df = fe_df[1:15,] #select 15 most significant pathway enriched
  
  fe_df$Term = factor(fe_df$Term,levels = rev(fe_df$Term))
  fe_df$label <- str_wrap(fe_df$Term, width = 40, exdent = 5,whitespace_only=FALSE) 
  
  plot = ggplot() +
    labs(title="",
         subtitle = paste(subtitle,'enrichment vs mock'), #Change paste statement for regina's data
         tag="") +
    geom_col(data=fe_df,aes(y=Term,
                            x= pmin(pmax(gene_count, 0), 30), #Squishes out of bounds Values
                            fill=Adjusted.P.value)) +
    geom_point(
      data = fe_df %>% filter(gene_count < 0 | gene_count > 30),  # Filter out-of-bound values
      aes(x = pmin(pmax(gene_count, 0), 30), y = Term),  # Plot at squished positions
      shape = 4,  # Cross marker for special ticks
      size = 3,
      color = "gray40"
    )+
    scale_fill_gradient(low = "red", high="blue",
                        oob = scales::squish,  # Include out-of-bounds values
                        limits = c(0.001, 0.01),
                        name = "Adj. P-value\n") +
    scale_y_discrete(labels = fe_df$label,
                     expand = expansion(mult = c(0.05, 0.05))) +
    scale_x_continuous(limits = c(0,30)) +
    xlab("Gene Count")+
    ylab("") +
    theme_bw() +
    theme(plot.subtitle = element_text(face = "italic"),
          legend.position = "right",
          legend.title = element_text(size = 10),
          legend.text = element_text(size = 8),
          axis.text.y = element_text(size = 6,
                                     lineheight = 0.8),
          plot.margin = margin(t = 5,
                               r = 30,
                               b = 10,
                               l = 10),
          axis.ticks.length = unit(0.2, "cm"))
}

#list which containns every plot generated
plt_list = list()

#Open pdf file to write to
pdf(file = "./03_plots/fe_BioPlanet_2019_VSmock_countsorted.pdf", #Change Filename based on function used
    width = 12,
    height = 7)

#Execution plot loop###
for(condition in colnames(func_enr_data)){
  print(condition)
  
  values = func_enr_data[,condition]
  idx = values==1
  genenames = rownames(func_enr_data)[idx]
  
  plot = fe_plot_countsorted(genenames,condition,'BioPlanet_2019')
  plot(plot)
  
  plt_list[[condition]] = plot
}
dev.off()


#Creating grid plots###

comparaisons <- c(12, 24, 36, 48)

pdf(file = "./03_plots/grid_fe_BioPlanet_2019_VSmock_countsorted.pdf", #Change Filename based on function used
    width = 12,
    height = 7)

for(time in comparaisons){
  subset = plt_list[grepl(paste0("_",time, "h"), plt_list)]
  print(plot_grid(plotlist = subset,label_size = 12,greedy=TRUE))
}
dev.off()
