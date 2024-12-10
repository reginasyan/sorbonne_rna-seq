# Load libraries
library("ggrepel")
library("enrichR")
library("ggplot2")
library("gplots")
library("DESeq2")
library("stringr")
library("dplyr")
library('cowplot')
library("reshape2")

# Set working directory
setwd("F:/Documents/Rstudio/Project-SARSHuman/")

#Creation Dataframe
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
func_enr_data <- data.frame(row.names = rownames(dds))

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
    
    ###GATING PARAMETERS###
    temp_data$enriched = 0
    temp_data$enriched[(temp_data$padj<0.05) & (temp_data$log2FoldChange > 0)] = 1 # Append column in res dataframe based on p. value adjusted
    temp_data$enriched[(temp_data$padj<0.05) & (temp_data$log2FoldChange < 0)] = -1
    
    func_enr_data[paste0(conds,'_',time,'h')] = temp_data$enriched # Append results to the global func_enr_data dataframe
  }
}

#Saving results
write.table(func_enr_data,file = "./01_data/func_enr_data_test.txt",sep = "\t",quote = FALSE)


##########################
###Function DEF section###
##########################

#Function creating a graph for each enriched gene name list for a given condition in an enrichR database
fe_plot_bar_ratiosorted = function(genenames,subtitle="",database_slc=''){
  
  #Vector of selectable databases
  dbs <- c("GO_Molecular_Function_2023",
           "GO_Cellular_Component_2023",
           "GO_Biological_Process_2023",
           "Human_Gene_Atlas",
           "BioPlanet_2019")
  
  enriched <- enrichr(genenames, dbs)
  
  fe_df = enriched[[database_slc]]
  
  fe_df = fe_df[fe_df$Adjusted.P.value < 0.05,] # Subset the vector for only the significant patwhays
  fe_df = fe_df[order(fe_df$Adjusted.P.value,decreasing = FALSE),] #Sort by Adjusted Pval
 
  fe_df = fe_df[1:15,] #select 15 most significant pathway enriched -> Chose number of pathways to compare
  
  fe_df$overlap_int = as.numeric(sapply(strsplit(fe_df$Overlap,split='/'), `[`, 1))
  fe_df$totalpath_int = as.numeric(sapply(strsplit(fe_df$Overlap,split='/'), `[`, 2))
  fe_df$gene_ratio = fe_df$overlap_int / fe_df$totalpath_int #Gene ratio is calculated by number of overlapping gene in a pathway divided by number in pathway
  
  
  fe_df$Term = factor(fe_df$Term,levels = rev(fe_df$Term))
  fe_df$label <- str_wrap(fe_df$Term, width = 40, exdent = 5,whitespace_only=FALSE) #Wraps y labels correctly
  
  #Create Plot object for the function call
  plot = ggplot() +     
    labs(title="",
         subtitle = paste(subtitle,'enrichment vs mock'), #Change paste statement for regina's data
         tag="") +
    geom_col(data=fe_df,aes(y=Term,x=gene_ratio,fill=Adjusted.P.value)) +
    scale_fill_gradient(low = "red", high="blue",
                        oob = scales::squish,  #Include out-of-bounds values
                        limits = c(0.0001, 0.01),
                        name = "Adj. P-value\n") +
    scale_y_discrete(labels = fe_df$label,
                     expand = expansion(mult = c(0.05, 0.05))) +
    xlab("Gene Ratio")+
    ylab("") +
    theme_bw() +
    theme(plot.subtitle = element_text(face = "italic"),
          legend.position = "right",
          legend.title = element_text(size = 10),
          legend.text = element_text(size = 8),
          axis.text.y = element_text(size = 6, #Change this value to 8 for single and 6 for grid
                                     lineheight = 0.8),
          plot.margin = margin(t = 5,
                               r = 30,
                               b = 10,
                               l = 10),
          axis.ticks.length = unit(0.2, "cm"))
    return(plot)
}

buble_plot = function(data_enriched,subtitle="",database_slc='') {
  """Needs to be fed func_enr_data object in data_enr argument"""
  
  #Vector of selectable databases
  dbs <- c("GO_Molecular_Function_2023",
           "GO_Cellular_Component_2023",
           "GO_Biological_Process_2023",
           "Human_Gene_Atlas",
           "BioPlanet_2019")
  
  enr_data = data_enriched[,subtitle]
  idx = enr_data==1
  genenames_up = rownames(data_enriched)[idx]
  
  idx = enr_data==-1 
  genenames_down = rownames(data_enriched)[idx]
  
  enriched_up <- enrichr(genenames_up, dbs)
  enriched_down <- enrichr(genenames_down, dbs)
  
  up = enriched_up[[database_slc]]
  up = up[up$Adjusted.P.value < 0.05,] # Subset the vector for only the significant patwhays
  
  down = enriched_down[[database_slc]]
  down = down[down$Adjusted.P.value < 0.05,] # Subset the vector for only the significant patwhays
  
  #Checks if the subset exist
  if(nrow(up) > 0) {
    up$regulation = 'Up-Regulated'
    up$gene_count = str_count(up$Genes, ";") + 1
    up$overlap_int = as.numeric(sapply(strsplit(up$Overlap,split='/'), `[`, 1))
    up$totalpath_int = as.numeric(sapply(strsplit(up$Overlap,split='/'), `[`, 2))
    up$gene_ratio = up$overlap_int / up$totalpath_int #Gene ratio is calculated by number of overlapping gene in a pathway divided by overall number in a given path
    up = up[order(up$Adjusted.P.value,decreasing = FALSE),] #Sort By Combined.Score
    ifelse(nrow(down)>0, up <- up[1:10,], up <- up[1:20,]) #select 15 most significant pathway enriched -> Chose number of pathways to compare
  } else { 
    #Placeholder Data
    up <- data.frame(
      Term = character(0),
      regulation = character(0),
      gene_count = numeric(0),
      gene_ratio = numeric(0),
      Adjusted.P.value = numeric(0))
  }
  
  if(nrow(down) > 0) {
    down$regulation = 'Down-Regulated'
    down$gene_count = str_count(down$Genes, ";") + 1
    down$overlap_int = as.numeric(sapply(strsplit(down$Overlap,split='/'), `[`, 1))
    down$totalpath_int = as.numeric(sapply(strsplit(down$Overlap,split='/'), `[`, 2))
    down$gene_ratio = down$overlap_int / down$totalpath_int #Gene ratio is calculated by number of overlapping gene in a pathway divided by overall number in a given path
    down = down[order(down$Adjusted.P.value,decreasing = FALSE),] #Sort By Combined.Score
    ifelse(nrow(up)>0, down <- down[1:10,], down <- down[1:20,]) #select 15 most significant pathway enriched -> Chose number of pathways to compare
  } else {
    #Placeholder Data
    down <- data.frame(
      Term = character(0),
      regulation = character(0),
      gene_count = numeric(0),
      gene_ratio = numeric(0),
      Adjusted.P.value = numeric(0))
  }

  #Creates buble data frame with both up and down regulated
  buble <- rbind(up,down)
  
  #Omit NA rows if less than 10 pathways enriched for each
  buble <- na.omit(buble)
  
  # In case Both Up and Down are empty and not significant
  if (nrow(buble) == 0) {
    message("No significant pathways found for either up- or down-regulated genes.")
    return(ggplot() + ggtitle("No significant pathways to display"))
  }
  
  #We Ensures Buble is also sorted by gene_ratio
  buble <- buble[order(buble$gene_ratio,decreasing = FALSE),]
  
  # Order by gene_ratio
  buble <- buble %>% mutate(Term_wrapped = str_wrap(Term, width = 40),
                            Term_wrapped = factor(Term_wrapped, levels = unique(Term_wrapped))
                            )

  #We create a plot object for the Bubble plot
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
                 
plt_buble_list = list() #List of buble plots

#Open pdf file to write to
pdf(file = "./03_plots/fe_BioPlanet_2019_VSmock_bubblesPlot.pdf", #Change Filename based on function used
    width = 12,
    height = 7)

#Looping
for(condition in colnames(func_enr_data)) {
  print(condition)
  
  plot = buble_plot(func_enr_data,condition,'BioPlanet_2019')
  plot(plot)
  
  plt_buble_list[[condition]] = plot #Saves plot to the environement
  
}
dev.off()

buble_plot(func_enr_data,"SARS.Cov_24h",'BioPlanet_2019')


########################
###BAR PLOT CALL LOOP###
########################

#list which containns every plot generated
plt_list = list()

#Open pdf file to write to
pdf(file = "./03_plots/fe_BioPlanet_2019_VSmock_combsorted_test.pdf", #Change Filename based on function used
    width = 12,
    height = 7)

#Execution plot loop for barpplot###
for(condition in colnames(func_enr_data)){
  print(condition)
  
  values = func_enr_data[,condition]
  idx = values==1
  genenames = rownames(func_enr_data)[idx]
  
  plot = fe_plot_bar_ratiosorted(genenames,condition,'BioPlanet_2019')
  plot(plot) 
  
  plt_list[[condition]] = plot #Saves plot to the environement
}
dev.off()

#Creating grid plots for bars from plt_list object###
comparaisons <- c(12, 24, 36, 48)

pdf(file = "./03_plots/grid_fe_BioPlanet_2019_VSmock_combsorted_test.pdf", #Change Filename based on function used
    width = 12,
    height = 7)

for(time in comparaisons){
  subset = plt_list[grepl(paste0("_",time, "h"), plt_list)]
  print(plot_grid(plotlist = subset,label_size = 12,greedy=TRUE))
}
dev.off()
