# R version 3.4.2 (2017-09-28)
# Platform: x86_64-w64-mingw32/x64 (64-bit)
# Running under: Windows 7 x64 (build 7601) Service Pack 1

# Packages - Bioconductor
if (!require(flowCore)) {source("https://bioconductor.org/biocLite.R")
  biocLite("flowCore")}
if (!require(FlowSOM)) {source("https://bioconductor.org/biocLite.R")
  biocLite("FlowSOM")}
if (!require(Biobase)) {source("https://bioconductor.org/biocLite.R")
  biocLite("Biobase")}

# Packages - Cran
if(!require(gtools)) install.packages("gtools")
if(!require(ggplot2)) install.packages("ggplot2")
if(!require(colorRamps)) install.packages("colorRamps")
if(!require(dplyr)) install.packages("dplyr")
if(!require(grid)) install.packages("grid")
if(!require(gridExtra)) install.packages("gridExtra")
if(!require(png)) install.packages("png")
if(!require(reshape2)) install.packages("reshape2")

# Packages - Github
if(!require(devtools)) install.packages("devtools") # If not already installed
if(!require(Rtsne.multicore)) devtools::install_github("RGLab/Rtsne.multicore")
if(!require(Rphenograph)) devtools::install_github("JinmiaoChenLab/Rphenograph")


# load all the librarys
library(gtools)
library(flowCore)         
library(Rtsne.multicore)
library(ggplot2)
library(Rphenograph)
library(FlowSOM)
library(colorRamps)
library(MEM)
library(dplyr)
library(grid)
library(gridExtra)
library(png)
library(Biobase)
library(reshape2)


set.seed(42)


## functions

getCyTOFData <- function(filename){
  
  fs <- read.FCS(filename)
  
  #params <- paste(fs@parameters@data$name, "_", fs@parameters@data$desc, sep = "")
  
  raw_data <- as.data.frame(fs@exprs)
  
  #params <- gsub("^.*?_", "", params)
  
  #names(raw_data) <- params 
  
  raw_data$Filename <- filename
  
  raw_data
}

getFlowData <- function(filename){
  
  fs <- read.FCS(filename)
  
  comp_matrix <- keyword(fs)$SPILL                                            ## generates the compensation matrix
  
  fs <- compensate(fs, comp_matrix) 
  
  raw_data <- as.data.frame(fs@exprs)
  
  gr <- grep("P(\\dS$|\\d\\dS$)", names(keyword(fs)))
  
  a <- unlist(fs@description[gr])[order(names(unlist(fs@description[gr])))]
  
  params <- a[mixedorder(names(a))]
  
  names(params) <- NULL
  
  names(raw_data) <- params
  
  raw_data$Filename <- filename
  
  raw_data
  
}

# 
# getData ----
#
# Function to read the FCS files and concatenate them into one large data.frame.
# Can use this to downsample the number of events for running a quick tsne. 
# 3 variables need to be entered, the files you want to read, True or False 
# for if you want to downsample, and the number you want to downsample to.
#

getData <- function(filenames, CyTOF = F, downsample = F, nsamp = NULL){
  
  combined_data <- NULL
  
  for (i in filenames){
  
    if (CyTOF) {
      raw_data <- getCyTOFData(i)
    } else {
      raw_data <- getFlowData(i)
    }                                                                            ## will be specific for each experiment
   
    
    if (downsample == T){
      raw_data <- raw_data[sample(nrow(raw_data), nsamp),]                                  ## if you want to downsample, this randomly
      raw_data <- raw_data[order(as.numeric(rownames(raw_data))),, drop = F]
    }                                                                           ## takes the number of rows in nsamp
    
    combined_data <- rbind(combined_data, raw_data)                                 ## creates the final table that will be returned
                                                                                ## at the end of the function
  }
  
  combined_data
}

#
# plotting_results ----
#
# Generates a tsne plot for each parameter listed (columns in the data.frame of 
# interest), also generates the flowSOM plot and colours the flowSOM onto the
# tsne plot. These are all saved as png images in the root directory

plotting_results <- function(parameters, fSOM = NULL, all_params = F, stats = F){
  
  if (!dir.exists("data_output")){
    dir.create("data_output")
  }
  
  path <- "./data_output/"
  
  if (all_params) {
    for (i in names(parameters)){
      if (i %in% c("tsne1", "tsne2", "Filename", "name", "FlowSOM")) {} else {
        ggplot(parameters, aes(tsne1, tsne2, col = parameters[, i])) +
        geom_point(size=0.01) +
        scale_colour_gradientn(colours=matlab.like2(10)) +
        labs(colour = i) +
        theme_bw() + 
        theme(panel.background = element_rect(colour = "white"), 
              panel.grid = element_blank(), legend.position = c(0.95,0.8))
    
        ggsave(paste0(path, "tsne of ", i,".png"), height = 8, width = 10)
      }
    }
  }
  
  if (!is.null(fSOM)) {
  png(filename = paste0(path, "01 FlowSOM spanning tree.png"), 
      width = 10, height = 8, units = "in", res = 300)
  PlotStars(fSOM[[1]],backgroundValues = as.factor(fSOM[[2]]))
  dev.off()
  
  png(filename = paste0(path, "02 FlowSOM matrix.png"), 
      width = 10, height = 8, units = "in", res = 300)
  PlotStars(flow.res[[1]], view = "grid", 
            backgroundValues = as.factor(flow.res[[2]]))
  dev.off()
  }
  #SOMplot
  
  ggplot(parameters, aes(tsne1, tsne2, col = FlowSOM)) +
    geom_point(size=0.01, show.legend = T) +
    guides(colour = guide_legend(override.aes = list(size = 5))) +
    theme_bw() + 
    theme(panel.background = element_rect(colour = "white"), panel.grid = element_blank())
  ggsave(paste0(path, "03 tsne of FlowSOM.png"), height = 8, width = 10)
  
  ggplot(parameters, aes(tsne1, tsne2, col = FlowSOM)) + 
    geom_point(size=0.01, show.legend = T) + 
    facet_wrap(~FlowSOM) +
    guides(colour = guide_legend(override.aes = list(size = 5))) +
    theme_bw() +
    theme(panel.background = element_rect(colour = "white"), panel.grid = element_blank())
  ggsave(paste0(path, "04 FlowSOM individual plots.png"), height = 8, width = 10)

  ggplot(parameters, aes(tsne1, tsne2, col = FlowSOM)) + 
    geom_point(size=0.01, show.legend = T) + 
    facet_wrap(~name) +
    guides(colour = guide_legend(override.aes = list(size = 5))) +
    theme_bw() +
    theme(panel.background = element_rect(colour = "white"), panel.grid = element_blank())
  ggsave(paste0(path, "05 tsne per Filename.png"), height = 8, width = 10)
  
  box_graph <- reshape2::melt(parameters, 
                              id.vars = c("Filename", "tsne1", "tsne2", 
                                          "FlowSOM", "name"))
  
  box_graph <- box_graph %>% group_by(variable) %>%
    mutate(med = median(value))
  
  ggplot(box_graph, aes(x=FlowSOM, y=value, fill = FlowSOM)) +
    geom_hline(aes(yintercept = med, group = variable)) +
    facet_wrap(~variable, scales = "free") +
    geom_boxplot(outlier.size = NULL) +
    guides(fill = F) +
    theme_bw() + 
    theme(panel.background = element_rect(colour = "white"), panel.grid = element_blank())
  ggsave(paste0(path, "06 boxplot of all parameters by clusters.png"), 
         height = 8, width = 10)
  
  if (stats) {
  
    info <- parameters %>% group_by(FlowSOM, name) %>% 
      count()
  
    interesting_stats <- as.data.frame(acast(info, FlowSOM~name, value.var = "n"))
  
    interesting_stats[is.na(interesting_stats)] <- 0
  
    interesting_stats$FlowSOM <- c(1:nrow(interesting_stats))
  
    interesting_perc <- round((interesting_stats[1:(length(names(interesting_stats))-1)] / 
                               sum(interesting_stats[,1]))*100, 2)
  
    interesting_perc$FlowSOM <- c(1:nrow(interesting_perc))
  
    bar_data <- melt(interesting_perc, id.vars = "FlowSOM")
  
    ggplot(bar_data, aes(x = as.factor(FlowSOM), y = value, fill = as.factor(FlowSOM))) + 
      geom_bar(stat="identity") + 
      facet_wrap(~variable, scales = "free_x") + 
      labs(x = "Cluster", y = "Percentage", fill = "Cluster")
  
    ggsave(paste0(path, "07 barchart of percentage populations by clusters.png"), 
         height = 8, width = 10)
  }
  }

#
# min_events ----
#
# Returns a vector with a single number which is the minimum number of events 
# in all the files, I use this so I can workout what I want to subset my data to
#
min_events <- function(filename){                       ### finds the minimum number of events in the FCS files
  events <- NULL
  
  for (i in filename){
    
    fs <- read.FCS(i)
    events <- c(events, dim(fs@exprs)[1])
  }
  return(min(events))
}

#
# merge.png.pdf ----
#
# because there are so many data points in the plots, you need to convert the 
# plots into PNGs before putting them into a PDF, the following function, 
# merge.png.pdf does this. It requires an input of the name of the
# resulting pdf, the png files you wish to add to the pdf and if you want to 
# delete the pngs.
#

merge.png.pdf <- function(pdfFile, pngFiles, deletePngFiles=FALSE) {   
  

    #### Package Install ####
  pngPackageExists <- require ("png")
  
  if ( !pngPackageExists ) {
    install.packages ("png")
    library ("png")
    
  }
  #########################
  
  theImages <- function(x) {        # generates a massive list of all the PNG files
    rasterGrob(readPNG(x))  
  }
  
  thePlots <- lapply(pngFiles, theImages) # read the pngs into one large list
  
  pdf(pdfFile, paper = "a4r", width = 20, height = 10)   # open the pdf with a size of a4 in the landscape orientation
  
  for (i in thePlots){                                   # puts a png on each page of the PDF
    grid.arrange(i, newpage = T)
    
  }
  
  dev.off()
  
  if (deletePngFiles) {
    
    unlink(pngFiles)
  }
  
}

#
# displayClusters ----
#
# This functions output is a graph via ggplot2. 
# Variables - 
# cluster - which cluster you want to look at, cluster 0 is the raw data
# x and y - the parameters you want to look at, ie "FSC-A" "SSC-A" "CD11b"
# filename - which files do you want to look at, default is fileNames which are
#            currently all the FCS files in the working directory
# binned - displays the data either binned or in raw format.
#

displayClusters <- function(variable, cluster = 0, x = "FSC-A", y = "SSC-A", 
                            binned = T) {
  
  
  if (binned == T) {
    bin <- geom_hex(bins = 200)
  } else {
    bin <- geom_point(size = 0.25)
  }
  
  
  tCluster <- NULL
  
  x <- grep(x, names(variable))
  y <- grep(y, names(variable))
  
  
  if (length(cluster) == 1 & cluster[1] == 0) {
    
    ggplot(working, aes(x = variable[x], y = variable[y])) + 
      bin + 
      labs(y = names(variable)[y], x = names(variable)[x]) +
      guides(colour = guide_legend(override.aes = list(size = 5))) +
      theme(panel.background = element_rect(colour = "white"), panel.grid = 
              element_blank()) +
      theme_bw()
    #facet_wrap(~Filename)
    
    # ggsave(paste0(variable[x], " x ", variable[y], " with  no Clusters ", ".png"), 
    #         units = "in")
    
  } else {
    
    
    for (i in cluster){
      tCluster <- rbind(tCluster, variable[variable$FlowSOM == i,])
    }
    
    
    # you need to order the clusters by their filenames and their rownames for 
    # the clusters to be displayed correctly
    Cluster <- tCluster[order(tCluster$Filename, 
                              as.numeric(rownames(tCluster))),, drop = F]  
    
    no_Clust <- paste(unique(Cluster$FlowSOM), collapse = ", ")                                        
    
    ggplot(variable, aes(x = variable[x], y = variable[y])) + 
      bin + 
      geom_point(size = 0.25, alpha = 0.2, data =Cluster, 
                 aes(x = Cluster[x], y = Cluster[y], col = as.factor(FlowSOM))) +
      scale_colour_discrete(name = "Cluster") +
      labs(y = names(variable)[y], x = names(variable)[x]) +
      guides(colour = guide_legend(override.aes = list(size = 5))) +
      theme(panel.background = element_rect(colour = "white"), 
            panel.grid = element_blank()) +
      theme_bw()
    #facet_wrap(~Filename)
    
    #ggsave(paste0(variable[x], " x ", variable[y], " with Clusters ", 
    #              no_Clust, ".png"), units = "in")
  }
}


#### variables ----

fileNames <- dir(pattern=".fcs$")

events <- min_events(fileNames)

CyTOF = T

fs <- read.FCS(fileNames[3])

params <- c(fs@parameters@data$desc, "Filename")
# Gets the data into a dataframe and downsamples, in this case to the minmum 
# number of events in the smallest FCS file
raw <- getData(fileNames, CyTOF = CyTOF, downsample = T, nsamp = events)


names(raw) <- params

names(raw)

# for CyTOF all channels are opened on the machine, this allows you to pick and
# choose which channels you want to use by editing the interested vector.
interested <- c(11:14,16,18,19,23,26,28,36,38,40,41,42,52,62, 63,65,66, 68,69,70,72, 73)   

working <- raw[interested]

# relabel your columns, very dependent on how you have labeled your data
names(working) <- gsub("^.*?_", "", names(working))
names(working)

# transform the data with an arcsin transformation or logicle depending on 

if (CyTOF) {
  working <- asinh(working[, names(working)] / 5)
} else {
  working <- flowFrame(exprs = as.matrix(raw[interested]))
  lgcl <- estimateLogicle(working, channels = names(raw[interested]))
  working <- transform(working, lgcl)
  working <- as.data.frame(working@exprs)
  }

p <- c(1:length(names(working))-1) # this is the columns that you will be used in the working variable
p <- 1:25 # channels that you want to look at

working$Filename <- raw$Filename


### Check the axis of all your data plots

# Using reshape to change the dataframe into something useful 
a <- melt(working, id.vars = c("Filename", "tsne1", "tsne2"))

# Plot the histograms of all your channels
ggplot(a, aes(value)) +
  geom_histogram(bins = 2000) +
  coord_cartesian(ylim = c(0, 100)) +
  theme_bw() + 
  theme(panel.background = element_rect(colour = "white"), panel.grid = element_blank()) +
  facet_wrap(~variable)


# Generating a tsne ----
# generates a tsne plot for the parameters you are interested in, I've turned 
# theta down to 0.25 to increase the similarity to the original tsne. 
# Default is 0.5

# test is a variable you can use to choose which channels to tSNE on.
test <- c(1,3, 7:10, 14, 18:19, 21, 22, 24)

# multicore tSNE
tsne <- Rtsne.multicore(working[p], theta = 0.5, perplexity = 30, 
                        max_iter = 2500, verbose = T, num_threads = 21, check_duplicates = F)

# puts the tsne into the flowData dataframe
working$tsne1 <- tsne$Y[,1]
working$tsne2 <- tsne$Y[,2]

# Phenograph algorithm ----
# runs the Phenograph clustering algorithm, used k =70 as I wanted similarity 
# with perplexity, default is 30
Rpheno <- Rphenograph(working[p], k = 30)

# binds the results into the dataframe
working$Phenograph <- as.factor(membership(Rpheno[[2]]))


# FlowSOM algorithm ----
#
# put the data you are interested into a new data.frame
# give the data.frame columns new names based on the parameters
SOMData <- working[test]
colnames(SOMData) <- params[interested]

# to run the FlowSOM need to convert the dataframe into a flowframe
fSOM <- flowFrame( exprs = as.matrix(SOMData))


set.seed(42)
# running the FlowSOM algorithm, You can change the nClus to maxMeta, and
# the algorithm will try to determine the best number of clusters. Generally
# i Start with maxMeta = 40, then rerun and increment my nClus up or down
# by 2 and examine the MEM output and map the clusters back onto 2D plots
# to verify my clustering.
flow.res <- FlowSOM(fSOM, colsToUse = c(1:length(SOMData)), seed = 42, xdim = 7,
                    ydim = 7, silent = F, maxMeta = 44, scale = F)

# Plot flowSOM graph ----
# generates spanning tree plot with the number of clusters it finds or place
# the clusters into grid display
PlotStars(flow.res[[1]], view = "grid", backgroundValues = as.factor(flow.res[[2]]))
PlotStars(flow.res[[1]],backgroundValues = as.factor(flow.res[[2]]))

# get the cluster information out of the flowSOM, change k to the number of 
# clusters in the plot above
metaClustering <- metaClustering_consensus(flow.res[[1]]$map$codes, k = 12)

metaClustering_perCell <- metaClustering[flow.res$FlowSOM$map$mapping[,1]]

working$FlowSOM <- as.factor(metaClustering_perCell)

#### building MEM plots ----

# you need to choose which channels you want to examine using MEM
fSOM_mem <- working[p]
fSOM_mem$cluster <- as.numeric(working$FlowSOM)

# the cluster column needs to be the very last column and needs to be labelled
# "cluster"
z <- c(getParam(fileNames)[c(-(1:6),-(16:17))], "cluster")

# changes the names of the fSOM_mem data.frame
names(fSOM_mem) <- z

# runs MEM
fSOM_mem <- MEM(fSOM_mem)


# generates a MEM heatmap, with hierarchical cluster on both the flurophores
# and the clusters. You can change the threshold for change ysing display.thresh
# If you are using this in a brand new project, you need to create a directory
# called "output files" in your working directory, otherwise it throws an error
# Readout are 2 PDFs and 2 csv files
build.heatmaps(fSOM_mem, cluster.MEM = "both", cluster.medians = "both", 
               display.thresh = 3, newWindow.heatmaps=F, output.files = T)


# plots the datas as png files
plotting_results(working, flow.res, all_params = T)

# turns the png files into a pdf
pngfiles <- dir(path = "./data_output/", pattern=".png$", full.names = T)
merge.png.pdf(pdfFile = "Data.pdf", pngFiles = pngfiles, deletePngFiles = F)

#### save relevant data to a .Rdata file type

save(raw, working, tsne, flow.res, fSOM_mem, file = "CyTOF Data.Rdata")

#
# write fcs ----
#
# This can be used to write the files into an FCS format.
# you will need to create a "Changed FCS" folder in the working directory.
# the FCS files have very strange scaling when you look at them in flowJo

a <- unique(flowData$Filename)

for (i in a){
  
  test <- flowData[flowData$Filename == i, ]
  
  file <- flowFrame(exprs = as.matrix(test[-18]))
  
  write.FCS(file, filename = (paste0("./Changed FCS/R_modified_", i)))
}

#
# Graphing cluster data onto 2D plots ----
#
# Run the displayClusters function, parameters to pass into function are
# clusters you want to run - can use c(1, 2, 3, 4) to look at multiple
# clusters in one 2D plot, assign your x and y parameter based on their name
# You can add just the files you want to display, it defaults to the global
# variable of fileNames (which is all the fcs files in the working directory)

displayClusters(working, cluster = c(2,5), x = "CD45", y = "CD54", bin = F)

#
# basic readout statistics ----
#
# Generate very basic statistics on each cluster, using number of events in
# cluster / total events to create percentage of Cluster. This then writes into
# a csv file.

test <- working %>% group_by(FlowSOM, name) %>% 
  count()

# acast is useful for changing the shape of the matrix, so it will be FlowSOM
# on the rows and the filenames on the columns. NA values are populated to
# positions you don't have data
interesting_stats <- as.data.frame(acast(test, FlowSOM~name, value.var = "n"))

# if there is an NA, return a 0 value
interesting_stats[is.na(interesting_stats)] <- 0

# adds a FlowSOM column, just so you can track the cluster
interesting_stats$FlowSOM <- c(1:nrow(interesting_stats))

# generate percentages
interesting_perc <- round((interesting_stats[1:(length(names(interesting_stats))-1)] / 
                             sum(interesting_stats[,1]))*100, 2)


interesting_perc$FlowSOM <- c(1:nrow(interesting_perc))

bar_data <- melt(interesting_perc, id.vars = "FlowSOM")

ggplot(bar_data, aes(x = as.factor(FlowSOM), y = value, fill = as.factor(FlowSOM))) + 
  geom_bar(stat="identity") + 
  facet_wrap(~variable, scales = "free_x") + 
  labs(x = "Cluster", y = "Percentage", fill = "Cluster")

write.csv(interesting_stats, "Events per cluster.csv", row.names = F)
write.csv(interesting_perc, "Percentage of population per cluster.csv", row.names = F)


# FIN ----

