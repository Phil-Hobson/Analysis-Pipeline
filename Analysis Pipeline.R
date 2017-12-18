# R version 3.4.2 (2017-09-28)
# Platform: x86_64-w64-mingw32/x64 (64-bit)
# Running under: Windows 7 x64 (build 7601) Service Pack 1

# load all the librarys
library(flowCore)         
library(Rtsne.multicore)
library(ggplot2)
library(Rphenograph)
library(FlowSOM)
library(colorRamps)
library(MEM)
library(gplots)
library(dplyr)
library(grid)
library(gridExtra)
library(png)


set.seed(42)


## functions

#
# getParam ----
#
# generates a vector that returns a all of the parameters

getParam <- function(filename){
  
  fs <- read.FCS(filename[1])
  
  gr <- grep("P(\\dS$|\\d\\dS$)", names(keyword(fs)))                        # finds the P*S variable for parameter names
  
  a <- unlist(fs@description[gr])[order(names(unlist(fs@description[gr])))]  # creates an unsorted vector of these names
  
  pos <- grep("P1S", names(a))                                               # finds P1S in the list
  
  a <- c(a[pos:(pos+8)], a[1:(pos-1)])                                       # based on P1S position sorts the names

  a                                                                          # returns the parameter names
}

# 
# getData ----
#
# Function to read the FCS files and concatenate them into one large data.frame.
# Can use this to downsample the number of events for running a quick tsne. 
# 3 variables need to be entered, the files you want to read, True or False 
# for if you want to downsample, and the number you want to downsample to.
#

getData <- function(filenames, downsample = F, nsamp = NULL){
  
  combined_data <- NULL
  
  for (i in 1:length(filenames)){
  
    fs <- read.FCS(filenames[[i]], package="flowCore", name.keyword="TUBE NAME", ## reads the FCS
                 phenoData=list(name="TUBE NAME", Filename="$FIL"))
  
    comp_matrix <- keyword(fs)$SPILL                                            ## generates the compensation matrix
  
    fs <- compensate(fs, comp_matrix)                                           ## applies the compensation matrix

    fs <- transform(fs, transformList(colnames(fs)[interested], logicleTransform() )) ## just looking at parameters to transform
                                                                                      ## will be specific for each experiment
  
    test <- NULL
    test <- as.data.frame(rbind(test, fs@exprs))                                ## create a dataframe containing the data
    test$Filename <- rep(fs@description$FILENAME, times = dim(fs)[1])           ## generate the filename for each entry in a new column called Filename
    
    if (downsample == T){
      test <- test[sample(nrow(test), nsamp),]                                  ## if you want to downsample, this randomly
      test <- test[order(as.numeric(rownames(test))),, drop = F]
    }                                                                           ## takes the number of rows in nsamp
    
    combined_data <- rbind(combined_data, test)                                 ## creates the final table that will be returned
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

plotting_results <- function(parameters){
  
  param <- getParam(fileNames[1])
    
  for (i in parameters){
    
    ggplot(flowData, aes(tsne1, tsne2, col = flowData[, i])) +
    geom_point(size=0.01) +
    facet_wrap(~Filename) +
    scale_colour_gradientn(colours=matlab.like2(10)) +
    labs(colour = param[i]) +
    theme_bw() + 
    theme(panel.background = element_rect(colour = "white"), panel.grid = element_blank())
    
    ggsave(paste0("tsne of ",param[i],".png"), units = "in")
  }
  
  png(filename = "01 FlowSOM spanning tree.png", width = 8, height = 6, units = "in", res = 300)
  PlotStars(flow.res[[1]],backgroundValues = as.factor(flow.res[[2]]))
  dev.off()
  #SOMplot
  
  
  ggplot(flowData, aes(tsne1, tsne2, col = FlowSOM)) + 
    geom_point(size=0.01) + 
    facet_wrap(~Filename) +
    scale_colour_gradientn(colours=primary.colors(max(flowData$FlowSOM))) + 
    theme_bw() +
    theme(panel.background = element_rect(colour = "white"), panel.grid = element_blank())
  
  ggsave("02 FlowSOM on tSNE.png", units = "in")
  
  ggplot(flowData, aes(tsne1, tsne2, col = FlowSOM)) + geom_point(size=0.01) + facet_wrap(~FlowSOM) +
    scale_colour_gradientn(colours=primary.colors(max(flowData$FlowSOM))) + theme_bw() +
    theme(panel.background = element_rect(colour = "white"), panel.grid = element_blank())
  ggsave("03 FlowSOM individual plots.png", units = "in")

}

#
# min_events ----
#
# Returns a vector with a single number which is the minimum number of events 
# in all the files, I use this so I can workout what I want to subset my data to
#
min_events <- function(){                       ### finds the minimum number of events in the FCS files
  events <- NULL
  
  for (i in fileNames){
    
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

displayClusters <- function(cluster = 0, x = "FSC-A", y = "SSC-A", 
                            filename = fileNames, binned = T) {
  
  
  if (binned == T) {
    bin <- geom_hex(bins = 200)
  } else {
    bin <- geom_point(size = 0.25)
  }
  
  
  tCluster <- NULL
  params <- getParam(filename[1])
  
  x <- grep(x, params)
  y <- grep(y, params)
  
  
  if (length(cluster) == 1 & cluster[1] == 0) {
    
    ggplot(flowData, aes(x = flowData[x], y = flowData[y])) + 
      bin + 
      labs(y = params[y], x = params[x]) +
      guides(colour = guide_legend(override.aes = list(size = 5))) +
      theme(panel.background = element_rect(colour = "white"), panel.grid = 
              element_blank()) +
      theme_bw() +
      facet_wrap(~Filename)
    
    ggsave(paste0(params[x], " x ", params[y], " with  no Clusters ", ".png"), 
           units = "in")
  
    } else {
    
  
    for (i in cluster){
      tCluster <- rbind(tCluster, flowData[flowData$FlowSOM == i,])
    }
  

    # you need to order the clusters by their filenames and their rownames for 
    # the clusters to be displayed correctly
    Cluster <- tCluster[order(tCluster$Filename, 
                              as.numeric(rownames(tCluster))),, drop = F]  
                                                                                          
    no_Clust <- paste(unique(Cluster$FlowSOM), collapse = ", ")                                        
  
  ggplot(flowData, aes(x = flowData[x], y = flowData[y])) + 
    bin + 
    geom_point(size = 0.25, alpha = 0.2, data =Cluster, 
               aes(x = Cluster[x], y = Cluster[y], col = as.factor(FlowSOM))) +
    scale_colour_discrete(name = "Cluster") +
    labs(y = params[y], x = params[x]) +
    guides(colour = guide_legend(override.aes = list(size = 5))) +
    theme(panel.background = element_rect(colour = "white"), 
          panel.grid = element_blank()) +
    theme_bw() +
    facet_wrap(~Filename)
  
  ggsave(paste0(params[x], " x ", params[y], " with Clusters ", 
                no_Clust, ".png"), units = "in")
  }
}


#### variables ----

fileNames <- dir(pattern=".fcs$")

fs <- read.FCS(fileNames[[1]], package="flowCore", name.keyword="TUBE NAME",
               phenoData=list(name="TUBE NAME", Filename="$FIL"))

params <- getParam(fileNames[1])
names(params) <- NULL

# these are the colours you are interested in based on the parameters 
# variable positions
interested <- 7:15      


# Gets the data into a dataframe and downsamples, in this case to the minmum 
# number of events in the smallest FCS file
flowData <- getData(fileNames, T, min_events())

# Generating a tsne ----
# generates a tsne plot for the parameters you are interested in, I've turned 
# theta down to 0.25 to increase the similarity to the original tsne. 
# Default is 0.5
tsne <- Rtsne.multicore(flowData[interested], theta = 0.5, perplexity = 70, 
                        max_iter = 5000, verbose = T, num_threads = 21)

# puts the tsne into the flowData dataframe
flowData <- cbind(flowData, tsne1 = tsne$Y[,1], tsne2 = tsne$Y[,2])

# Phenograph algorithm ----
# runs the Phenograph clustering algorithm, used k =70 as I wanted similarity 
# with perplexity, default is 30
Rpheno <- Rphenograph(flowData[interested], k = 30)

# binds the results into the dataframe
flowData$Rphenograph <- as.vector(membership(Rpheno[[2]]))


# FlowSOM algorithm ----
#
# put the data you are interested into a new data.frame
# give the data.frame columns new names based on the parameters
SOMData <- flowData[interested]
colnames(SOMData) <- params[interested]

# to run the FlowSOM need to convert the dataframe into a flowframe
test <- flowFrame( exprs = as.matrix(SOMData))


set.seed(42)
# running the FlowSOM algorithm, You can change the nClus to maxMeta, and
# the algorithm will try to determine the best number of clusters. Generally
# i Start with maxMeta = 40, then rerun and increment my nClus up or down
# by 2 and examine the MEM output and map the clusters back onto 2D plots
# to verify my clustering.
flow.res <- FlowSOM(test, colsToUse = c(1:length(SOMData)), seed = 42, xdim = 7,
                    ydim = 7, silent = F, maxMeta = 44, scale = T)

# Plot flowSOM graph ----
# generates spanning tree plot with the number of clusters it finds
PlotStars(flow.res[[1]],backgroundValues = as.factor(flow.res[[2]]))

# get the cluster information out of the flowSOM, change k to the number of 
# clusters in the plot above
metaClustering <- metaClustering_consensus(flow.res[[1]]$map$codes, k = 13)

metaClustering_perCell <- metaClustering[flow.res$FlowSOM$map$mapping[,1]]

flowData$FlowSOM <- as.vector(metaClustering_perCell)

#### building MEM plots ----

# you need to choose which channels you want to examine using MEM
fSOM_mem <- flowData[,c(-(1:6),-(16:20))]

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
build.heatmaps(fSOM_mem, cluster.MEM = "both", cluster.medians = "none", 
               display.thresh = 3, newWindow.heatmaps=F, output.files = T)


# plots the datas as png files
plotting_results(interested)

# turns the png files into a pdf
pngfiles <- dir(pattern=".png$")
merge.png.pdf(pdfFile = "Data.pdf", pngFiles = pngfiles, deletePngFiles = T)

#### save relevant data to a .Rdata file type

save(flowData, tsne, flow.res, fSOM_mem, file = "13 Cluster Data.Rdata")

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

displayClusters(c(2), x = "SiglecH", y = "CD11b")

#
# basic readout statistics ----
#
# Generate very basic statistics on each cluster, using number of events in
# cluster / total events to create percentage of Cluster. This then writes into
# a csv file.

stats <- data.frame(flowData %>% 
  group_by(FlowSOM, Filename) %>%
  summarise("Cluster" = n())) 

# note this is based on using the min_events(), but will only work if that is
# what you downsampled too. Otherwise you will need to change that value for
# this code to work.
stats$Total <- rep(min_events(), time = dim(stats)[1])

stats$"% Cluster" <- round((stats$Cluster / stats$Total) * 100, 2)

write.csv(stats, "Basic Statiscs.csv")


# FIN ----
