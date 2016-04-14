

setwd("~/Desktop/2016Spring/Phylometh/Final Project Wooliver")
source("makeIndividualTree.R")

myTree <- read.tree("Euc.final.dated.noNit.tre")
myTree$tip.label <- gsub("_"," ", myTree$tip.label)
species <- myTree$tip.label

ala<-read.csv("ALA_Tassie Euc full data.csv") # environmental factors for over 70,000 eucalypt individuals across Tasmania, obtained from Atlas of Living Australia (see descriptions in "ALA and Tenure descriptions.R")

alaDat <- na.omit(data.frame(alaInTree[,c("Species...matched","Temperature...annual.mean..Bio01.")]))


makeIndividualTree(myTree, alaDat)