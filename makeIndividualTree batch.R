##### final project: More than a few studies looking at variation in traits (functional traits, habitats) across species across species have failed to incorporate within-species variation in their analyses. Rather, these studies use only means of individual data to represent species values (Smith & Beaulieu 2009). However, not incorporating within species variation in these situations can lead to increased type I error rates of phylogenetic independent contrasts (up to 17%; Harmon & Losos 2005) and have considerable effect on parameter estimation (reviewed by Garamszegi 20--). Moreover, within-species variation can sometimes be greater than among-species variation, and is often biologically meaningful. Thus, it is important to incorporate within-species variation into phylogenetic comparative analyses. Although there are methods that can incorporate either standard errors or variances around species' means (Ives et al. 2007; Hansen and Bartoszek 2012), there is no method available (in R) to allow the incorporation of each individual-level datum into species-level phylogenetic analyses. I have developed a method to create individual-level phylogenies to do so, by adding tips to an input species-level phylogeny to represent each indivudual in an input dataset with multiple individuals per species. This approach can be used in any phylogenetic comparative method.


# This function is called 'makeIndividualTree' and requires a phylogeny of class phylo and a dataset, in which the first column contains species' names (that match the phylogeny) and the second column that contains trait values. The dataset can (but doesn't have to) include additional columns of data.






# Example 1: test whether this method yields OUwie estimates that differ from species-level tree plus standard error 

install.packages("OUwie"); library(OUwie)
set.seed(2)

nSpp <- 20
repSpp <- 3
exTree <- rtree(nSpp)
exData <- data.frame(species=rep(exTree$tip.label, each=repSpp), trait1=rep(c(1,2,3,4), each=15), trait2=rnorm(nSpp*repSpp, 10, 2))

exIndTree <- makeIndividualTree(exTree, exData)
par(mfrow=c(1,2))
plot.phylo(exTree, no.margin=TRUE) # species-level tree
plot.phylo(exIndTree$phy, cex=0.7, no.margin=TRUE) # individual-level tree



# individual-level OUwie analysis:

ouwieDat <- data.frame(species=rownames(exIndTree$data), regime=as.numeric(exIndTree$data[,2]), continuousTrait=as.numeric(exIndTree$data[,3]))
ppTree <- rayDISC(exIndTree$phy, ouwieDat, model="ER", node.states="marginal")

indBM1est <- OUwie(ppTree$phy, ouwieDat, model="BM1")
indBMSest <- OUwie(ppTree$phy, ouwieDat, model="BMS")
indOU1est <- OUwie(ppTree$phy, ouwieDat, model="OU1")
indOUMest <- OUwie(ppTree$phy, ouwieDat, model="OUM")
indOUMVest <- OUwie(ppTree$phy, ouwieDat, model="OUMV")
indOUMAest <- OUwie(ppTree$phy, ouwieDat, model="OUMA")
indOUMVAest <- OUwie(ppTree$phy, ouwieDat, model="OUMVA")
mods <- c("BM1","BMS","OU1","OUM","OUMV","OUMA","OUMVA")
fits <- round(c(indBM1est$AICc, indBMSest$AICc, indOU1est$AICc, indOUMest$AICc, indOUMVest$AICc, indOUMVest$AICc, indOUMVAest$AICc),2); fits[7] <- NA 
# AICc for OUMVA model is wacky; we'll exclude it
deltaAICc <- round(fits-min(na.omit(fits)),2)
cbind(mods, fits, deltaAICc) 
print(paste("For the individual-level method, OU1 is the best model. The alpha estimate is ", round(indOU1est$solution[1,1], 2), ", the sigma estimate is ", round(indOU1est$solution[2,1], 2), ", and the theta estimate is ", round(indOU1est$theta[1,1], 2), " with a standard error of ", round(indOU1est$theta[1,2], 2), sep=""))

# species-level OUwie analysis:

	spMeans <- aggregate(trait2~species, data= exData, FUN="mean")
	rownames(spMeans) <- spMeans[,1]; spMeans <- spMeans[exTree$tip.label,]
	
	spSD <- aggregate(trait2~species, data= exData, FUN="sd")
	rownames(spSD) <- spSD[,1]; spSD <- spSD[exTree$tip.label,]
	
	spL <- aggregate(trait2~species, data= exData, FUN="length") 
	rownames(spL) <- spL[,1]; spL <- spL[exTree$tip.label,]

ouwieDat2 <- data.frame(species=exTree$tip.label, regime=rep(c(1,2,3,4), each=5), continuousTrait= spMeans[,2], se= spSD[,2]/sqrt(spL[,2]))
ppTree2 <- rayDISC(exTree, ouwieDat2, model="ER", node.states="marginal")
spOUMest <- OUwie(ppTree2$phy, ouwieDat2, model="OUM")

spBM1est <- OUwie(ppTree2 $phy, ouwieDat2, model="BM1", mserr="known")
spBMSest <- OUwie(ppTree2 $phy, ouwieDat2, model="BMS", mserr="known")
spOU1est <- OUwie(ppTree2 $phy, ouwieDat2, model="OU1", mserr="known")
spOUMest <- OUwie(ppTree2 $phy, ouwieDat2, model="OUM", mserr="known")
spOUMVest <- OUwie(ppTree2 $phy, ouwieDat2, model="OUMV", mserr="known")
spOUMAest <- OUwie(ppTree2 $phy, ouwieDat2, model="OUMA", mserr="known")
spOUMVAest <- OUwie(ppTree2 $phy, ouwieDat2, model="OUMVA", mserr="known")
fits2 <- round(c(spBM1est$AICc, spBMSest$AICc, spOU1est$AICc, spOUMest$AICc, spOUMVest$AICc, spOUMVest$AICc, spOUMVAest$AICc),2)
deltaAICc2 <- round(fits2-min(na.omit(fits2)),2)
cbind(mods, fits2, deltaAICc2) 
print(paste("For the species-level method, BM1 is the best model. The sigma estimate is ", round(spBM1est $solution[2,1], 2), ", and the theta estimate is ", round(spBM1est $theta[1,1], 2), " with a standard error of ", round(spBM1est $theta[1,2], 2), sep=""))

# The methods favor different OUwie models and estimates of evolutionary rates (sigma). In this example we see that sigma is estimated to be orders of magnitude higher with the individual-level method compared to the species-level method. This may be because it allows more room for variation within species, with individuals stemming from shorter branch lengths compared to branch lengths farther back in the phylogeny. Thus, we may be underestimating evolutionary rates when using species-level phylogenies in phylogenetic comparative analyses.






# Example 2: test the function using a phylogeny of Eucalyptus species that are native to Tasmania, Australia, and a dataset (derived from the Atlas of Living Australia) of environmental variables for over 70,000 naturally occuring individuals.

setwd("~/Desktop/2016Spring/Phylometh/Final Project Wooliver")
source("makeIndividualTree.R")

myTree <- read.tree("Euc.final.dated.noNit.tre") # phylogeny of 28 Eucalyptus species
myTree$tip.label <- gsub("_"," ", myTree$tip.label)

ala<-read.csv("ALA_Tassie Euc full data.csv") # raw data
alaDat <- na.omit(data.frame(ala[,c("Species...matched","Temperature...annual.mean..Bio01.", "Moisture.Index...annual.mean..Bio28.")])) # create the input dataset


myIndTree <- makeIndividualTree(myTree, alaDat)










