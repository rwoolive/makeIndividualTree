##### final project: More than a few studies looking at variation in traits (functional traits, habitats) across species across species have failed to incorporate within-species variation in their analyses. Rather, these studies use only means of individual data to represent species values (Smith & Beaulieu 2009). However, not incorporating within species variation in these situations can lead to increased type I error rates of phylogenetic independent contrasts (up to 17%; Harmon & Losos 2005) and have considerable effect on parameter estimation (reviewed by Garamszegi 20--). Moreover, within-species variation can sometimes be greater than among-species variation, and is often biologically meaningful. Thus, it is important to incorporate within-species variation into phylogenetic comparative analyses. Although there are methods that can incorporate either standard errors or variances around species' means (Ives et al. 2007; Hansen and Bartoszek 2012), there is no method available (in R) to allow the incorporation of each individual-level datum into species-level phylogenetic analyses. I have developed a method to create individual-level phylogenies to do so, by adding tips to an inputted species-level phylogeny to represent each indivudual. This approach can be used in any phylogenetic comparative method.


# this function requires a phylogeny of class 'phylo' and a dataset, in which the first column contains species' names that match the phylogeny and the second column contains trait values

makeIndividualTree <- function(phy,data){
	# first, prune the dataset to species that are in the tree
	install.packages("ape"); library(ape)
	install.packages("geiger"); library(geiger)
	data <- data[data[,1] %in% phy$tip.label,] 
	speciesInData <- unique(data[,1])
	# then, prune the tree to species that are in the dataset
	toKeep <- which(phy$tip.label %in% speciesInData) 
	phy <- drop.tip(phy, seq(1,length(speciesInData))[-toKeep])
	originalN <- length(phy$tip.label)
	originalSpecies <- phy$tip.label
		
	# determine how many individuals of each species are in the dataset (ultimately how many tips to add to the original tree + 1):
	numberInds <- data.frame(number=table(data[,1]))
	# order this data frame by the species in the phylogeny:
	rownames(numberInds) <- names(table(data[,1]))
	numberInds <- numberInds[phy$tip.label,]
	
	# make star phylogenies for individuals of each species - 1:
	indStarPhys <- list()
	minLength <- min(phy$edge.length)/2

	for(s in 1:dim(numberInds)[1]){
		# make edge matrix:
		n <- (numberInds[s,2])
		edge<-matrix(NA,n-1,2)
		edge[,1]<-n; edge[,2]<-1:(n-1)
		# and assign branch lengths
		edge.length= rep(minLength, (n-1))
		tip.labels <- NA
		reps <- seq(2,n)
		for(p in 1:(n-1)){
			tip.labels[p] <- paste(phy$tip.label[s], "_rep", reps[p], sep="")
		}

		indStarPhys [[s]]<-list(edge=edge,edge.length=edge.length,Nnode=1,tip.label=as.vector(tip.labels))
		class(indStarPhys [[s]])<-"phylo"
		phy <- bind.tree(phy, indStarPhys[[s]], s, minLength)
	}
		phy$tip.label[c(1:originalN)] <- paste(phy$tip.label[c(1:originalN)], "_rep1", sep="")
		
		# reorder dataset according to the order of species on the phylogeny:
		numberInds <- numberInds[order(rownames(numberInds)),]
		dataOrdered <- data[order(data[,1]),]
		dataReps <- NA
		for(i in 1:dim(numberInds)[1]){
			dataReps <- c(dataReps, seq(1, numberInds[i,2]))
		}
		dataReps <- dataReps[-1]
		rownames(dataOrdered) <- paste(as.character(dataOrdered[,1]), "_rep", dataReps, sep="")
		treeDat <- treedata(phy, dataOrdered, sort=TRUE)
		return(treeDat)
}





