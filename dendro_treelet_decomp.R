dendro_treelet_decomp = function(grm_name = NA, num_clusters = NA){	
	
	#--------------------------------------------------
	# get full treelet object
	#--------------------------------------------------
	message("getting the treelet_decomp object")
	
	tree = treelet_decomp(grm_name)
	
	#--------------------------------------------------
	# unpack treelet_decomp object
	#--------------------------------------------------

	merges = tree$merges 
	treetop = tree$level
	
	#--------------------------------------------------
	# set num_clusters
	#--------------------------------------------------

	if(is.na(num_clusters)){
		message("num_clusters is NA: setting it to single cluster")
		num_clusters = 1 
	}
	if(num_clusters>treetop){
		message("num_clusters greater than n - 1 : seeting it to maximum value (n-1)")
		num_clusters = 1 
	}

	#--------------------------------------------------
	# Get ordering of leaf nodes  
	#--------------------------------------------------
	message("building the plot")

	
	p = nrow(merges)+1
	ordering = rep(NA,p)
	ordering[1:2] = merges[p-1,]
	for(l in 2:(p-1)){
		which.order = which(ordering == merges[p-l,1])
		ordering[(which.order+2):p] = ordering[(which.order+1):(p-1)]
		ordering[(which.order+1)] = merges[(p-l),2]
	}
	
	#--------------------------------------------------
	# Build the hclust object 
	#--------------------------------------------------
	
	merges.kept.orig = merges[(nrow(merges)-(treetop)):nrow(merges),]

	merges.kept = NA
	for(i in 1:length(as.vector(merges.kept.orig))){
		merges.kept[i] = which(sort(unique(as.vector(merges.kept.orig))) == as.vector(merges.kept.orig)[i])	
	}

	merges.kept = matrix(merges.kept,ncol=2)
	id.kept = unique(as.vector(merges.kept.orig))

	order.id = NA
	for(o in 1:length(id.kept)){
		order.id[o] = which(ordering == id.kept[o])
	}

	order.id.orig = id.kept[order(order.id)]
	order.id.kept = rank(order.id.orig)


	merge.new = matrix(NA,nrow=nrow(merges.kept),ncol=2)
	for(i in 1:nrow(merges.kept)){	
	
		alpha.tmp = merges.kept[i,1]
		beta.tmp = merges.kept[i,2]

		n.used.al = which(matrix(merges.kept[1:i,],ncol=2,byrow=F)==alpha.tmp,arr.ind=TRUE)
		n.used.bet = which(matrix(merges.kept[1:i,],ncol=2,byrow=F)==beta.tmp,arr.ind=TRUE)

		if(nrow(n.used.al)==1 & nrow(n.used.bet)==1){
			merge.new[i,] = -1*merges.kept[i,]
		}
	
		if(nrow(n.used.al)>1 & nrow(n.used.bet)==1){
			merge.new[i,2] = -1*merges.kept[i,2]
			merge.new[i,1] = n.used.al[(nrow(n.used.al)-1),1]
		}
		if(nrow(n.used.al)>1 & nrow(n.used.bet)>1){
			merge.new[i,1] = n.used.al[(nrow(n.used.al)-1),1]
			merge.new[i,2] = n.used.bet[(nrow(n.used.bet)-1),1]
		}
		if(nrow(n.used.al)==1 & nrow(n.used.bet)>1){
			merge.new[i,1] = -1*merges.kept[i,1]
			merge.new[i,2] = n.used.bet[(nrow(n.used.bet)-1),1]
		}
	}



	dend.tree = list()
	dend.tree$merge = merge.new
	dend.tree$height = seq(1,nrow(merges.kept),1)
	dend.tree$order = order.id.kept
	dend.tree$labels = as.character(sort(order.id.orig))

	class(dend.tree) = "hclust"


	#--------------------------------------------------
	# Format the plot and make it look pretty
	#--------------------------------------------------

	source("http://addictedtor.free.fr/packages/A2R/lastVersion/R/code.R")
	library(ggdendro)

	title = paste("Clustering Path:", grm_name)

	if(num_clusters < 2){
		ggdendrogram(dend.tree, size = 4, theme_dendro = FALSE) + 
		theme_classic() + 
		labs(x = "", y = "", title = title)+ 
		theme(axis.text.x = element_text(angle = 90, hjust = 1),
			line = element_blank(), 
			axis.text.y = element_blank(), 
			axis.ticks.y = element_blank()) 
		
		}else{
		A2Rplot(dend.tree, k = num_clusters, boxes = FALSE, col.up = "gray50", main = "")
		mtext(title, side=3, adj=0, line=1.2)
		
	}
	
	
}
