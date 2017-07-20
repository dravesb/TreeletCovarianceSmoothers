plot_treelet_decomp = function(grm_name = NA, tree_cut = NA){	
	
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
	rel = tree$rel_value
	
	
	#--------------------------------------------------
	# set tree_cut
	#--------------------------------------------------

	if(is.na(tree_cut)){
		message("tree_cut is NA: seeting it to maximum value (n-1)")
		tree_cut = nrow(merges) 
	}
	if(tree_cut>treetop){
		message("tree_cut greater than n - 1 : seeting it to maximum value (n-1)")
		tree_cut = nrow(merges) 
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

	rel.kept = rel[(nrow(merges)-(treetop)):nrow(merges),1]

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


	d = 1-rel.kept
	ep = .0001
	for(i in 2:length(d)){
		if(d[i]<=d[i-1]){
			d[i] = d[i-1]+ep
		} else{ d[i] = d[i] }
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

	if(treetop - tree_cut + 1 < 2){
		ggdendrogram(dend.tree, size = 4, theme_dendro = FALSE) + 
		theme_classic() + 
		labs(x = "", y = "", title = title)+ 
		theme(axis.text.x = element_text(angle = 90, hjust = 1),
			line = element_blank(), 
			axis.text.y = element_blank(), 
			axis.ticks.y = element_blank()) 
		
		}else{
		A2Rplot(dend.tree, k = treetop - tree_cut + 1, boxes = FALSE, col.up = "gray50", main = "")
		mtext(title, side=3, adj=0, line=1.2)
		
	}
	
	
}
