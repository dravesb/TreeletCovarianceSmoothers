treelet_covariance_smoother = function(tree, lambda = 0){
	
	#-----------------------
	# Unpack treelet object 
	#-----------------------

	G = tree$G 
	V = tree$basis
	n = dim(G)[1]
	sum_indices = tree$sum_indices
	treetop = tree$level	


	#---------------------------------
	# Check for errors in lambda input 
	#---------------------------------	
	
	lambda = sort(lambda)
	if(any(lambda < 0) || any(lambda) >= 1){
		message("lambda can only be in [0,1)\n using only these values")
		lambda = lambda[which(0 <= lambda & lambda < 1)]
	}

	if(length(lambda) == 0){
		message("there were no applicable values of lambda\n using the defult: lamdba = 0")
		lambda = 0 
	}
	#---------------------------------
	# Treshold G based on lambda input 
	#---------------------------------

	G_thres = list()

	for(i in 1:length(lambda)){

		#Threshold matrix
		G[which( G < lambda[i], arr.ind = TRUE)] = 0
		G_thres[[i]] = G 
	}

	#--------------------------------
	# Package Smoothed A - estimates  
	#--------------------------------


	full_A = lapply(G_thres, function(x) V %*% x %*% t(V))

	if(treetop == n-1){
		sum_A = lapply(G_thres, function(x) V[,sum_indices] %*% as.matrix( x[sum_indices, sum_indices] * t(V[,sum_indices])))
	}else{
		sum_A = lapply(G_thres, function(x) V[,sum_indices] %*% x[sum_indices, sum_indices] %*% t(V[,sum_indices]))
	}
	
	#-------------------------------------------
	# Package object
	#-------------------------------------------
	
	treelet_covariance_smoother = function(full, sum,lev, lam){
		#Set the value
		value = list(full_basis_matrices = full, sum_basis_matrices = sum,level = lev, lambda = lam)
		
		#Set the class
		attr(value, "class") = "treelet_covariance_smoother"
		value 
	}	

	obj = treelet_covariance_smoother(full_A, sum_A, treetop, lambda)

	#-------------------------------------------
	# Return the object 
	#-------------------------------------------
	
	return(obj)
	
		

}
