treelet_decomp = function(grm_name = NA,treetop = NA){
	
	#--------------------------------
	#Check for errors in the grm file
	#--------------------------------
	
	if(is.na(grm_name)){
		message("please specify grm file name \n this should be in the current directory")
		stop()
	}	

	#-------------------------------------------
	#Format the GRM file -> relationship matrix  
	#-------------------------------------------
	
	file = paste(grm_name, ".grm.gz", sep = "")
	stacked_data = as.matrix(read.table(file), header = F)
	n = max(stacked_data[,1])
		
	A = matrix(NA, nrow = n, ncol = n)
	A[stacked_data[,1:2]] = stacked_data[,4]
	A[stacked_data[,2:1]] = stacked_data[,4]

	#----------------------------------
	#Check errors in treetop 
	#----------------------------------

	if(is.na(treetop)){
		message("setting the tree level to the maximum possible value (n - 1)")
		treetop = n - 1
	}
	
	if(treetop > n-1){
		message("treetop is greater than the dimension of the relationship matrix. \n setting the tree level to the maximum possible value (n - 1)")
		treetop = n - 1
	}
		

	#-------------------------------------------
	# Initial correlation values
	#-------------------------------------------
	
	message("computing the correlation....")
	M = cov2cor(A)
	
	#-------------------------------------------
	# Initial looping variables
	#-------------------------------------------
	
	message("initializing looping variables....")
	clusters = 1:n
	merged_indices = matrix()
	pc_ratio = numeric(treetop)
	clusters_by_iteration = matrix(NA, nrow = treetop, ncol = n )
	basis = diag(rep(1,n))
	merges = matrix(NA, ncol = 2, nrow = treetop)	
		
	#-------------------------------------------
	# Begin iterative procedure 
	#-------------------------------------------
	for(lev in 1:treetop){
		
		#print out progress 
		if(lev %% 100 == 0 ){
    			message(paste("Computing level:",lev))
    		}
		
		#mask_M <- matrix to find hightest pairwise relatedness
		mask_M = M
		
		#set non - applicable entries to -1  
		#mask_M[which(mask_M == 0, arr.ind = TRUE)] = mask_M[merged_indices, ] = mask_M[, merged_indices] = -1 
		mask_M[upper.tri(M, diag = TRUE)] = mask_M[merged_indices, ] = mask_M[, merged_indices] = -1 

		#find highest related pair 
		alpha_beta = which(mask_M == max(mask_M), arr.ind = TRUE)[1,]
		
		#get relatedness submatrix 
		sub_A = A[alpha_beta, alpha_beta]	
		

		#if they're unrelated, rotation angle = 0 
		#so no need to do local PCA  
		if(sub_A[1,2] == 0){
			
			#set up new variables
			rotated_A = A 
			updated_M = M
			rotation = diag(c(1,1))
			theta = 0 
			index = c(1,2) 
			
		#if they are related, find rotation angle 
		#and complete local PCA 	
		}else{
			#get the angle 
			theta = 1/2 * atan(2*sub_A[1,2]/(sub_A[1,1] - sub_A[2,2]))
			
			#find rotation matrix 
			cs = cos(theta)
			sn = sin(theta)
			rotation = rbind(c(cs, -sn), c(sn, cs)) 
			
			#rotate 
			tempA = A 
			tempA[alpha_beta,] = t(rotation) %*% A[alpha_beta,]
			A = tempA 
			A[,alpha_beta] = tempA[,alpha_beta] %*% rotation
			
			#replace and update variables after rotation
			sub_A = A[alpha_beta, alpha_beta]
			index = c(sub_A[1,1], sub_A[2,2])
			index = sort.list(index, decreasing = TRUE)
			
			new_diag = diag(A)
			temp = sqrt(matrix(new_diag[alpha_beta], ncol = 1) %*% new_diag)
			temp = A[alpha_beta,]/temp
			
			#update similarity matrix 
			M[alpha_beta,] = temp
			M[,alpha_beta] = t(temp)
				
		}
		
		#record which individuals were merged 
		merges[lev, ] = alpha_beta[index]
		
		#order the PC indicies & store the cluster names
		princ_component_index = alpha_beta[index]
		clusters[princ_component_index] = cbind(n + lev, -(n+lev)) 
		
		#update merged_indicies
		if(lev == 1){
			merged_indices = as.matrix(princ_component_index[2])
		} else{
			merged_indices = cbind(merged_indices,princ_component_index[2])
		}

		
		#record PC ratio
		pc_ratio[lev] =  A[princ_component_index[2], princ_component_index[2]]/A[princ_component_index[1], princ_component_index[1]]
		
		#record clusters total 
		clusters_by_iteration[lev,] = clusters
		
		#update basis set
		tmp1 = basis[alpha_beta,]
		tmp2 = t(rotation) %*% tmp1 
		basis[alpha_beta,] = tmp2
		
	}
	
	
	#--------------------------------------------------
	# Package basis set and projected covariance matrix
	#--------------------------------------------------


	message("finishing up...")
	
	V = t(basis)
	G = as.matrix(A)
	sum_indices = which(clusters>0)
	
	
	#-------------------------------------------
	# Package object
	#-------------------------------------------
	
	treelet_decomp = function(base, gamma, s_indices, merge, cp, treetop){
		#Set the value
		value = list(basis = base, G = gamma, sum_indices = s_indices, merges = merge, cluster_paths = cp, level = treetop)
		
		#Set the class
		attr(value, "class") = "treelet_decomp"
		value 
	}	
	
	obj = treelet_decomp(V, G, sum_indices, merges, clusters_by_iteration, treetop)

	#-------------------------------------------
	# Return the object 
	#-------------------------------------------
	
	return(obj)
	
}



