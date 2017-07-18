treelet_decomposition = function(grm_name = NA,treetop = NA){
	
	#--------------------------------
	#Check for erros in the grm file
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
	
	#-------------------------------------------
	# Initial correlation values
	#-------------------------------------------
	
	print("Computing the correlation....")
	M = cov2cor(A)
	
	#-------------------------------------------
	# Initial looping variables
	#-------------------------------------------
	
	print("initializing looping variables....")
	merges = matrix(NA, ncol = 2, nrow = treetop)
	clusters = 1:n
	merged_indices = matrix()
	pc_ratio = numeric(treetop)
	clusters_by_iteration = matrix(NA, nrow = treetop, ncol = n )
	basis = diag(rep(1,n))
	
	
	#-------------------------------------------
	# Begin iterative procedure 
	#-------------------------------------------
	for(lev in 1:treetop){
		
		#print out progress 
		if(lev %% 100 == 0 ){
    		print(paste("Computing level:",lev))
    	}
		
		#mask_M <- matrix to find hightest pairwise relatedness
		mask_M = upper.tri(M) * M
		
		#set non - applicable entries to -1  
		mask_M[lower.tri(mask_M, diag = TRUE)] = mask_M[merged_indices, ] = mask_M[, merged_indices] = -1 
		
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
		merges[lev, ] = alpha_beta
		
		#order the PC indicies & store the cluster names
		princ_component_index = alpha_beta[index]
		clusters[alpha_beta] = cbind(n + lev, 0) 
		
		#update merged_indicies
		merged_indices = ifelse(lev == 1, as.matrix(princ_component_index[2]), cbind(merged_indices,princ_component_index[2]))
		
		#record PC ratio
		pc_ratio[lev] =  A[princ_component_index[2], princ_component_index[2]]/A[princ_component_index[1], princ_component_index[1]]
		
		#record clusters total 
		clusters_by_iteration[lev,] = clusters
		
		#update basis set
		tmp1 = basis[alpha_beta,]
		tmp2 = t(rotation) %*% tmp1 
		basis[alpha_beta,] = tmp2
		
	}
	
	#-------------------------------------------
	# Repackage smoothed relationship matrices
	#-------------------------------------------
	
	V = t(basis)
	G = as.matrix(A)
	sum_indices = which(clusters != 0)

	full_A = V %*% G %*% t(V) 
	sum_A = V[,sum_indices] %*% G[sum_indices, sum_indices] %*% t(V[,sum_indices])
	
	#-------------------------------------------
	# Package object
	#-------------------------------------------
		
	
	
	
	
	
}



