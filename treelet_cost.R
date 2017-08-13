treelet_cost = function(lev_set = NA, lam_set = seq(0,.1, .01), num_test = 50, no_cores = NA){
	
	#-------------------------------------------
	#Format the GRM file -> relationship matrix  
	#-------------------------------------------
	
	file = paste(getwd(), "/cv_grms/train.grm.gz", sep = "")
	stacked_data = as.matrix(read.table(file), header = F)
	n = max(stacked_data[,1])
	
	
	A = matrix(NA, nrow = n, ncol = n)
	A[stacked_data[,1:2]] = stacked_data[,4]
	A[stacked_data[,2:1]] = stacked_data[,4]

	#-------------------------------------------
	#Format lev set  
	#-------------------------------------------
	
	if(is.na(lev_set)){
		message("setting lev_set to max height")
		lev_set = c(n-1)
	}

	#set treetop
	treetop = max(lev_set)

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
	
	#--------------------------------
	# Set number of cores to use
	#--------------------------------
	
	if(is.na(no_cores)){
		message("number of cores not specified \n setting to no_cores - 1")
		no_cores = detectCores() - 1
	}	
	
	if(no_cores > detectCores()){
		message("number of cores specified greater than total number of cores \n setting to no_cores - 1")
		no_cores = detectCores() - 1
	}
	
	#set up a cluster
	cl = makeCluster(no_cores)
	registerDoParallel(cl)
	
	#-----------------------------------
	# Read in testing sets
	#-----------------------------------
	
	#object to hold testing matrices
	testers = list()

	for(i in 1:num_test){
		
		#Get the filename 		
		file = paste(getwd(),"/cv_grms/test", i, ".grm.gz", sep = "")
		
		#Read in and format as matrix 
		temp = as.matrix(read.table(file, header = F))	
		A_l = matrix(NA,n,n)
		A_l[temp[,1:2]] = A_l[temp[,2:1]] = temp[,4]
		
		#add to the list
		testers[[i]] = A_l
		
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
	# Set cost matrix
	#-------------------------------------------
		
	cost = matrix(NA, nrow = length(lev_set), ncol = length(lam_set))
	rownames(cost) = as.character(lev_set)
	colnames(cost) = as.character(lam_set)	
		
		
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
		
		
		#-----------------------------------------------
		# Get cost if we're in the lev_set
		#-----------------------------------------------
		
		if(lev %in% lev_set){
			
			#format the basis sets/sum_indices
			V = t(basis)
			G = as.matrix(A)
			sum_indices = which(clusters>0)
			
			#loop over lam_set
			for(lam in lam_set){
				
				#thresholding...
				G_thres = G
				G_thres[which(G_thres < lam)] = 0 
				
				
				#get the smoothed matrix 
				if(lev == n-1){
					A_til = V[,sum_indices]%*%(as.matrix(G_thres[sum_indices,sum_indices]*t(V[,sum_indices])))
				}else{
	    			A_til = V[,sum_indices]%*%G_thres[sum_indices,sum_indices]%*%t(G_thres[,sum_indices])
				}
				
				#multiprocessing to get cost at this level
				result = foreach(i=1:num_test, .combine = "+")%dopar%{
					
					#get cost 
					test = testers[[i]]
					cost.here = sum(abs(G[lower.tri(G)]) * (test[lower.tri(test)] - A_til[lower.tri(A_til)])^2)
					
					#return it here
					return(cost.here)
				}
				
				#fill the cost matrix 
				cost[which(as.numeric(rownames(cost)) == lev), which(as.numeric(colnames(cost)) == lam)] = result
				
			}
			
			
		}
		
		
	}
	
	#close the cluster
	stopCluster(cl)
	

	#-----------------------------------------------
	# return cost
	#-----------------------------------------------
		
	return(cost)
	
	
}