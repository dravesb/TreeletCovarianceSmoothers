treelet_cv = function(grm_name = NA, num_test = 50, snp_set_size = NA,lev_set = NA, lam_set = seq(0,.1, .01), no_cores = NA){
	
	#-------------------------------------------
	#Prepare the files and testing sets
	#-------------------------------------------
	
	treelet_prepare_cv(grm_name, num_test, snp_set_size)	

	#-------------------------------------------
	#Get the cost 
	#-------------------------------------------
		
	cost = treelet_cost(lev_set = c(500), lam_set, num_test, no_cores)
	
	#-------------------------------------------
	# make the plots 
	#-------------------------------------------
	
	if(dim(cost)[1] == 1){
		#make a line plot

		df = data.frame(X = lam_set, Y = as.vector(cost))
		p = ggplot(df, aes(x = X, y = Y)) + geom_line(col = "orange")+ geom_point(col = "blue", shape = 1) + theme_gdocs() + xlab("Lambda Set") + ylab("Cost")
		
	}else{
		#make a surface plot 
		
		
		
	}
	
	
	
}
