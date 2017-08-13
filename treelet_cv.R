treelet_cv = function(grm_name = NA, num_test = 50, snp_set_size = NA,lev_set = NA, lam_set = seq(0,.1, .01), no_cores = NA, remove = T){
	
	#-------------------------------------------
	#Prepare the files and testing sets
	#-------------------------------------------
	
	treelet_prepare_cv(grm_name, num_test, snp_set_size)	

	#-------------------------------------------
	#Get the cost 
	#-------------------------------------------
		
	cost = treelet_cost(lev_set, lam_set, num_test, no_cores)
	
	#-------------------------------------------
	# make the plots 
	#-------------------------------------------
	
	if(dim(cost)[1] == 1){
		#make a line plot

		df = data.frame(X = lam_set, Y = as.vector(cost))
		p = ggplot(df, aes(x = X, y = Y)) + geom_line(col = "orange")+ geom_point(col = "blue", shape = 1) + theme_gdocs() + xlab("Lambda Set") + ylab("Cost")
		
	}else{
		#make a surface plot 
		p = plot_ly(x = lam_set, y = lev_set, z = cost) %>% add_surface() %>% layout(title = paste(grm_name, "Cost"), scene = list(xaxis = list(title = "Lambda Set"), yaxis = list(title = "Level Set"), zaxis = list(title = "Cost")))
		
	}
	
	#-------------------------------------------
	# Remove the cost grms/snp_sets
	#-------------------------------------------
	
	if(remove = TRUE){
		unlink("./cv_grms", recursive = TRUE)
		unlink("./snp_sets", recursive = TRUE)	
	}
	
	
	
	
	
	
	#-------------------------------------------
	# Package object
	#-------------------------------------------
	
	cost_obj = list(p = p, cost = cost)
	class(cost_obj) = "Cost"
	
	#return the object
	return(cost_obj) 
		
}
