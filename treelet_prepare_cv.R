treelet_prep_cv = funtion(grm_name, no_cores = NA, num_test = 50){
	
	#--------------------------------
	#Check for errors in the grm file
	#--------------------------------
	
	if(is.na(grm_name)){
		stop("please specify grm file name \n this should be in the current directory")
	}	
	
	#--------------------------------
	# Set number of cores to use
	#--------------------------------
	
	if(is.na(no_cores)){
		message("number of cores not specified \n setting to no_cores - 1")
		no_cores = detectCores() - 1
	}	
	
	#-------------------------------------------
	# read in SNP file 
	#-------------------------------------------
	
	snp_name = paste(grm_name, ".bim.gz", sep = "")
	snps = as.matrix(read.table(snp_name))[,2]
	
	if(!exists("snps")){
		stop("couldn't find the file")
	}
	
	
	#--------------------------------
	#Make training and testing sets
	#--------------------------------
	message("writting out testing/training sets")
	
	
	#get train/test indices
	indices = 1:length(snps)
	train_index = sample(indices, round(length(snps)/2))
	test_index = indices[-train_index]
	
	#create directory for snp sets
	dir.create("./snp_sets", showWarnings = FALSE)
	setwd("./snp_sets")
	
	#write out training snps
	write.table(snps[train_index], "train", col.names = F, row.names = F, quote = F)
	
	#split the testing sets once more
	size = floor(length(test_index)/50)
	
	for(i in 1:num_test){
		
		#get the smaller snp_set
		small_snp_set = sample(test_index, size)
		
		#remove small_snp_set from test_index 
		test_index = test_index[!test_index %in% small_snp_set]
		
		#write out small_snp_set
		file = paste("test",i, sep = "")
		write.table(small_snp_set, file, col.names = F, row.names = F, quote = F)
	}
	
	
	#go back a directory
	setwd("..")
	
	#--------------------------------
	#Make training grms
	#--------------------------------
	
	#create directory for the new grms
	dir.create("./cv_grms", showWarnings = FALSE)

	#make the call to gcta 

	result = foreach(i=1:num_test)%dopar%{
		
		#format snp file names3
		snp_file = paste(getwd(), "/snp_sets/test",i, sep = "")
		out = paste(getwd(), "/cv_grms/train", sep = "")	
	}

	
	
	
	system(paste("./gcta_mac --bfile ",grm_name," --extract ",snp_file," --make-grm --out "sep = ""))
	
	
	#--------------------------------
	#Make testing grms
	#--------------------------------

	
	
}