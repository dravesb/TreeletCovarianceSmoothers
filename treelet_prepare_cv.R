treelet_prepare_cv = function(grm_name, num_test = 50, snp_set_size = NA){
	
	#--------------------------------
	#Check for errors in the grm file
	#--------------------------------
	
	if(is.na(grm_name)){
		stop("please specify grm file name \n this should be in the current directory")
	}	
	
	#-------------------------------------------
	# read in SNP file 
	#-------------------------------------------
	
	snp_name = paste(grm_name, ".bim", sep = "")
	snps = as.matrix(read.table(snp_name))[,2]
	
	if(!exists("snps")){
		stop("couldn't find the file")
	}
	
	
	if(is.na(snp_set_size)){
		
		#include entire snp set 
		snp_set_size = length(snps)
	}
	
	
	if(length(snps)<snp_set_size){
		message("snp_set_size is greater than the number of snps \n setting it to the max value")
		snp_set_size = length(snps)
	}
	
	snps = sample(snps, snp_set_size)
	
	
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
		write.table(snps[small_snp_set], file, col.names = F, row.names = F, quote = F)
	}
	
	
	#go back a directory
	setwd("..")
	
	
	message("writting out grms - this may take a very long time...")
	
	#--------------------------------
	#Make training grms
	#--------------------------------
	
	#create directory for the new grms
	dir.create("./cv_grms", showWarnings = FALSE)
		
	#format snp file names
	snp_file = paste(getwd(), "/snp_sets/train", sep = "")
	out = paste(getwd(), "/cv_grms/train", sep = "")	

	#make the call to gcta 
	system(paste("./gcta_mac --bfile ",grm_name," --extract ",snp_file," --make-grm --out ", out, sep = ""))
	
	
	#--------------------------------
	#Make testing grms
	#--------------------------------
	for(i in 1:num_test){
		
		#format snp file names
		snp_file = paste(getwd(), "/snp_sets/test",i, sep = "")
		out = paste(getwd(), "/cv_grms/test",i, sep = "")	
		
		#make the call to gcta
		system(paste("./gcta_mac --bfile ",grm_name," --extract ",snp_file," --make-grm --out ",out,sep = ""))
	}

	
	
}