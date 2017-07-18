#Read in grm files
tbl = read.table("no_smooth.grm.gz", header = F)
tbl2 = read.table("no_smooth.grm.id", header =F) 

#Make TCS GRM 
list_tree = plot_matrix(500,"no_smooth")

basis = list_tree[[1]]
C = list_tree[[2]]

#Treshold matrix
tmp_C = C
tmp_C[which(abs(tmp_C)<tcs_opt)] = 0 

#Compute A 
A = basis%*%tmp_C%*%t(basis)

#Write out results
newRel = A[upper.tri(A, diag = TRUE)]
tbl[,4] = newRel		
write.table(tbl,"tcs_opt_lam.grm.gz", quote=F, col.names =F, row.names =F)
write.table(tbl2,"tcs_opt_lam.grm.id", quote=F, col.names =F, row.names =F)