    tbl = read.table("no_smooth.grm.gz", header = F)
    tbl2 = read.table("no_smooth.grm.id", header =F)

	#format GRM
	rel.tmp=as.matrix(read.table("no_smooth.grm.gz",header=F))
	rel1 = matrix(NA,nrow=n,ncol=n)		
	rel1[rel.tmp[,1:2]] = rel.tmp[,4]
	rel1[rel.tmp[,2:1]] = rel.tmp[,4]
	
	#Set tree height
	J = tcb_opt
	
	
	print("computing the Correlation....")
	C = rel1
	cc = cov2cor(C)
	
	#Initialize All variables Here 
	print("initializing variables....")
	
	#Build J Tree Variables
	myCs=list()
    dim_C = dim(C)[1]
    maxlevel = nrow(C)-1
    #J = maxlevel 
    Z = matrix(rep(0, J * 2), ncol = 2)
    T = list()
    theta = rep(0, J)
    PCidx = matrix(rep(0, J * 2), ncol = 2)
    L = 1
    maskno = matrix()
    nodes = seq(1, dim_C, by = 1)
    dlabels = rep(0, dim_C)
    PC_ratio = rep(0, dim_C - 1)
    Zpos = matrix(rep(0, J * 2), ncol = 2)
    all_d = matrix(rep(0, J * dim_C), ncol = dim_C)
    all_nodes = matrix(rep(0, J * dim_C), ncol = dim_C)
    cc.out = rep(NA,maxlevel)
    C.out = rep(NA,maxlevel)
     
    #J Tree Basis Variables 
    m = dim(all_nodes)[2]
    tmpfilts = diag(rep(1, m))
    ind = list()
    sums = matrix(rep(0, m * J), ncol = m)
    difs = matrix(rep(0, m * J), ncol = m)
    myout = list()
    
    
    for (lev in 1:J) {
    	
    	if(lev %%100 == 0 ){
    		print(paste("computing the level",lev, "............"))
    	}
  
    	
    	#Build J Tree Here
    	
        mask_C = upper.tri(cc) * cc
        k = (mask_C == 0)
        mask_C[k] = -1
        mask_C[maskno, ] = -1
        mask_C[, maskno] = -1
        compno = which(mask_C == max(mask_C), arr.ind = TRUE)[1,]
        Cred = C[compno, compno]
		cc.out[lev] = cc[compno,compno][1,2]
		C.out[lev] = C[compno,compno][1,2]
 
        if(Cred[1, 2] == 0){
            Cnew = C
            ccnew = cc
            R = diag(c(1, 1))
            theta = 0
            idx = c(1, 2)
        }else{
            C11 = Cred[1, 1]
            C22 = Cred[2, 2]
            C12 = Cred[1, 2]
            th = 1/2 * atan(2 * C12/(C11 - C22))
            cs = cos(th)
            sn = sin(th)
            R = rbind(c(cs, -sn), c(sn, cs))
            M = C
            M[compno, ] = t(R) %*% C[compno, ]
            C = M
            C[, compno] = M[, compno] %*% R
            Cred = C[compno, compno]
            idx = c(Cred[1, 1], Cred[2, 2])
            idx = sort.list(idx, decreasing = TRUE)
            dnew = diag(C)
            temp = sqrt(matrix(dnew[compno], ncol = 1) %*% dnew)
            temp = C[compno, ]/temp
            cc[compno, ] = temp
            cc[, compno] = t(temp)
        }
 
        #PCidx[lev, ] = idx
        theta[lev] = th
        #T[[lev]] = R
        Z[lev, ] = nodes[compno]
        pind = compno[idx]
        p1 = pind[1]
        p2 = pind[2]
        nodes[pind] = cbind(dim_C + lev, 0)
        
        dlabels[p2] = lev
        if (lev == 1) {
            maskno = p2
            maskno = as.matrix(maskno)
        }else {
            maskno = cbind(maskno, p2)
        }
        PC_ratio[lev] = C[p2, p2]/C[p1, p1]
        #Zpos[lev, ] = compno
        all_d[lev, ] = t(dlabels)
        all_nodes[lev, ] = nodes
		
		
		#J Tree Basis Here
		#s = tmpfilts[Zpos[lev, ], ]
        s = tmpfilts[compno, ]
        #R = T[[lev]]
	    y = t(R)%*%s
        #tmpfilts[Zpos[lev, ], ] = y
        tmpfilts[compno, ] = y
        #y = y[PCidx[lev, ], ]
        y = y[idx, ]
        
        
        sums[lev, ] = y[1, ]
        difs[lev, ] = y[2, ]
	
		#Make A hat
	


		
	
    if(lev == J){
    	
    	train_basis = t(tmpfilts)
		train_C = C 
    	
    	indices = which(nodes != 0)
		if(tcb_opt == 499){
			A = train_basis[,indices]%*%(as.matrix(train_C[indices,indices]*t(train_basis[,indices])))
		}else{
	    	A = train_basis[,indices]%*%train_C[indices,indices]%*%t(train_basis[,indices])
		}
    	
        #Write out results
        newRel = A[upper.tri(A, diag = TRUE)]
        tbl[,4] = newRel        
        #write.table(tbl,"tcb_opt_lev.grm.gz", quote=F, col.names =F, row.names =F)
        #write.table(tbl2,"tcb_opt_lev.grm.id", quote=F, col.names =F, row.names =F)    
    }	
	
	
	
}
