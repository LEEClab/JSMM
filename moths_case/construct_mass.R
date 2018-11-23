construct_mass =function(FEMbasis,Agrid,sparse,k,parameter)
{
  nodes = FEMbasis$mesh$nodes
  triangles = FEMbasis$mesh$triangles
  detJ = FEMbasis$detJ
  order = FEMbasis$order
  
  nele  = dim(triangles)[1]
  nnod  = dim(nodes)[1]
  
  if (order < 1 || order > 2){
    stop("ORDER not 1 or 2")
  }
  
  K0M = NULL
  
  if (order ==2)
  {   
    #  the integrals of products of basis functions for master element:
    
    K0M = matrix(c( 6, -1, -1, -4,  0,  0,
                    -1,  6, -1,  0, -4,  0,
                    -1, -1,  6,  0,  0, -4,
                    -4,  0,  0, 32, 16, 16,
                    0, -4,  0, 16, 32, 16,
                    0,  0, -4, 16, 16, 32), ncol=6, nrow=6, byrow=T)/360
  }else if (order == 1)
  {  
    #  the integrals of products of basis functions for master element:
    K0M = matrix( c( 2,  1,  1,
                     1,  2,  1,
                     1 , 1,  2), ncol=3, nrow=3, byrow=T) / 24
  }
  
  # assemble the mass matrix
  
  if (sparse){
    nele_local = ncol(triangles)^2
    K0_triplet = matrix(0,nrow=nele*nele_local,ncol=3)
    
    vec = detJ*k*parameter
    lK = as.vector(as.vector(K0M) %o% vec) 
    K0 = sparseMatrix(i=Agrid[,1], j=Agrid[,2], x=lK, dims=c(nnod,nnod))
    
    # for (el in 1:nele) {  
    #   lK = K0M * detJ[el] * k[el] * parameter[el]
    #   grid = Agrid[nele_local*(el-1)+(1:nele_local),]
    #   vals = cbind(grid, c(lK))
    #   K0_triplet[(el-1)*nele_local+(1:nele_local),] = vals
    # }
    # K0<-sparseMatrix(i=K0_triplet[,1], j=K0_triplet[,2], x = K0_triplet[,3], dims = c(nnod,nnod))
  }
  else
  {
    K0 = matrix(0,nrow=nnod,ncol=nnod)
    for (el in 1:nele)
    {  
      ind = triangles[el,]
      K0[ind,ind] = K0[ind,ind] + K0M * detJ[el] * k[el] * parameter[el]
    }
    K0<-Matrix(K0)
  }
  K0
}
