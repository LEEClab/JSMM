construct_stiff = function(FEMbasis,Agrid,sparse,k,parameter)
{
  nodes = FEMbasis$mesh$nodes
  triangles = FEMbasis$mesh$triangles
  
  nele  = dim(triangles)[[1]]
  nnod  = dim(nodes)[[1]]
  detJ     = FEMbasis$detJ
  order = FEMbasis$order
  metric = FEMbasis$metric

  KXX = NULL
  KXY = NULL
  KYY = NULL
  
  if (order < 1 || order > 2){
    stop("ORDER not 1 or 2")
  }
  
  if (order == 2)
  {   
    #  values of K1 for master elements
    KXX = matrix( c( 3, 1, 0, 0, 0,-4, 
                     1, 3, 0, 0, 0,-4,
                     0, 0, 0, 0, 0, 0,
                     0, 0, 0, 8,-8, 0,
                     0, 0, 0,-8, 8, 0,
                     -4,-4, 0, 0, 0, 8), ncol=6, nrow=6, byrow=T)/6
    
    KXY = matrix(c( 3, 0, 1, 0,-4, 0, 
                    1, 0,-1, 4, 0,-4, 
                    0, 0, 0, 0, 0, 0,
                    0, 0, 4, 4,-4,-4,
                    0, 0,-4,-4, 4, 4,
                    -4, 0, 0,-4, 4, 4), ncol=6, nrow=6, byrow=T)/6
    
    KYY = matrix( c( 3, 0, 1, 0,-4, 0,
                     0, 0, 0, 0, 0, 0,
                     1, 0, 3, 0,-4, 0,
                     0, 0, 0, 8, 0,-8,
                     -4, 0,-4, 0, 8, 0,
                     0, 0, 0,-8, 0, 8), ncol=6, nrow=6, byrow=T)/6
  }else if (order == 1)
  {   
    KXX = matrix( c(  1, -1,  0,
                      -1,  1,  0,
                      0,  0,  0), ncol=3, nrow=3, byrow=T) /2
    
    KXY = matrix( c(  1,  0, -1,
                      -1,  0,  1,
                      0,  0,  0), ncol=3, nrow=3, byrow=T) /2
    
    KYY = matrix( c(  1,  0, -1,
                      0,  0,  0,
                      -1,  0,  1), ncol=3, nrow=3, byrow=T) /2
  }
  
  #  assemble the stiffness matrix
  
  if (sparse){
     
     mat = metric*(as.vector(detJ)*k*parameter)
     vals = as.vector(aperm(mat[,1,1]%o%KXX + mat[,1,2]%o%KXY + mat[,2,1]%o%t(KXY) + mat[,2,2]%o%KYY,c(2,3,1)))
     K1 = sparseMatrix(i=Agrid[,1], j=Agrid[,2], x=vals, dims = c(nnod,nnod))

    # nele_local = ncol(triangles)^2
    # K1_triplet = matrix(0,nrow=nele*nele_local,ncol=3)
    # for (el in 1:nele)
    # {  
    #   K1M = (metric[el,1,1]*KXX    + metric[el,1,2]*KXY +
    #            metric[el,2,1]*t(KXY) + metric[el,2,2]*KYY)
    #   lK1<-K1M*detJ[el] * k[el] * parameter[el]
    #   grid <- Agrid[(nele_local*(el-1)+1):(nele_local*(el)),]
    #   vals<-cbind(grid, c(lK1))
    #   K1_triplet[((el-1)*nele_local+1):(el*nele_local),] = vals
    # }
    # K1<-sparseMatrix(i = K1_triplet[,1], j=K1_triplet[,2], x = K1_triplet[,3], dims = c(nnod,nnod))
  }
  else
  {
    K1   = matrix(0,nrow=nnod,ncol=nnod)
    for (el in 1:nele)
    {
      ind    = triangles[el,]
      K1M = (metric[el,1,1]*KXX    + metric[el,1,2]*KXY +
               metric[el,2,1]*t(KXY) + metric[el,2,2]*KYY)
      K1[ind,ind] = K1[ind,ind] + K1M*detJ[el] * k[el] * parameter[el]
    }
    K1<-Matrix(K1)
  }
  K1
}
