construct_int = function(FEMbasis,k){
  nodes = FEMbasis$mesh$nodes
  triangles = FEMbasis$mesh$triangles

  nele  = nrow(triangles)
  nnod  = nrow(nodes)
  int0 = rep(0,nnod)
  detJ = FEMbasis$detJ #twice the area of each triangle
  for (i in 1:nele){
     int0[triangles[i,]] <- int0[triangles[i,]] + k[i]*detJ[i]/6
   }
 int0
}
