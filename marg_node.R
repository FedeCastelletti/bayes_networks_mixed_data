###########################
## Preliminary functions ##
###########################

# To find parents and family of a node in a DAG

pa = function(node, DAG) {
    pa = which(DAG[,node] != 0)
    return(pa)
}

fa = function(node, DAG) {
  pa = which(DAG[,node] != 0)
  fa = c(node, pa)
  return(fa)
}

###################
## Main function ##
###################

marg_j = function(node, DAG, tXX, n, a, U){
  
  ###########
  ## Input ##
  ###########
  
  # node : numerical label of the node
  # DAG  : (q,q) adjacency matrix of the DAG
  # tXX  : (q,q) design matrix t(X)%*%X with X the (n,q) data matrix
  # n    : number of observations (rows) in the data matrix X
  # a    : shape hyperparameter of the DAG Wishart prior
  # U    : position hyperparameter of the DAG Wishart prior
  
  ############
  ## Output ##
  ############
  
  # marg_node : logarithm of the marginal likelihood of the node
  
  j  = node
  pa = pa(j, DAG)
  q  = ncol(tXX)
  
  aj = a + length(pa) - q + 1
  
  Upost = U + tXX
  
  if(length(pa) == 0){
    
    U_jj = U[j,j]
    Upost_jj = Upost[j,j]
    
    prior.normcost = -lgamma(aj/2) + aj/2*log(U_jj/2)
    post.normcost  = -lgamma(aj/2 + n/2) + (aj/2 + n/2)*log(Upost_jj/2)
    
  }else{
    
    U_paj.j = U[pa,j]
    U_jj    = U[j,j] - t(U_paj.j)%*%solve(U[pa,pa])%*%U_paj.j
    Upost_paj.j = Upost[pa,j]
    Upost_jj    = Upost[j,j] - t(Upost_paj.j)%*%solve(Upost[pa,pa])%*%Upost_paj.j
     
    prior.normcost = -lgamma(aj/2) + aj/2*log(U_jj/2) + 0.5*log(det(as.matrix(U[pa,pa])))
    post.normcost  = -lgamma(aj/2 + n/2) + (aj/2 + n/2)*log(Upost_jj/2) + 0.5*log(det(as.matrix(Upost[pa,pa])))
    
  }
  
  marg_node = -n/2*log(2*pi) + prior.normcost - post.normcost
  
  return(marg_node)
}
