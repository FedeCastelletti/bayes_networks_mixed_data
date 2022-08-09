####################
## Main functions ##
####################

# Sample from Multivariate-Normal-Inverse-Gamma (prior/posterior) distribution relative to node j in DAG
# i.e. DAG-Wishart distribution relative to node j

# It corresponds to the distribution of the parameters of a Normal linear regression model
# Response variable is "node" and covariates are given by the parents of node in DAG

sample_DL_j = function(node, DAG, aj, U){
  
  ###########
  ## Input ##
  ###########
  
  # node : numerical label of the node in the DAG
  # DAG  : (q,q) adjacency matrix of the DAG
  # aj   : node-shape hyperparameter of the DAG-Wishart
  # U    : position hyperparameter of the DAG-Wishart, a (q, q) s.p.d. matrix
  
  ############
  ## Output ##
  ############
  
  # out : a list with two elements; a vector with one draw for the (vector) regression coefficient and a scalar with one draw for the conditional variance
  
  j  = node
  pa = pa(j, DAG)

  out = list(Djj = 0, Lj = 0)

  if(length(pa) == 0){
    U_jj = U[j,j]
    out$Djj = rgamma(1, shape = aj/2, rate = U_jj/2)^-1
  }else{
    U_paj.j   = U[pa,j]
    invU_papa = solve(U[pa,pa])
    U_jj      = U[j,j] - t(U_paj.j)%*%invU_papa%*%U_paj.j

    out$Djj = rgamma(1, shape = aj/2, rate = U_jj/2)^-1
    out$Lj  = mvtnorm::rmvnorm(1, -invU_papa%*%U_paj.j, out$Djj*invU_papa)
  }

  return(out)
}

# Random sample from a DAG-Wishart distribution

rDAGWishart = function(DAG, a, U){
  
  ###########
  ## Input ##
  ###########
  
  # node : numerical label of the node in the DAG
  # DAG : (q,q) adjacency matrix of the DAG
  # a   : common shape hyperparameter of the DAG-Wishart (a > q - 1)
  
  ############
  ## Output ##
  ############
  
  # D_out : sampled (q,q) diagonal matrix collecting node-conditional variances
  # L_out : sampled (q,q) matrix collecting regression coefficients
  
  q   = ncol(DAG)
  ajs = sapply(1:q, function(j) a + sum(DAG[,j]==1) - q + 1)

  L_out = matrix(0, q, q)
  D_out = matrix(0, q, q)

  params = lapply(1:q, function(j) sample_DL_j(j, DAG, ajs[j], U))
  sigmas = sapply(1:q, function(x) params[[x]]$Djj)

  L = lapply(1:q, function(x) params[[x]]$Lj)

  D_out = diag(sigmas)

    for(j in 1:q){
      whc = which(DAG[,j] == 1)
      L_out[whc,j] = as.numeric(L[[j]])
    }

  diag(L_out) = 1

  return(list(D_out = D_out, L_out = L_out))

}
