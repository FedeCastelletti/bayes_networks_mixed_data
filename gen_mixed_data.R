#################################################
## Function for DAG and data random generation ##
#################################################

gen_mixed_data = function(i, n, q, prob, type.var, ind){
  
  library(gRbase)
  library(mvtnorm)
  
  # i    : seed to be set
  # n    : number of observations
  # q    : number of variables
  # prob : probability of edge inclusion for the true DAG
  # K    : number of levels for each categorical variable
  
  # type.var : type of variables to be generated in (gaussian, binary, count, mixed)
  
  # ind : (q,1) of indicators for variable types: "g" Gaussian, "b" binary, "c" count
  
  if(type.var == "mixed" & is.null(ind)){
    
    ind = rep("NA", q); ind[1:round(q/3)] = "o"; ind[(round(q/3)+1):round(2*q/3)] = "c"; ind[ind == "NA"] = "b"
    
  }
  
  set.seed(i)
  
  # Function to generate a topologically-ordered DAG
  
  rDAG = function(q, w){
    DAG = matrix(0, q, q); colnames(DAG) = rownames(DAG) = 1:q
    DAG[lower.tri(DAG)] = rbinom(n = q*(q-1)/2, size = 1, prob = w)
    return(DAG)
  }
  
  A = rDAG(q = q, w = prob)
  
  B = A*matrix(runif(q*q, 0.1, 1), q, q)*sample(c(-1, 1), size = q*q, replace = TRUE); diag(B) = 1
  
  Sigma_cond = diag(rep(1, q))
  
  Sigma = solve(t(B))%*%Sigma_cond%*%solve(B); mu = c(rep(0, q))
  
  library(mvtnorm)
  
  # Generate the latent (Gaussian) data
  
  Z = data.frame(rmvnorm(n, mu, Sigma))
  Z = t((t(Z) - apply(Z, 2, mean))/apply(Z, 2, sd))
  
  X_gaussian = Z
  
  # Discretize the latent and obtain the dataset X_binary, X_count, X_ordinal
  
  set.seed(i)
  
  probs = runif(q, 0.2, 0.8)
  
  X_binary = sapply(1:q, function(j) qbinom(p = pnorm(Z[,j]), prob = probs[j], size = 1))
  
  set.seed(i)
  
  lambda = runif(q, 5, 10)
  
  X_count = sapply(1:q, function(j) qpois(p = pnorm(Z[,j]), lambda = lambda[j]))
  
  set.seed(i)
  
  probs = runif(q, 0.2, 0.8)
  
  ns = rep(5, q)
  
  X_ordinal = sapply(1:q, function(j) qbinom(p = pnorm(Z[,j]), prob = probs[j], size = ns[j]))
  
  
  if(type.var == "binary"){
    
    return(list(dag = A, X = X_binary, types = rep("b", q), Sigma = Sigma))
    
  }
  
  if(type.var == "count"){
    
    return(list(dag = A, X = X_count, types = rep("c", q), Sigma = Sigma))
    
  }
  
  if(type.var == "ordinal"){
    
    return(list(dag = A, X = X_ordinal, types = ind, Sigma = Sigma))
    
  }
  
  if(type.var == "mixed"){
    
    X_mixed = matrix(NA, n, q)
    
    X_mixed[, ind == "b"] = X_binary[, ind == "b"]
    X_mixed[, ind == "c"] = X_count[, ind == "c"]
    X_mixed[, ind == "o"] = X_ordinal[, ind == "o"]
    
    
    return(list(dag = A, X = X_mixed, types = ind, Sigma = Sigma))
    
  }
  
}

# Example

# gen_mixed_data(i = 1, n = 20, q = 10, prob = 0.2, type.var = "gaussian", ind = NULL)
# gen_mixed_data(i = 1, n = 20, q = 10, prob = 0.2, type.var = "binary", ind = NULL)
# gen_mixed_data(i = 1, n = 20, q = 10, prob = 0.2, type.var = "count", ind = NULL)
# gen_mixed_data(i = 1, n = 20, q = 10, prob = 0.2, type.var = "mixed", ind = NULL)

# ind = c("b", "b", "b", "c", "c", "g", "g", "b", "c", "b")

# gen_mixed_data(i = 1, n = 20, q = 10, prob = 0.2, type.var = "mixed", ind = ind)
