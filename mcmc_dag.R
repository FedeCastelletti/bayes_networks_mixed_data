######################
## Main MCMC scheme ##
######################

mcmc_dag = function(X, S, burn, a_pi, b_pi, A_constr, output = "full"){
  
  ###########
  ## Input ##
  ###########
  
  # X        : the (n,q) data matrix
  # S        : the (integer) number of MCMC iterations
  # burn     : the (integer) burn-in period
  # a_pi     : first shape hyperparameter of the Beta prior for the probability of edge inclusion
  # b_pi     : second shape hyperparameter of the Beta prior for the probability of edge inclusion
  # A_constr : (q,q) matrix representing the constraints imposed to the DAG space
  
  # Obs : A_constr[u,v] = NA means edge u -> v is not allowed
  
  ############
  ## Output ##
  ############
  
  # If output = "full" :
  
  # X          : the input (n,q) data matrix
  # Graph_post : an array collecting the (S - burn) (q,q) adjacency matrices of the DAGs visited by the MCMC
  # Sigma_post : an array collecting the (S - burn) (q,q) covariance matrices sampled by the MCMC
  
  # If output = "bma" :
  
  # X         : the input (n,q) data matrix
  # Probs_hat : a (q,q) matrix collecting posterior probabilities of edge inclusion for each edge (u,v) in the DAG space
  # Sigma_hat : a (q,q) Bayesian Model Averaging (BMA) estimate of the covariance matrix
  
  
  if(is.null(a_pi)){a_pi = 1}
  if(is.null(b_pi)){b_pi = 1}
  
  source("move_dag.r")
  source("marg_node.r")
  source("rDAGWishart.r")
  
  library(mvtnorm)
  library(gRbase)
  
  q = ncol(X)
  n = nrow(X)
  
  Sigma_post = array(NA, c(q, q, S))
  Graph_post = array(NA, c(q, q, S))
  
  # Initialize the chain
  
  Graph = matrix(0, q, q)
  
  a = q
  g = 1/n
  U = diag(g, q)
  
  set.seed(4)
  
  Z = qnorm(apply(X, 2, rank, ties.method = "random")/(n + 1))
  Zfill = matrix(rnorm(n*q), n, q)
  
  Z[is.na(X)] = Zfill[is.na(X)]
  
  # Initialize Sigma

  tZZ = t(Z)%*%(Z)
  
  Upost = U + tZZ
  
  DL_inits = rDAGWishart(DAG = Graph, a + n, U = Upost)
  
  D = DL_inits$D_out
  L = DL_inits$L_out
  
  Sigma = solve(t(L))%*%D%*%solve(L);
  
  Sigma_post[,,1] = Sigma
  Graph_post[,,1] = Graph
  
  R = NULL
  
  for(j in 1:q){
    
    R = cbind(R, match(X[,j], sort(unique(X[,j]))))
    
  }
  
  set = 1:q
  
  Z = t((t(Z) - apply(Z, 2, mean))/apply(Z, 2, sd))
  
  cat("MCMC sampling")
  pb = utils::txtProgressBar(min = 2, max = S, style = 3)
  
  for(s in 1:S){
    
    # update Z[,j]
    
    for(j in sample(set)){
      
      Beta_0 = L; diag(Beta_0) = 0
      
      cond_mean = -as.matrix(Z[,-j])%*%Beta_0[-j,j]
      
      for(r in sort(unique(R[,j]))){
        
        ir = (1:n)[R[,j] == r & !is.na(R[,j])]
        
        lb = suppressWarnings(max(Z[R[,j] < r, j], na.rm = T))
        ub = suppressWarnings(min(Z[R[,j] > r, j], na.rm = T))
        
        Z[ir,j] = qnorm(runif(length(ir), pnorm(lb, cond_mean[ir], sqrt(D[j,j])), 
                              pnorm(ub, cond_mean[ir], sqrt(D[j,j]))), cond_mean[ir], sqrt(D[j,j]))
        
      }
      
      # to handle NAs
      
      ir      = (1:n)[is.na(R[,j])]
      Z[ir,j] = rnorm(length(ir), cond_mean[ir], sqrt(D[j,j]))
      
    }
    
    Z = t((t(Z))/apply(Z, 2, sd))
    
    tZZ = t(Z)%*%(Z)
    
    ## Update the graph given the latents
    
    Graph_move = move(A = Graph, q = q, A_constr = A_constr)
    
    Graph_prop = Graph_move$A_new
    nodes_prop = Graph_move$nodes
    
    type.operator = Graph_move$type.operator
    
    
    # Multiplicity correction (log)prior ratio
    
    logprior.new = lgamma(n.edge(Graph_prop) + a_pi) + 
      lgamma(q*(q-1)/2 - n.edge(Graph_prop) + b_pi - 1)
    
    logprior.old = lgamma(n.edge(Graph) + a_pi) + 
      lgamma(q*(q-1)/2 - n.edge(Graph) + b_pi - 1)
    
    logprior = logprior.new - logprior.old
    
    
    # Marginal likelihood ratio (distinguish 3 cases):
    
    if(type.operator == 1){
      
      # (1) Insert a directed edge
      
      marg_star = marg_j(node = nodes_prop[2], DAG = Graph_prop, tXX = tZZ, n = n, a = a, U = U)
      marg      = marg_j(node = nodes_prop[2], DAG = Graph, tXX = tZZ, n = n, a = a, U = U)
      
    }else{
      
      if(type.operator == 2){
        
        # (2) Delete a directed edge
        
        marg_star = marg_j(node = nodes_prop[2], DAG = Graph_prop, tXX = tZZ, n = n, a = a, U = U)
        marg      = marg_j(node = nodes_prop[2], DAG = Graph, tXX = tZZ, n = n, a = a, U = U)
        
      }else{
        
        # (3) Reverse a directed edge
        
        marg_star = marg_j(node = nodes_prop[1], DAG = Graph_prop, tXX = tZZ, n = n, a = a, U = U) +
          marg_j(node = nodes_prop[2], DAG = Graph_prop, tXX = tZZ, n = n, a = a, U = U)
        
        marg = marg_j(node = nodes_prop[1], DAG = Graph, tXX = tZZ, n = n, a = a, U = U) +
          marg_j(node = nodes_prop[2], DAG = Graph, tXX = tZZ, n = n, a = a, U = U)
        
      }
      
    }
    
    # acceptance ratio
    
    ratio_D = min(0, marg_star - marg + logprior)
    
    # accept move
    
    if(log(runif(1)) < ratio_D){
      
      Graph = Graph_prop
      
    }
    
    Graph_post[,,s] = Graph
    
    
    ## Sample Sigma
    
    Upost = U + tZZ
    
    DL = rDAGWishart(DAG = Graph, a + n, U = Upost)
    
    L  = DL$L_out
    D  = DL$D_out
    
    Sigma = solve(t(L))%*%D%*%solve(L);
    
    Sigma_post[,,s] = Sigma
    
    
    utils::setTxtProgressBar(pb, s)
    close(pb)
    
  }
  
  probs_hat = round(apply(Graph_post[,,(burn + 1):S], MARGIN = c(1,2), mean), 2)
  sigma_hat = round(apply(Sigma_post[,,(burn + 1):S], MARGIN = c(1,2), mean), 2)
  
  if(output == "full"){
    
    return(out = list(X = X,
                      Graph_post = Graph_post,
                      Sigma_post = Sigma_post))
    
  }
  
  if(output == "bma"){
    
    return(out = list(X = X,
                      Probs_hat = probs_hat,
                      Sigma_hat = sigma_hat))
    
    
  }
  
}