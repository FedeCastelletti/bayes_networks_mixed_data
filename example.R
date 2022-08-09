#################################
## Example with simulated data ##
#################################

q = 6
n = 500

# Generate a dataset with mixed variables:

# X1, X2 are ordinal variables with 5 levels (Binomial)
# X3, X4 are count variables (Poisson)
# X5, X6 are binary variables (Bernoulli)

type.var = "mixed"

source("gen_mixed_data.r")

out_data = gen_mixed_data(i = 4, n = n, q = q, prob = 2/q, type.var = "mixed", ind = NULL)

X   = out_data$X
DAG = out_data$dag

source("mcmc_dag.r")

S    = 2500
burn = 500

# Suppose nodes 1 and 2 are response variables
# no edges are allowe from 1 and 2 to other nodes

response_nodes = 1:2
A_constr = matrix(0, q, q)
A_constr[upper.tri(A_constr)] = NA; A_constr[-response_nodes,] = 0

A_constr

# Run the MCMC

out_mcmc = mcmc_dag(X = X, S = S, burn = burn, a_pi = 1, b_pi = 5, A_constr = A_constr)

# Produce posterior summaries:
# Posterior probabilities of edge inclusion
# Bayesian Model Averaging (BMA) estimate of the covariance matrix Sigma

probs_hat = round(apply(out_mcmc$Graph_post[,,(burn + 1):S], MARGIN = c(1,2), mean), 2)
sigma_hat = round(apply(out_mcmc$Sigma_post[,,(burn + 1):S], MARGIN = c(1,2), mean), 2)

probs_hat
sigma_hat

# Correlation matrix:

cov2cor(out_data$Sigma)

# A DAG estimate (Median Probability DAG model) obtained by including edges whose posterior probability of inclusion is higher than 0.5

est_dag = (probs_hat > 0.5)*1

plot(as(est_dag, "graphNEL"))
