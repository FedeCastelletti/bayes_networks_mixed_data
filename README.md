# bayes_networks_mixed_data

These R codes implement our Bayesian method for network structure learning given mixed data.
Specifically:

gen_mixed_data.R : randomly generates a DAG structure on q nodes with a given probability of edge inclusion and generates a dataset containing mixed (binary, ordinal/binomial, count/poisson) data.

mcmc_dag.R : contains the main MCMC algorithm for posterior inference on DAG structures and parameters

move_dag.R    : performs one move from a DAG to an adjacent DAG (implements the proposal distribution over the space of DAGs)

marg_node.R   : computes the (log)marginal likelihood of a DAG model relative to a node-component of the DAG

rDAGWishart.R : samples from a DAG-Wishart (posterior) distribution

example.R : implements mcmc_dag.R on a simulated dataset with mixed data generated using gen_mixed_data.R
