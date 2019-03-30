1. MarkovSimulate

	Given a probability transition matrix, PTM, and a number of periods, T, MarkovSimulate() will generate a time series of state variables based on the probability transition based with the first 5% of options burned to remove dependence on initial conditions. If no argument for an initial state, S0, is given, the default argument begin simulation from state 1.

2. StationaryDistribution

	This function takes a probability transition matrix, PTM, and iterates on an equispaced initial guess until machine epsilon convergence. If no maximum number of iterations, maxIter, is specified, StationaryDistribution() will conclude that no stationary distribution exists if $norm(pi_{2000}-pi_{1999}) > 10^(-16)$

3. Rouwenhorst

	This function applies the Rouwenhorst method as suggested by Karen Kopecky (2010) to discretize a stationary AR(1) process with normally distributed errors of the form: $z_t = rho z_{t-1} + e_t$

4. MarkovMoments

	Given a probability transition matrix, PTM, and a vector of values for each state, X, MarkovMoments() will produce a table output that displays the mean, variance, covariance, and autocorrelation associated with the PTM. Note that these are not based on simulation, but on the explicit formulas suggested by Paul Klein (2019). 

  
