# NumericalMethods

MarkovProcesses

1. Rouwenhorst
2. StationaryDistribution
3. Moments

		Given a probability transition matrix, PTM, and a vector of values for each state, X, Moments() will produce a table output that      displays the mean, variance, covariance, and autocorrelation associated with the PTM. Note that these are not based on simulation, but on the explicit formulas suggested by Paul Klein (2019). 

4. Simulate

		Given a probability transition matrix, PTM, and a number of periods, T, Simulate() will generate a time series of state variables based on the probability transition based with the first 5% of options burned to remove dependence on initial conditions. If no argument for an initial state, S0, is given, the default argument begin simulation from state 1.
