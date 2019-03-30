function F          = MarkovMoments(PTM,X)
    
% Written by Curtis Aquino, 2019
% Requires: StationaryDistribution()

%##########################################################################
% Given a probability transition matrix, PTM, and a vector of values for
% each state, X, Moments() will produce a table output that displays the
% mean, variance, covariance, and autocorrelation associated with the PTM.
% Note that these are not based on simulation, but on the explicit formulas
% suggested by Paul Klein (2019). 
%##########################################################################

% ********************************************************
% Fixes user input
% ********************************************************

if isrow(X)
    X = X';
end
if istable(PTM)
   PTM = PTM{:,:}; 
end

% ********************************************************
% Derives the stationary distribution, if it exists
% ********************************************************
    
StDs            = StationaryDistribution(PTM);

% ********************************************************
% Computes moments following Klein (2019)
% ********************************************************

Mean            = sum(StDs.*X);
Variance        = sum(((X-Mean).^2).*StDs);
Covariance      = sum(StDs'.*sum(X(:)'*X.*PTM,2))-Mean^2;
Autocorrelation = Covariance/Variance;

% ********************************************************
% Produces output as a table
% ********************************************************

F               = array2table([Mean,Variance],'VariableNames',{'Mean','Variance'});

end