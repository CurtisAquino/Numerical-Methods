function F          = StationaryDistribution(PTM,maxIter)
    
% Written by Curtis Aquino, 2019

%##########################################################################
% This function takes a probability transition matrix, PTM, and iterates on
% an equispaced initial guess until machine epsilon convergence. If no
% maximum number of iterations, maxIter, is specified,
% StationaryDistribution() will conclude that no stationary distribution
% exists if norm(pi_{2000}-pi_{1999}) > 10^(-16)
%##########################################################################  

% ********************************************************
% Default argument
% ********************************************************    
    
if nargin < 2
    maxIter = 2000;
end
    
% ********************************************************
% Ensures that input is a PTM
% ********************************************************

if range(size(PTM)) ~= 0
    error('This probability transition matrix is not a square')
end

% *********************************************************
% Finds the stationary distribution depending on data class
% *********************************************************

type                = class(PTM);
switch type 
    case 'table'
        pi0         = ones(1,size(PTM,1))/size(PTM,1);
        pi1         = Inf;
        itr         = 1;
        thr         = Inf;
        while thr > eps
            pi1     = pi0 * PTM{:,:};
            pi0     = pi1;
            itr     = itr + 1;
            thr     = norm(pi1-pi0);
            if itr > maxIter
                error('No stationary distribution exists')
            end
        end
    case 'double'
        pi0         = ones(1,size(PTM,1))/size(PTM,1);
        pi1         = Inf;
        itr         = 1;
        thr         = Inf;
        while thr > eps
            pi1     = pi0 * PTM(:,:);
            thr     = norm(pi1-pi0);
            pi0     = pi1;
            itr     = itr + 1;
            if itr > maxIter
                error('No stationary distribution exists')
            end
        end
    otherwise
        error('Incorrect input type')
end

% ********************************************************
% Results
% ********************************************************

F                   = pi0';

end