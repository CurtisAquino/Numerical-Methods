function F          = MarkovSimulate(PTM,T,S0)
    
% Written by Curtis Aquino, 2019

%##########################################################################
% Given a probability transition matrix, PTM, and a number of periods, T,
% Simulate() will generate a time series of state variables based on the
% probability transition based with the first 5% of options burned to
% remove dependence on initial conditions. If no argument for an initial
% state, S0, is given, the default argument begin simulation from state 1.
%##########################################################################
    
% ********************************************************
% Default argument
% ********************************************************    
    
if nargin < 3
    S(1)    = 1;
else
    S(1)    = S0;
end
    
% ********************************************************
% Checks whether the PTM is valid
% ********************************************************

if range(size(PTM)) ~= 0
    error('This is not a probability transition matrix')
end

% ********************************************************
% Builds bounds for uniform draws
% ********************************************************
    
TT              = T + round(0.05*T);
for i = 1:size(PTM,1)
    Bounds      = cumsum(PTM(i,:));
    Bd(:,:,i)   = [[0;Bounds(1:(end-1))'],Bounds(:)];
end
UniformDraws    = rand(TT,1);
    
% ********************************************************
% Simulation
% ********************************************************

for i = 2:TT
    S(i)         = find(UniformDraws(i-1)>Bd(:,1,S(i-1)) & UniformDraws(i-1)<Bd(:,2,S(i-1)));
end

% ********************************************************
% Burns 5% to mitigate dependence on initial conditions
% ********************************************************

F = S(round(0.05*T)+1:end)';
    
    
end