function [Y,PTM]    = Rouwenhorst(N,shockvar,p,q)
    
% Written by Curtis Aquino, 2019

%##########################################################################
% This function applies the Rouwenhorst method as suggested by Karen
% Kopecky (2010) to discretize a stationary AR(1) process with normally
% distributed errors of the form: zt = rho z_{t-1} + e_t
%##########################################################################    
    
% ********************************************************
% Default is symmetry of p and q
% ********************************************************

if nargin < 5
    q           = p; 
end
            
% ********************************************************
% Recursive formulation
% ********************************************************

Pi = cell(1,N);
for i = 2:N
    if i == 2
        Pi{i}   = [p,1-p;1-p,p];
    else
        Ze      = zeros(i-1,1);
        Pi{i}   = p*[Pi{i-1},Ze;Ze',0]+(1-p)*[Ze,Pi{i-1};0,Ze']+(1-q)*[Ze',0;Pi{i-1},Ze]+q*[0,Ze';Ze,Pi{i-1}];  
        Pi{i}(2:(end-1),:) = Pi{i}(2:(end-1),:)/2;
    end
end
PTM             = Pi{end};
            
% ********************************************************
% Generate the state space
% ********************************************************

Psi             = sqrt(N-1)*sqrt(shockvar);
Y               = fliplr(linspace(-Psi,Psi,N))';
            
end