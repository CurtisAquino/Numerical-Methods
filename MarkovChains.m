classdef MarkovChains
    methods(Static)
        %%
function [Pi,Y] = Rouwenhorst(N,shockvar,p,q)
            % Given zt = rho z_{t-1} + e_t
            
            %N = 4; shockvar = 1; p = (1+0.975)/2;

            % Symmetry is default
            if nargin < 4; q = p; end
            
            % Recursion
            Pi = cell(1,N);
            for i = 2:N
                if i == 2
                    Pi{i}   = [p,1-p;1-p,p];
                else
                    Ze = zeros(i-1,1);
                    Pi{i} = p*[Pi{i-1},Ze;Ze',0]+(1-p)*[Ze,Pi{i-1};0,Ze']+(1-q)*[Ze',0;Pi{i-1},Ze]+q*[0,Ze';Ze,Pi{i-1}];  
                    Pi{i}(2:(end-1),:) = Pi{i}(2:(end-1),:)/2;
                end
            end
            Pi = Pi{end};
            
            % State space
            Psi= sqrt(N-1)*sqrt(shockvar);
            Y = linspace(-Psi,Psi,N);
            
end
    
function Lam = InvariantRouwenhorst(N,shockvar,p,q)
    
    % Symmetry is default
    if nargin < 4; q = p; end
    
    s = (1-q)/(2-(p+q));
    Lam     = zeros(N,1);
    for i = 1:N
        Lam(i) = nchoosek(N-1,i-1).*s.^(N-i).*(1-s).^(i-1);
    end
    
end
        
    end
end