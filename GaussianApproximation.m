classdef GaussianApproximation
    methods(Static)
        %% GAUSSIANAPPROXIMATION TOOLKIT
        % Written by Curtis Aquino (2019)
        % Contains:
        %   1. Initialization()
        %   2. Polynomial()
        %   3. NodesWeights()
        %   4. Interpolate()
        
%% 1. INITIALIZATION      
% This function is required for all functions to run.
%==========================================================================
function [T,Spprt]      = Initialization(GaussianPolynomial)
    Types   = {
            'Gauss-Legendre';...
            'Gauss-Chebyshev';...
            'Gauss-Hermite';...
            'Gauss-Laguerre'
            };
    T       = find(strcmp(GaussianPolynomial,Types));
    Spprt   = length(Types);
end

%% 2. POLYNOMIAL
% Creates the Nth order polynomials for a class of Gaussian polynomials.
%==========================================================================
function Polynomial     = Polynomial(PolynomialOrder,GaussianPolynomial)
    [T,Spprt]   = GaussianApproximation.Initialization(GaussianPolynomial);
    N           = PolynomialOrder+1;
    X           = sym('X','real');
    Coef        = sym(zeros(2,N+2,Spprt));
    Coef(:,:,1) = [[0,0,X.*(2.*(1:N)+1)./(2:(N+1))];[0,0,(1:N)./(2:(N+1))]];
    Coef(:,:,2) = [[0,0,repmat(2*X,1,N)]; [0,0,ones(1,N)]];
    Coef(:,:,3) = [[0,0,repmat(2*X,1,N)]; [0,0,2*(1:N)]];
    Coef(:,:,4) = [[0,0,(2.*(1:N)+1-X)./(2:(N+1))];[0,0,(1:N)./(2:(N+1))]];
    F2          = [X,X,2*X,1-X];
    Polynomial           = [1,F2(T)];
    for i = 3:N
        Polynomial(:,i)  = Coef(1,i,T).*Polynomial(:,i-1)-Coef(2,i,T).*Polynomial(:,i-2); 
    end
end

%% 3. NODESWEIGHTS
% Finds the nodes and weights for a variety of Gaussian polynomials.
%==========================================================================
function [Nodes,Weights]= NodesWeights(PolynomialOrder,GaussianPolynomial)
    N       = PolynomialOrder;
    [T,~]   = GaussianApproximation.Initialization(GaussianPolynomial);
    if T == 2
        Nodes       = -cos((2*(1:N)-1)*pi/(2*N))';
        Weights       = repmat(pi/N,N,1);
    else
        Var     = GaussianApproximation.Polynomial(N+1,GaussianPolynomial);
        Var     = Var(end-1);
        Jac     = jacobian(Var);
        Nodes       = real(double(vpasolve(Var)));    
        if T == 1 
            Weights       = 2./double((1-Nodes.^2).*subs(Jac,Nodes).^2);
        elseif T ==3 
            Weights       = double((2^(N+1)*factorial(N)*sqrt(pi))./(subs(Jac,Nodes).^2));
        elseif T == 4
            Weights       = double(1./(Nodes.*(subs(Jac,Nodes).^2)));
        end
    end
end

%% 4. INTERPOLATE
% Given Y, uses M nodes to interpolate an Nth order polynomial.
%==========================================================================
function InterpolatedFunction = Interpolate(YData,PolynomialOrder,MNodes,GaussianPolynomial)
    if length(YData) ~= MNodes
        fprintf('Error: X and Y data do not match')
        return
    end
    N       = PolynomialOrder;
    [X,]    = GaussianApproximation.NodesWeights(MNodes,GaussianPolynomial);
    poly    = GaussianApproximation.Polynomial(N,GaussianPolynomial);
    for i = 1:(N+1)
        Coef(i) = double(sum(YData.*subs(poly(i),X))/sum(subs(poly(i),X).^2));
    end
    InterpolatedFunction = sum(Coef.*poly);
end

    end
end