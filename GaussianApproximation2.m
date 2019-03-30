classdef GaussianApproximation2
    methods(Static)
        %% GAUSSIANAPPROXIMATION TOOLKIT
        %
        % Written by Curtis Aquino (2019). Contains:
        %
        % # Initializiation() is required for all functions to run.
        % # SymbolicPolynomial() creates symbolic expressions for the first N polynomials for a class of Gaussian polynomials. 
        % # NumericPolynomial() returns numerical solutions for the first N polynomials for a class of Gaussian polynomials evaluated at a set of points. 
        % # SymbolicNodesWeights() finds the nodes and weights for a variety of Gaussian polynomials using MATLAB's symbolic solver.
        % # NumericNodesWeights() finds the nodes and weights for a variety of Gaussian polynomials numerically.
        % # SymbolicInterpolate() interpolates the Y data using M nodes by an Nth order polynomial for a variety of Gaussian polynomials.
        % # NumericInterpolate() interpolates the Y data using M nodes by an Nth order polynomial for a variety of Gaussian polynomials.
   
%% 1. INITIALIZATION      
function [T,Sp]     = Initialization(GaussianPolynomial)
    Types   = {
            'Gauss-Legendre';...
            'Gauss-Chebyshev';...
            'Gauss-Hermite';...
            'Gauss-Laguerre'
            };
    T       = find(strcmp(GaussianPolynomial,Types));
    Sp   = length(Types);
end
%% 2. SYMBOLICPOLYNOMIAL
function P          = SymbolicPolynomial(PolynomialOrder,GaussianPolynomial)
    
    [T,~]           = GaussianApproximation2.Initialization(GaussianPolynomial);
    N               = PolynomialOrder+1;
    X               = sym('X','real');
    P      = sym(zeros(N,1));
    P(1)   = 1;
    
    if T == 1
       Coef     = [[0,0,X.*(2.*(1:N)+1)./(2:(N+1))];[0,0,(1:N)./(2:(N+1))]];
       P(2) = X;
    end 
    
    if T == 2
       Coef     = [[0,0,repmat(2*X,1,N)]; [0,0,ones(1,N)]];
       P(2) = X;
    end
        
    if T == 3
       Coef     = [[0,0,repmat(2*X,1,N)]; [0,0,2*(1:N)]];
       P(2) = 2*X;
    end
    
    if T == 4
       Coef     = [[0,0,(2.*(1:N)+1-X)./(2:(N+1))];[0,0,(1:N)./(2:(N+1))]];
       P(2) = 1-X;
    end
    
    for i = 3:N
       P(i) =  Coef(1,i)*P(i-1)-Coef(2,i)*P(i-2);
    end
end
%% 3. NUMERICPOLYNOMIAL
function P          = NumericPolynomial(PolynomialOrder,GaussianPolynomial,X)
    
    % Initialization
    [T,~]           = GaussianApproximation2.Initialization(GaussianPolynomial);
    N               = PolynomialOrder+1;
    XLen            = length(X);
    P(:,1) = ones(XLen,1);
    
    % Fixes input size
    if isrow(X); X = X'; end
    
    % Initialize each Gaussian polynomial
    if T == 1
       Coef(:,:,1) = [zeros(XLen,2),(X.*(2.*(1:N)+1))./(2:(N+1))]; 
       Coef(:,:,2) = repmat([0,0,(1:N)./(2:(N+1))],XLen,1);
       P(:,2) = X;
    end 
    if T == 2
       Coef(:,:,1) = [zeros(XLen,2),repmat(2*X,1,N)];
       Coef(:,:,2) = repmat([0,0,ones(1,N)],XLen,1);
       P(:,2) = X;
    end
    if T == 3
       Coef(:,:,1) = [zeros(XLen,2),repmat(2*X,1,N)];
       Coef(:,:,2) = repmat([0,0,2*(1:N)],XLen,1);
       P(:,2) = 2*X;
    end
    if T == 4
       Coef(:,:,1) = [zeros(XLen,2),(2.*(1:N)+1-X)./(2:(N+1))]; 
       Coef(:,:,2) = repmat([0,0,(1:N)./(2:(N+1))],XLen,1);
       P(:,2) = 1-X;
    end
    
    % Recursion relation
    for i = 3:N
       P(:,i) = Coef(:,i,1).*P(:,i-1)-Coef(:,i,2).*P(:,i-2);
    end
    
    % Output table
    for i = 1:N
        ColName{i} = sprintf('Degree%i',i-1); %#ok<AGROW>
    end
    P  = array2table(P,'VariableNames',ColName);
end
%% 4. SYMBOLICNODESWEIGHTS
function [Nd,Wt]    = SymbolicNodesWeights(PolynomialOrder,GaussianPolynomial)
    N       = PolynomialOrder;
    [T,~]   = GaussianApproximation2.Initialization(GaussianPolynomial);
    if T == 2
        Nd       = -cos((2*(1:N)-1)*pi/(2*N))';
        Wt       = repmat(pi/N,N,1);
    else
        Var     = GaussianApproximation2.SymbolicPolynomial(N+1,GaussianPolynomial);
        Var     = Var(end-1);
        Jac     = jacobian(Var);
        Nd       = real(double(vpasolve(Var)));    
        if T == 1 
            Wt       = 2./double((1-Nd.^2).*subs(Jac,Nd).^2);
        elseif T ==3 
            Wt       = double((2^(N+1)*factorial(N)*sqrt(pi))./(subs(Jac,Nd).^2));
        elseif T == 4
            Wt       = double(1./(Nd.*(subs(Jac,Nd).^2)));
        end
    end
end
%% 5. NUMERICNODESWEIGHTS
function [Nd,Wt]    = NumericNodesWeights(PolynomialOrder,GaussianPolynomial)
    
    N       = PolynomialOrder;
    [T,~]   = GaussianApproximation2.Initialization(GaussianPolynomial);
    Fn      = @(X) GaussianApproximation2.NumericPolynomial(N,GaussianPolynomial,X);
    
    if T == 2
        Nd      = -cos((2*(1:N)-1)*pi/(2*N))';
        Wt      = repmat(pi/N,N,1);
    else
        if T == 1
            Range       = [-1,1];
        end
        if T == 3
            Range       = [-sqrt(4*N+1),sqrt(4*N+1)];
        end
        if T == 4
            Range       = [0,N+(N-1)*sqrt(N)];
        end
        
        % Find root grid
        GridPrecision = 100000;
        IntBi   = linspace(Range(1),Range(2),GridPrecision);
        Roots   = Fn(IntBi);
        Roots   = Roots{:,end};
        SgnChg  = [diff(sign(Roots)) ~= 0];
        Grid    = [[Range(1),IntBi(SgnChg)]',[IntBi(SgnChg),Range(2)]'];
        
        % Find root subgrid
        SubGridPrec     = 5000;
        SGrid           = [];
        for j = 1:size(Grid,1)
            IntBi   = linspace(Grid(j,1),Grid(j,2),SubGridPrec);
            Roots   = Fn(IntBi);
            Roots   = Roots{:,end};
            SgnChg  = [diff(sign(Roots)) ~= 0;false];
            Idx     = find(SgnChg == 1);
            
            if Idx == 1
                SGrid   = [SGrid;IntBi([find(SgnChg == 1),find(SgnChg == 1)+1])];
            elseif Idx == SubGridPrec
                SGrid   = [SGrid;IntBi([find(SgnChg == 1)-1,find(SgnChg == 1)])];
            else 
                SGrid   = [SGrid;IntBi([find(SgnChg == 1)-1,find(SgnChg == 1)+1])];
            end
        end
        
        % Bisection
        Bisection   = table(SGrid(:,1),SGrid(:,2),mean(SGrid,2),zeros(size(SGrid,1),1),'VariableNames',{'LoBi','HiBi','Guess','F'});
        itr = 1; thr1 = Inf; thr2 = Inf;
        while thr1 > 10^(-16) && thr2 > 10^(-16) && itr < 5000
            Roots       = Fn(Bisection.Guess);
            Bisection.F = Roots{:,N+1};
            Bisection(mod(1:size(Bisection,1),2)==1,:).F = Bisection(mod(1:size(Bisection,1),2)==1,:).F*-1;
            Bisection.LoBi(Bisection.F<0) = Bisection.Guess(Bisection.F<0);
            Bisection.HiBi(Bisection.F>0) = Bisection.Guess(Bisection.F>0);
            thr2        = norm(Bisection.Guess - mean([Bisection.LoBi,Bisection.HiBi],2));
            itr         = itr + 1;
            thr1        = norm(Bisection.F);
            Bisection.Guess = mean([Bisection.LoBi,Bisection.HiBi],2);
        end
        Nd       = Bisection.Guess;
        
        % Finds unique values within machine precision 10^(-16)
        for i = 16:-1:1
            Ndtest  = uniquetol(Nd,10^(-i));
            if length(Ndtest) == N
               Nd = sort(Ndtest);
               i = 1;
            end
        end
        
        % Calculate weights
        if T == 1 
            temp = GaussianApproximation2.NumericPolynomial(N+1,GaussianPolynomial,Nd);
            Wt = (2*(1-Nd.^2))./((N+1)^2*temp{:,end}.^2);
        end
        if T ==3 
            temp = GaussianApproximation2.NumericPolynomial(N+1,GaussianPolynomial,Nd);
            Wt = (2^(N+1)*factorial(N)*sqrt(pi))./(temp{:,end}.^2);
        end
        if T == 4
            temp = GaussianApproximation2.NumericPolynomial(N+1,GaussianPolynomial,Nd);
            Wt = Nd./((N+1)^2*temp{:,end}.^2);
        end
    end 
end    
%% 6. SYMBOLICINTERPOLATE
function I          = SymbolicInterpolate(YData,PolynomialOrder,GaussianPolynomial)
    
%     YData = 1./(linspace(-1,1)+1.1);
%     PolynomialOrder = 4;
%     GaussianPolynomial = 'Gauss-Chebyshev';
    
    MNodes  = length(YData);
    N       = PolynomialOrder;
    [X,~]   = GaussianApproximation2.SymbolicNodesWeights(MNodes,GaussianPolynomial);
    poly    = GaussianApproximation2.SymbolicPolynomial(N,GaussianPolynomial);
    
    
    
    for i = 1:(N+1)
        Coef(i) = double(sum(YData.*subs(poly(i),X))/sum(subs(poly(i),X).^2));
    end
    I = sum(Coef'.*poly);
end
%% 7. NUMERICINTERPOLATE
function I          = NumericInterpolate(YData,PolynomialOrder,GaussianPolynomial)
    if isrow(YData); YData = YData'; end
    MNodes  = length(YData);
    N       = PolynomialOrder;
    [X,~]   = GaussianApproximation2.NumericNodesWeights(MNodes,GaussianPolynomial);
    poly    = GaussianApproximation2.NumericPolynomial(N,GaussianPolynomial,X);
    for i = 1:(N+1)
        Coef(i) = double(sum(YData.*poly{:,i})/sum(poly{:,i}.^2));
    end    
    I   = @(X) sum(table2array(GaussianApproximation2.NumericPolynomial(N,GaussianPolynomial,X)).*Coef);
end
    end
end