classdef Interpolation
    methods(Static)
        %% INTERPOLATION TOOLKIT
        %
        % Written by Curtis Aquino (2019). Contains:
        % 
        % # PiecewiseLinearInterpolation() performs linear piecewise interpolation given x- and y-data.
        % # ButlandSchumakerSpline() implements the shape-preserving spline described in Iqbar (1993).
        
%% 1. PIECEWISELINEARINTERPOLATION
function F = PiecewiseLinearInterpolation(x,y,x0)
    [~,idx] = sort(x);
    x       = x(idx);
    y       = y(idx);
    bds     = [x(1:(end-1)),x(2:end)];
    n       = length(x);
    m       = length(x0);
    a       = (y(2:n) - y(1:(n-1)))./(x(2:n)-x(1:(n-1)));
    b       = (y(1:(n-1)).*x(2:n)-y(2:n).*x(1:(n-1)))./(x(2:n)-x(1:(n-1)));
    for i = 1:m
        idx     = find(bds(:,1) <= x0(i) & x0(i) <= bds(:,2));
        if length(idx) > 1
            idx     = find(bds(:,1) <= x0(i) & x0(i) < bds(:,2));
        end
        if isempty(idx)
            if x0(i) < x(1)
                F(i) = x(1);
            elseif x0(i) > x(end)
                F(i) = x(end);
            end
        else
            F(i)    = a(idx)*x0(i)+b(idx);
        end
    end
end
%% 2. BUTLANDSCHUMAKERSPLINE
function F = ButlandSchumakerSpline(X,Y,X0)
    % Fixes user inputs
    if isrow(X0); X0 = X0'; end
    if isrow(X); X = X'; end
    if isrow(Y); Y = Y'; end
    [~,idx]     = sort(X);
    X           = X(idx);
    Y           = Y(idx);
    
    % Build what Iqbal (1992) calls the Butland formulas
    del     = (Y(2:end)-Y(1:(end-1)))./(X(2:end)-X(1:(end-1)));
    temp1   = del(2:end).*del(1:(end-1));
    d       = ((2*temp1)./(del(2:end)+del(1:(end-1)))).*(temp1>0);
    d1      = (2*del(1)-d(1)).*(del(1)*(2*del(1)-d(1))>0);
    dn      = (2*del(end)-d(end)).*(del(end)*(2*del(end)-d(end))>0);
    d       = [d1;d;dn];
    
    % Builds what Iqbal (1992) calls xi and dbar
    xi1c    = ((d(2:end)-del).*(d(1:(end-1))-del))<0;
    xi2c    = ((d(2:end)-del).*(d(1:(end-1))-del))>=0;
    xi1     = X(2:end)+(d(1:(end-1))-del).*(X(2:end)-X(1:(end-1)))./(d(2:end)-d(1:(end-1)));
    xi2     = (X(1:(end-1))+X(2:end))/2;
    xi      = xi1.*xi1c+xi2.*xi2c;
    dbar    = (2*del-d(2:end))+(d(2:end)-d(1:(end-1))).*(xi-X(1:(end-1)))./(X(2:end)-X(1:(end-1)));
    
    % Builds all of the coefficients
    Coef.A  = [Y(1:(end-1)),Y(1:(end-1)),Y(1:(end-1))+d(1:(end-1)).*(xi-X(1:(end-1)))+(dbar-d(1:(end-1))).*(xi-X(1:(end-1)))/2];
    Coef.B  = [d(1:(end-1)),d(1:(end-1)),dbar];
    Coef.C  = [0.5*(d(2:end)-d(1:(end-1)))./(X(2:end)-X(1:(end-1))),0.5*(dbar-d(1:(end-1)))./(xi-X(1:(end-1))),0.5*(d(2:end)-dbar)./(X(2:end)-xi)];
    
    % Checks which points need interior knots
    intknot = d(1:(end-1))+d(2:end) == 2*del;
    
    % Builds coefficients and ranges
    Range = []; Coefs = [];
    for i = 1:length(intknot)
        if intknot(i) == 1
            Range   = [Range;[X(i),X(i+1)]];
            Coefs   = [Coefs;Coef.A(i,1),Coef.B(i,1),Coef.C(i,1)];
        else
            Range   = [Range;[X(i),xi(i)];[xi(i),X(i+1)]];
            Coefs   = [Coefs;Coef.A(i,2),Coef.B(i,2),Coef.C(i,2)];
            Coefs   = [Coefs;Coef.A(i,3),Coef.B(i,3),Coef.C(i,3)];
        end
    end
    
    % Interpolates
    for i = 1:size(X0,1)
        temp        = Range(:,1) <= X0(i) & X0(i) <= Range(:,2);
        F(i)        = max(sum(Coefs(temp,:).*[ones(size(Coefs(temp,:),1),1),(X0(i)-Range(temp,1)),(X0(i)-Range(temp,1)).^2],2));
    end
    
end
    end
end