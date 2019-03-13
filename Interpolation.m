classdef Interpolation
    methods(Static)
        %% INTERPOLATION TOOLKIT
        % Written by Curtis Aquino (2019). Contains:
        % 
        % # PiecewiseLinearInterpolation() performs linear piecewise interpolation given x- and y-data.
        % # CubicSplinePlus() uses MATLAB's spline() function but endogenously turns it into a piecewise symbolic function.
        % # SchumakerSpline() implements the shape-preserving spline proposed by Schumaker.
        
%% PIECEWISELINEARINTERPOLATION
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
%% CUBICSPLINEPLUS
function CubicSplinePlus(X,Y)
    temp        = spline(X,Y);
    F           = sum(temp.coefs.*(sym('X')-X(1:(end-1))').^(3:-1:0),2); %#ok<NASGU>
    % Endogenous Piecewise
    syms X
    Text = sprintf('piecewise(');
    for i = 1:(length(temp.breaks)-1)
        if i == (length(temp.breaks)-1)
            Text = [Text, sprintf('temp.breaks(%i) <= X & X <= temp.breaks(%i),F(%i))',i,i+1,i)];       %#ok<*AGROW>
        elseif i == 1
            Text = [Text, sprintf('temp.breaks(%i) <= X & X <= temp.breaks(%i),F(%i),',i,i+1,i)];       %#ok<*AGROW>
        else
            Text = [Text, sprintf('temp.breaks(%i) <= X & X <= temp.breaks(%i),F(%i),',i,i+1,i)];
        end
    end
    F = eval(Text); %#ok<NASGU>
end
%% SCHUMAKERSPLINE
function F = SchumakerSpline(x)
    x = 1;
    F = x;
end
    end
end