classdef Interpolation
    methods(Static)
        %% INTERPOLATION TOOLKIT
        % Written by Curtis Aquino (2019)
        % Contains:
        %   1. PiecewiseLinearInterpolation()
        
%% 1. PIECEWISELINEARINTERPOLATION
% This function performs linear interpolation given x- and y-data.
%==========================================================================
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

    end
end