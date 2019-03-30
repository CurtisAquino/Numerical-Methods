function F  = sample(X,N)

% Fix input sizes if required
if size(X,2) > size(X,1)
    X = X';
end

% Default size is to draw 
if nargin < 2
   N = size(X,1); 
end

% Do not allow sizes above the number of observations
if N > size(X,1)
    N = size(X,1);
end
    
F = X(unidrnd(size(X,1),[N,1]),:);

end