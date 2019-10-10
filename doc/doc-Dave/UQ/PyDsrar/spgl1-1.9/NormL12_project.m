function x = NormL12_project(g,x,weights,tau)
% Projection with number of groups equal to g

% Convert to matrix
m = round(length(x) / g); n = g;
x = reshape(x,m,n);

% Compute two-norms of rows
if isreal(x)
   xa  = sqrt(sum(x.^2,2));
else
   xa  = sqrt(sum(abs(x).^2,2));
end

% Project one one-norm ball
idx = xa < eps;
xc  = oneProjector(xa,weights,tau);

% Scale original
xc  = xc ./ xa; xc(idx) = 0;
x   = spdiags(xc,0,m,m)*x;

% Vectorize result
x = x(:);
