function x = NormGroupL2_project(groups,x,weights,tau)
% Projection binary group matrix

% Compute two-norms of rows
if isreal(x)
   xa  = sqrt(sum(groups * x.^2,2));
else
   xa  = sqrt(sum(groups * abs(x).^2,2));
end

% Project one one-norm ball
idx = xa < eps;
xc  = oneProjector(xa,weights,tau);

% Scale original
xc  = xc ./ xa; xc(idx) = 0;
x   = full(groups' * xc).*x;
