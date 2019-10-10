function d = NormL1NN_dual(x,weights)
% Dual of non-negative L1 gauge function

x(x < 0) = 0;
d = norm(x./weights,inf);
