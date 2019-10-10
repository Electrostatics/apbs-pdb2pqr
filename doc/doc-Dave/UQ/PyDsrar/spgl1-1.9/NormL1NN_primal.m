function p = NormL1NN_primal(x,weights)
% Non-negative L1 gauge function

p = norm(x.*weights,1);
if any(x < 0) p = Inf; end;

