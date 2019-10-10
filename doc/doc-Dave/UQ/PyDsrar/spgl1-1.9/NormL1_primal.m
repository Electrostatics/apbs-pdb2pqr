function p = NormL1_primal(x,weights)

p = norm(x.*weights,1);
