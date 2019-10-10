function d = NormL1_dual(x,weights)

d = norm(x./weights,inf);
