function d = NormGroupL2_dual(groups,x,weights)

if isreal(x)
   d = norm(sqrt(sum(groups * x.^2,2))./weights,inf);
else
   d = norm(sqrt(sum(groups * abs(x).^2,2))./weights,inf);
end
