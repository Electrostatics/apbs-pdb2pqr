function p = NormGroupL2_primal(groups,x,weights)

if isreal(x)
   p = sum(weights.*sqrt(sum(groups * x.^2,2)));
else
   p = sum(weights.*sqrt(sum(groups * abs(x).^2,2)));
end
