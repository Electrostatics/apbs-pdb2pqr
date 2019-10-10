function d = NormL12_dual(g,x,weights)

m = round(length(x) / g); n = g;

if isreal(x)
   d = norm(sqrt(sum(reshape(x,m,n).^2,2))./weights,inf);
else
   d = norm(sqrt(sum(abs(reshape(x,m,n)).^2,2))./weights,inf);
end

  
