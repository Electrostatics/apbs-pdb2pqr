function p = NormL12_primal(g,x,weights)

m = round(length(x) / g); n = g;

if isreal(x)
   p = sum(weights.*sqrt(sum(reshape(x,m,n).^2,2)));
else
   p = sum(weights.*sqrt(sum(abs(reshape(x,m,n)).^2,2)));
end
