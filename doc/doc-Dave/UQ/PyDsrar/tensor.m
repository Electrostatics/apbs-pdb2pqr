% Give the index of a full tensor product basis:
%
%   \phi_a = \phi_a1 \phi_a2 \phi_a3 ... \phi_ad,
%
% where \sum_i^d ai = order.
% 

function [matrix] = tensor(dim, order)
matrix = [];
if (order == 0)
  matrix = zeros(1,dim);
elseif (order == 1)
  matrix = eye(dim);
else
  % when one of x_i = order
  for ind = 1:dim-1
    tmp = zeros(1,dim);
    tmp(ind) = order;
    matrix = [matrix; tmp];
    for p = order-1 : -1 : 1
      submatrix = tensor(dim-ind, order-p);
      subrow = size(submatrix, 1);
      tmp = [zeros(subrow, ind-1)  p*ones(subrow, 1) submatrix];
      matrix = [matrix; tmp];
    end
  end
  matrix = [matrix; [zeros(1, dim-1) order]];
end
