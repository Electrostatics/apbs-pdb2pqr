% Give the index of a full tensor product basis:
%
%   \phi_a = \phi_a1 \phi_a2 \phi_a3 ... \phi_ad,
%
% where \sum_i^d ai <= order.
%

function [matrix] = full_tensor(tensor, dim, order)
matrix = zeros(1, dim);
for p = 1:order
  submatrix = tensor(dim, p);
  matrix = [matrix; submatrix];
end
