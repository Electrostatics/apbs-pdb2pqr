function coeff_analysis_single(poly_order, fileName)

  base = 'PyDsrar/output/';

  %data = load([base 'd_20_coeff_p3_scale_1.dat']);
  data = load([base fileName]);
  %poly_order = 3;
  poly_order = double(poly_order);

  [row col] = size(data);
  num_basis = row;
  d_nominal_coeff = data;
  d_nominal_coeff(1,1) = 1.0;

  d_nominal_coeff = zeros(num_basis,num_basis);
  d_nominal_coeff(1,1) = 1.0;

  for ll = 2:num_basis
    d_nominal_coeff(ll,ll) = data(ll,ll);
    for mm = 1:ll-1
      d_nominal_coeff(ll,:) = d_nominal_coeff(ll,:) - data(ll,mm)*d_nominal_coeff(mm,:);
    end
  end

  %dlmwrite([base 'd_nominal_20_coeff_p3_scale_1.dat'], d_nominal_coeff, 'delimiter', ' ', 'precision', 12);
  %dlmwrite([base 'd_nominal_' Nrandomdim 'coeff_p' poly_order '_scale_1.dat'], d_nominal_coeff, 'delimiter', ' ', 'precision', 12);
  dlmwrite([base 'nominal_coeff_scale_1.dat'], d_nominal_coeff, 'delimiter', ' ', 'precision', 12);
