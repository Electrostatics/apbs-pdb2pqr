
function orthogonal_basis_construct(poly_order, dim, fileName)

  %dim = 20; %SAME AS Nrandomdim, USER CHOSEN
  %poly_order = 3; %ALSO USER CHOSEN

  dim = double(dim);
  poly_order = double(poly_order);

  base = 'PyDsrar/output/';
  %xfull = load('../../generate_points/rand_sample_d4_N_200000.dat');
  %xfull = load([base 'random_MD_Traj_dim_20_N_125000.dat']);
  xfull = load([base fileName]);
  [row column] = size(xfull);
  % num_sample is the first half of the data set
  num_sample = floor(row/2);
  %num_sample = floor(row/4);

  %xfull = 2.0*xfull - 1;

  %[row col] = size(xfull);
  %mean_xfull = mean(xfull);
  %std_xfull = std(xfull);
  %xfull = (xfull - ones(row,1)*mean_xfull);
  %xfull = (xfull - ones(row,1)*mean_xfull)./(ones(row,1)*std_xfull);

  train_sample = xfull(1:num_sample,1:dim);
  test_sample = xfull(num_sample+1:end,1:dim);

  [num_sample ncol] = size(train_sample);

  num_basis = nchoosek(dim+poly_order, poly_order);
  val_basis_sample = ones(num_sample, num_basis);
  indx_mat = full_tensor(@tensor, dim, poly_order);
  d_coeff = zeros(num_basis, num_basis);

  %num_basis = 3;
  for ii = 2:num_basis
    for pp = 1:dim
      val_basis_sample(:,ii) = val_basis_sample(:,ii).*train_sample(:,pp).^indx_mat(ii,pp);
    end
    %% compute the d_coeff from 1 to ii-1
    for jj = 1:ii-1
      d_coeff(ii,jj) = val_basis_sample(:,jj)'*val_basis_sample(:,ii)/num_sample;
    end
    %% put all the terms together
    for jj = 1:ii-1
      val_basis_sample(:,ii) = val_basis_sample(:,ii) - d_coeff(ii,jj)*val_basis_sample(:,jj);
    end
    %% normalization
    sum_tmp = val_basis_sample(:,ii)'*val_basis_sample(:,ii)/num_sample;
    %val_basis_sample(:,ii)'*val_basis_sample(:,ii)
    d_coeff(ii,1:ii-1) = d_coeff(ii,1:ii-1)/sqrt(sum_tmp);
    d_coeff(ii,ii) = 1.0/sqrt(sum_tmp);
    val_basis_sample(:,ii) = val_basis_sample(:,ii)/sqrt(sum_tmp);
    %% update the measurement matrix
    M_full_basis(ii,ii) = 1.0;
    for jj = 1:ii-1
      M_full_basis(ii,jj) = val_basis_sample(:,jj)'*val_basis_sample(:,ii)/num_sample;
      M_full_basis(jj,ii) = M_full_basis(ii,jj);
    end
  end

  % VALIDATION
  %%% compute the measurematrix on train
  measure_train = ones(num_sample, num_basis);
  for ii = 2:num_basis
    for pp = 1:dim
      measure_train(:,ii) = measure_train(:,ii).*train_sample(:,pp).^indx_mat(ii,pp);
    end
    measure_train(:,ii) = measure_train(:,ii)*d_coeff(ii,ii);

    for jj = 1:ii-1
      measure_train(:,ii) = measure_train(:,ii) - d_coeff(ii,jj)*measure_train(:,jj);
    end
  end

  %%% compute the measurematrix on test
  [num_sample ncol] = size(test_sample);
  measure_test = ones(num_sample, num_basis);
  for ii = 2:num_basis
    for pp = 1:dim
      measure_test(:,ii) = measure_test(:,ii).*test_sample(:,pp).^indx_mat(ii,pp);
    end
    measure_test(:,ii) = measure_test(:,ii)*d_coeff(ii,ii);
    for jj = 1:ii-1
      measure_test(:,ii) = measure_test(:,ii) - d_coeff(ii,jj)*measure_test(:,jj);
    end
  end

  %% max value of the measurement on train
  max_train_list = max(abs(measure_train));
  max_test_list = max(abs(measure_test));

  %dlmwrite([base 'd_' num2str(dim) '_coeff_p' num2str(poly_order) '_scale_1' '.dat'], d_coeff, 'delimiter', ' ', 'precision', 10);
  dlmwrite([base 'coeff_scale_1.dat'], d_coeff, 'delimiter', ' ', 'precision', 10);

  %dlmwrite([base 'maxlist_' num2str(dim) '_coeff_p' num2str(poly_order) '_scale_1' '.dat'], [(1:num_basis)' max_train_list' max_test_list'], 'delimiter', ' ', 'precision', 10);
  dlmwrite([base 'maxlist_coeff__scale_1.dat'], [(1:num_basis)' max_train_list' max_test_list'], 'delimiter', ' ', 'precision', 10);

  %save([base 'M_first_half_full_basis_d_' num2str(dim) '_coeff_p' num2str(poly_order) '_scale_1' '.mat'], 'M_full_basis');
  save([base 'M_first_half_full_basis_coeff_scale_1.mat'], 'M_full_basis');

  %save([base 'Measurematrix_train_d_' num2str(dim) '_coeff_p' num2str(poly_order) '_scale_1' '.mat'], 'measure_train');
  save([base 'Measurematrix_train_coeff_scale_1.mat'], 'measure_train');

  %save([base 'Measurematrix_test_d_' num2str(dim) '_coeff_p' num2str(poly_order) '_scale_1' '.mat'], 'measure_test');
  save([base 'Measurematrix_test_coeff_scale_1.mat'], 'measure_test');
