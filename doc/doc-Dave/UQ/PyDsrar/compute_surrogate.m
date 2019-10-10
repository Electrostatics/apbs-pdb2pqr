%function compute_residue_p3_fun(sparse_OMP_delta)
% addpath('~/matlab_package/spgl1-1.7/');
function compute_surrogate(poly_order, Ndim, mdName, enName, startN, stopN, stepS)

  %dim = 20;
  %dim_orthogonal = 20;
  %poly_order = 3;
  dim=double(Ndim);
  dim_orthogonal=double(Ndim);
  poly_order = double(poly_order);

  %num_all_sample = 100000;
  nevery0 = 10;
  relative_tol = 1.0e-04;
  %relative_tol = 0;
  %num_all_sample = 100000;
  num_basis = nchoosek(dim+poly_order, poly_order);
  %load(['kernel_dim' num2str(dim) '_P' num2str(poly_order) '.mat']);
  seed = 20000;
  rng(seed);

  base = 'PyDsrar/output/';

  % num_sample_list  = 200:200:1400;
  num_sample_list = startN:stepS:stopN;

  %fname = [base 'random_MD_Traj_dim_20_N_125000.dat'];
  fname = [base mdName];
  f_read = sprintf('%s', fname);
  all_sample = load(f_read);
  [row col] = size(all_sample);

  mean_all_sample = mean(all_sample);
  std_all_sample = std(all_sample);
  all_sample = (all_sample - ones(row,1)*mean_all_sample);

  num_train_sample = floor(row/2);
  training_index = 1:num_train_sample;
  test_index = num_train_sample+1:row;
  %test_index = training_index;

  training_sample = all_sample(training_index,1:dim);
  test_sample = all_sample(test_index,1:dim);
  [num_test ncol] = size(test_sample);
  fprintf('******\n');
  size(test_sample);

  %data_full = load([base 'global_energy_dim_20_ALA_num_125000']);
  data_full = load([base enName]);
  all_obser = data_full(test_index);
  normalization = norm(all_obser);

  display('finish reading data');
  tstart = tic;

  full_measure_mat = zeros(num_test, num_basis);
  %%%%%%%%%measure_mat_name = ['full_measure_mat_dim_' num2str(dim) '_num_sample_' num2str(num_all_sample) '.mat'];
  fpath = 'PyDsrar/output/';
  measure_file_name  = 'Measurematrix_test_coeff_scale_1.mat';
  measure_mat_name = sprintf('%s%s',fpath, measure_file_name);

  if(exist(measure_mat_name, 'file'))
    load(measure_mat_name);
    full_measure_mat = measure_test;
  else
    parfor k = 1:num_test
      orthogonal_vec = eval_tensor_2D_poly(@orthogonal_2D, all_sample_orthogonal(k,:), ...
                                             poly_order, indx_mat_orthogonal, d_coeff, indx_orthogonal_map);
      hermite_vec = eval_tensor_poly(@hermite_norm, all_sample_hermite(k,:)', dim-dim_orthogonal, ...
                                             poly_order, indx_mat_hermite);
      full_measure_mat(k,:) = (orthogonal_vec.*hermite_vec)';
    end
    save(measure_mat_name, 'full_measure_mat', '-v7.3');
  end

  %%% add normalization to the measurematrix %%%
  train_measure_file_name  = 'Measurematrix_train_coeff_scale_1.mat';
  train_measure_mat_name = sprintf('%s%s',fpath, train_measure_file_name);
  if(exist(train_measure_mat_name, 'file'))
    load(train_measure_mat_name);
  else
    display('error, matrix not found');
  end
  maxlist = max(measure_train);
  %full_measure_mat = full_measure_mat.*(ones(num_test,1)*(maxlist.^(-1)));
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  num_iteration = 0;
  %recon_coef = zeros(num_basis, num_iteration+1);
  num_trial = 8;
  recon_coef = zeros(num_basis, num_trial);
  threshold = 1e-1;
  rmse = zeros(num_trial, num_iteration+1);
  surrogate_val = zeros(num_test, num_trial);

  %randn('seed', 1);
  %all_sample = randn(dim, 400*num_trial)';

  num_sets = 3;
  % reconstruction index
  recon_indx = cell(num_sets,1);
  % validation index
  valid_indx = cell(num_sets,1);
  % reconstruction measure matrix
  recon_matrix = cell(num_sets,1);
  % validation measure matrix
  valid_matrix = cell(num_sets,1);
  % reconstruction observation
  recon_obser = cell(num_sets,1);
  % validation measure matrix
  valid_obser = cell(num_sets,1);

  %rotate_U = eye(dim);
  rotate_U = eye(dim-2);

  %for num_sample = 60:20:60
  for ll = 1:length(num_sample_list)

    num_sample = num_sample_list(ll);
    delta = zeros(num_trial, 1);
    delta2= zeros(num_trial, 1);
    rmse(:) = 0;

    Nr = floor((num_sets-1)*num_sample/num_sets);
    sample_pnt = zeros(num_sample, dim);
    training = zeros(num_sample, 1);
    measure_mat = zeros(num_sample, num_basis);
    sample_indx = zeros(num_sample, num_trial);

    for i = 1:num_sets
      valid_indx{i} = i:num_sets:num_sample;
      recon_indx{i} = setxor(1:num_sample, valid_indx{i});
    end

    for t = 1:num_trial
      %shift = (t-1)*num_sample*floor(nevery/25);
      %shift = floor((t-1)*num_sample);
      %nevery = floor(1200/num_sample*nevery0);
      %sample_pnt = all_sample(shift+1:nevery:shift+num_sample*nevery,:);
      %sample_pnt = test_sample(shift+1:nevery:shift+num_sample*nevery,:);
      %measure_mat = full_measure_mat(shift+1:nevery:shift+num_sample*nevery,:);
      %training = all_obser(shift+1:nevery:shift+num_sample*nevery);
      rand_idx = randperm(num_test,num_sample);
      sample_pnt = test_sample(rand_idx,:);
      measure_mat = full_measure_mat(rand_idx,:);
      training = all_obser(rand_idx);
      sample_indx(:,t) = rand_idx';

      for i = 1:num_sets
        recon_matrix{i} = measure_mat(recon_indx{i}, :);
        valid_matrix{i} = measure_mat(valid_indx{i}, :);
        recon_obser{i} = training(recon_indx{i}, :);
        valid_obser{i} = training(valid_indx{i}, :);
      end

      %delta0 = delta_initial*(1600/num_sample);
      %delta0 = delta_initial*(1600/num_sample)^1.5;
      delta0 = relative_tol*norm(training)*3000/num_sample; %will be user specified
      %delta0 = relative_tol*norm(training);
      %delta0 = delta_initial;
      %num_delta = 16;
      %d_delta = 40^(1/(num_delta-1));
      num_delta = 5;
      d_delta = 10^(1/(num_delta-1));
      delta_r = zeros(num_delta, 1);
      delta_v = zeros(num_delta, 1);
      for i = 1:num_delta
        delta_r(i) = delta0;
        tmp = 0;
        for j = 1:num_sets
          %opts = spgSetParms('verbosity', 0, 'iteration', 10*num_basis);
          %opts = spgSetParms('verbosity', 0, 'iteration', max(80000,10*num_basis), 'bpTol', 1e-03, 'optTol', 1e-03, 'decTol', 1e-05);
          opts = spgSetParms('verbosity', 0, 'iteration', max(80000,10*num_basis));
          c = spg_bpdn(recon_matrix{j}, recon_obser{j}, delta_r(i), opts);
          % check the error
          tmp = tmp + norm(valid_matrix{j}*c-valid_obser{j}, 2.0);
        end
        delta_v(i) = tmp/num_sets;
        delta0 = delta0/d_delta;
      end
      [tmp, min_ind] = min(delta_v);
      delta_r(min_ind)
      delta(t) = delta_r(min_ind)*sqrt(num_sets/(num_sets-1));

      %opts = spgSetParms('verbosity', 0);
      [num_sample t delta(t)]
      %opts = spgSetParms('verbosity', 1, 'iteration', max(80000,10*num_basis), 'bpTol', 1e-07, 'optTol', 1e-07, 'decTol', 1e-08);
      opts = spgSetParms('verbosity', 1, 'iteration', max(80000,10*num_basis));
      %opts = spgSetParms('verbosity', 1, 'iteration', max(400000,10*num_basis), 'bpTol', 1e-14, 'optTol', 1e-14, 'decTol', 2e-15);
      c = spg_bpdn(measure_mat, training, delta(t), opts);
      %c = SolveOMP(measure_mat, training, num_basis, 500, 0, 0, 0, delta(t));
      %recon_coef(:,1) = c;
      recon_coef(:,t) = c;
      display('l1  validation error')
      l1_val = full_measure_mat*c;

      surrogate_val(:,t) = l1_val;
      rmse(t,1) = norm(full_measure_mat*c-all_obser)/normalization

      % Rotate
      %display('start rotating');

      % turn off rotate %for comparing with previous work (legendre hermite)
      if 1
      rotate_pnt = sample_pnt;
      all_rotate_pnt = all_sample;
      tmp_delta = delta(t);

      tmp_delta = delta(t)/(sqrt(num_sets/(num_sets-1)));
      for  iter = 1:num_iteration
        A = zeros(dim);
        for i = 1:dim
          %A(i,i) = c'*kernel{i,i}*c/2.0;
          A(i,i) = c'*ktest{i,i}*c/2.0;
          for j = i+1:dim
            A(i,j) = c'*ktest{i,j}*c;
            %A(i,j) = c'*kernel{i,j}*c;
          end
        end
        A = A+A';
        A_hermite = A(3:dim,3:dim);

        [U,S,V]=svd(A_hermite);
        if (norm(S) < 1e-6)
          %display('no rotation');
          S_stop(t) = iter-1;
          %break;
        end

        rotate_pnt_orthogonal = rotate_pnt(:,1:2);
        rotate_pnt_hermite = rotate_pnt(:,3:end);
        rotate_pnt_hermite = rotate_pnt_hermite*U;
        rotate_pnt(:,3:end) = rotate_pnt_hermite;

        all_rotate_pnt_orthogonal = all_rotate_pnt(:,1:2);
        all_rotate_pnt_hermite = all_rotate_pnt(:,3:end);
        all_rotate_pnt_hermite = all_rotate_pnt_hermite*U;
        all_rotate_pnt(:,3:end) = all_rotate_pnt_hermite;

        %rotate_pnt = rotate_pnt*U;
        %all_rotate_pnt = all_rotate_pnt*U;

        rotate_U = rotate_U*U;

        %parfor k = 1:num_sample
        for k = 1:num_sample
          orthogonal_vec = eval_tensor_2D_poly(@orthogonal_2D, rotate_pnt_orthogonal(k,:), ...
                                             poly_order, indx_mat_orthogonal, d_coeff);
          hermite_vec = eval_tensor_poly(@hermite_norm, rotate_pnt_hermite(k,:)', dim-2, ...
                                             poly_order, indx_mat_hermite);
          measure_mat(k,:) = (orthogonal_vec.*hermite_vec)';
          %measure_mat(k,:) = eval_tensor_poly(@hermite_norm, rotate_pnt(k,:)', dim, ...
          %                                    poly_order, indx_mat)';
        end

        for i = 1:num_sets
          recon_matrix{i} = measure_mat(recon_indx{i}, :);
          valid_matrix{i} = measure_mat(valid_indx{i}, :);
          recon_obser{i} = training(recon_indx{i}, :);
          valid_obser{i} = training(valid_indx{i}, :);
        end
        delta0 = tmp_delta;
        num_delta2 = 2;
        d_delta2 = 2^(1/(num_delta2-1));
        delta2_r = zeros(num_delta2, 1);
        delta2_v = zeros(num_delta2, 1);
        for i = 1:num_delta2
          delta2_r(i) = delta0;
          tmp = 0;
          for j = 1:num_sets
            opts = spgSetParms('verbosity', 0);
            c = spg_bpdn(recon_matrix{j}, recon_obser{j}, delta2_r(i), opts);
            % check the error
            tmp = tmp + norm(valid_matrix{j}*c-valid_obser{j}, 2.0);
          end
          delta2_v(i) = tmp/num_sets;
          delta0 = delta0/d_delta2;
        end
        [tmp, min_ind] = min(delta2_v);
        %tmp_delta = delta2_r(min_ind);
        delta2(t) = delta2_r(min_ind)*sqrt(num_sets/(num_sets-1));

        opts = spgSetParms('verbosity', 0);
        c = spg_bpdn(measure_mat, training, delta2(t), opts);
        recon_coef(:,iter+1) = c;
        display('rotation error')
        for k = 1:num_all_sample
        %parfor k = 1:num_all_sample
          orthogonal_vec = eval_tensor_2D_poly(@orthogonal_2D, all_rotate_pnt_orthogonal(k,:), ...
                                             poly_order, indx_mat_orthogonal, d_coeff);
          hermite_vec = eval_tensor_poly(@hermite_norm, all_rotate_pnt_hermite(k,:)', dim-2, ...
                                             poly_order, indx_mat_hermite);
          full_rotate_measure_mat(k,:) = (orthogonal_vec.*hermite_vec)';
          %full_rotate_measure_mat(k,:) = eval_tensor_poly(@hermite_norm, all_rotate_pnt(k,:)', dim, ...
          %                                                poly_order, indx_mat)';
        end
        rotate_val = full_rotate_measure_mat*c;
        rmse(t,1+iter) = norm(full_rotate_measure_mat*c-all_obser)/normalization
        end
      end
    end
    %save([base 'circle_orth_basis_test_rand_poly_basis_unscale_set_4_dim' num2str(dim) '_p' num2str(poly_order) '_sam' num2str(num_sample) '.mat'], ...
    %      'rmse', 'recon_coef', 'num_trial','l1_val', 'relative_tol',  'sample_indx');
    save([base 'circle_orth_basis_test_rand_poly_basis_unscale_set_4_dim_sam' num2str(num_sample) '.mat'], ...
          'rmse', 'recon_coef', 'num_trial','l1_val', 'relative_tol',  'sample_indx');
  end
  tellapsed = toc(tstart)
