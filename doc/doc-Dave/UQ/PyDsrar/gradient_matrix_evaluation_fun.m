function gradient_matrix_evaluation_fun(num_train, poly_order, dim, fileName)
tic
%dim = 20;
%poly_order = 3;
%num_train = 200;
poly_order = double(poly_order);

base = 'PyDsrar/output/';

%xfull = load('../../generate_points/rand_sample_dim_25_N_200000.dat');
%xfull = load([base 'random_MD_Traj_dim_20_N_125000.dat']);
xfull = load([base fileName]);
[row column] = size(xfull);
% num_sample is the first half of the data set
%num_sample = floor(row/2);
num_sample = floor(row/2);
%num_sample = floor(row/4);

[row col] = size(xfull);
mean_xfull = mean(xfull);
std_xfull = std(xfull);

train_sample = xfull(1:num_sample,1:dim);
test_sample = xfull(num_sample+1:end,1:dim);

[num_sample ncol] = size(train_sample);

num_basis = nchoosek(dim+poly_order, poly_order);
indx_mat = full_tensor(@tensor, dim, poly_order);
%val_basis_sample = ones(num_sample, num_basis);
%der_basis_sample = ones(num_sample, num_basis, dim);
%der_nominal_basis = ones(num_sample, num_basis, dim);
der_nominal_basis = ones(num_sample, num_basis);
der_QoI_sample = ones(num_sample, dim);
der_QoI_ll = ones(num_sample, 1);

M_basis_transfer = zeros(num_basis, num_basis);

%%%% read the d_coeff and d_nominal_coeff of the constructed  basis
%d_coeff = zeros(num_basis, num_basis);
%d_nominal_coeff = zeros(num_basis,num_basis);
%fpath = '../d_25_p_3/';
%d_coeff_name = 'd_25_coeff_p3_scale_1.dat';
%d_nominal_coeff_name = 'd_nominal_25_coeff_p3_scale_1.dat';
%d_coeff = load([fpath d_coeff_name']);
d_coeff = load([base 'coeff_scale_1.dat']);
d_nominal_coeff = load([base 'nominal_coeff_scale_1.dat']); %load from coeff analysis single mat

%%% read the gPC coeff %%%%%%%
%load(['/people/leih646/work/uq/biomolecule/blb/dim_12/test_orth_basis/dim_12_p_4_set_dim_12/circle_orth_basis_test_rand_poly_basis_unscale_set_2_dim12_p4_sam' num2str(num_train) '.mat']);
load([base 'circle_orth_basis_test_rand_poly_basis_unscale_set_4_dim_sam' num2str(num_train) '.mat']);
d_nominal_gPC_coeff = d_nominal_coeff'*recon_coef(:,1);

%%% pre-computed the polynomial basis derivative
for ll = 1:dim
  der_nominal_basis(:) = 1.0;
  der_nominal_basis(:,1) = 0;
  for ii = 2:num_basis
    indx_ll_val = indx_mat(ii,ll);
    indx_ll_mat = indx_mat(ii,:);
    indx_ll_mat(1,ll) = indx_ll_mat(1,ll)-1;
    for pp = 1:dim
      der_nominal_basis(:,ii) = der_nominal_basis(:,ii).*train_sample(:,pp).^indx_ll_mat(1,pp);
    end
    der_nominal_basis(:,ii) = indx_ll_val*der_nominal_basis(:,ii);
  end
  %der_basis_ll = der_nominal_basis*d_nominal_coeff';
  der_QoI_ll = der_nominal_basis*d_nominal_gPC_coeff;
  der_QoI_sample(:,ll) = der_QoI_ll;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Gradient_mat = der_QoI_sample'*der_QoI_sample/num_sample;
%save([base 'Gradient_mat_orth_basis_train_d_' num2str(dim) '_coeff_p' num2str(poly_order) '_sam_' num2str(num_train) '.mat'], 'Gradient_mat');
save([base 'Gradient_mat_orth_basis_train_coeff_sam' num2str(num_train) '.mat'], 'Gradient_mat');

[V D] = eig(Gradient_mat);
[dummy, sort_ind]=sort(diag(D), 'descend');
D=D(sort_ind,sort_ind);
V=V(:,sort_ind);
xfull_rotate = xfull * V;

%dlmwrite([base 'rotate_orth_basis_by_train_sam_' num2str(num_train) '_rand_sample_dim_20_N_125000.dat'], xfull_rotate, 'precision', 15, 'delimiter', ' ');
dlmwrite([base 'rotate_orth_basis_by_train_sam' num2str(num_train) '_rand_sample.dat'], xfull_rotate, 'precision', 15, 'delimiter', ' ');

%save(['Derivative_vec_orth_basis_train_d_' num2str(dim) '_coeff_p' num2str(poly_order) '_sam_' num2str(num_train) '_der_' num2str(ll) '.mat'], 'der_QoI_ll');
%save(['Derivative_basis_mat_orth_basis_train_d_' num2str(dim) '_coeff_p' num2str(poly_order) '_sam_' num2str(num_train) '_der_' num2str(ll) '.mat'], 'der_basis_ll');
toc
