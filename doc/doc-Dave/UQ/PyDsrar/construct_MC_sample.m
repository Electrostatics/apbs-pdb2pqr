function construct_MC_sample(fileName, Nrandomdim, Npart, Nperpart, Natom)

  % clear all;
  % close all;

  tic
  %Natom = 22;
  Natom = double(Natom);
  Ndim = Natom*3;
  % Nrandomdim = 20;
  % Nsample = 125000*Natom;
  Nsample = Nperpart * Npart * Natom;

  base = 'PyDsrar/output/';

  %Nburn = 500*Natom+1;
  Nburn = 0*Natom+1;
  % Nreplica = 125000;
  Nreplica = Nperpart * Npart;
  nevery = 1;
  seed = 6767;
  rng(seed);

  % Npart = 1000;
  Nreplica_each = Nreplica/Npart;


  % MD_file_name = 'mol_ALA.txt';
  MD_file_name = fileName;

  atom_pos_file_name = 'PyDsrar/output/atom_pos.mat';

  if(~exist(atom_pos_file_name, 'file'))
    % rawdata = load('mol_ALA.txt');
    rawdata = load(fileName);
    rawdata = 10*rawdata;
    %xpos = rawdata(:,4:6);
    xpos_full = rawdata(:,4:6);
    fprintf('-- Saving file: %s\n', atom_pos_file_name);
    save(atom_pos_file_name, 'xpos_full');
  else
    load(atom_pos_file_name);
  end

  %xpos = rawdata(Nburn:Nsample,4:6);
  xpos = xpos_full(Nburn:Nsample,:);
  [nrow ncol] = size(xpos);

  nsize = nrow*ncol;
  xlist = reshape(xpos',[1 nsize]);
  xarray_full = reshape(xlist, [Ndim nsize/Ndim]);
  xarray_full = xarray_full';
  xarray = xarray_full(1:nevery:end,:);
  nloop = nsize/Ndim/nevery;

  xcen = mean(xarray(:,1:Ndim));
  %%%% construct the covariance matrix %%%%%
  xcov = zeros(Ndim, Ndim);
  xcov_tmp = xcov;

  for ll = 1:nloop
      xtmp_vec = xarray(ll,:) - xcen(1,:);
      xcov_tmp = xtmp_vec'*xtmp_vec;
      xcov = xcov + xcov_tmp;
  end

  xcov = xcov/nloop;
  %[V D] = eigs(xcov, 24);
  [V D] = eig(xcov);
  %% sort the matrix
  [D_temp, ind]=sort(diag(D),'descend');
  D = diag(sort(diag(D), 'descend'));
  V = V(:, ind);
  %dlmwrite('xcov.dat', xcov);

  D_reduce = zeros(Ndim, Ndim);
  Deigen_reduce = zeros(Ndim, Ndim);

  for i = 1:Nrandomdim
    D_reduce(i,i) = 1.0;
    Deigen_reduce(i,i) = 1.0/D(i,i);
  end

  MT = V*D_reduce*V';
  Deigen_reduce_half = Deigen_reduce^0.5;

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %% write out the equilibrium set of data files
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  file_base = 'PyDsrar/output/testmol_eq';
  file_tail = '.pqr';
  file_name = sprintf('%s%s', file_base,file_tail);
  fid = fopen(file_name,'w');
  fprintf('-- Writing to: %s\n', file_name);

  for i = 1:Natom
      %dx = perturb_full(3*i-2, loop);
      %dy = perturb_full(3*i-1, loop);
      %dz = perturb_full(3*i, loop);
    fprintf(fid, '%4.9f %4.9f %4.9f\n', xcen(1,3*i-2), xcen(1,3*i-1), xcen(1,3*i));
  end

  fclose(fid);

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %% read the position from pqr data file
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %xyz = load('xyz_eq.dat');
  %xeq = reshape(xyz', [1 Ndim]);

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %% generate MC input
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  xcoord = zeros(Ndim, Nreplica);
  npoint3D = randn(Nrandomdim,Nreplica);

  for ll = 1:Nreplica
      %xtmp_vec = xarray(2*ll,:) - xcen(1,:);
      xtmp_vec = xarray(ll,:) - xcen(1,:);
      xtransform = MT * xtmp_vec';
      %xcoord(:,ll) = xeq' + xtransform;
      xcoord(:,ll) = xcen' + xtransform;
      rtmp_vec = Deigen_reduce_half*V'*xtmp_vec';
      npoint3D(:,ll) = rtmp_vec(1:Nrandomdim);
  end
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


  if 1
    %tmpstr = [base 'random_MD_Traj_dim_' num2str(Nrandomdim) '_N_' num2str(Nreplica) '.dat'];
    tmpstr = [base 'random_MD_Traj.dat'];

    fprintf('-- Writing file: %s\n', tmpstr);
    dlmwrite(tmpstr, npoint3D', 'delimiter', ' ', 'precision', 15);

    %tmpstr = [base 'ala_MD_Traj_dim_' num2str(Nrandomdim) '_N_' num2str(Nreplica) '.dat'];
    tmpstr = [base 'mol_MD_Traj.dat'];
    fprintf('-- Writing file: %s\n', tmpstr);
    dlmwrite(tmpstr, xcoord', 'delimiter', ' ', 'precision', 15);

    fprintf('-- Writing files to: PyDsrar/output/APBS\n');
    for ii = 1:Npart
      %tmpstr = [base 'APBS/x_sample_dim_' num2str(Nrandomdim) '_part_' num2str(ii) '_list.dat'];
      tmpstr = [base 'APBS/x_sample_' num2str(ii) '_list.dat'];
      %fprintf('--Writing file: %s\n', tmpstr);
      dlmwrite(tmpstr, ...
            xcoord(1:3:Ndim-2,Nreplica_each*(ii-1)+1:Nreplica_each*ii), 'delimiter', ' ', 'precision', 15);

      %tmpstr = [base 'APBS/y_sample_dim_' num2str(Nrandomdim) '_part_' num2str(ii) '_list.dat'];
      tmpstr = [base 'APBS/y_sample_' num2str(ii) '_list.dat'];
      %fprintf('-- Writing file: %s\n', tmpstr);
      dlmwrite(tmpstr, ...
            xcoord(2:3:Ndim-1,Nreplica_each*(ii-1)+1:Nreplica_each*ii), 'delimiter', ' ', 'precision', 15);

      %tmpstr = [base 'APBS/z_sample_dim_' num2str(Nrandomdim) '_part_' num2str(ii) '_list.dat'];
      tmpstr = [base 'APBS/z_sample_' num2str(ii) '_list.dat'];
      %fprintf('-- Writing file: %s\n', tmpstr);
      dlmwrite(tmpstr, ...
            xcoord(3:3:Ndim,Nreplica_each*(ii-1)+1:Nreplica_each*ii), 'delimiter', ' ', 'precision', 15);
    end
  end

  toc

end
