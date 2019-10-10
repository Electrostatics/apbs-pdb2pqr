function spgdemo(interactive)
%DEMO  Demonstrates the use of the SPGL1 solver
%
% See also SPGL1.
    
%   demo.m
%   $Id: spgdemo.m 1079 2008-08-20 21:34:15Z ewout78 $
%
%   ----------------------------------------------------------------------
%   This file is part of SPGL1 (Spectral Projected Gradient for L1).
%
%   Copyright (C) 2007 Ewout van den Berg and Michael P. Friedlander,
%   Department of Computer Science, University of British Columbia, Canada.
%   All rights reserved. E-mail: <{ewout78,mpf}@cs.ubc.ca>.
%
%   SPGL1 is free software; you can redistribute it and/or modify it
%   under the terms of the GNU Lesser General Public License as
%   published by the Free Software Foundation; either version 2.1 of the
%   License, or (at your option) any later version.
%
%   SPGL1 is distributed in the hope that it will be useful, but WITHOUT
%   ANY WARRANTY; without even the implied warranty of MERCHANTABILITY
%   or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser General
%   Public License for more details.
%
%   You should have received a copy of the GNU Lesser General Public
%   License along with SPGL1; if not, write to the Free Software
%   Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301
%   USA
%   ----------------------------------------------------------------------
    if nargin < 1 || isempty(interactive), interactive = true; end
    
    % Initialize random number generators 
    rand('state',0);
    randn('state',0);
    
    % Create random m-by-n encoding matrix and sparse vector
    m = 50; n = 128; k = 14;
    [A,Rtmp] = qr(randn(n,m),0);
    A  = A';
    p  = randperm(n); p = p(1:k);
    x0 = zeros(n,1); x0(p) = randn(k,1);
 
    
    % -----------------------------------------------------------
    % Solve the underdetermined LASSO problem for ||x||_1 <= pi:
    %
    %    minimize ||Ax-b||_2 subject to ||x||_1 <= 3.14159...
    %
    % -----------------------------------------------------------
    fprintf(['%% ', repmat('-',1,78), '\n']);
    fprintf('%% Solve the underdetermined LASSO problem for   \n');
    fprintf('%%                                               \n');
    fprintf('%%   minimize ||Ax-b||_2 subject to ||x||_1 <= 3.14159...\n');
    fprintf('%%                                               \n');
    fprintf(['%% ', repmat('-',1,78), '\n']);

    fprintf('\nPress <return> to continue ... \n');
    if interactive, pause; end

    % Set up vector b, and run solver    
    b = A * x0;
    tau = pi;
    x = spg_lasso(A, b, tau);

    fprintf([repmat('-',1,35), ' Solution ', repmat('-',1,35), '\n']);
    fprintf(['nonzeros(x) = %i,   ' ...
             '||x||_1 = %12.6e,   ' ...
             '||x||_1 - pi = %13.6e\n'], ...
            length(find(abs(x)>1e-5)), norm(x,1), norm(x,1)-pi);
    fprintf([repmat('-',1,80), '\n']);
    fprintf('\n\n');

    
    % -----------------------------------------------------------
    % Solve the basis pursuit (BP) problem:
    %
    %    minimize ||x||_1 subject to Ax = b
    %
    % -----------------------------------------------------------
    fprintf(['%% ', repmat('-',1,78), '\n']);
    fprintf('%% Solve the basis pursuit (BP) problem:\n');
    fprintf('%%                                      \n');
    fprintf('%%   minimize ||x||_1 subject to Ax = b \n');
    fprintf('%%                                      \n');
    fprintf(['%% ', repmat('-',1,78), '\n']);

    fprintf('\nPress <return> to continue ... \n');
    if interactive, pause; end

    % Set up vector b, and run solver 
    b = A * x0;  % Signal
    opts = spgSetParms('verbosity',1);
    x = spg_bp(A, b, opts);
    
    figure(1); subplot(2,4,1);
    plot(1:n,x,'b', 1:n,x0,'ro');
    legend('Recovered coefficients','Original coefficients');
    title('(a) Basis Pursuit');

    fprintf([repmat('-',1,35), ' Solution ', repmat('-',1,35), '\n']);
    fprintf('See figure 1(a).\n');
    fprintf([repmat('-',1,80), '\n']);
    fprintf('\n\n');

    
    % -----------------------------------------------------------
    % Solve the basis pursuit denoise (BPDN) problem:
    %
    %    minimize ||x||_1 subject to ||Ax - b||_2 <= 0.1
    %
    % -----------------------------------------------------------
    fprintf(['%% ', repmat('-',1,78), '\n']);
    fprintf('%% Solve the basis pursuit denoise (BPDN) problem:  \n');
    fprintf('%%                                                  \n');
    fprintf('%%   minimize ||x||_1 subject to ||Ax - b||_2 <= 0.1\n');
    fprintf('%%                                                  \n');
    fprintf(['%% ', repmat('-',1,78), '\n']);

    fprintf('\nPress <return> to continue ... \n');
    if interactive, pause; end

    % Set up vector b, and run solver
    b = A * x0 + randn(m,1) * 0.075;
    sigma = 0.10;       % Desired ||Ax - b||_2
    opts = spgSetParms('verbosity',1);
    x = spg_bpdn(A, b, sigma, opts);
    
    figure(1); subplot(2,4,2);
    plot(1:n,x,'b', 1:n,x0,'ro');
    legend('Recovered coefficients','Original coefficients');
    title('(b) Basis Pursuit Denoise');

    fprintf([repmat('-',1,35), ' Solution ', repmat('-',1,35), '\n']);
    fprintf('See figure 1(b).\n');
    fprintf([repmat('-',1,80), '\n']);
    fprintf('\n\n');


    % -----------------------------------------------------------
    % Solve the basis pursuit (BP) problem in COMPLEX variables:
    %
    %    minimize ||z||_1 subject to Az = b
    %
    % -----------------------------------------------------------
    fprintf(['%% ', repmat('-',1,78), '\n']);
    fprintf('%% Solve the basis pursuit (BP) problem in COMPLEX variables:\n');
    fprintf('%%                                                    \n');
    fprintf('%%   minimize ||z||_1 subject to Az = b               \n');
    fprintf('%%                                                    \n');
    fprintf(['%% ', repmat('-',1,78), '\n']);

    fprintf('\nPress <return> to continue ... \n');
    if interactive, pause; end

    % Create partial Fourier operator with rows idx
    idx = randperm(n); idx = idx(1:m);
    opA = @(x,mode) partialFourier(idx,n,x,mode);

    % Create sparse coefficients and b = 'A' * z0;
    z0 = zeros(n,1);
    z0(p) = randn(k,1) + sqrt(-1) * randn(k,1);
    b = opA(z0,1);
    
    opts = spgSetParms('verbosity',1);
    z = spg_bp(opA,b,opts);
    
    figure(1); subplot(2,4,3);
    plot(1:n,real(z),'b+',1:n,real(z0),'bo', ...
         1:n,imag(z),'r+',1:n,imag(z0),'ro');
    legend('Recovered (real)', 'Original (real)', ...
           'Recovered (imag)', 'Original (imag)');
    title('(c) Complex Basis Pursuit');
    
    fprintf([repmat('-',1,35), ' Solution ', repmat('-',1,35), '\n']);
    fprintf('See figure 1(c).\n');
    fprintf([repmat('-',1,80), '\n']);
    fprintf('\n\n');
    
    
    % -----------------------------------------------------------
    % Sample the Pareto frontier at 100 points:
    %
    %    phi(tau) = minimize ||Ax-b||_2 subject to ||x|| <= tau
    %
    % -----------------------------------------------------------
    fprintf(['%% ', repmat('-',1,78), '\n']);
    fprintf('%% Sample the Pareto frontier at 100 points:\n');
    fprintf('%%                                              \n');
    fprintf('%%   phi(tau) = minimize ||Ax-b||_2 subject to ||x|| <= tau\n');
    fprintf('%%                                              \n');
    fprintf(['%% ', repmat('-',1,78), '\n']);

    fprintf('\nPress <return> to continue ... \n');
    if interactive, pause; end
    
    fprintf('\nComputing sample');
    

    % Set up vector b, and run solver    
    b = A*x0;
    x = zeros(n,1);
    tau = linspace(0,1.05 * norm(x0,1),100);
    phi = zeros(size(tau));
    
    opts = spgSetParms('iterations',1000,'verbosity',0);
    for i=1:length(tau)
        [x,r] = spgl1(A,b,tau(i),[],x,opts);
        phi(i) = norm(r);
        if ~mod(i,10), fprintf('...%i',i); end
    end
    fprintf('\n');
    
    figure(1); subplot(2,4,4);
    plot(tau,phi);
    title('(d) Pareto frontier');
    xlabel('||x||_1'); ylabel('||Ax-b||_2');
    
    fprintf('\n');
    fprintf([repmat('-',1,35), ' Solution ', repmat('-',1,35), '\n']);
    fprintf('See figure 1(d).\n');
    fprintf([repmat('-',1,80), '\n']);
    fprintf('\n\n');


    % -----------------------------------------------------------
    % Solve
    %
    %    minimize ||y||_1 subject to AW^{-1}y = b
    %    
    % and the weighted basis pursuit (BP) problem:
    %
    %    minimize ||Wx||_1 subject to Ax = b
    %
    % followed by setting y = Wx.
    % -----------------------------------------------------------
    fprintf(['%% ', repmat('-',1,78), '\n']);
    fprintf('%% Solve                                          \n');
    fprintf('%%                                                \n');
    fprintf('%% (1) minimize ||y||_1 subject to AW^{-1}y = b   \n');
    fprintf('%%                                                \n');
    fprintf('%% and the weighted basis pursuit (BP) problem:   \n');
    fprintf('%%                                                \n');
    fprintf('%% (2) minimize ||Wx||_1 subject to Ax = b        \n');
    fprintf('%%                                                \n');
    fprintf('%% followed by setting y = Wx.                    \n');
    fprintf(['%% ', repmat('-',1,78), '\n']);

    fprintf('\nPress <return> to continue ... \n');
    if interactive, pause; end

    % Sparsify vector x0 a bit more to get exact recovery
    k = 9;
    x0 = zeros(n,1); x0(p(1:k)) = randn(k,1);

    % Set up weights w and vector b 
    w     = rand(n,1) + 0.1;       % Weights
    b     = A * (x0 ./ w);         % Signal
    
    % Run solver for both variants
    opts = spgSetParms('iterations',1000,'verbosity',1);
    AW   = A * spdiags(1./w,0,n,n);
    x = spg_bp(AW, b, opts);
    x1 = x;                        % Reconstruct solution, no weighting
    
    opts = spgSetParms('iterations',1000,'verbosity',1,'weights',w);
    x = spg_bp(A, b, opts);
    x2 = x .* w;                   % Reconstructed solution, with weighting
    
    figure(1); subplot(2,4,5);
    plot(1:n,x1,'m*',1:n,x2,'b', 1:n,x0,'ro');
    legend('Coefficients (1)','Coefficients (2)','Original coefficients');
    title('(e) Weighted Basis Pursuit');
    
    fprintf('\n');
    fprintf([repmat('-',1,35), ' Solution ', repmat('-',1,35), '\n']);
    fprintf('See figure 1(e).\n');
    fprintf([repmat('-',1,80), '\n']);
    fprintf('\n\n');

    
    % -----------------------------------------------------------
    % Solve the multiple measurement vector (MMV) problem
    %
    %    minimize ||Y||_1,2 subject to AW^{-1}Y = B
    %    
    % and the weighted MMV problem (weights on the rows of X):
    %
    %    minimize ||WX||_1,2 subject to AX = B
    %
    % followed by setting Y = WX.
    % -----------------------------------------------------------
    fprintf(['%% ', repmat('-',1,78), '\n']);
    fprintf('%% Solve the multiple measurement vector (MMV) problem      \n');
    fprintf('%%                                                          \n');
    fprintf('%% (1) minimize ||Y||_1,2 subject to AW^{-1}Y = B           \n');
    fprintf('%%                                                          \n');
    fprintf('%% and the weighted MMV problem (weights on the rows of X): \n');
    fprintf('%%                                                          \n');
    fprintf('%% (2) minimize ||WX||_1,2 subject to AX = B                \n');
    fprintf('%%                                                          \n');
    fprintf('%% followed by setting Y = WX.                              \n');
    fprintf(['%% ', repmat('-',1,78), '\n']);

    fprintf('\nPress <return> to continue ... \n');
    if interactive, pause; end

    % Initialize random number generator
    randn('state',0); rand('state',0);
    
    % Create problem
    m = 100; n = 150; k = 12; l = 6;
    A = randn(m,n);
    p = randperm(n); p = p(1:k);
    X0= zeros(n,l); X0(p,:) = randn(k,l);
    
    weights = 3 * rand(n,1) + 0.1;
    W = spdiags(1./weights,0,n,n);
    
    B = A*W*X0;
    
    % Solve unweighted version
    opts = spgSetParms('verbosity',1);
    x    = spg_mmv(A*W,B,0,opts);
    x1   = x;
    
    % Solve weighted version
    opts = spgSetParms('verbosity',1,'weights',weights);
    x    = spg_mmv(A,B,0,opts);
    x2   = spdiags(weights,0,n,n) * x;
    
    % Plot results
    figure(1); subplot(2,4,6);
    plot(x1(:,1),'b-'); hold on;
    plot(x2(:,1),'b.');
    plot(X0,'ro');
    plot(x1(:,2:end),'-');
    plot(x2(:,2:end),'b.');
    legend('Coefficients (1)','Coefficients (2)','Original coefficients');
    title('(f) Weighted Basis Pursuit with Multiple Measurement Vectors');

    fprintf('\n');
    fprintf([repmat('-',1,35), ' Solution ', repmat('-',1,35), '\n']);
    fprintf('See figure 1(f).\n');
    fprintf([repmat('-',1,80), '\n']);
    fprintf('\n\n');

    
    % -----------------------------------------------------------
    % Solve the group-sparse Basis Pursuit problem
    %
    %    minimize    sum_i ||y(group == i)||_2
    %    subject to  AW^{-1}y = b,
    %    
    % with W(i,i) = w(group(i)), and the weighted group-sparse
    % problem
    %
    %    minimize    sum_i w(i)*||x(group == i)||_2
    %    subject to  Ax = b,
    %
    % followed by setting y = Wx.
    % -----------------------------------------------------------
    fprintf(['%% ', repmat('-',1,78), '\n']);
    fprintf('%% Solve the group-sparse Basis Pursuit problem            \n');
    fprintf('%%                                                         \n');
    fprintf('%% (1) minimize    sum_i ||y(group == i)||_2               \n');
    fprintf('%%     subject to  AW^{-1}y = b,                           \n');
    fprintf('%%                                                         \n');
    fprintf('%% with W(i,i) = w(group(i)), and the weighted group-sparse\n');
    fprintf('%% problem                                                 \n');
    fprintf('%%                                                         \n');
    fprintf('%% (2) minimize    sum_i w(i)*||x(group == i)||_2          \n');
    fprintf('%%     subject to  Ax = b,                                 \n');
    fprintf('%%                                                         \n');
    fprintf('%% followed by setting y = Wx.                             \n');
    fprintf(['%% ', repmat('-',1,78), '\n']);

    fprintf('\nPress <return> to continue ... \n');
    if interactive, pause; end

    % Initialize random number generator
    randn('state',0); rand('state',2); % 2
    
    % Set problem size and number of groups
    m = 100; n = 150; nGroups = 25; groups = [];
    
    % Generate groups with desired number of unique groups
    while (length(unique(groups)) ~= nGroups)
       groups  = sort(ceil(rand(n,1) * nGroups)); % Sort for display purpose
    end

    % Determine weight for each group
    weights = 3*rand(nGroups,1) + 0.1;
    W       = spdiags(1./weights(groups),0,n,n);

    % Create sparse vector x0 and observation vector b
    p   = randperm(nGroups); p = p(1:3);
    idx = ismember(groups,p);
    x0  = zeros(n,1); x0(idx) = randn(sum(idx),1);
    b   = A*W*x0;
    
    % Solve unweighted version
    opts = spgSetParms('verbosity',1);
    x    = spg_group(A*W,b,groups,0,opts);
    x1   = x;

    % Solve weighted version
    opts = spgSetParms('verbosity',1,'weights',weights);
    x    = spg_group(A,b,groups,0,opts);
    x2   = spdiags(weights(groups),0,n,n) * x;
    
    % Plot results
    figure(1); subplot(2,4,7);
    plot(x1); hold on;
    plot(x2,'b+');
    plot(x0,'ro'); hold off;
    legend('Coefficients (1)','Coefficients (2)','Original coefficients');
    title('(g) Weighted Group-sparse Basis Pursuit');

    fprintf('\n');
    fprintf([repmat('-',1,35), ' Solution ', repmat('-',1,35), '\n']);
    fprintf('See figure 1(g).\n');
    fprintf([repmat('-',1,80), '\n']);
    
    
end % function demo    

function y = partialFourier(idx,n,x,mode)
    if mode==1
       % y = P(idx) * FFT(x)
       z = fft(x) / sqrt(n);
       y = z(idx);
    else
       z = zeros(n,1);
       z(idx) = x;
       y = ifft(z) * sqrt(n);
    end
 end % function partialFourier
