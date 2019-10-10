function [x,r,g,info] = spg_group( A, b, groups, sigma, options )
%SPG_GROUP  Solve jointly-sparse basis pursuit denoise (BPDN)
%
%   SPG_GROUP is designed to solve the basis pursuit denoise problem
%
%   (BPDN)  minimize    sum_k ||X_{i : GROUPS(i) = k}||_2
%           subject to  ||A X - B||_2,2 <= SIGMA,
%
%   where A is an M-by-N matrix, B is a vector, GROUPS is vector
%   containing the group number of the corresponding index in X, and
%   SIGMA is a nonnegative scalar.  In all cases below, A can be an
%   explicit M-by-N matrix or matrix-like object for which the
%   operations A*x and A'*y are defined (i.e., matrix-vector
%   multiplication with A and its adjoint.)
%
%   Also, A can be a function handle that points to a function with the
%   signature
%
%   v = A(w,mode)   which returns  v = A *w  if mode == 1;
%                                  v = A'*w  if mode == 2. 
%   
%   X = SPG_GROUP(A,B,G,SIGMA) solves the BPDN problem.  If SIGMA=0,
%   SIGMA=[] or SIGMA is omitted, then the jointly-sparse basis
%   pursuit (BP) problem is solved; i.e., the constraints in the BPDN
%   problem are taken as AX=B.
%
%   X = SPG_GROUP(A,B,G,SIGMA,OPTIONS) specifies options that are set
%   using SPGSETPARMS.
%
%   [X,R,G,INFO] = SPG_GROUP(A,B,GROUPS,SIGMA,OPTIONS) additionally
%   returns the residual R = B - A*X, the objective gradient G = A'*R,
%   and an INFO structure.  (See SPGL1 for a description of this last
%   output argument.)
%
%   See also spgl1, spgSetParms, spg_bp, spg_lasso.

%   Thanks to Aswin Sankaranarayanan for pointing out performance
%   issues with an earlier version of the group preprocessing code.

%   Copyright 2008, Ewout van den Berg and Michael P. Friedlander
%   http://www.cs.ubc.ca/labs/scl/spgl1
%   $Id$

if ~exist('options','var'), options = []; end
if ~exist('sigma','var') || isempty(sigma), sigma = 0; end
if ~exist('groups','var') || isempty(groups)
    error('Third argument cannot be empty.');
end
if ~exist('b','var') || isempty(b)
    error('Second argument cannot be empty.');
end
if ~exist('A','var') || isempty(A)
    error('First argument cannot be empty.');
end

% Preprocess the groups, normalize numbering
g = groups(:);
n = length(g);
[gidx,idx1,idx2] = unique(g);
groups = sparse(idx2,1:n,ones(1,n),length(gidx),n);

% Set projection specific functions
options.project     = @(x,weight,tau) NormGroupL2_project(groups,x,weight,tau);
options.primal_norm = @(x,weight    ) NormGroupL2_primal(groups,x,weight);
options.dual_norm   = @(x,weight    ) NormGroupL2_dual(groups,x,weight);

tau = 0;
x0  = [];
[x,r,g,info] = spgl1(A,b,tau,sigma,x0,options);
