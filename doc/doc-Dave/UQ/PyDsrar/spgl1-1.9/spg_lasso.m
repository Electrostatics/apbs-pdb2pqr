function [x,r,g,info] = spg_lasso(A,b,tau,options )
%SPG_LASSO  Solve the LASSO problem
%
%   SPG_LASSO is designed to solve the LASSO problem
%
%   (LASSO)  minimize  ||AX - B||_2  subject to  ||X||_1 <= tau,
%
%   where A is an M-by-N matrix, B is an M-vector, and TAU is a
%   nonnegative scalar.  In all cases below, A can be an explicit M-by-N
%   matrix or matrix-like object for which the operations  A*x  and  A'*y
%   are defined (i.e., matrix-vector multiplication with A and its
%   adjoint.)
%
%   Also, A can be a function handle that points to a function with the
%   signature
%
%   v = A(w,mode)   which returns  v = A *w  if mode == 1;
%                                  v = A'*w  if mode == 2. 
%   
%   X = SPG_LASSO(A,B,TAU) solves the LASSO problem.
%
%   X = SPG_LASSO(A,B,TAU,OPTIONS) specifies options that are set using
%   SPGSETPARMS.
%
%   [X,R,G,INFO] = SPG_LASSO(A,B,TAU,OPTIONS) additionally returns the
%   residual R = B - A*X, the objective gradient G = A'*R, and an INFO
%   structure.  (See SPGL1 for a description of this last output argument.)
%
%   See also spgl1, spgSetParms, spg_bp, spg_bpdn.

%   Copyright 2008, Ewout van den Berg and Michael P. Friedlander
%   http://www.cs.ubc.ca/labs/scl/spgl1
%   $Id: spg_lasso.m 1074 2008-08-19 05:24:28Z ewout78 $

if ~exist('options','var'), options = []; end
if ~exist('tau','var'), tau = []; end
if ~exist('b','var') || isempty(b)
    error('Second argument cannot be empty.');
end
if ~exist('A','var') || isempty(A)
    error('First argument cannot be empty.');
end

sigma = [];
x0  = [];
[x,r,g,info] = spgl1(A,b,tau,sigma,x0,options);
