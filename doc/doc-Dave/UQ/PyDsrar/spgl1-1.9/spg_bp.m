function [x,r,g,info] = spg_bp(A,b,options )
%SPG_BP  Solve the basis pursuit (BP) problem
%
%   SPG_BP is designed to solve the basis pursuit problem
%
%   (BP)  minimize  ||X||_1  subject to  AX = B,
%
%   where A is an M-by-N matrix, B is an M-vector, and SIGMA is a
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
%   X = SPG_BP(A,B) solves the BP problem.
%
%   X = SPG_BP(A,B,OPTIONS) specifies options that are set using
%   SPGSETPARMS.
%
%   [X,R,G,INFO] = SPG_BP(A,B,OPTIONS) additionally returns the
%   residual R = B - A*X (which should be small), the objective gradient G
%   = A'*R, and an INFO structure.  (See SPGL1 for a description of this
%   last output argument.)
%
%   See also spgl1, spgSetParms, spg_bpdn, spg_lasso.

%   Copyright 2008, Ewout van den Berg and Michael P. Friedlander
%   http://www.cs.ubc.ca/labs/scl/spgl1
%   $Id: spg_bp.m 1074 2008-08-19 05:24:28Z ewout78 $

if ~exist('options','var'), options = []; end
if ~exist('b','var') || isempty(b)
    error('Second argument cannot be empty.');
end
if ~exist('A','var') || isempty(A)
    error('First argument cannot be empty.');
end

sigma = 0;
tau = 0;
x0  = [];
[x,r,g,info] = spgl1(A,b,tau,sigma,x0,options);
