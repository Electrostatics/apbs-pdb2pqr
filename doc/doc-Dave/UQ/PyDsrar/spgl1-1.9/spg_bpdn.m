function [x,r,g,info] = spg_bpdn( A, b, sigma, options )
%SPG_BPDN  Solve the basis pursuit denoise (BPDN) problem
%
%   SPG_BPDN is designed to solve the basis pursuit denoise problem
%
%   (BPDN)  minimize  ||X||_1  subject to  ||A X - B|| <= SIGMA,
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
%   X = SPG_BPDN(A,B,SIGMA) solves the BPDN problem.  If SIGMA=0 or
%   SIGMA=[], then the basis pursuit (BP) problem is solved; i.e., the
%   constraints in the BPDN problem are taken as AX=B.
%
%   X = SPG_BPDN(A,B,SIGMA,OPTIONS) specifies options that are set using
%   SPGSETPARMS.
%
%   [X,R,G,INFO] = SPG_BPDN(A,B,SIGMA,OPTIONS) additionally returns the
%   residual R = B - A*X, the objective gradient G = A'*R, and an INFO
%   structure.  (See SPGL1 for a description of this last output argument.)
%
%   See also spgl1, spgSetParms, spg_bp, spg_lasso.

%   Copyright 2008, Ewout van den Berg and Michael P. Friedlander
%   http://www.cs.ubc.ca/labs/scl/spgl1
%   $Id: spg_bpdn.m 1389 2009-05-29 18:32:33Z mpf $

if ~exist('options','var'), options = []; end
if ~exist('sigma','var'), sigma = 0; end
if ~exist('b','var') || isempty(b)
    error('Second argument cannot be empty.');
end
if ~exist('A','var') || isempty(A)
    error('First argument cannot be empty.');
end

tau = 0;
x0  = [];
[x,r,g,info] = spgl1(A,b,tau,sigma,x0,options);
