SPGL1: Spectral Projected Gradient for L1 minimization
------------------------------------------------------

1. Introduction
===============

Thank you for downloading the SPGL1 solver!  SPGL1 is a Matlab solver
for large-scale one-norm regularized least squares.  It is designed to
solve any of the following three problems:

1. Basis pursuit denoise (BPDN):
   minimize  ||x||_1  subject to  ||Ax - b||_2 <= sigma,

2. Basis pursuit (BP):
   minimize   ||x||_1  subject to  Ax = b
 
3. Lasso:
   minimize  ||Ax - b||_2  subject to  ||x||_1 <= tau,

The matrix A can be defined explicily, or as an operator (i.e., a
function) that return both both Ax and A'y.  SPGL1 can solve these
three problems in both the real and complex domains.


Home page: https://www.math.ucdavis.edu/~mpf/spgl1/


2. Quick start
==============

Start Matlab and make sure the working directory is set to the
directory containing the SPGL1 source files. When this is done, run

 >> spgdemo

at the Matlab prompt.  This script illustrates various uses of SPGL1:

- Solve (BPDN) for some sigma > 0
- Solve (Lasso)
- Solve (BP)
- Solve a (BP) problem in complex variables
- Sample the entire Pareto frontier (i.e., ||Ax-b||_2 vs ||x||_1)
  for a small test problem.


3. Installation
===============

3.1  MEX interface
------------------

A vital component of SPGL1 is a routine (oneProjector.m) for
projecting vectors onto the one-norm ball.  The default distribution
includes a pure Matlab version of oneProjector which should work on
all platforms, and also a C-version of this routine that is more
efficient on large problems.  Precompiled MEX interfaces to the C
implementation of oneProjector are included for Windows
(oneProjector.dll), Linux/x86 (oneProjector.mexglx) and MacOSX/Intel
(oneProjector.mexmaci).  If you need to compile the MEX interface on
your own machine, run the following command at the Matlab prompt:

 >> spgsetup

or, equivalently, change to the "private" directory and issue the
command

 >> mex oneProjector.c oneProjector_core.c -output oneProjector -DNDEBUG

If the MEX interface cannot be found, SPGL1 falls back to the slower
Matlab implementation of oneProjector.

3.2  Path
---------

In order to use SPGL1 from any directory other than the one
containing the main spgl1 routine, add the SPGL1 package to your
default path:

 >> addpath <dir-name>

where <dir-name> is the location of spgl1.m.  You can also add this
command to your startup.m file.

4. References
=============

The algorithm implemented by SPGL1 is described in the paper

- E. van den Berg and M. P. Friedlander, "Probing the Pareto frontier
  for basis pursuit solutions", SIAM J. on Scientific Computing,
  31(2):890-912, November 2008

- Sparse optimization with least-squares constraints E. van den Berg
  and M. P. Friedlander, Tech. Rep. TR-2010-02, Dept of Computer
  Science, Univ of British Columbia, January 2010

