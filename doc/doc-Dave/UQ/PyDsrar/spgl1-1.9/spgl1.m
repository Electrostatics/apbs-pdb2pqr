function [x,r,g,info] = spgl1( A, b, tau, sigma, x, options )
%SPGL1  Solve basis pursuit, basis pursuit denoise, and LASSO
%
% [x, r, g, info] = spgl1(A, b, tau, sigma, x0, options)
%
% ---------------------------------------------------------------------
% Solve the basis pursuit denoise (BPDN) problem
%
% (BPDN)   minimize  ||x||_1  subj to  ||Ax-b||_2 <= sigma,
%
% or the l1-regularized least-squares problem
%
% (LASSO)  minimize  ||Ax-b||_2  subj to  ||x||_1 <= tau.
% ---------------------------------------------------------------------
%
% INPUTS
% ======
% A        is an m-by-n matrix, explicit or an operator.
%          If A is a function, then it must have the signature
%
%          y = A(x,mode)   if mode == 1 then y = A x  (y is m-by-1);
%                          if mode == 2 then y = A'x  (y is n-by-1).
%
% b        is an m-vector.
% tau      is a nonnegative scalar; see (LASSO).
% sigma    if sigma != inf or != [], then spgl1 will launch into a
%          root-finding mode to find the tau above that solves (BPDN).
%          In this case, it's STRONGLY recommended that tau = 0.
% x0       is an n-vector estimate of the solution (possibly all
%          zeros). If x0 = [], then SPGL1 determines the length n via
%          n = length( A'b ) and sets  x0 = zeros(n,1).
% options  is a structure of options from spgSetParms. Any unset options
%          are set to their default value; set options=[] to use all
%          default values.
%
% OUTPUTS
% =======
% x        is a solution of the problem
% r        is the residual, r = b - Ax
% g        is the gradient, g = -A'r
% info     is a structure with the following information:
%          .tau     final value of tau (see sigma above)
%          .rNorm   two-norm of the optimal residual
%          .rGap    relative duality gap (an optimality measure)
%          .gNorm   Lagrange multiplier of (LASSO)
%          .stat    = 1 found a BPDN solution
%                   = 2 found a BP sol'n; exit based on small gradient
%                   = 3 found a BP sol'n; exit based on small residual
%                   = 4 found a LASSO solution
%                   = 5 error: too many iterations
%                   = 6 error: linesearch failed
%                   = 7 error: found suboptimal BP solution
%                   = 8 error: too many matrix-vector products
%          .time    total solution time (seconds)
%          .nProdA  number of multiplications with A
%          .nProdAt number of multiplications with A'
%
% OPTIONS
% =======
% Use the options structure to control various aspects of the algorithm:
%
% options.fid         File ID to direct log output
%        .verbosity   0=quiet, 1=some output, 2=more output.
%        .iterations  Max. number of iterations (default if 10*m).
%        .bpTol       Tolerance for identifying a basis pursuit solution.
%        .optTol      Optimality tolerance (default is 1e-4).
%        .decTol      Larger decTol means more frequent Newton updates.
%        .subspaceMin 0=no subspace minimization, 1=subspace minimization.
%
% EXAMPLE
% =======
%   m = 120; n = 512; k = 20; % m rows, n cols, k nonzeros.
%   p = randperm(n); x0 = zeros(n,1); x0(p(1:k)) = sign(randn(k,1));
%   A  = randn(m,n); [Q,R] = qr(A',0);  A = Q';
%   b  = A*x0 + 0.005 * randn(m,1);
%   opts = spgSetParms('optTol',1e-4);
%   [x,r,g,info] = spgl1(A, b, 0, 1e-3, [], opts); % Find BP sol'n.
%
% AUTHORS
% =======
%  Ewout van den Berg (vandenberg.ewout@gmail.com)
%  Michael P. Friedlander (mpf@math.ucdavis.edu)
%
% BUGS
% ====
% Please send bug reports or comments to
%            Michael P. Friedlander (mpf@math.ucdavis.edu)
%            Ewout van den Berg (vandenberg.ewout@gmail.com)

% 15 Apr 07: First version derived from spg.m.
% 17 Apr 07: Added root-finding code.
% 18 Apr 07: sigma was being compared to 1/2 r'r, rather than
%            norm(r), as advertised.  Now immediately change sigma to
%            (1/2)sigma^2, and changed log output accordingly.
% 24 Apr 07: Added quadratic root-finding code as an option.
% 24 Apr 07: Exit conditions need to guard against small ||r||
%            (ie, a BP solution).  Added test1,test2,test3 below.
% 15 May 07: Trigger to update tau is now based on relative difference
%            in objective between consecutive iterations.
% 15 Jul 07: Added code to allow a limited number of line-search
%            errors. 
% 23 Feb 08: Fixed bug in one-norm projection using weights. Thanks
%            to Xiangrui Meng for reporting this bug.
% 26 May 08: The simple call spgl1(A,b) now solves (BPDN) with sigma=0.
% 18 Mar 13: Reset f = fOld if curvilinear line-search fails.
%            Avoid computing the Barzilai-Borwein scaling parameter
%            when both line-search algorithms failed.
% 09 Sep 13: Recompute function information at new x when tau decreases.
%            Fixed bug in subspace minimization. Thanks to Yang Lei
%            for reporting this bug.
    
%   ----------------------------------------------------------------------
%   This file is part of SPGL1 (Spectral Projected-Gradient for L1).
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
REVISION = '1.9';
DATE     = '29 Apr 2015';

tic;              % Start your watches!
m = length(b);

%----------------------------------------------------------------------
% Check arguments. 
%----------------------------------------------------------------------
if ~exist('options','var'), options = []; end
if ~exist('x','var'), x = []; end
if ~exist('sigma','var'), sigma = []; end
if ~exist('tau','var'), tau = []; end

if nargin < 2 || isempty(b) || isempty(A)
   error('At least two arguments are required');
elseif isempty(tau) && isempty(sigma)
   tau = 0;
   sigma = 0;
   singleTau = false;
elseif isempty(sigma) % && ~isempty(tau)  <-- implied
   singleTau = true;
else
   if isempty(tau)
      tau = 0;
   end
   singleTau = false;
end

%----------------------------------------------------------------------
% Grab input options and set defaults where needed. 
%----------------------------------------------------------------------
defaultopts = spgSetParms(...
'fid'        ,      1 , ... % File ID for output
'verbosity'  ,      2 , ... % Verbosity level
'iterations' ,   10*m , ... % Max number of iterations
'nPrevVals'  ,      3 , ... % Number previous func values for linesearch
'bpTol'      ,  1e-06 , ... % Tolerance for basis pursuit solution 
'lsTol'      ,  1e-06 , ... % Least-squares optimality tolerance
'optTol'     ,  1e-04 , ... % Optimality tolerance
'decTol'     ,  1e-04 , ... % Req'd rel. change in primal obj. for Newton
'stepMin'    ,  1e-16 , ... % Minimum spectral step
'stepMax'    ,  1e+05 , ... % Maximum spectral step
'rootMethod' ,      2 , ... % Root finding method: 2=quad,1=linear (not used).
'activeSetIt',    Inf , ... % Exit with EXIT_ACTIVE_SET if nnz same for # its.
'subspaceMin',      0 , ... % Use subspace minimization
'iscomplex'  ,    NaN , ... % Flag set to indicate complex problem
'maxMatvec'  ,    Inf , ... % Maximum matrix-vector multiplies allowed
'weights'    ,      1 , ... % Weights W in ||Wx||_1
'project'    , @NormL1_project , ...
'primal_norm', @NormL1_primal  , ...
'dual_norm'  , @NormL1_dual      ...
   );
options = spgSetParms(defaultopts, options);

fid           = options.fid;
logLevel      = options.verbosity;
maxIts        = options.iterations;
nPrevVals     = options.nPrevVals;
bpTol         = options.bpTol;
lsTol         = options.lsTol;
optTol        = options.optTol;
decTol        = options.decTol;
stepMin       = options.stepMin;
stepMax       = options.stepMax;
activeSetIt   = options.activeSetIt;
subspaceMin   = options.subspaceMin;
maxMatvec     = max(3,options.maxMatvec);
weights       = options.weights;

maxLineErrors = 10;     % Maximum number of line-search failures.
pivTol        = 1e-12;  % Threshold for significant Newton step.

%----------------------------------------------------------------------
% Initialize local variables.
%----------------------------------------------------------------------
iter          = 0;  itnTotLSQR = 0; % Total SPGL1 and LSQR iterations.
nProdA        = 0;  nProdAt    = 0;
lastFv        = -inf(nPrevVals,1);  % Last m function values.
nLineTot      = 0;                  % Total no. of linesearch steps.
printTau      = false;
nNewton       = 0;
bNorm         = norm(b,2);
stat          = false;
timeProject   = 0;
timeMatProd   = 0;
nnzIter       = 0;                  % No. of its with fixed pattern.
nnzIdx        = [];                 % Active-set indicator.
subspace      = false;              % Flag if did subspace min in current itn.
stepG         = 1;                  % Step length for projected gradient.
testUpdateTau = 0;                  % Previous step did not update tau

% Determine initial x, vector length n, and see if problem is complex
explicit = ~(isa(A,'function_handle'));
if isempty(x)
   if isnumeric(A)
      n = size(A,2);
      realx = isreal(A) && isreal(b);
   else
      x = Aprod(b,2);
      n = length(x);
      realx = isreal(x) && isreal(b);
   end
   x = zeros(n,1);
else
   n     = length(x);
   realx = isreal(x) && isreal(b);
end
if isnumeric(A), realx = realx && isreal(A); end;

% Override options when options.iscomplex flag is set
if (~isnan(options.iscomplex)), realx = (options.iscomplex == 0); end

% Check if all weights (if any) are strictly positive. In previous
% versions we also checked if the number of weights was equal to
% n. In the case of multiple measurement vectors, this no longer
% needs to apply, so the check was removed.
if ~isempty(weights)
  if any(~isfinite(weights))
     error('Entries in options.weights must be finite');
  end
  if any(weights <= 0)
     error('Entries in options.weights must be strictly positive');
  end
else
  weights = 1;
end

% Quick exit if sigma >= ||b||.  Set tau = 0 to short-circuit the loop.
if bNorm <= sigma
   printf('W: sigma >= ||b||.  Exact solution is x = 0.\n');
   tau = 0;  singleTau = true;
end 
  
% Don't do subspace minimization if x is complex.
if ~realx && subspaceMin
   printf('W: Subspace minimization disabled when variables are complex.\n');
   subspaceMin = false;
end

% Pre-allocate iteration info vectors
xNorm1 = zeros(min(maxIts,10000),1);
rNorm2 = zeros(min(maxIts,10000),1);
lambda = zeros(min(maxIts,10000),1);

% Exit conditions (constants).
EXIT_ROOT_FOUND    = 1;
EXIT_BPSOL_FOUND   = 2;
EXIT_LEAST_SQUARES = 3;
EXIT_OPTIMAL       = 4;
EXIT_ITERATIONS    = 5;
EXIT_LINE_ERROR    = 6;
EXIT_SUBOPTIMAL_BP = 7;
EXIT_MATVEC_LIMIT  = 8;
EXIT_ACTIVE_SET    = 9; % [sic]

%----------------------------------------------------------------------
% Log header.
%----------------------------------------------------------------------
printf('\n');
printf(' %s\n',repmat('=',1,80));
printf(' SPGL1  v.%s (%s)\n', REVISION, DATE);
printf(' %s\n',repmat('=',1,80));
printf(' %-22s: %8i %4s'   ,'No. rows'          ,m       ,'');
printf(' %-22s: %8i\n'     ,'No. columns'       ,n          );
printf(' %-22s: %8.2e %4s' ,'Initial tau'       ,tau     ,'');
printf(' %-22s: %8.2e\n'   ,'Two-norm of b'     ,bNorm      );
printf(' %-22s: %8.2e %4s' ,'Optimality tol'    ,optTol  ,'');
if singleTau
   printf(' %-22s: %8.2e\n'  ,'Target one-norm of x'  ,tau       );
else
   printf(' %-22s: %8.2e\n','Target objective'  ,sigma      );
end
printf(' %-22s: %8.2e %4s' ,'Basis pursuit tol' ,bpTol   ,'');
printf(' %-22s: %8i\n'     ,'Maximum iterations',maxIts     );
printf('\n');
if singleTau
   logB = ' %5i  %13.7e  %13.7e  %9.2e  %6.1f  %6i  %6i';
   logH = ' %5s  %13s  %13s  %9s  %6s  %6s  %6s\n';
   printf(logH,'Iter','Objective','Relative Gap','gNorm','stepG','nnzX','nnzG');
else
   logB = ' %5i  %13.7e  %13.7e  %9.2e  %9.3e  %6.1f  %6i  %6i';
   logH = ' %5s  %13s  %13s  %9s  %9s  %6s  %6s  %6s  %13s\n';
   printf(logH,'Iter','Objective','Relative Gap','Rel Error',...
          'gNorm','stepG','nnzX','nnzG','tau');
end    
    
% Project the starting point and evaluate function and gradient.
x         = project(x,tau);
r         = b - Aprod(x,1);  % r = b - Ax
g         =   - Aprod(r,2);  % g = -A'r
f         = r'*r / 2; 

% Required for nonmonotone strategy.
lastFv(1) = f;
fBest     = f;
xBest     = x;
fOld      = f;

% Compute projected gradient direction and initial steplength.
dx     = project(x - g, tau) - x;
dxNorm = norm(dx,inf);
if dxNorm < (1 / stepMax)
   gStep = stepMax;
else
   gStep = min( stepMax, max(stepMin, 1/dxNorm) );
end

%----------------------------------------------------------------------
% MAIN LOOP.
%----------------------------------------------------------------------
while 1
    
    %------------------------------------------------------------------
    % Test exit conditions.
    %------------------------------------------------------------------

    % Compute quantities needed for log and exit conditions.
    gNorm   = options.dual_norm(-g,weights);
    rNorm   = norm(r, 2);
    gap     = r'*(r-b) + tau*gNorm;
    rGap    = abs(gap) / max(1,f);
    aError1 = rNorm - sigma;
    aError2 = f - sigma^2 / 2;
    rError1 = abs(aError1) / max(1,rNorm);
    rError2 = abs(aError2) / max(1,f);

    % Count number of consecutive iterations with identical support.
    nnzOld = nnzIdx;
    [nnzX,nnzG,nnzIdx,nnzDiff] = activeVars(x,g,nnzIdx,options);

    if nnzDiff
       nnzIter = 0;
    else
       nnzIter = nnzIter + 1;
       if nnzIter >= activeSetIt, stat=EXIT_ACTIVE_SET; end
    end
    
    % Single tau: Check if we're optimal.
    % The 2nd condition is there to guard against large tau.
    if singleTau
       if rGap <= optTol || rNorm < optTol*bNorm
          stat  = EXIT_OPTIMAL;
       end
 
    % Multiple tau: Check if found root and/or if tau needs updating.
    else
       
       % Test if a least-squares solution has been found
       if gNorm <= lsTol * rNorm
          stat = EXIT_LEAST_SQUARES;
       end
       
       if rGap <= max(optTol, rError2) || rError1 <= optTol
          % The problem is nearly optimal for the current tau.
          % Check optimality of the current root.
          test1 = rNorm       <=   bpTol * bNorm;
       %  test2 = gNorm       <=   bpTol * rNorm;
          test3 = rError1     <=  optTol;
          test4 = rNorm       <=  sigma;
          
          if test4, stat=EXIT_SUBOPTIMAL_BP;end  % Found suboptimal BP sol.
          if test3, stat=EXIT_ROOT_FOUND;   end  % Found approx root.
          if test1, stat=EXIT_BPSOL_FOUND;  end  % Resid minim'zd -> BP sol.
        % 30 Jun 09: Large tau could mean large rGap even near LS sol.
        %            Move LS check out of this if statement.
        % if test2, stat=EXIT_LEAST_SQUARES; end % Gradient zero -> BP sol.
       end
       
       testRelChange1 = (abs(f - fOld) <= decTol * f);
       testRelChange2 = (abs(f - fOld) <= 1e-1 * f * (abs(rNorm - sigma)));
       testUpdateTau  = ((testRelChange1 && rNorm >  2 * sigma) || ...
                         (testRelChange2 && rNorm <= 2 * sigma)) && ...
                         ~stat && ~testUpdateTau;
       
       if testUpdateTau
          % Update tau.
          tauOld   = tau;
          tau      = max(0,tau + (rNorm * aError1) / gNorm);
          nNewton  = nNewton + 1;
          printTau = abs(tauOld - tau) >= 1e-6 * tau; % For log only.
          if tau < tauOld
             % The one-norm ball has decreased.  Need to make sure that the
             % next iterate if feasible, which we do by projecting it.
             x = project(x,tau);
             
             % Update the residual, gradient, and function value.
             r = b - Aprod(x,1);  % r = b - Ax
             g =   - Aprod(r,2);  % g = -A'r
             f = r'*r / 2;
             
             % Reset the function value history.
             lastFv = -inf(nPrevVals,1);  % Last m function values.
             lastFv(1) = f;
          end
       end
    end

    % Too many its and not converged.
    if ~stat  &&  iter >= maxIts
        stat = EXIT_ITERATIONS;
    end

    %------------------------------------------------------------------
    % Print log, update history and act on exit conditions.
    %------------------------------------------------------------------
    if logLevel >= 2 || singleTau || printTau || iter == 0 || stat
       tauFlag = '              '; subFlag = '';
       if printTau, tauFlag = sprintf(' %13.7e',tau);   end
       if subspace, subFlag = sprintf(' S %2i',itnLSQR); end
       if singleTau
          printf(logB,iter,rNorm,rGap,gNorm,log10(stepG),nnzX,nnzG);
          if subspace
             printf('  %s',subFlag);
          end
       else
          printf(logB,iter,rNorm,rGap,rError1,gNorm,log10(stepG),nnzX,nnzG);
          if printTau || subspace
             printf(' %s',[tauFlag subFlag]);
          end
       end
       printf('\n');
    end
    printTau = false;
    subspace = false;
    
    % Update history info
    xNorm1(iter+1) = options.primal_norm(x,weights);
    rNorm2(iter+1) = rNorm;
    lambda(iter+1) = gNorm;
    
    if stat, break; end % Act on exit conditions.
        
    %==================================================================
    % Iterations begin here.
    %==================================================================
    iter = iter + 1;
    xOld = x;  fOld = f;  gOld = g;  rOld = r;

    try
       %---------------------------------------------------------------
       % Projected gradient step and linesearch.
       %---------------------------------------------------------------
       [f,x,r,nLine,stepG,lnErr] = ...
           spgLineCurvy(x,gStep*g,max(lastFv),@Aprod,b,@project,tau);
       nLineTot = nLineTot + nLine;
       if lnErr
          % Projected backtrack failed. Retry with feasible dir'n linesearch.
          x    = xOld;
          f    = fOld;
          r    = rOld;
          dx   = project(x - gStep*g, tau) - x;
          gtd  = g'*dx;
          [f,x,r,nLine,lnErr] = spgLine(f,x,dx,gtd,max(lastFv),@Aprod,b);
          nLineTot = nLineTot + nLine;
       end
       if lnErr
          % Failed again. Revert to previous iterates and damp max BB step.
          x = xOld;
          f = fOld;
          if maxLineErrors <= 0
             stat = EXIT_LINE_ERROR;
          else
             stepMax = stepMax / 10;
             printf(['W: Linesearch failed with error %i. '...
                     'Damping max BB scaling to %6.1e.\n'],lnErr,stepMax);
             maxLineErrors = maxLineErrors - 1;
          end
       end

       %---------------------------------------------------------------
       % Subspace minimization (only if active-set change is small).
       %---------------------------------------------------------------
       doSubspaceMin = false;
       if subspaceMin
          g = - Aprod(r,2);
          [nnzX,nnzG,nnzIdx,nnzDiff] = activeVars(x,g,nnzOld,options);
          if ~nnzDiff
              if nnzX == nnzG, itnMaxLSQR = 20;
              else             itnMaxLSQR = 5;
              end
              nnzIdx = abs(x) >= optTol; 
              doSubspaceMin = true;
          end
       end

       if doSubspaceMin
     
          % LSQR parameters
          damp       = 1e-5;
          aTol       = 1e-1;
          bTol       = 1e-1;
          conLim     = 1e12;
          showLSQR   = 0;
       
          ebar   = sign(x(nnzIdx));
          nebar  = length(ebar);
          Sprod  = @(y,mode)LSQRprod(@Aprod,nnzIdx,ebar,n,y,mode);
       
          [dxbar, istop, itnLSQR] = ...
             lsqr(m,nebar,Sprod,r,damp,aTol,bTol,conLim,itnMaxLSQR,showLSQR);
              
          itnTotLSQR = itnTotLSQR + itnLSQR;
       
          if istop ~= 4  % LSQR iterations successful. Take the subspace step.
             % Push dx back into full space:  dx = Z dx.
             dx = zeros(n,1);
             dx(nnzIdx) = dxbar - (1/nebar)*(ebar'*dxbar)*ebar;

             % Find largest step to a change in sign.
             block1 = nnzIdx  &  x < 0  &  dx > +pivTol;
             block2 = nnzIdx  &  x > 0  &  dx < -pivTol;
             alpha1 = Inf; alpha2 = Inf;
             if any(block1), alpha1 = min(-x(block1) ./ dx(block1)); end
             if any(block2), alpha2 = min(-x(block2) ./ dx(block2)); end
             alpha = min([1  alpha1  alpha2]);
             ensure(alpha >= 0);
             ensure(ebar'*dx(nnzIdx) <= optTol);          
          
             % Update variables.
             x    = x + alpha*dx;
             r    = b - Aprod(x,1);
             f    = r'*r / 2;
             subspace = true;
          end
       end
       
       ensure(options.primal_norm(x,weights) <= tau+optTol);

       %---------------------------------------------------------------
       % Update gradient and compute new Barzilai-Borwein scaling.
       %---------------------------------------------------------------
       if (~lnErr)
          g    = - Aprod(r,2);
          s    = x - xOld;
          y    = g - gOld;
          sts  = s'*s;
          sty  = s'*y;
          if   sty <= 0,  gStep = stepMax;
          else            gStep = min( stepMax, max(stepMin, sts/sty) );
          end
       else
          gStep = min( stepMax, gStep );
       end
       
    catch err % Detect matrix-vector multiply limit error
       if strcmp(err.identifier,'SPGL1:MaximumMatvec')
         stat = EXIT_MATVEC_LIMIT;
         iter = iter - 1;
         x = xOld;  f = fOld;  g = gOld;  r = rOld;
         break;
       else
         rethrow(err);
       end
    end

    %------------------------------------------------------------------
    % Update function history.
    %------------------------------------------------------------------
    if singleTau || f > sigma^2 / 2 % Don't update if superoptimal.
       lastFv(mod(iter,nPrevVals)+1) = f;
       if fBest > f
          fBest = f;
          xBest = x;
       end
    end

end % while 1
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Restore best solution (only if solving single problem).
if singleTau && f > fBest
   rNorm = sqrt(2*fBest);
   printf('\n Restoring best iterate to objective %13.7e\n',rNorm);
   x = xBest;
   r = b - Aprod(x,1);
   g =   - Aprod(r,2);
   gNorm = options.dual_norm(g,weights);
   rNorm = norm(r,  2);
end

% Final cleanup before exit.
info.tau         = tau;
info.rNorm       = rNorm;
info.rGap        = rGap;
info.gNorm       = gNorm;
info.rGap        = rGap;
info.stat        = stat;
info.iter        = iter;
info.nProdA      = nProdA;
info.nProdAt     = nProdAt;
info.nNewton     = nNewton;
info.timeProject = timeProject;
info.timeMatProd = timeMatProd;
info.itnLSQR     = itnTotLSQR;
info.options     = options;
info.timeTotal   = toc;

info.xNorm1      = xNorm1(1:iter);
info.rNorm2      = rNorm2(1:iter);
info.lambda      = lambda(1:iter);

% Print final output.
switch (stat)
   case EXIT_OPTIMAL
      printf('\n EXIT -- Optimal solution found\n')
   case EXIT_ITERATIONS
      printf('\n ERROR EXIT -- Too many iterations\n');
   case EXIT_ROOT_FOUND
      printf('\n EXIT -- Found a root\n');
   case {EXIT_BPSOL_FOUND}
      printf('\n EXIT -- Found a BP solution\n');
   case {EXIT_LEAST_SQUARES}
      printf('\n EXIT -- Found a least-squares solution\n');
   case EXIT_LINE_ERROR
      printf('\n ERROR EXIT -- Linesearch error (%i)\n',lnErr);
   case EXIT_SUBOPTIMAL_BP
      printf('\n EXIT -- Found a suboptimal BP solution\n');
   case EXIT_MATVEC_LIMIT
      printf('\n EXIT -- Maximum matrix-vector operations reached\n');
   case EXIT_ACTIVE_SET
      printf('\n EXIT -- Found a possible active set\n');
   otherwise
      error('Unknown termination condition\n');
end
printf('\n');
printf(' %-20s:  %6i %6s %-20s:  %6.1f\n',...
   'Products with A',nProdA,'','Total time   (secs)',info.timeTotal);
printf(' %-20s:  %6i %6s %-20s:  %6.1f\n',...
   'Products with A''',nProdAt,'','Project time (secs)',timeProject);
printf(' %-20s:  %6i %6s %-20s:  %6.1f\n',...
   'Newton iterations',nNewton,'','Mat-vec time (secs)',timeMatProd);
printf(' %-20s:  %6i %6s %-20s:  %6i\n', ...
   'Line search its',nLineTot,'','Subspace iterations',itnTotLSQR);
printf('\n');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% NESTED FUNCTIONS.  These share some vars with workspace above.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
function z = Aprod(x,mode)
   if (nProdA + nProdAt >= maxMatvec)
     error('SPGL1:MaximumMatvec','');
   end
     
   tStart = toc;
   if mode == 1
      nProdA = nProdA + 1;
      if   explicit, z = A*x;
      else           z = A(x,1);
      end
   elseif mode == 2
      nProdAt = nProdAt + 1;
      if   explicit, z = A'*x;
      else           z = A(x,2);
      end
   else
      error('Wrong mode!');
   end
   timeMatProd = timeMatProd + (toc - tStart);
end % function Aprod

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function printf(varargin)
  if logLevel > 0
     fprintf(fid,varargin{:});
  end
end % function printf

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function x = project(x, tau)
   tStart      = toc;

   x = options.project(x,weights,tau);
   
   timeProject = timeProject + (toc - tStart);
end % function project

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% End of nested functions.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end % function spg

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PRIVATE FUNCTIONS.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [nnzX,nnzG,nnzIdx,nnzDiff] = activeVars(x,g,nnzIdx,options)
% Find the current active set.
% nnzX    is the number of nonzero x.
% nnzG    is the number of elements in nnzIdx.
% nnzIdx  is a vector of primal/dual indicators.
% nnzDiff is the no. of elements that changed in the support.
  xTol    = min(.1,10*options.optTol);
  gTol    = min(.1,10*options.optTol);
  gNorm   = options.dual_norm(g,options.weights);
  nnzOld  = nnzIdx;

  % Reduced costs for postive & negative parts of x.
  z1 = gNorm + g;
  z2 = gNorm - g;

  % Primal/dual based indicators.
  xPos    = x >  xTol  &  z1 < gTol; %g < gTol;%
  xNeg    = x < -xTol  &  z2 < gTol; %g > gTol;%
  nnzIdx  = xPos | xNeg;

  % Count is based on simple primal indicator.
  nnzX    = sum(abs(x) >= xTol);
  nnzG    = sum(nnzIdx);
  
  if isempty(nnzOld)
     nnzDiff = inf;
  else
     nnzDiff = sum(nnzIdx ~= nnzOld);
  end
  
end % function spgActiveVars

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function z = LSQRprod(Aprod,nnzIdx,ebar,n,dx,mode)
% Matrix multiplication for subspace minimization.
% Only called by LSQR.
  nbar = length(ebar);
   if mode == 1
      y = zeros(n,1);
      y(nnzIdx) = dx - (1/nbar)*(ebar'*dx)*ebar; % y(nnzIdx) = Z*dx
      z = Aprod(y,1);                            % z = S Z dx
   else
      y = Aprod(dx,2);
      z = y(nnzIdx) - (1/nbar)*(ebar'*y(nnzIdx))*ebar;
   end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [fNew,xNew,rNew,iter,err] = spgLine(f,x,d,gtd,fMax,Aprod,b)
% Nonmonotone linesearch.

EXIT_CONVERGED  = 0;
EXIT_ITERATIONS = 1;
maxIts = 10;
step   = 1;
iter   = 0;
gamma  = 1e-4;
gtd    = -abs(gtd); % 03 Aug 07: If gtd is complex,
                    % then should be looking at -abs(gtd).
while 1

    % Evaluate trial point and function value.
    xNew = x + step*d;
    rNew = b - Aprod(xNew,1);
    fNew = rNew'*rNew / 2;

    % Check exit conditions.
    if fNew < fMax + gamma*step*gtd  % Sufficient descent condition.
       err = EXIT_CONVERGED;
       break
    elseif  iter >= maxIts           % Too many linesearch iterations.
       err = EXIT_ITERATIONS;
       break
    end
    
    % New linesearch iteration.
    iter = iter + 1;
    
    % Safeguarded quadratic interpolation.
    if step <= 0.1
       step  = step / 2;
    else
       tmp = (-gtd*step^2) / (2*(fNew-f-step*gtd));
       if tmp < 0.1 || tmp > 0.9*step || isnan(tmp)
          tmp = step / 2;
       end
       step = tmp;
    end
    
end % while 1

end % function spgLine

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [fNew,xNew,rNew,iter,step,err] = ...
    spgLineCurvy(x,g,fMax,Aprod,b,project,tau)
% Projected backtracking linesearch.
% On entry,
% g  is the (possibly scaled) steepest descent direction.

EXIT_CONVERGED  = 0;
EXIT_ITERATIONS = 1;
EXIT_NODESCENT  = 2;
gamma  = 1e-4;
maxIts = 10;
step   =  1;
sNorm  =  0;
scale  =  1;      % Safeguard scaling.  (See below.)
nSafe  =  0;      % No. of safeguarding steps.
iter   =  0;
debug  =  false;  % Set to true to enable log.
n      =  length(x);

if debug
   fprintf(' %5s  %13s  %13s  %13s  %8s\n',...
           'LSits','fNew','step','gts','scale');   %#ok<UNRCH>
end
   
while 1

    % Evaluate trial point and function value.
    xNew     = project(x - step*scale*g, tau);
    rNew     = b - Aprod(xNew,1);
    fNew     = rNew'*rNew / 2;
    s        = xNew - x;
    gts      = scale * real(g' * s);
%   gts      = scale * (g' * s);
    if gts >= 0
       err = EXIT_NODESCENT;
       break
    end

    if debug
       fprintf(' LS %2i  %13.7e  %13.7e  %13.6e  %8.1e\n',...
               iter,fNew,step,gts,scale);   %#ok<UNRCH>
    end
    
    % 03 Aug 07: If gts is complex, then should be looking at -abs(gts).
    % 13 Jul 11: It's enough to use real part of g's (see above).
    if fNew < fMax + gamma*step*gts
%   if fNew < fMax - gamma*step*abs(gts)  % Sufficient descent condition.
       err = EXIT_CONVERGED;
       break
    elseif iter >= maxIts                 % Too many linesearch iterations.
       err = EXIT_ITERATIONS;
       break
    end
    
    % New linesearch iteration.
    iter = iter + 1;
    step = step / 2;

    % Safeguard: If stepMax is huge, then even damped search
    % directions can give exactly the same point after projection.  If
    % we observe this in adjacent iterations, we drastically damp the
    % next search direction.
    % 31 May 07: Damp consecutive safeguarding steps.
    sNormOld  = sNorm;
    sNorm     = norm(s) / sqrt(n);
    %   if sNorm >= sNormOld
    if abs(sNorm - sNormOld) <= 1e-6 * sNorm
       gNorm = norm(g) / sqrt(n);
       scale = sNorm / gNorm / (2^nSafe);
       nSafe = nSafe + 1;
    end
    
end % while 1

end % function spgLineCurvy
