function [x, itn] = oneProjectorMex(b,d,tau)
% [x, itn] = oneProjectorMex(b,d,tau) 
% Return the orthogonal projection of the vector b >=0 onto the
% (weighted) L1 ball. In case vector d is specified, matrix D is
% defined as diag(d), otherwise the identity matrix is used.
%
% On exit,
% x      solves   minimize  ||b-x||_2  st  ||Dx||_1 <= tau.
% itn    is the number of elements of b that were thresholded.
%
% See also spgl1, oneProjector.

%   oneProjectorMex.m
%   $Id: oneProjectorMex.m 1200 2008-11-21 19:58:28Z mpf $
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

if nargin < 3
   tau = d;
   d   = 1;
end

if isscalar(d)
  [x,itn] = oneProjectorMex_I(b,tau/abs(d));
else
  [x,itn] = oneProjectorMex_D(b,d,tau);
end

end % function oneProjectorMex


% ----------------------------------------------------------------------
function [x,itn] = oneProjectorMex_I(b,tau)
% ----------------------------------------------------------------------

   % Initialization
   n     = length(b);
   x     = zeros(n,1);
   bNorm = norm(b,1);

   % Check for quick exit.
   if (tau >= bNorm), x = b; itn = 0; return; end
   if (tau <  eps  ),        itn = 0; return; end

   % Preprocessing (b is assumed to be >= 0)
   [b,idx] = sort(b,'descend'); % Descending.

   csb       = -tau;
   alphaPrev = 0;
   for j= 1:n
      csb       = csb + b(j);
      alpha     = csb / j;
   
      % We are done as soon as the constraint can be satisfied
      % without exceeding the current minimum value of b
      if alpha >= b(j)
         break;
      end
   
      alphaPrev = alpha;
   end

   % Set the solution by applying soft-thresholding with
   % the previous value of alpha
   x(idx) = max(0,b - alphaPrev);

   % Set number of iterations
   itn = j;
end


% ----------------------------------------------------------------------
function [x,itn] = oneProjectorMex_D(b,d,tau)
% ----------------------------------------------------------------------

   % Initialization
   n = length(b);
   x = zeros(n,1);

   % Check for quick exit.
   if (tau >= norm(d.*b,1)), x = b; itn = 0; return; end
   if (tau <  eps         ),        itn = 0; return; end

   % Preprocessing (b is assumed to be >= 0)
   [bd,idx] = sort(b ./ d,'descend'); % Descending.
   b  = b(idx);
   d  = d(idx);

   % Optimize
   csdb = 0; csd2 = 0;
   soft = 0; alpha1 = 0; i = 1;
   while (i <= n)
      csdb = csdb + d(i).*b(i);
      csd2 = csd2 + d(i).*d(i);
  
      alpha1 = (csdb - tau) / csd2;
      alpha2 = bd(i);

      if alpha1 >= alpha2
         break;
      end
    
      soft = alpha1;  i = i + 1;
   end
   x(idx(1:i-1)) = b(1:i-1) - d(1:i-1) * max(0,soft);

   % Set number of iterations
   itn = i;
end
