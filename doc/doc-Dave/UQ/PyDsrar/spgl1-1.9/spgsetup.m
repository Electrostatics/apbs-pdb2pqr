%SETUP  Compile the SPGL1 MEX interfaces

%   setup.m
%   $Id: spgsetup.m 1019 2008-06-17 19:33:48Z mpf $
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

root = pwd;
try
    cd('private')
    mex oneProjectorMex.c oneProjectorCore.c heap.c -output oneProjectorMex -DNDEBUG
    fprintf('Successfully compiled oneProjector.\n');
    cd(root)
catch
    cd(root)
    fprintf('Could not compile oneProjector.');
    fprintf('You can still use the slower ".m" version.');
    rethrow(lasterr);
end
    
