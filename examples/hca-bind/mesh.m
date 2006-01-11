%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Dimension:     3
% Domain:        Unit cube divided into six tetrahedra.
% Boundary:      Four dirichlet bndry nodes (one cube face)
%
% rcsid="$Id$"
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% APBS -- Adaptive Poisson-Boltzmann Solver
% 
%   Nathan A. Baker (baker@biochem.wustl.edu)
%   Dept. Biochemistry and Molecular Biophysics
%   Center for Computational Biology
%   Washington University in St. Louis
% 
%   Additional contributing authors listed in the code documentation.
% 
% Copyright (c) 2002-2005.  Washington University in St. Louis.
% All Rights Reserved.
% Portions Copyright (c) 1999-2002.  The Regents of the University of 
% California.
% Portions Copyright (c) 1995.  Michael Holst.
% 
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 2 of the License, or
% (at your option) any later version.
% 
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
% 
% You should have received a copy of the GNU General Public License
% along with this program; if not, write to the Free Software
% Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
%
% Linking APBS statically or dynamically with other modules is making a
% combined work based on APBS. Thus, the terms and conditions of the GNU
% General Public License cover the whole combination.
% 
% SPECIAL GPL EXCEPTION
% In addition, as a special exception, the copyright holders of APBS
% give you permission to combine the APBS program with free software
% programs or libraries that are released under the GNU LGPL and with
% code included in releases of ISIM, PMV, PyMOL, SMOL, VMD.  This
% special exception permission is also extended to any software listed
% in the SPECIAL GPL EXCEPTION clauses by the PMG, FEtk, MC, or MALOC
% libraries.
% 
% Note that people who make modified versions of APBS are not obligated
% to grant this special exception for their modified versions; it is
% their choice whether to do so. The GNU General Public License gives
% permission to release a modified version without this exception; this
% exception also makes it possible to release a modified version which
% carries forward this exception.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

mcsf_begin=1;

      dim=3;         % intrinsic manifold dimension
    dimii=3;         % imbedding manifold dimension
 vertices=8;         % number of vertices
simplices=6;         % number of simplices

vert=[
%-------- ---- ----------------- ----------------- -----------------
% Vert-ID Chrt X-Coordinate      Y-Coordinate      Z-Coordinate
%-------- ---- ----------------- ----------------- -----------------
0         0    -4.5000000000e+01 -4.5000000000e+01 -4.5000000000e+01
1         0     4.5000000000e+01 -4.5000000000e+01 -4.5000000000e+01
2         0    -4.5000000000e+01  4.5000000000e+01 -4.5000000000e+01
3         0     4.5000000000e+01  4.5000000000e+01 -4.5000000000e+01
4         0    -4.5000000000e+01 -4.5000000000e+01  4.5000000000e+01
5         0     4.5000000000e+01 -4.5000000000e+01  4.5000000000e+01
6         0    -4.5000000000e+01  4.5000000000e+01  4.5000000000e+01
7         0     4.5000000000e+01  4.5000000000e+01  4.5000000000e+01
];

simp=[
%-------- ---- ---- ------------------- ---------------------------------------
% Simp-ID Grp  Mat  Face-Types          Vertex-Numbers
%-------- ---- ---- ------------------- ---------------------------------------
0         0    0    0    1    0    1    0 5 1 2
1         0    0    0    1    1    0    0 5 2 4
2         0    0    0    1    0    1    1 5 3 2
3         0    0    0    1    0    1    3 5 7 2
4         0    0    1    1    0    0    2 5 7 6
5         0    0    1    1    0    0    2 5 6 4
];

mcsf_end=1;

