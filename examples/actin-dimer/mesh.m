%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Dimension:     3
% Domain:        Unit cube divided into six tetrahedra.
% Boundary:      Four dirichlet bndry nodes (one cube face)
%
% rcsid="$Id$"
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% APBS -- Adaptive Poisson-Boltzmann Solver
%
%  Nathan A. Baker (nbaker@wasabi.ucsd.edu)
%  Dept. of Chemistry and Biochemistry
%  Dept. of Mathematics, Scientific Computing Group
%  University of California, San Diego 
%
%  Additional contributing authors listed in the code documentation.
%
% Copyright © 1999. The Regents of the University of California (Regents).
% All Rights Reserved. 
% 
% Permission to use, copy, modify, and distribute this software and its
% documentation for educational, research, and not-for-profit purposes,
% without fee and without a signed licensing agreement, is hereby granted,
% provided that the above copyright notice, this paragraph and the
% following two paragraphs appear in all copies, modifications, and
% distributions.
% 
% IN NO EVENT SHALL REGENTS BE LIABLE TO ANY PARTY FOR DIRECT, INDIRECT,
% SPECIAL, INCIDENTAL, OR CONSEQUENTIAL DAMAGES, INCLUDING LOST PROFITS,
% ARISING OUT OF THE USE OF THIS SOFTWARE AND ITS DOCUMENTATION, EVEN IF
% REGENTS HAS BEEN ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.  
% 
% REGENTS SPECIFICALLY DISCLAIMS ANY WARRANTIES, INCLUDING, BUT NOT
% LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A
% PARTICULAR PURPOSE.  THE SOFTWARE AND ACCOMPANYING DOCUMENTATION, IF
% ANY, PROVIDED HEREUNDER IS PROVIDED "AS IS".  REGENTS HAS NO OBLIGATION
% TO PROVIDE MAINTENANCE, SUPPORT, UPDATES, ENHANCEMENTS, OR
% MODIFICATIONS. 
%%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

mcsf_begin=1;

      dim=3;         % intrinsic manifold dimension
    dimii=3;         % imbedding manifold dimension
 vertices=8;         % number of vertices
simplices=6;         % number of simplices

vert=[
%-------- ---- ----------------- ----------------- -----------------
% Vert-ID Chrt X-Coordinate      Y-Coordinate      Z-Coordinate
%-------- ---- ----------------- ----------------- -----------------
0         0    -1.0000000000e+02 -1.0000000000e+02 -1.0000000000e+02
1         0     1.0000000000e+02 -1.0000000000e+02 -1.0000000000e+02
2         0    -1.0000000000e+02  1.0000000000e+02 -1.0000000000e+02
3         0     1.0000000000e+02  1.0000000000e+02 -1.0000000000e+02
4         0    -1.0000000000e+02 -1.0000000000e+02  1.0000000000e+02
5         0     1.0000000000e+02 -1.0000000000e+02  1.0000000000e+02
6         0    -1.0000000000e+02  1.0000000000e+02  1.0000000000e+02
7         0     1.0000000000e+02  1.0000000000e+02  1.0000000000e+02
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

