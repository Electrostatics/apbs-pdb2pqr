/**
 *  @file    mgmesh.c
 *  @author  Nathan Baker
 *  @brief   Small program to write out acceptable combinations of grid
 *           dimensions, multigrid levels, etc. for PMG
 *  @version $Id$
 *
 *  @attention
 *  @verbatim
 *
 * APBS -- Adaptive Poisson-Boltzmann Solver
 *
 * Nathan A. Baker (baker@biochem.wustl.edu)
 * Dept. of Biochemistry and Molecular Biophysics
 * Center for Computational Biology
 * Washington University in St. Louis
 *
 * Additional contributing authors listed in the code documentation.
 *
 * Copyright (c) 2003.  Washington University in St. Louis.
 * All Rights Reserved.
 * Portions Copyright (c) 1999-2003.  The Regents of the University of
 * California.  
 * Portions Copyright (c) 1995.  Michael Holst.
 *
 * This file is part of APBS.
 *
 * APBS is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * APBS is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with APBS; if not, write to the Free Software
 * Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307  USA
 *
 * @endverbatim
 */

#include "apbscfg.h"
#include "maloc/maloc.h"


int main(int argc, char **argv) {

    int i, lev;
    int maxvert = 700;
    int minlev = 3;
    double newval, oldval;

    Vio_start();

    Vnm_print(1, "\n\nThis program determines the acceptable meshpoint number\n"
                 "and level combinations for the PMG multigrid libraries and\n"
                 "%d or more levels in the mesh (because you typically use\n"
                 "one less than the max number of levels)\n\n\n", minlev);

    for (i=2; i<maxvert; i++) { 
        /* the number of times it's divisible. */
        lev = 0;
        newval = (double)(i-1);
        oldval = (double)(i-1);
        while (1) {
           oldval = newval;
           newval = newval/2.0;
           if ((floor(newval) != newval) || (ceil(newval) != newval)) break;
           lev++;
        } 
        if (lev >= minlev) {
            Vnm_print(1, "%4d verts/direction => %d levels\n", i, lev);
            Vnm_print(1, "                        %d verts on coarsest level\n",
              (int)oldval); 
            Vnm_print(1, "                        ~%g MB memory (for %d^3 mesh)\n",
              (double)(i*i*i)*160.0/1024.0/1024.0, i);
        }
    }

#if 0
    int i, maxmult = 30;
    int j, maxlev = 5;
    double log2 = log(2.0);
    double x;

    for (i=0; i<maxlev; i++) {
        for (j=0; j<maxmult; j++) {
            printf("%g ", j*pow(2,i+1) + 1);
        }
        printf("\n");
    }
#endif
   
    return 0;
}
