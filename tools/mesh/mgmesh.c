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
 * Copyright (c) 2002-2006.  Washington University in St. Louis.
 * All Rights Reserved.
 * Portions Copyright (c) 1999-2002.  The Regents of the University of
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
 * Linking APBS statically or dynamically with other modules is making a
 * combined work based on APBS. Thus, the terms and conditions of the GNU
 * General Public License cover the whole combination.
 * 
 * SPECIAL GPL EXCEPTION
 * In addition, as a special exception, the copyright holders of APBS
 * give you permission to combine the APBS program with free software
 * programs and libraries that are released under the GNU LGPL or with
 * code included in releases of ISIM, Ion Simulator Interface, PMV, PyMOL
 * SMOL, VMD, and Vision. Such combined software may be linked with APBS and 
 * redistributed together in original or modified form as mere aggregation
 * without requirement that the entire work be under the scope of the GNU 
 * General Public License. This special exception permission is also extended
 * to any software listed in the SPECIAL GPL EXCEPTION clauses by the PMG,
 * FEtk, MC, or MALOC libraries.
 * 
 * Note that people who make modified versions of APBS are not obligated
 * to grant this special exception for their modified versions; it is
 * their choice whether to do so. The GNU General Public License gives
 * permission to release a modified version without this exception; this
 * exception also makes it possible to release a modified version which
 * carries forward this exception.
 *
 * @endverbatim
 */

#include "apbscfg.h"
#include "maloc/maloc.h"
#include "apbs/apbs.h"


int main(int argc, char **argv) {

    int i, lev;
    int maxvert = 700;
    int minlev = VMGNLEV;
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
