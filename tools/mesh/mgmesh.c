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
 * Copyright (c) 2002-2009, Washington University in St. Louis.
 * Portions Copyright (c) 2002-2009.  Nathan A. Baker
 * Portions Copyright (c) 1999-2002.  The Regents of the University of California.
 * Portions Copyright (c) 1995.  Michael Holst
 *
 * All rights reserved.
 * 
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are met: 
 *
 * -  Redistributions of source code must retain the above copyright notice, this
 * list of conditions and the following disclaimer.  
 * 
 * - Redistributions in binary form must reproduce the above copyright notice,
 * this list of conditions and the following disclaimer in the documentation
 * and/or other materials provided with the distribution.
 * 
 * - Neither the name of Washington University in St. Louis nor the names of its
 * contributors may be used to endorse or promote products derived from this
 * software without specific prior written permission.
 * 
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
 * "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
 * LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
 * A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR
 * CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
 * EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
 * PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
 * PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
 * LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
 * NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
 * SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE. *
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
