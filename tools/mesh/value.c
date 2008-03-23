/**
*  @file    value.c
 *  @author  Nathan Baker
 *  @brief   Get information about solution at a point
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
 * Copyright (c) 2002-2007, Washington University in St. Louis.
 * Portions Copyright (c) 2002-2007.  Nathan A. Baker
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
#include "apbs/apbs.h"  

VEMBED(rcsid="$Id$")

int usage(int rc) {
	
    char *usage = "\n\n\
----------------------------------------------------------------------\n\
    This driver program reads in data and prints solution information at a\n\
    point.  It is invoked as:\n\
	value <x> <y> <z> <file.dx>\n\n\
    where <x>, <y>, and <z> are points and <file.dx> is an OpenDX-format\n\
    file.\n\
----------------------------------------------------------------------\n\n";

    Vnm_print(2, usage);

    exit(rc);

    return 0;
}

int main(int argc, char **argv) {
	
    Vgrid *grid;
    int inorm;
    char *path;
    double pt[3], val, grad[3];
	
    /* *************** CHECK INVOCATION ******************* */
    Vio_start();
    Vnm_redirect(1);
    Vnm_print(1, "\n");
    if (argc != 5) {
        Vnm_print(2, "Error -- got %d arguments, expected 5.\n", argc);
        usage(2);
    }
    sscanf(argv[1], "%lf", &(pt[0]));
    sscanf(argv[2], "%lf", &(pt[1]));
    sscanf(argv[3], "%lf", &(pt[2]));
    path = argv[4];
	
    /* *************** READ DATA ******************* */
    Vnm_print(1, "Reading data from %s...\n", path);
    grid = Vgrid_ctor(0, 0, 0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, VNULL);
    if (!Vgrid_readDX(grid, "FILE", "ASC", VNULL, path)) {
        Vnm_print(2, "main:  Problem reading OpenDX-format grid from %s\n",
				  path);
        return 2;
    }
	
    /* *************** READ DATA ******************* */
    Vnm_print(1, "\nData at (%g, %g, %g):\n", pt[0], pt[1], pt[2]);
    if (Vgrid_value(grid, pt, &val)) {
        Vnm_print(1, "Value = %1.12E kT/e\n", val);
    } else  Vnm_print(1, "Unable to get value.\n");
    if (Vgrid_gradient(grid, pt, grad)) {
        Vnm_print(1, "Gradient = (%1.12E, %1.12E, %1.12E) kT/e/A\n", 
				  grad[0], grad[1], grad[2]);
    } else  Vnm_print(1, "Unable to get gradient.\n");
    /* 
	if (Vgrid_curvature(grid, pt, 0, &val)) {
		Vnm_print(1, "Reduced maximal curvature = %1.12E kT/e/A/A\n", val);
	} else Vnm_print(1, "Unable to get curvature.\n");
	if (Vgrid_curvature(grid, pt, 1, &val)) {
		Vnm_print(1, "Mean curvature (Laplace) = %1.12E kT/e/A/A\n", val);
	} else Vnm_print(1, "Unable to get curvature.\n");
	if (Vgrid_curvature(grid, pt, 2, &val)) {
		Vnm_print(1, "Gauss curvature = %1.12E kT/e/A/A\n", val);
	} else Vnm_print(1, "Unable to get curvature.\n");
	if (Vgrid_curvature(grid, pt, 3, &val)) {
		Vnm_print(1, "True maximal curvature = %1.12E kT/e/A/A\n", val);
	} else Vnm_print(1, "Unable to get curvature.\n");
	*/
	
    Vnm_print(1, "\n");
    return 0;
	
}
