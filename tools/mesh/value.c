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
	
	double x[81];
	double y[5][81];
	
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
	int i,k,s;
	double j,m;
	char *files[] = {"f1f94_para0.dx","f1f94_para1.dx","f1f94_para2.dx","f1f94_para3.dx"};
	for(s=0;s<4;s++){
		grid = Vgrid_ctor(0, 0, 0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, VNULL);
		if (!Vgrid_readDX(grid, "FILE", "ASC", VNULL, files[s])) {
			Vnm_print(2, "main:  Problem reading OpenDX-format grid from %s\n",
					  path);
			return 2;
		}
		double offset[] = {-30.45781,-15.01276,0.4322926,15.87734};
		for(k=0,m=0;k<5;k++,m+=0.17385){
			for(i=0,j=0.0;i<81;i++,j+=0.19029){
				pt[0] = offset[s] + j;
				pt[1] = -17.33536 + m;
				pt[2] = -22.732;
				Vgrid_value(grid, pt, &val);
				x[i] = pt[0];
				y[k][i] = val;
			}
		}
		for(i=0;i<81;i++){
			printf("%f %f %f %f %f %f\n",x[i],y[0][i],y[1][i],y[2][i],y[3][i],y[4][i]);
		}
	}
	

    return 0;

}
