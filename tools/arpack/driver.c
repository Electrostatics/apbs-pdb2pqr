/**
 *  @file    driver.c
 *  @author  Nathan Baker
 *  @brief   Drives ARPACK eigenvalue calculations
 *  @version $Id$
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
#include "apbs/vstring.h"
#include "maloc/maloc.h"

#define ERRRC 13

#define F77READHBHEAD VF77_MANGLE(readhbhead, READHBHEAD)
VEXTERNC void F77READHBHEAD(int *nrow, int *ncol, int *nnzero, char *path);
#define F77READHB VF77_MANGLE(readhb, READHB)
VEXTERNC void F77READHB(int *nrow, int *ncol, int *nnzero, double *values, 
                        int *rowind, int *colptr, char *path, char *title,
                        char *key, char *mxtype);
#define F77GETEIGS VF77_MANGLE(geteigs, GETEIGS)
VEXTERNC void F77GETEIGS(int *nrow, int *ncol, int *nnzero, double *values, 
                         int *colptr, int *rowind, int *nev, int *ncv, 
                         int *key, double *v, double *workl, double *workd,
                         double *d, double *resid, double *ax, int *select);

int main(int argc, char **argv) {

    int key, nev, ncv, ncol, nrow, nnzero;
    int *colptr, *rowind, *select;
    double *values, *v, *workl, *workd, *d, *resid, *ax;
    char *keystr, *numstr, path[256], title[72], mxtype[3], hbkey[8];

    char *usage = "\n\n"
                  "This is a program which uses ARPACK to compute\n"
                  "portions of the eigenvalue spectruc of operators.\n\n"
                  "Usage: driver <key> <number> <path>\n\n"
                  "Where <key>    is the portion of the spec:\n"
                  "               lm => largest magnitude\n"
                  "               sm => smallest magnitude\n"
                  "               be => from both ends of spectrum\n"
                  "      <number> is the number of eigenvalues to compute\n"
                  "      <path>   is the path to the matrix in Harwell-\n"
                  "               Boeing format.\n\n";

    Vio_start();

    if (argc != 4) {
        Vnm_print(2, "ERROR:  Got %d arguments, expected 4!", argc);
        Vnm_print(1, "%s", usage);
        return ERRRC;
    } 

    /* Figure out which eigenvalues to compute */
    keystr = argv[1];
    if (Vstring_strcasecmp(keystr, "lm") == 0) {
        key = 0;
        Vnm_print(1, "Computing largest magnitude portion of spectrum.\n");
    } else if (Vstring_strcasecmp(keystr, "sm") == 0) {
        key = 1;
        Vnm_print(1, "Computing smallest magnitude portion of spectrum.\n");
    } else if (Vstring_strcasecmp(keystr, "be") == 0) {
        key = 2;
        Vnm_print(1, "Computing portions from both ends of spectrum.\n");
    } else {
        Vnm_print(2, "ERROR:  Unrecognized key = %s\n", keystr);
        Vnm_print(2, "%s", usage);
        return ERRRC;
    }

    /* Figure out how many eigenvalues to compute */
    numstr = argv[2];
    if (sscanf(numstr, "%d", &nev) == 0) {
        Vnm_print(2, "ERROR:  Unrecognized number of eigenvalues = %s\n", 
          numstr);
        Vnm_print(2, "%s", usage);
        return ERRRC;
    }
    Vnm_print(1, "Computing %d eigenvalue(s).\n", nev);

    /* Read in matrix */
    strncpy(path, argv[3], 255);
    Vnm_print(1, "Reading header from matrix file at %s...\n", path);
    F77READHBHEAD(&nrow, &ncol, &nnzero, path);
    Vnm_print(1, "Matrix %s contains %d rows, %d columns, and %d nonzeros.\n",
      path, nrow, ncol, nnzero);
    Vnm_print(1, "Allocating space for %d integers (colptr)...\n", (ncol+1));
    colptr = Vmem_malloc(VNULL, (ncol+1), sizeof(int));
    Vnm_print(1, "Allocating space for %d integers (rowind)...\n", nnzero);
    rowind = Vmem_malloc(VNULL, nnzero, sizeof(int));
    Vnm_print(1, "Allocating space for %d doubles (values)...\n", nnzero);
    values = Vmem_malloc(VNULL, nnzero, sizeof(double));
    Vnm_print(1, "Reading matrix from matrix file at %s...\n", path);
    F77READHB(&nrow, &ncol, &nnzero, values, rowind, colptr, path, title, 
      hbkey, mxtype);
    Vnm_print(1, "Matrix %s (still) contains %d rows, %d columns, and %d \
nonzeros.\n", path, nrow, ncol, nnzero);
    Vnm_print(1, "Harwell-Boeing matrix MX type = %c%c%c\n", 
      mxtype[0], mxtype[1], mxtype[2]);

    /* Set up work space */
    ncv = VMIN2(nev*2, nrow/2);
    Vnm_print(1, "Using %d Lanczos vectors in calculation.\n", ncv);
    Vnm_print(1, "Allocating space for %d doubles (eigvecs)...\n", 
      nrow*ncv);
    v = Vmem_malloc(VNULL, nrow*ncv, sizeof(double));
    Vnm_print(1, "Allocating space for %d doubles (work)...\n", 
      ncv*(ncv+8));
    workl = Vmem_malloc(VNULL, ncv*(ncv+8), sizeof(double));
    Vnm_print(1, "Allocating space for %d doubles (work)...\n", 
      3*ncol);
    workd = Vmem_malloc(VNULL, 3*ncol, sizeof(double));
    Vnm_print(1, "Allocating space for %d doubles (eigvals)...\n", 
      2*ncv);
    d = Vmem_malloc(VNULL, 2*ncv, sizeof(double));
    Vnm_print(1, "Allocating space for %d doubles (residuals)...\n", 
      nrow);
    resid = Vmem_malloc(VNULL, nrow, sizeof(double));
    Vnm_print(1, "Allocating space for %d doubles (work)...\n", 
      nrow);
    ax = Vmem_malloc(VNULL, nrow, sizeof(double));
    Vnm_print(1, "Allocating space for %d integers (work)...\n", 
      nrow);
    select = Vmem_malloc(VNULL, ncv, sizeof(int));
    Vnm_print(1, "CURRENTLY USING %4.3f MB MEMORY.\n", 
      Vmem_bytesTotal()/1024./1024.);
   
    /* Calculate */ 
    Vnm_print(1, "Calculating %d eigenvalues with ARPACK...\n", nev);
	F77GETEIGS(&nrow, &ncol, &nnzero, values, colptr, rowind, &nev, &ncv, &key,
      v, workl, workd , d, resid, ax, select);

    

    return 0;

}
