/**
 *  @file    acc.c
 *  @author  Nathan Baker
 *  @brief   Small program to calculate volumes, areas, etc. of molecules
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

#include <sys/times.h>
#include <unistd.h>

#include "apbscfg.h"
#include "apbs/vatom.h"
#include "apbs/valist.h"
#include "apbs/vacc.h"

int main(int argc, char **argv) {

    /* OBJECTS */
    Valist *alist;
    Vacc  *acc;
    Vatom *atom;

    /* VACC PARAMETERS */
    double probe_radius = 1.4;
    int nsphere = 200;
    int nx = 50;
    int ny = 50;
    int nz = 50;
    int method;

    /* SYSTEM PARAMETERS */
    char *path;

    /* QUADRATURE VARIABLES */
    double vdwVol = 0.0;
    double ivdwVol = 0.0;
    double molVol = 0.0;
    double sasa = 0.0;
    double atom_sasa = 0.0;
    /* Quadrature steps */
    double dx, dy, dz;
    double x_len, y_len, z_len;
    double x_cen, y_cen, z_cen;
    double x, y, z;
    double w = 1.0;
    double vec[3];
    int i;

    /* TIMING VARIABLES */
    double t;
    clock_t start_t;
    clock_t stop_t;
    struct tms tmst;
    long clktck = sysconf(_SC_CLK_TCK);
 
    char *usage = "\n\nUsage: acc <-vdw|-ivdw|-mol|-sasa> <molecule.pqr>"
                  "\tWhere -vdw   => volume inside VdW surface\n"
                  "\t\t-ivdw  => volume inside VdW surface inflated by probe\
 radius\n"
                  "\t\t-mol   => volume inside molecular surface\n"
                  "\t\t-sasa  => unweighted solvent-accessible surface area\n"
                  "\t\t-wsasa => weighted solvent-accessible surface area\n"
                  "\t\t          (uses \"charge\" column of PQR format for\
 weights)\n"
                  "\t\t-vsasa => verbose unweighted solvent-accessible surface\
 area (prints\n"
                  "\t\t          per-atom contributions for custom-weighted\
 calculations)\n"
                  "\tand <molecule.pqr> is the path to the molecule structure\
 in PQR format\n\n";

    if (argc != 3) {
        Vnm_print(1, "\n*** Syntax error: got %d arguments, expected 2.\n",
          argc);
        Vnm_print(1, "%s", usage);
        return 13;
    };
    if (strcmp(argv[1], "-vdw") == 0) method = 0;
    else if (strcmp(argv[1], "-ivdw") == 0) method = 1;
    else if (strcmp(argv[1], "-mol") == 0) method = 2;
    else if (strcmp(argv[1], "-sasa") == 0) method = 3;
    else if (strcmp(argv[1], "-wsasa") == 0) method = 4;
    else if (strcmp(argv[1], "-vsasa") == 0) method = 5;
    else {
        Vnm_print(1, "\n*** Syntax error: invalid argument '%s'.\n", argv[1]);
        Vnm_print(1, "%s", usage);
        exit(666);
    }
    path = argv[2];

    Vnm_print(1, "Setting up atom list from %s\n", path);
    alist = Valist_ctor();
    Valist_readPQR(alist,"FILE","ASC",VNULL,path);

    Vnm_print(1, "Setting up VACC object\n");
    Vnm_print(1, "\tProbe radius = %4.3f\n", probe_radius);
    Vnm_print(1, "\tQuad sphere has %d pts\n", nsphere);
    acc = Vacc_ctor(alist,probe_radius,nx,ny,nz,nsphere);
    VASSERT(acc != VNULL);
    dx = acc->hx; dy = acc->hy; dz = acc->hzed;
    x_len = (acc->nx * acc->hx);
    y_len = (acc->ny * acc->hy);
    z_len = (acc->nz * acc->hzed);
    x_cen = (acc->grid_lower_corner[0] + x_len/2);
    y_cen = (acc->grid_lower_corner[1] + y_len/2);
    z_cen = (acc->grid_lower_corner[2] + z_len/2);

    switch (method) {
        case 0:
            Vnm_print(1, "\nVAN DER WAALS VOLUME QUADRATURE\n");
            /* Start the work steps */
            start_t = times(&tmst);
            for (x=x_cen-x_len; x<=x_cen+x_len; x=x+dx) {
                for (y=y_cen-y_len; y<=y_cen+y_len; y=y+dy) {
                    for (z=z_cen-z_len; z<=z_cen+z_len; z=z+dz) {
        
                        vec[0] = x; vec[1] = y; vec[2] = z;
                        w = 1.0;
                        if ((x==(x_cen-x_len)) || (x==(x_cen+x_len))) w=w*0.5;
                        if ((y==(y_cen-y_len)) || (y==(y_cen+y_len))) w=w*0.5;
                        if ((z==(z_cen-z_len)) || (z==(z_cen+z_len))) w=w*0.5;
                        if (Vacc_vdwAcc(acc,vec) == 0) vdwVol += dx*dy*dz*w;
                    }
                }
            }
            stop_t = times(&tmst);
            t = (double)(stop_t - start_t)/((double)clktck);
            Vnm_print(1, "\tApprox. volume = %4.3f A^3\n",vdwVol);
            Vnm_print(1, "\t%g sec for %d quad pts implies %g sec per pt.\n",
              t, VRINT(8*x_len*y_len*z_len/(dx*dy*dz)), 
              t/(8*x_len*y_len*z_len/(dx*dy*dz)));
            fflush(stdout);
            break;

        case 1:
            Vnm_print(1, "\nINFLATED VAN DER WAALS VOLUME QUADRATURE\n");
            /* Start the work steps */
            start_t = times(&tmst);
            for (x=x_cen-x_len; x<=x_cen+x_len; x=x+dx) {
                for (y=y_cen-y_len; y<=y_cen+y_len; y=y+dy) {
                    for (z=z_cen-z_len; z<=z_cen+z_len; z=z+dz) {
        
                        vec[0] = x; vec[1] = y; vec[2] = z;
                        w = 1.0;
                        if ((x==(x_cen-x_len)) || (x==(x_cen+x_len))) w=w*0.5;
                        if ((y==(y_cen-y_len)) || (y==(y_cen+y_len))) w=w*0.5;
                        if ((z==(z_cen-z_len)) || (z==(z_cen+z_len))) w=w*0.5;
                        if (Vacc_ivdwAcc(acc,vec,probe_radius) == 0) 
                          ivdwVol += dx*dy*dz*w;
                    }
                }
            }
            stop_t = times(&tmst);
            t = (double)(stop_t - start_t)/((double)clktck);
            Vnm_print(1, "\tApprox. volume = %4.3f A^3\n",ivdwVol);
            Vnm_print(1, "\t%g sec for %d quad pts implies %g sec per pt.\n",
              t, VRINT(8*x_len*y_len*z_len/(dx*dy*dz)), 
              t/(8*x_len*y_len*z_len/(dx*dy*dz)));
            fflush(stdout);
            break;

        case 2:
            Vnm_print(1, "\nMOLECULAR (CONNOLLY) VOLUME QUADRATURE\n");
            /* Start the work steps */
            start_t = times(&tmst);
            for (x=x_cen-x_len; x<=x_cen+x_len; x=x+dx) {
                for (y=y_cen-y_len; y<=y_cen+y_len; y=y+dy) {
                    for (z=z_cen-z_len; z<=z_cen+z_len; z=z+dz) {
        
                        vec[0] = x; vec[1] = y; vec[2] = z;
                        w = 1.0;
                        if ((x==(x_cen-x_len)) || (x==(x_cen+x_len))) w=w*0.5;
                        if ((y==(y_cen-y_len)) || (y==(y_cen+y_len))) w=w*0.5;
                        if ((z==(z_cen-z_len)) || (z==(z_cen+z_len))) w=w*0.5;
                        if (Vacc_molAcc(acc,vec,probe_radius) == 0) 
                          molVol += dx*dy*dz*w;
                    }
                }
            }
            stop_t = times(&tmst);
            t = (double)(stop_t - start_t)/((double)clktck);
            Vnm_print(1, "\tApprox. volume = %4.3f A^3\n",molVol);
            Vnm_print(1, "\t%g sec for %d quad pts implies %g sec per pt.\n",
              t, VRINT(8*x_len*y_len*z_len/(dx*dy*dz)), 
              t/(8*x_len*y_len*z_len/(dx*dy*dz)));
            fflush(stdout);
            break;

        case 3:
            Vnm_print(1, "\nUNWEIGHTED SOLVENT-ACCESSIBLE SURFACE AREA\n");
            sasa = Vacc_totalSASA(acc, probe_radius);
            Vnm_print(1, "\tApprox. SASA = %1.12e A^2\n", sasa);
        
            Valist_dtor(&alist);
            Vacc_dtor(&acc);
            break;

        case 4:
            Vnm_print(1, "\nWEIGHTED SOLVENT-ACCESSIBLE SURFACE AREA\n");
            Vnm_print(1, "\tNOTE: Using the \"charge\" column of the PQR file\
 as weights...\n");
            sasa = 0;
            for (i=0; i<Valist_getNumberAtoms(alist); i++) {
                atom = Valist_getAtom(alist, i);
                sasa += (Vatom_getCharge(atom)*Vacc_atomSASA(acc, 
                         probe_radius, i));
            }
            Vnm_print(1, "\tApprox. weighted SASA = %1.12e A^2\n", sasa);
            Valist_dtor(&alist);
            Vacc_dtor(&acc);
            break;

        case 5:
            Vnm_print(1, "\nVERBOSE UNWEIGHTED SOLVENT-ACCESSIBLE SURFACE\
 AREA\n");
            sasa = 0;
            for (i=0; i<Valist_getNumberAtoms(alist); i++) {
                atom = Valist_getAtom(alist, i);
                atom_sasa = Vacc_atomSASA(acc, probe_radius, i);
                Vnm_print(1, "\tAtom %d:  %1.12E A^2\n", i, atom_sasa);
                sasa += atom_sasa;
            }
            Vnm_print(1, "\t--------------------------------------\n", sasa);
            Vnm_print(1, "\tTOTAL:  %1.12e A^2\n", sasa);
            Valist_dtor(&alist);
            Vacc_dtor(&acc);
            break;


        default:
            break;
    }
        
    return 0;
}
