/* ///////////////////////////////////////////////////////////////////////////
/// APBS -- Adaptive Poisson-Boltzmann Solver
///
///  Nathan A. Baker (nbaker@wasabi.ucsd.edu)
///  Dept. of Chemistry and Biochemistry
///  Dept. of Mathematics, Scientific Computing Group
///  University of California, San Diego 
///
///  Additional contributing authors listed in the code documentation.
///
/// Copyright © 1999. The Regents of the University of California (Regents).
/// All Rights Reserved. 
/// 
/// Permission to use, copy, modify, and distribute this software and its
/// documentation for educational, research, and not-for-profit purposes,
/// without fee and without a signed licensing agreement, is hereby granted,
/// provided that the above copyright notice, this paragraph and the
/// following two paragraphs appear in all copies, modifications, and
/// distributions.
/// 
/// IN NO EVENT SHALL REGENTS BE LIABLE TO ANY PARTY FOR DIRECT, INDIRECT,
/// SPECIAL, INCIDENTAL, OR CONSEQUENTIAL DAMAGES, INCLUDING LOST PROFITS,
/// ARISING OUT OF THE USE OF THIS SOFTWARE AND ITS DOCUMENTATION, EVEN IF
/// REGENTS HAS BEEN ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.  
/// 
/// REGENTS SPECIFICALLY DISCLAIMS ANY WARRANTIES, INCLUDING, BUT NOT
/// LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A
/// PARTICULAR PURPOSE.  THE SOFTWARE AND ACCOMPANYING DOCUMENTATION, IF
/// ANY, PROVIDED HEREUNDER IS PROVIDED "AS IS".  REGENTS HAS NO OBLIGATION
/// TO PROVIDE MAINTENANCE, SUPPORT, UPDATES, ENHANCEMENTS, OR
/// MODIFICATIONS. 
////////////////////////////////////////////////////////////////////////////
/// rcsid="$Id$"
//////////////////////////////////////////////////////////////////////////// */

/* ///////////////////////////////////////////////////////////////////////////
// File:     acc.c
//
// Purpose:  Calculate solvent/ion accessible volumes, areas, etc.
//
// Author:   Nathan Baker
/////////////////////////////////////////////////////////////////////////// */

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
    /* Quadrature steps */
    double dx, dy, dz;
    double x_len, y_len, z_len;
    double x_cen, y_cen, z_cen;
    double x, y, z;
    double w = 1.0;
    double vec[3];

    /* TIMING VARIABLES */
    double t;
    clock_t start_t;
    clock_t stop_t;
    struct tms tmst;
    long clktck = sysconf(_SC_CLK_TCK);
 
    char *usage = "\n\nUsage: acc <-vdw|-ivdw|-mol|-sasa> <molecule.pqr>"
                  "\tWhere -vdw   => volume inside VdW surface\n"
                  "\t\t-ivdw  => volume inside VdW surface inflated by probe radius\n"
                  "\t\t-mol   => volume inside molecular surface\n"
                  "\t\t-sasa  => solvent-accessible surface area\n"
                  "\tand <molecule.pqr> is the path to the molecule structure in PQR format\n\n";

    if (argc != 3) {
        printf("\n*** Syntax error: got %d arguments, expected 2.\n",argc);
        printf("%s", usage);
        exit(666);
    };
    if (strcmp(argv[1], "-vdw") == 0) method = 0;
    else if (strcmp(argv[1], "-ivdw") == 0) method = 1;
    else if (strcmp(argv[1], "-mol") == 0) method = 2;
    else if (strcmp(argv[1], "-sasa") == 0) method = 3;
    else {
        printf("\n*** Syntax error: invalid argument '%s'.\n", argv[1]);
        printf("%s", usage);
        exit(666);
    }
    path = argv[2];

    printf("Setting up atom list from %s\n", path);
    alist = Valist_ctor();
    Valist_readPQR(alist,"FILE","ASC",VNULL,path);

    printf("Setting up VACC object\n");
    printf("\tProbe radius = %4.3f\n", probe_radius);
    printf("\tQuad sphere has %d pts\n", nsphere);
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
            printf("\nVAN DER WAALS VOLUME QUADRATURE\n");
            /* Start the work steps */
            start_t = times(&tmst);
            for (x=x_cen-x_len; x<=x_cen+x_len; x=x+dx) {
                for (y=y_cen-y_len; y<=y_cen+y_len; y=y+dy) {
                    for (z=z_cen-z_len; z<=z_cen+z_len; z=z+dz) {
        
                        vec[0] = x; vec[1] = y; vec[2] = z;
                        w = 1.0;
                        if ((x==(x_cen - x_len)) || (x==(x_cen + x_len))) w=w*0.5;
                        if ((y==(y_cen - y_len)) || (y==(y_cen + y_len))) w=w*0.5;
                        if ((z==(z_cen - z_len)) || (z==(z_cen + z_len))) w=w*0.5;
                        if (Vacc_vdwAcc(acc,vec) == 0) vdwVol += dx*dy*dz*w;
                    }
                }
            }
            stop_t = times(&tmst);
            t = (double)(stop_t - start_t)/((double)clktck);
            printf("\tApprox. volume = %4.3f A^3\n",vdwVol);
            printf("\t%g sec for %d quad pts implies %g sec per pt.\n",
              t, VRINT(8*x_len*y_len*z_len/(dx*dy*dz)), 
              t/(8*x_len*y_len*z_len/(dx*dy*dz)));
            fflush(stdout);
            break;

        case 1:
            printf("\nINFLATED VAN DER WAALS VOLUME QUADRATURE\n");
            /* Start the work steps */
            start_t = times(&tmst);
            for (x=x_cen-x_len; x<=x_cen+x_len; x=x+dx) {
                for (y=y_cen-y_len; y<=y_cen+y_len; y=y+dy) {
                    for (z=z_cen-z_len; z<=z_cen+z_len; z=z+dz) {
        
                        vec[0] = x; vec[1] = y; vec[2] = z;
                        w = 1.0;
                        if ((x==(x_cen - x_len)) || (x==(x_cen + x_len))) w=w*0.5;
                        if ((y==(y_cen - y_len)) || (y==(y_cen + y_len))) w=w*0.5;
                        if ((z==(z_cen - z_len)) || (z==(z_cen + z_len))) w=w*0.5;
                        if (Vacc_ivdwAcc(acc,vec,probe_radius) == 0) ivdwVol += dx*dy*dz*w;
                    }
                }
            }
            stop_t = times(&tmst);
            t = (double)(stop_t - start_t)/((double)clktck);
            printf("\tApprox. volume = %4.3f A^3\n",ivdwVol);
            printf("\t%g sec for %d quad pts implies %g sec per pt.\n",
              t, VRINT(8*x_len*y_len*z_len/(dx*dy*dz)), 
              t/(8*x_len*y_len*z_len/(dx*dy*dz)));
            fflush(stdout);
            break;

        case 2:
            printf("\nMOLECULAR (CONNOLLY) VOLUME QUADRATURE\n");
            /* Start the work steps */
            start_t = times(&tmst);
            for (x=x_cen-x_len; x<=x_cen+x_len; x=x+dx) {
                for (y=y_cen-y_len; y<=y_cen+y_len; y=y+dy) {
                    for (z=z_cen-z_len; z<=z_cen+z_len; z=z+dz) {
        
                        vec[0] = x; vec[1] = y; vec[2] = z;
                        w = 1.0;
                        if ((x==(x_cen - x_len)) || (x==(x_cen + x_len))) w=w*0.5;
                        if ((y==(y_cen - y_len)) || (y==(y_cen + y_len))) w=w*0.5;
                        if ((z==(z_cen - z_len)) || (z==(z_cen + z_len))) w=w*0.5;
                        if (Vacc_molAcc(acc,vec,probe_radius) == 0) molVol += dx*dy*dz*w;
                    }
                }
            }
            stop_t = times(&tmst);
            t = (double)(stop_t - start_t)/((double)clktck);
            printf("\tApprox. volume = %4.3f A^3\n",molVol);
            printf("\t%g sec for %d quad pts implies %g sec per pt.\n",
              t, VRINT(8*x_len*y_len*z_len/(dx*dy*dz)), 
              t/(8*x_len*y_len*z_len/(dx*dy*dz)));
            fflush(stdout);
            break;

        case 3:
            printf("\nSOLVENT-ACCESSIBLE SURFACE AREA\n");
            sasa = Vacc_totalSASA(acc, probe_radius);
            printf("\tApprox. SASA = %4.3f A^2\n", sasa);
        
            Valist_dtor(&alist);
            Vacc_dtor(&acc);
            break;

        default:
            break;
    }
        
    return 0;
}
