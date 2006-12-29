/**
 *  @file    smooth.c
 *  @author  Nathan Baker
 *  @brief   Convolve grid data with various filters
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
 * Copyright (c) 2002-2007.  Washington University in St. Louis.
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

#define IJK(i,j,k)  (((k)*(nx)*(ny))+((j)*(nx))+(i))
#define ERRRC 2

typedef enum Smooth_Filter {
    SM_GAUSSIAN  /**< Gaussian filter */
} Smooth_Filter;


VEMBED(rcsid="$Id$")

int gaussian(Vgrid *grid, double stddev, double bandwidth);

int usage(int rc) {

    char *usage = "\n\n\
    ----------------------------------------------------------------------\n\
    This driver program convolves grid data with various filters.  It is\n\
    invoked as:\n\
      smooth <args>\n\n\
    where <args> is a list of the following arguments:\n\
      REQUIRED GENERAL ARGUMENTS:\n\
      --format=<format>  where <format> specifies the data format and is one\n\
                          of the following: dx (OpenDX)\n\
      --input=<file>     where <file> is the input mesh data in the\n\
                         specified format\n\
      --output=<file>    where <file> is the output mesh data in the\n\
                         specified format\n\
      --filter=<filter>  where <filter> is the filter with which the data\n\
                         will be convolved and is one of the following:\n\
                         gaussian (Gaussian filter)\n\
      REQUIRED FILTER-SPECIFIC ARGUMENTS:\n\
        Gaussian filter:\n\
        --stddev=<n>     the standard deviation of the filter (in A)\n\
        --bandwidth=<n>  the bandwith of the filter (in units of stddev)\n\
    ----------------------------------------------------------------------\n\n";

    Vnm_print(2, usage);

    exit(rc);
 
    return 0;
}

int main(int argc, char **argv) {

    /* *************** VARIABLES ******************* */
    int i;
    Vgrid *grid = VNULL;
    /* Input parameters */
    Vdata_Format format; int gotFormat = 0;
    char inPath[VMAX_BUFSIZE]; int gotInPath = 0;
    char outPath[VMAX_BUFSIZE]; int gotOutPath = 0;
    Smooth_Filter filter; int gotFilter = 0;
    double stddev; int gotStddev = 0;
    double bandwidth; int gotBandwidth = 0;

    char *header = "\n\n\
    ----------------------------------------------------------------------\n\
    ----------------------------------------------------------------------\n\
    \n\n";

    /* *************** CHECK INVOCATION ******************* */
    Vio_start();
    Vnm_redirect(1);
    Vnm_print(1, "%s", header);
    for (i=1; i<argc; i++) {
        Vnm_print(1, "Parsing: %s...\n", argv[i]);
        if (strstr(argv[i], "--format") != NULL) {
            if (strstr(argv[i], "dx") != NULL) {
                gotFormat = 1;
                format = VDF_DX;
            } else {
                Vnm_print(2, "Error:  %s\n", argv[i]);
                usage(2);
            }
        } else if (strstr(argv[i], "--input") != NULL) {
            if (sscanf(argv[i], "--input=%s", inPath) == 1) gotInPath = 1;
            else {
                Vnm_print(2, "Error:  %s\n", argv[i]);
                usage(2);
            }
        } else if (strstr(argv[i], "--output") != NULL) {
            if (sscanf(argv[i], "--output=%s", outPath) == 1) gotOutPath = 1;
            else {
                Vnm_print(2, "Error:  %s\n", argv[i]);
                usage(2);
            }
        } else if (strstr(argv[i], "--filter") != NULL) {
            if (strstr(argv[i], "gaussian") != NULL) {
                gotFilter = 1;
                filter = SM_GAUSSIAN;
            } else {
                Vnm_print(2, "Error:  %s\n", argv[i]);
                usage(2);
            }
        } else if (strstr(argv[i], "--stddev") != NULL) {
            if (sscanf(argv[i], "--stddev=%lf", &stddev) == 1) gotStddev = 1;
            else {
                Vnm_print(2, "Error:  %s\n", argv[i]);
                usage(2);
            }
        } else if (strstr(argv[i], "--bandwidth") != NULL) {
            if (sscanf(argv[i], "--bandwidth=%lf", &bandwidth) == 1) 
                gotBandwidth = 1;
            else {
                Vnm_print(2, "Error:  %s\n", argv[i]);
                usage(2);
            }
        } else {
            Vnm_print(2, "Error:  %s\n", argv[i]);
            usage(2);
        }
    }
    if (!gotFormat) {
        Vnm_print(2, "Error:  --format not specified!\n");
        usage(2);
    } 
    if (!gotInPath) {
        Vnm_print(2, "Error:  --input not specified!\n");
        usage(2);
    }
    if (!gotOutPath) {
        Vnm_print(2, "Error:  --output not specified!\n");
        usage(2);
    }
    if (!gotFilter) {
        Vnm_print(2, "Error:  --filter not specified!\n");
        usage(2);
    }
    if (filter == SM_GAUSSIAN) {
        if (!gotStddev) {
            Vnm_print(2, "Error:  --stddev not specified!\n");
            usage(2);
        }
        if (!gotBandwidth) {
            Vnm_print(2, "Error:  --bandwidth not specified!\n");
            usage(2);
        }
    }

    /* *************** READ DATA ******************* */
    Vnm_print(1, "main:  Reading data from %s...\n", inPath);
    grid = Vgrid_ctor(0, 0, 0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, VNULL);
    if (format == VDF_DX) {
        if (!Vgrid_readDX(grid, "FILE", "ASC", VNULL, inPath)) {
            Vnm_print(2, "main:  Problem reading OpenDX-format grid from %s\n",
              inPath);
            return ERRRC;
        }
    }

    /* *************** SMOOTH ******************* */
    switch(filter) {
        case SM_GAUSSIAN:
           Vnm_print(1, "Smoothing data with Gaussian filter...\n");
           gaussian(grid, stddev, bandwidth);
           break;
        default:
           Vnm_print(2, "Invalid format (%d)!\n", format);
           usage(2);
    }

    /* *************** READ DATA ******************* */
    Vnm_print(1, "main:  Writing data to %s...\n", outPath);
    if (format == VDF_DX) 
      Vgrid_writeDX(grid, "FILE", "ASC", VNULL, outPath, "Smoothed data",
        VNULL);

    return 0;
}

int gaussian(Vgrid *grid, double stddev, double bandwidth) {

    int nx, ny, nz, iband, jband, kband, i, j, k, ii, jj, kk;
    int kkmin, jjmin, iimin, kkmax, jjmax, iimax;
    double hx, hy, hzed, xmin, ymin, zmin, *newData, *oldData, scal, norm;
    double u, ga, dist2;

    Vnm_print(1, "Gaussian filter:  std. dev. = %g A, bandwidth = %g A.\n",
      stddev, bandwidth*stddev);

    nx = grid->nx; ny = grid->ny; nz = grid->nz;
    hx = grid->hx; hy = grid->hy; hzed = grid->hzed;
    xmin = grid->xmin; ymin = grid->ymin; zmin = grid->zmin;
    Vnm_print(1, "Grid:  %d x %d x %d points\n", nx, ny, nz);
    Vnm_print(1, "Grid:  %g, %g, %g A spacing\n", hx, hy, hzed);
    Vnm_print(1, "Grid:  (%g, %g, %g) A origin\n", xmin, ymin, zmin);

    /* Convert HALF bandwidth to grid units */
    iband = (int)(stddev*bandwidth/hx);
    jband = (int)(stddev*bandwidth/hy);
    kband = (int)(stddev*bandwidth/hzed);
    Vnm_print(1, "Bandwidth converted to %d x %d x %d grid units.\n",
      iband, jband, kband);
    Vnm_print(1, "This means any non-zero data within (%g, %g, %g) of the\n", 
      (iband+1)*hx, (jband+1)*hy, (kband+1)*hzed);
    Vnm_print(1, "domain boundary will be convolved differently.\n");

    /* Get exponent scaling factor */
    scal = 2.0 * stddev * stddev;
    VASSERT(scal > 0);
    scal = 1.0/scal;

    /* Get data */
    oldData = grid->data;
    newData = Vmem_malloc(VNULL, (nx*ny*nz), sizeof(double));

    /* Apply filter */
    for (k=1; k<(nz-1); k++) {
        kkmin = VMAX2(1, (k - kband));
        kkmax = VMIN2((nz-1), (k + kband));
        for (j=1; j<(ny-1); j++) {
            jjmin = VMAX2(1, (j - jband));
            jjmax = VMIN2((ny-1), (j + jband));
            for (i=1; i<(nx-1); i++) {
                iimin = VMAX2(1, (i - iband));
                iimax = VMIN2((nx-1), (i + iband));
                u = 0;
                /* We normalize within the loop to conserve densities */
                norm = 0;
                for (kk=kkmin; kk<kkmax; kk++) {
                    for (jj=jjmin; jj<jjmax; jj++) {
                        for (ii=iimin; ii<iimax; ii++) {
                            dist2 = VSQR(hx*(i-ii)) + VSQR(hy*(j-jj)) +
                              VSQR(hzed*(k-kk));
                            ga = VEXP(-dist2*scal);
                            u += (ga*oldData[IJK(ii,jj,kk)]);
                            norm += ga;
                        } /* ii loop */
                    } /* jj loop */
                } /* kk loop */
                if (norm > 0) newData[IJK(i,j,k)] = u/norm;
                else newData[IJK(i,j,k)] = 0;
            } /* i loop */
        } /* j loop */
    } /* k loop */

    /* Replace data */
    for (i=0; i<(nx*ny*nz); i++) grid->data[i] = newData[i];
    Vmem_free(VNULL, (nx*ny*nz), sizeof(double), (void **)&newData);

    return 0;
}
