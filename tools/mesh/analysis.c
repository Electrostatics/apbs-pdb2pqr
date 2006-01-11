/**
 *  @file    analysis.c
 *  @author  Nathan Baker
 *  @brief   Program that computes various analyses for a scalar field
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
 * Copyright (c) 2002-2005.  Washington University in St. Louis.
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
 * programs or libraries that are released under the GNU LGPL and with
 * code included in releases of ISIM, PMV, PyMOL, SMOL, VMD.  This
 * special exception permission is also extended to any software listed
 * in the SPECIAL GPL EXCEPTION clauses by the PMG, FEtk, MC, or MALOC
 * libraries.
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

/**
 * @brief  Prints usage information and exits
 * @author  Nathan Baker
 * @param  rc  Exit code */
void usage(int rc) {
    char *usage = "\n\
A program to calculate various metrics for a scalar field\n\
  Usage:  analysis <req args> [opts]\n\
where <req args> are the required arguments:\n\
  --format=<format>  The input file format.  Acceptable values include\n\
       dx  OpenDX format\n\
  --scalar=<path>  The path to the scalar data file\n\
and where [opts] are the options:\n\
  --help  Print this message\n\
  --mask=<path>  A file with scalar values specifying a \"mask\" or\n\
    characteristic function for the similarity calculation.  This file\n\
    contains values between 1 and 0 which are multiplied against the\n\
    scalar data set analysis calculations are performed.\n\
    \n";

    Vnm_print(2, usage);
    exit(rc);
}

/**
 * @brief  Read a grid
 * @param  grid  Pointer to Vgrid file to be read
 * @param  path  Path to read from 
 * @param  format  Format to read
 * @return 1 if successful, 0 otherwise */
int readGrid(Vgrid **grid, char *path, Vdata_Format format) {

    *grid = Vgrid_ctor(0, 0, 0, 
            0.0, 0.0, 0.0,
            0.0, 0.0, 0.0,
            VNULL);

    switch (format) {
        case VDF_DX:
            return Vgrid_readDX(*grid, "FILE", "ASC", VNULL, path);
            break;
        case VDF_UHBD:
            Vnm_print(2, "Sorry, UHBD input not supported yet!\n");
            return 0;
            break;
        case VDF_AVS:
            Vnm_print(2, "Sorry, AVS input not supported yet!\n");
            return 0;
            break;
        default:
            Vnm_print(2, "Unknown data format (%d)!\n", format);
            return 0;
    }

    return 1;
}

int main(int argc, char **argv) {

    /* *************** VARIABLES ******************* */
    int i, j, k, onGridS, onGridV, nx, ny, nz;
    double hx, hy, hzed, xmin, ymin, zmin, dvol, svol, gvol;
    double norm_L1, norm_L2, snorm_H1, norm_H1, norm_Linf;
    double maxS, maxSpt[3], maxG2, maxG2pt[3];
    double minS, minSpt[3], minG2, minG2pt[3];
    int haveMaxS, haveMaxG2, haveMinS, haveMinG2;
    double val, sval, mval, pt[3];
    double gval2, gval[3];
    Vgrid *scalar, *mask;
    char scalarPath[VMAX_ARGLEN];
    int gotScalar = 0;
    char maskPath[VMAX_ARGLEN];
    int gotMask = 0;
    Vdata_Format format;
    int gotFormat = 0;
    char *tstr, *targ;
 
    /* *************** CHECK INVOCATION ******************* */
    Vio_start();
    /* Parse args */
    for (i=1; i<argc; i++) {
        targ = argv[i];
        /* FORMAT */
        tstr = strstr(targ, "format");
        if (tstr != NULL) {
            tstr = tstr + 7;
            if (strcmp(tstr, "dx") == 0) {
                format = VDF_DX;
                gotFormat = 1;
            } else {
                Vnm_print(2, "Error!  Unknown format (%s)!\n", tstr);
                usage(2);
            }
        }
        /* SCALAR */
        tstr = strstr(targ, "scalar");
        if (tstr != NULL) {
            tstr = tstr + 7;
            strncpy(scalarPath, tstr, VMAX_ARGLEN);
            gotScalar = 1;
        }
        /* HELP */
        tstr = strstr(targ, "help");
        if (tstr != NULL) usage(0);
        /* MASK */
        tstr = strstr(targ, "mask");
        if (tstr != NULL) {
            tstr = tstr + 5;
            strncpy(maskPath, tstr, VMAX_ARGLEN);
            gotMask = 1;
        }
    }
    /* Check and print args */
    if (!gotFormat) {
        Vnm_print(2, "Error!  --format not specified!\n");
        usage(2);
    } else {
        switch (format) {
            case VDF_DX:
                Vnm_print(1, "format:  OpenDX\n");
                break;
            case VDF_UHBD:
                Vnm_print(1, "format:  UHBD\n");
                break;
            case VDF_AVS:
                Vnm_print(1, "format:  AVS\n");
                break;
            default:
                Vnm_print(2, "Error!  Unknown format (%d)!\n", format);
                usage(2);
        }
    }
    if (!gotScalar) {
        Vnm_print(2, "Error!  --scalar not specified!\n");
        usage(2);
    } else {
        Vnm_print(1, "Data set:  %s\n", scalarPath);
    }
    if (gotMask) {
        Vnm_print(1, "Mask:  %s\n", maskPath);
    }

    /* Read scalar set */
    Vnm_print(1, "Reading scalar data set from %s...\n", scalarPath);
    if (!readGrid(&scalar, scalarPath, format)) {
        Vnm_print(2, "Error reading scalar data set!\n");
        return 2;
    }
    Vnm_print(1, "Read %d x %d x %d grid.\n", 
            scalar->nx, scalar->ny, scalar->nz);

    /* Read mask */
    if (gotMask) {
        Vnm_print(1, "Reading mask data set from %s...\n", maskPath);
        if (!readGrid(&mask, maskPath, format)) {
            Vnm_print(2, "Error reading mask data set!\n");
            return 2;
        }
        Vnm_print(1, "Read %d x %d x %d grid.\n", 
                mask->nx, mask->ny, mask->nz);
    }

    /* Calculate relative L2 norm of difference */
    Vnm_print(1, "Calculating metrics...\n");
    nx = scalar->nx; ny = scalar->ny; nz = scalar->nz; 
    hx = scalar->hx; hy = scalar->hy; hzed = scalar->hzed; 
    dvol = (hx*hy*hzed);
    xmin = scalar->xmin; ymin = scalar->ymin; zmin = scalar->zmin; 
    norm_L1 = 0; norm_L2 = 0; snorm_H1 = 0; norm_H1 = 0;
    haveMaxS = 0; haveMinS = 0; haveMaxG2 = 0; haveMinG2 = 0;
    svol = 0; gvol = 0;
    for (i=0; i<nx; i++) {
        pt[0] = i*hx + xmin;
        for (j=0; j<ny; j++) {
            pt[1] = j*hy + ymin;
            for (k=0; k<nz; k++) {

                /* Grid value */
                pt[2] = k*hzed + zmin;
                onGridS = Vgrid_value(scalar, pt, &sval);
                onGridV = Vgrid_gradient(scalar, pt, gval);
                if (onGridV) {
                    gval2 = 0.0;
                    gval2 = VSQR(gval[0]) + VSQR(gval[1]) + VSQR(gval[2]);
                } else gval2 = 0.0;
                if (gotMask) onGridS = Vgrid_value(mask, pt, &mval);
                else mval = 1.0; 

                /* Max/min */
                if (mval > 0) {
                    if ((!haveMaxS) || (sval > maxS) ) {
                        haveMaxS = 1;
                        maxS = sval;
                        maxSpt[0] = pt[0];
                        maxSpt[1] = pt[1];
                        maxSpt[2] = pt[2];
                    }
                    if ((!haveMinS) || (sval < minS) ) {
                        haveMinS = 1;
                        minS = sval;
                        minSpt[0] = pt[0];
                        minSpt[1] = pt[1];
                        minSpt[2] = pt[2];
                    }
                    if ((!haveMaxG2) || (gval2 > maxG2) ) {
                        haveMaxG2 = 1;
                        maxG2 = gval2;
                        maxG2pt[0] = pt[0];
                        maxG2pt[1] = pt[1];
                        maxG2pt[2] = pt[2];
                    }
                    if ((!haveMinG2) || (gval2 < minG2) ) {
                        haveMinG2 = 1;
                        minG2 = gval2;
                        minG2pt[0] = pt[0];
                        minG2pt[1] = pt[1];
                        minG2pt[2] = pt[2];
                    }
                }

                if (onGridS) {
                    val = sval*mval;

                    /* L2 */
                    norm_L2 += VSQR(val);
                    /* L1 */
                    norm_L1 += VABS(val);
                    /* Volume */
                    svol += dvol;

                }

                if (onGridV && onGridS) {
                    val = mval*(VSQR(gval[0]) + VSQR(gval[1]) + VSQR(gval[2]));
                    snorm_H1 += VSQR(val);
                    gvol += dvol;
                }
            }
        }
    }
    norm_Linf = VMAX2(VABS(maxS), VABS(minS));
    norm_L2 = VSQRT(norm_L2*dvol);
    norm_L1 = (norm_L1*dvol);
    snorm_H1 = VSQRT(snorm_H1*dvol);
    norm_H1 = VSQRT(VSQR(snorm_H1)+VSQR(norm_L2));

    Vnm_print(1, "Volume used to calculate L2 and L1 norms = %1.12E\n", 
            svol);
    Vnm_print(1, "Volume used to calculate H1 norms        = %1.12E\n", 
            gvol);
    Vnm_print(1, "Max scalar value                         = %1.12E\n", 
            maxS);
    Vnm_print(1, "Max scalar value location                = (%4.3f, %4.3f, %4.3f)\n", 
            maxSpt[0], maxSpt[1], maxSpt[2]);
    Vnm_print(1, "Min scalar value                         = %1.12E\n", 
            minS);
    Vnm_print(1, "Min scalar value location                = (%4.3f, %4.3f, %4.3f)\n", 
            minSpt[0], minSpt[1], minSpt[2]);
    Vnm_print(1, "Max gradient-squared value               = %1.12E\n", 
            maxG2);
    Vnm_print(1, "Max gradient-squared value location      = (%4.3f, %4.3f, %4.3f)\n", 
            maxG2pt[0], maxG2pt[1], maxG2pt[2]);
    Vnm_print(1, "Min gradient-squared value               = %1.12E\n", 
            minG2);
    Vnm_print(1, "Min gradient-squared value location      = (%4.3f, %4.3f, %4.3f)\n", 
            minG2pt[0], minG2pt[1], minG2pt[2]);
    Vnm_print(1, "L2 norm                                  = %1.12E\n", 
            norm_L2);
    Vnm_print(1, "L1 norm                                  = %1.12E\n", 
            norm_L1);
    Vnm_print(1, "Linf norm                                = %1.12E\n", 
            norm_Linf);
    Vnm_print(1, "H1 semi-norm                             = %1.12E\n", 
            snorm_H1);
    Vnm_print(1, "H1 norm                                  = %1.12E\n", 
            norm_H1);

    return 0;
}
