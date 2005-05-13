/**
 *  @file    vpmg.c
 *  @author  Nathan Baker
 *  @brief   Class Vpmg methods
 *  @ingroup Vpmg
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
 * Copyright (c) 2002-2004.  Washington University in St. Louis.
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
 * @endverbatim
 */

#include "apbscfg.h"
#include "vpmg-private.h"
#include "apbs/vpmg.h"
#include "apbs/vhal.h"

VEMBED(rcsid="$Id$")

#if !defined(VINLINE_VPMG)

VPUBLIC unsigned long int Vpmg_memChk(Vpmg *thee) {
    if (thee == VNULL) return 0;
    return Vmem_bytes(thee->vmem);
}

#endif /* if !defined(VINLINE_VPMG) */


VPUBLIC void Vpmg_printColComp(Vpmg *thee, char path[72], char title[72], 
  char mxtype[3], int flag) {

    int nn, nxm2, nym2, nzm2, ncol, nrow, nonz; 
    double *nzval;
    int *colptr, *rowind;

    /* Calculate the total number of unknowns */
    nxm2 = thee->pmgp->nx - 2;
    nym2 = thee->pmgp->ny - 2;
    nzm2 = thee->pmgp->nz - 2;
    nn = nxm2*nym2*nzm2;
    ncol = nn;
    nrow = nn;

    /* Calculate the number of non-zero matrix entries:
     *    nn       nonzeros on diagonal
     *    nn-1     nonzeros on first off-diagonal
     *    nn-nx    nonzeros on second off-diagonal
     *    nn-nx*ny nonzeros on third off-diagonal
     *
     *    7*nn-2*nx*ny-2*nx-2 TOTAL non-zeros
     */
    nonz = 7*nn - 2*nxm2*nym2 - 2*nxm2 - 2;
    nzval  = Vmem_malloc(thee->vmem, nonz, sizeof(double));
    rowind = Vmem_malloc(thee->vmem, nonz, sizeof(int));
    colptr = Vmem_malloc(thee->vmem, (ncol+1), sizeof(int));

#ifndef VAPBSQUIET
    Vnm_print(1, "Vpmg_printColComp:  Allocated space for %d nonzeros\n",
      nonz);
#endif

    F77BCOLCOMP(thee->iparm, thee->rparm, thee->iwork, thee->rwork,
      nzval, rowind, colptr, &flag);

#if 0
    for (i=0; i<nn; i++) {
        Vnm_print(1, "nnz(%d) = %g\n", i, nzval[i]);
    }
#endif

    /* I do not understand why I need to pass nzval in this way, but it
     * works... */
    F77PCOLCOMP(&nrow, &ncol, &nonz, &(nzval[0]), rowind, colptr, path, title, 
      mxtype);

    Vmem_free(thee->vmem, (ncol+1), sizeof(int), (void **)&colptr);
    Vmem_free(thee->vmem, nonz, sizeof(int), (void **)&rowind);
    Vmem_free(thee->vmem, nonz, sizeof(double), (void **)&nzval);

}

VPUBLIC Vpmg* Vpmg_ctor(Vpmgp *pmgp, Vpbe *pbe, int focusFlag, 
        Vpmg *pmgOLD, MGparm *mgparm, PBEparm_calcEnergy energyFlag) {

    Vpmg *thee = VNULL;

    thee = Vmem_malloc(VNULL, 1, sizeof(Vpmg) );
    VASSERT(thee != VNULL);
    VASSERT( Vpmg_ctor2(thee, pmgp, pbe, focusFlag, pmgOLD, mgparm, 
                energyFlag) );

    return thee;
}


VPUBLIC void Vpmg_solve(Vpmg *thee) {

    int i;
    double zkappa2;

    thee->filled = 0;

    /* This is a really disgusting hack, however, it preserves the relative
	 * clarity of the code elsewhere.  There are two paths the code can follow
	 * based on whether we're solving the NPBE or LPBE.  
	 * - For the NPBE, we need to keep track of (possibly asymmetric)
	 *   contributions from different ion species.  In this case, ccf should
	 *   only contain a characteristic function which describes ion
	 *   accessibility.  The appropriate scaling coefficients are included in
	 *   mypde.f (as set by F77MYPDEFINIT).
     * - For the LPBE, we don't need to keep track of individual ion
	 *   contributions.  In this case, ccf should contain a characteristic
	 *   function _ALREADY SCALED BY THE APPROPRIATE COEFFICIENTS_.
     * In all cases, the fillco functions in vpmg-setup.c only fill ccf with
	 * the values of the characteristic function.  This would be fine if all
	 * paths of execution used the functions in mypde.f.  Unforunately, the
	 * functions in mypde.f are not called by PMG in the case of the LPBE.
	 * Rather than modifying PMG (we want to maintain as much compatibility as
	 * possible), we will scale ccf here and then immediately unscale it at the
	 * end of this function. */
    zkappa2 = Vpbe_getZkappa2(thee->pbe);
    if (zkappa2 > VPMGSMALL) {
        for (i=0; i<(thee->pmgp->nx*thee->pmgp->ny*thee->pmgp->nz); i++) {
            thee->ccf[i] = zkappa2*thee->ccf[i];
        }
    }

    switch(thee->pmgp->meth) {
        /* CGMG (linear) */
        case 0:
            F77CGMGDRIV(thee->iparm, thee->rparm, thee->iwork, thee->rwork,
              thee->u, thee->xf, thee->yf, thee->zf, thee->gxcf, thee->gycf,
              thee->gzcf, thee->a1cf, thee->a2cf, thee->a3cf, thee->ccf,
              thee->fcf, thee->tcf);
            thee->filled = 0;
            break;
        /* Newton (nonlinear) */
        case 1:
            F77NEWDRIV(thee->iparm, thee->rparm, thee->iwork, thee->rwork, 
              thee->u, thee->xf, thee->yf, thee->zf, thee->gxcf, thee->gycf,
              thee->gzcf, thee->a1cf, thee->a2cf, thee->a3cf, thee->ccf, 
              thee->fcf, thee->tcf);
            thee->filled = 0;
            break;
        /* MG (linear/nonlinear) */
        case 2:
	    F77MGDRIV(thee->iparm, thee->rparm, thee->iwork, thee->rwork,
	      thee->u, thee->xf, thee->yf, thee->zf, thee->gxcf, thee->gycf,
	      thee->gzcf, thee->a1cf, thee->a2cf, thee->a3cf, thee->ccf,
              thee->fcf, thee->tcf);
            thee->filled = 0;
            break;
        /* CGHS (linear/nonlinear) */
        case 3: 
	    F77NCGHSDRIV(thee->iparm, thee->rparm, thee->iwork, thee->rwork,
	      thee->u, thee->xf, thee->yf, thee->zf, thee->gxcf, thee->gycf,
	      thee->gzcf, thee->a1cf, thee->a2cf, thee->a3cf, thee->ccf,
              thee->fcf, thee->tcf);
            thee->filled = 0;
            break;
        /* SOR (linear/nonlinear) */
        case 4:
	    F77NSORDRIV(thee->iparm, thee->rparm, thee->iwork, thee->rwork,
	      thee->u, thee->xf, thee->yf, thee->zf, thee->gxcf, thee->gycf,
	      thee->gzcf, thee->a1cf, thee->a2cf, thee->a3cf, thee->ccf,
              thee->fcf, thee->tcf);
            thee->filled = 0;
            break;
        /* GSRB (linear/nonlinear) */
        case 5:
	    F77NGSRBDRIV(thee->iparm, thee->rparm, thee->iwork, thee->rwork,
	      thee->u, thee->xf, thee->yf, thee->zf, thee->gxcf, thee->gycf,
	      thee->gzcf, thee->a1cf, thee->a2cf, thee->a3cf, thee->ccf,
              thee->fcf, thee->tcf); 
            thee->filled = 0;
            break;
        /* WJAC (linear/nonlinear) */
        case 6:
	    F77NWJACDRIV(thee->iparm, thee->rparm, thee->iwork, thee->rwork,
	      thee->u, thee->xf, thee->yf, thee->zf, thee->gxcf, thee->gycf,
	      thee->gzcf, thee->a1cf, thee->a2cf, thee->a3cf, thee->ccf,
              thee->fcf, thee->tcf);
            thee->filled = 0;
            break;
        /* RICH (linear/nonlinear) */
        case 7:
	    F77NRICHDRIV(thee->iparm, thee->rparm, thee->iwork, thee->rwork,
	      thee->u, thee->xf, thee->yf, thee->zf, thee->gxcf, thee->gycf,
	      thee->gzcf, thee->a1cf, thee->a2cf, thee->a3cf, thee->ccf,
              thee->fcf, thee->tcf);
            thee->filled = 0;
            break;
        /* Error handling */
        default: 
            Vnm_print(2, "Vpgm_solve: invalid solver method key (%d)\n",
              thee->pmgp->key);
            break;
    }

    /* Un-scale ccf (see long comment above) */
    if (zkappa2 > VPMGSMALL) {
        for (i=0; i<(thee->pmgp->nx*thee->pmgp->ny*thee->pmgp->nz); i++) {
            thee->ccf[i] = thee->ccf[i]/zkappa2;
        }
    }

}

    
VPUBLIC void Vpmg_dtor(Vpmg **thee) {
    
    if ((*thee) != VNULL) {
        Vpmg_dtor2(*thee);
        Vmem_free(VNULL, 1, sizeof(Vpmg), (void **)thee);
        (*thee) = VNULL;
    }

}

VPUBLIC void Vpmg_dtor2(Vpmg *thee) { 

    /* Clear out the FORTRAN arrays */
    F77MYPDEFCLEAR();

    /* Clean up the storage */
    Vmem_free(thee->vmem, 100, sizeof(int), (void **)&(thee->iparm));
    Vmem_free(thee->vmem, 100, sizeof(double), (void **)&(thee->rparm));
    Vmem_free(thee->vmem, thee->pmgp->niwk, sizeof(int), 
      (void **)&(thee->iwork));
    Vmem_free(thee->vmem, thee->pmgp->nrwk, sizeof(double), 
      (void **)&(thee->rwork));
    Vmem_free(thee->vmem, thee->pmgp->narr, sizeof(double),
      (void **)&(thee->a1cf));
    Vmem_free(thee->vmem, thee->pmgp->narr, sizeof(double), 
      (void **)&(thee->a2cf));
    Vmem_free(thee->vmem, thee->pmgp->narr, sizeof(double),
      (void **)&(thee->a3cf));
    Vmem_free(thee->vmem, thee->pmgp->narr, sizeof(double),
      (void **)&(thee->ccf));
    Vmem_free(thee->vmem, thee->pmgp->narr, sizeof(double), 
      (void **)&(thee->fcf));
    Vmem_free(thee->vmem, thee->pmgp->narr, sizeof(double), 
      (void **)&(thee->tcf));
    Vmem_free(thee->vmem, thee->pmgp->narr, sizeof(double), 
      (void **)&(thee->u));
    Vmem_free(thee->vmem, 5*(thee->pmgp->nx), sizeof(double),
      (void **)&(thee->xf));
    Vmem_free(thee->vmem, 5*(thee->pmgp->ny), sizeof(double),
      (void **)&(thee->yf));
    Vmem_free(thee->vmem, 5*(thee->pmgp->nz), sizeof(double),
      (void **)&(thee->zf));
    Vmem_free(thee->vmem, 10*(thee->pmgp->ny)*(thee->pmgp->nz), sizeof(double),
      (void **)&(thee->gxcf));
    Vmem_free(thee->vmem, 10*(thee->pmgp->nx)*(thee->pmgp->nz), sizeof(double),
      (void **)&(thee->gycf));
    Vmem_free(thee->vmem, 10*(thee->pmgp->nx)*(thee->pmgp->ny), sizeof(double),
      (void **)&(thee->gzcf));
    Vmem_free(thee->vmem, (thee->pmgp->nx)*(thee->pmgp->ny)*(thee->pmgp->nz), 
      sizeof(double), (void **)&(thee->pvec));

    Vmem_dtor(&(thee->vmem));
}

VPUBLIC void Vpmg_setPart(Vpmg *thee, double lowerCorner[3], 
        double upperCorner[3], int bflags[6]) {

    Valist *alist;
    Vatom *atom;
    int i, j, k, nx, ny, nz;
    double xmin, ymin, zmin, x, y, z, hx, hy, hzed, xok, yok, zok;
    double x0,x1,y0,y1,z0,z1;

    nx = thee->pmgp->nx;
    ny = thee->pmgp->ny;
    nz = thee->pmgp->nz;
    hx = thee->pmgp->hx;
    hy = thee->pmgp->hy;
    hzed = thee->pmgp->hzed;
    xmin = thee->pmgp->xcent - 0.5*hx*(nx-1);
    ymin = thee->pmgp->ycent - 0.5*hy*(ny-1);
    zmin = thee->pmgp->zcent - 0.5*hzed*(nz-1);

    xok = 0;
    yok = 0;
    zok = 0;

    /* We need have called Vpmg_fillco first */

    alist = thee->pbe->alist;

    Vnm_print(0, "Vpmg_setPart:  lower corner = (%g, %g, %g)\n",
      lowerCorner[0], lowerCorner[1], lowerCorner[2]);
    Vnm_print(0, "Vpmg_setPart:  upper corner = (%g, %g, %g)\n",
      upperCorner[0], upperCorner[1], upperCorner[2]);
    Vnm_print(0, "Vpmg_setPart:  actual minimums = (%g, %g, %g)\n",
      xmin, ymin, zmin);
    Vnm_print(0, "Vpmg_setPart:  actual maximums = (%g, %g, %g)\n",
      xmin+hx*(nx-1), ymin+hy*(ny-1), zmin+hzed*(nz-1));
    Vnm_print(0, "Vpmg_setPart:  bflag[FRONT] = %d\n", 
      bflags[VAPBS_FRONT]);
    Vnm_print(0, "Vpmg_setPart:  bflag[BACK] = %d\n", 
      bflags[VAPBS_BACK]);
    Vnm_print(0, "Vpmg_setPart:  bflag[LEFT] = %d\n", 
      bflags[VAPBS_LEFT]);
    Vnm_print(0, "Vpmg_setPart:  bflag[RIGHT] = %d\n", 
      bflags[VAPBS_RIGHT]);
    Vnm_print(0, "Vpmg_setPart:  bflag[UP] = %d\n", 
      bflags[VAPBS_UP]);
    Vnm_print(0, "Vpmg_setPart:  bflag[DOWN] = %d\n", 
      bflags[VAPBS_DOWN]);

    /* Identify atoms as inside, outside, or on the border
       If on the border, use the bflags to determine if there
       is an adjacent processor - if so, this atom should be equally
       shared. */

    for (i=0; i<Valist_getNumberAtoms(alist); i++) {
        atom = Valist_getAtom(alist, i);
        if ((atom->position[0] < upperCorner[0]) &&
            (atom->position[0] > lowerCorner[0])) xok = 1;
        else {
            if ((VABS(atom->position[0] - lowerCorner[0]) < VPMGSMALL) && 
                (bflags[VAPBS_LEFT] == 0)) xok = 1;
            else if ((VABS(atom->position[0] - lowerCorner[0]) < VPMGSMALL) && 
                (bflags[VAPBS_LEFT] == 1)) xok = 0.5;
            else if ((VABS(atom->position[0] - upperCorner[0]) < VPMGSMALL) &&
                (bflags[VAPBS_RIGHT] == 0)) xok = 1;
            else if ((VABS(atom->position[0] - upperCorner[0]) < VPMGSMALL) &&
                (bflags[VAPBS_RIGHT] == 1)) xok = 0.5;
            else xok = 0;
        }
        if ((atom->position[1] < upperCorner[1]) &&
            (atom->position[1] > lowerCorner[1])) yok = 1;
        else {
            if ((VABS(atom->position[1] - lowerCorner[1]) < VPMGSMALL) && 
                (bflags[VAPBS_BACK] == 0)) yok = 1;
            else if ((VABS(atom->position[1] - lowerCorner[1]) < VPMGSMALL) && 
                (bflags[VAPBS_BACK] == 1)) yok = 0.5;
            else if ((VABS(atom->position[1] - upperCorner[1]) < VPMGSMALL) &&
                (bflags[VAPBS_FRONT] == 0)) yok = 1;
            else if ((VABS(atom->position[1] - upperCorner[1]) < VPMGSMALL) &&
                (bflags[VAPBS_FRONT] == 1)) yok = 0.5;
            else yok = 0;
        }
        if ((atom->position[2] < upperCorner[2]) &&
            (atom->position[2] > lowerCorner[2])) zok = 1;
        else {
            if ((VABS(atom->position[2] - lowerCorner[2]) < VPMGSMALL) && 
                (bflags[VAPBS_DOWN] == 0)) zok = 1;
            else if ((VABS(atom->position[2] - lowerCorner[2]) < VPMGSMALL) && 
                (bflags[VAPBS_DOWN] == 1)) zok = 0.5;
            else if ((VABS(atom->position[2] - upperCorner[2]) < VPMGSMALL) &&
                (bflags[VAPBS_UP] == 0)) zok = 1;
            else if ((VABS(atom->position[2] - upperCorner[2]) < VPMGSMALL) &&
                (bflags[VAPBS_UP] == 1)) zok = 0.5;
            else zok = 0;
        }

        atom->partID = xok*yok*zok;     
    }

    /* Load up pvec -
       For all points within h{axis}/2 of a border - use a gradient
       to determine the pvec weight.
       Points on the boundary depend on the presence of an adjacent
       processor. */

    for (i=0; i<(nx*ny*nz); i++) thee->pvec[i] = 0.0;

    for (i=0; i<nx; i++) {
        xok = 0.0;
        x = i*hx + xmin;
        if ( (x < (upperCorner[0]-hx/2)) && 
             (x > (lowerCorner[0]+hx/2))
           ) xok = 1.0;
        else if ( (VABS(x - lowerCorner[0]) < VPMGSMALL) && 
                  (bflags[VAPBS_LEFT] == 0)) xok = 1.0;
        else if ((VABS(x - lowerCorner[0]) < VPMGSMALL) && 
                 (bflags[VAPBS_LEFT] == 1)) xok = 0.5;
        else if ((VABS(x - upperCorner[0]) < VPMGSMALL) &&
                 (bflags[VAPBS_RIGHT] == 0)) xok = 1.0;
        else if ((VABS(x - upperCorner[0]) < VPMGSMALL) &&
                 (bflags[VAPBS_RIGHT] == 1)) xok = 0.5;
        else if ((x > (upperCorner[0] + hx/2)) || (x < (lowerCorner[0] - hx/2))) xok = 0.0;
        else if ((x < (upperCorner[0] + hx/2)) || (x > (lowerCorner[0] - hx/2))) {
            x0 = VMAX2(x - hx/2, lowerCorner[0]);
            x1 = VMIN2(x + hx/2, upperCorner[0]);
            xok = VABS(x1-x0)/hx;

            if (xok < 0.0) {
                if (VABS(xok) < VPMGSMALL) xok = 0.0;
                else {
                    Vnm_print(2, "Vpmg_setPart:  fell off x-interval (%1.12E)!\n",
                            xok);
                    VASSERT(0);
                }
            } 
            if (xok > 1.0) {
                if (VABS(xok - 1.0) < VPMGSMALL) xok = 1.0;
                else {
                    Vnm_print(2, "Vpmg_setPart:  fell off x-interval (%1.12E)!\n",
                            xok);
                    VASSERT(0);
                }
            } 

        } else xok = 0.0;
     
        for (j=0; j<ny; j++) {
            yok = 0.0;
            y = j*hy + ymin;
            if ((y < (upperCorner[1]-hy/2)) && (y > (lowerCorner[1]+hy/2))) yok = 1.0;
            else if ((VABS(y - lowerCorner[1]) < VPMGSMALL) && 
                     (bflags[VAPBS_BACK] == 0)) yok = 1.0;
            else if ((VABS(y - lowerCorner[1]) < VPMGSMALL) && 
                     (bflags[VAPBS_BACK] == 1)) yok = 0.5;
            else if ((VABS(y - upperCorner[1]) < VPMGSMALL) &&
                     (bflags[VAPBS_FRONT] == 0)) yok = 1.0;
            else if ((VABS(y - upperCorner[1]) < VPMGSMALL) &&
                     (bflags[VAPBS_FRONT] == 1)) yok = 0.5;
            else if ((y > (upperCorner[1] + hy/2)) || (y < (lowerCorner[1] - hy/2))) yok=0.0;
            else if ((y < (upperCorner[1] + hy/2)) || (y > (lowerCorner[1] - hy/2))){
                y0 = VMAX2(y - hy/2, lowerCorner[1]);
                y1 = VMIN2(y + hy/2, upperCorner[1]);
                yok = VABS(y1-y0)/hy;

                if (yok < 0.0) {
                    if (VABS(yok) < VPMGSMALL) yok = 0.0;
                    else {
                        Vnm_print(2, "Vpmg_setPart:  fell off y-interval (%1.12E)!\n",
                                yok);
                        VASSERT(0);
                    }
                } 
                if (yok > 1.0) {
                    if (VABS(yok - 1.0) < VPMGSMALL) yok = 1.0;
                    else {
                        Vnm_print(2, "Vpmg_setPart:  fell off y-interval (%1.12E)!\n",
                                yok);
                        VASSERT(0);
                    }
                } 
            }
            else yok=0.0;

            for (k=0; k<nz; k++) {
                zok = 0.0; 
                z = k*hzed + zmin;
                if ((z < (upperCorner[2]-hzed/2)) && (z > (lowerCorner[2]+hzed/2))) zok = 1.0;
                else if ((VABS(z - lowerCorner[2]) < VPMGSMALL) && 
                         (bflags[VAPBS_DOWN] == 0)) zok = 1.0;
                else if ((VABS(z - lowerCorner[2]) < VPMGSMALL) && 
                         (bflags[VAPBS_DOWN] == 1)) zok = 0.5;
                else if ((VABS(z - upperCorner[2]) < VPMGSMALL) &&
                         (bflags[VAPBS_UP] == 0)) zok = 1.0;
                else if ((VABS(z - upperCorner[2]) < VPMGSMALL) &&
                         (bflags[VAPBS_UP] == 1)) zok = 0.5;
                else if ((z > (upperCorner[2] + hzed/2)) || (z < (lowerCorner[2] - hzed/2))) zok=0.0;
                else if ((z < (upperCorner[2] + hzed/2)) || (z > (lowerCorner[2] - hzed/2))){
                    z0 = VMAX2(z - hzed/2, lowerCorner[2]);
                    z1 = VMIN2(z + hzed/2, upperCorner[2]);
                    zok = VABS(z1-z0)/hzed;

                    if (zok < 0.0) {
                        if (VABS(zok) < VPMGSMALL) zok = 0.0;
                        else {
                            Vnm_print(2, "Vpmg_setPart:  fell off z-interval (%1.12E)!\n",
                                    zok);
                            VASSERT(0);
                        }
                    } 
                    if (zok > 1.0) {
                        if (VABS(zok - 1.0) < VPMGSMALL) zok = 1.0;
                        else {
                            Vnm_print(2, "Vpmg_setPart:  fell off z-interval (%1.12E)!\n",
                                    zok);
                            VASSERT(0);
                        }
                    } 
                }
                else zok = 0.0;
                
                if (VABS(xok*yok*zok) < VPMGSMALL) thee->pvec[IJK(i,j,k)] = 0.0;
                else thee->pvec[IJK(i,j,k)] = xok*yok*zok;
               
            }
        }
    } 
}

VPUBLIC void Vpmg_unsetPart(Vpmg *thee) {

    int i, nx, ny, nz;
    Vatom *atom;
    Valist *alist;

    VASSERT(thee != VNULL);

    nx = thee->pmgp->nx;
    ny = thee->pmgp->ny;
    nz = thee->pmgp->nz;
    alist = thee->pbe->alist;

    for (i=0; i<(nx*ny*nz); i++) thee->pvec[i] = 1;
    for (i=0; i<Valist_getNumberAtoms(alist); i++) {
        atom = Valist_getAtom(alist, i);
        atom->partID = 1;
    }
}

VPUBLIC void Vpmg_fillArray(Vpmg *thee, double *vec, Vdata_Type type, 
  double parm, Vhal_PBEType pbetype) {

    Vacc *acc = VNULL;
    Vpbe *pbe = VNULL;
    Vgrid *grid = VNULL;
    double position[3], hx, hy, hzed, xmin, ymin, zmin;
    double grad[3], eps, epsp, epss, zmagic;
    int i, j, k, l, nx, ny, nz, ichop;

    pbe = thee->pbe;
    acc = Vpbe_getVacc(pbe);
    nx = thee->pmgp->nx;
    ny = thee->pmgp->ny;
    nz = thee->pmgp->nz;
    hx = thee->pmgp->hx;
    hy = thee->pmgp->hy;
    hzed = thee->pmgp->hzed;
    xmin = thee->pmgp->xmin;
    ymin = thee->pmgp->ymin;
    zmin = thee->pmgp->zmin;
    epsp = Vpbe_getSoluteDiel(pbe);
    epss = Vpbe_getSolventDiel(pbe);
    zmagic = Vpbe_getZmagic(pbe);

    switch (type) {

        case VDT_CHARGE:

            /* Call the coefficient discretization routine */
            if (!thee->filled) 
              Vpmg_fillco(thee, thee->surfMeth, thee->splineWin,
                thee->chargeMeth,
                thee->useDielXMap, thee->dielXMap,
                thee->useDielYMap, thee->dielYMap,
                thee->useDielZMap, thee->dielZMap,
                thee->useKappaMap, thee->kappaMap,
                thee->useChargeMap, thee->chargeMap);
                
            /* Copy the charge array into the argument vector */
            for (i=0; i<nx*ny*nz; i++) vec[i] = thee->fcf[i]/zmagic;
            break;

        case VDT_DIELX:

            /* Call the coefficient discretization routine */
            if (!thee->filled)
              Vpmg_fillco(thee, thee->surfMeth, thee->splineWin,
                thee->chargeMeth,
                thee->useDielXMap, thee->dielXMap,
                thee->useDielYMap, thee->dielYMap,
                thee->useDielZMap, thee->dielZMap,
                thee->useKappaMap, thee->kappaMap,
                thee->useChargeMap, thee->chargeMap);

            /* Copy the x-shifted dielectric array into the argument vector */
            for (i=0; i<nx*ny*nz; i++) vec[i] = thee->a1cf[i];
            break;

        case VDT_DIELY:

            /* Call the coefficient discretization routine */
            if (!thee->filled)
              Vpmg_fillco(thee, thee->surfMeth, thee->splineWin,
                thee->chargeMeth,
                thee->useDielXMap, thee->dielXMap,
                thee->useDielYMap, thee->dielYMap,
                thee->useDielZMap, thee->dielZMap,
                thee->useKappaMap, thee->kappaMap,
                thee->useChargeMap, thee->chargeMap);

            /* Copy the y-shifted dielectric array into the argument vector */
            for (i=0; i<nx*ny*nz; i++) vec[i] = thee->a2cf[i];
            break;

        case VDT_DIELZ:

            /* Call the coefficient discretization routine */
            if (!thee->filled)
              Vpmg_fillco(thee, thee->surfMeth, thee->splineWin,
                thee->chargeMeth,
                thee->useDielXMap, thee->dielXMap,
                thee->useDielYMap, thee->dielYMap,
                thee->useDielZMap, thee->dielZMap,
                thee->useKappaMap, thee->kappaMap,
                thee->useChargeMap, thee->chargeMap);

            /* Copy the z-shifted dielectric array into the argument vector */
            for (i=0; i<nx*ny*nz; i++) vec[i] = thee->a3cf[i];
            break;

        case VDT_KAPPA:

            /* Call the coefficient discretization routine */
            if (!thee->filled)
              Vpmg_fillco(thee, thee->surfMeth, thee->splineWin,
                thee->chargeMeth,
                thee->useDielXMap, thee->dielXMap,
                thee->useDielYMap, thee->dielYMap,
                thee->useDielZMap, thee->dielZMap,
                thee->useKappaMap, thee->kappaMap,
                thee->useChargeMap, thee->chargeMap);

            /* Copy the kappa array into the argument vector */
            for (i=0; i<nx*ny*nz; i++) vec[i] = thee->ccf[i];
            break;

        case VDT_POT:

            for (i=0; i<nx*ny*nz; i++) vec[i] = thee->u[i];
            break;

        case VDT_SMOL:
 
            for (k=0; k<nz; k++) {
                for (j=0; j<ny; j++) {
                    for (i=0; i<nx; i++) {

                        position[0] = i*hx + xmin;
                        position[1] = j*hy + ymin;
                        position[2] = k*hzed + zmin;

                        vec[IJK(i,j,k)] = (Vacc_molAcc(acc,position,parm));
                    }
                }
            }
            break;

        case VDT_SSPL:

            for (k=0; k<nz; k++) {
                for (j=0; j<ny; j++) {
                    for (i=0; i<nx; i++) {

                        position[0] = i*hx + xmin;
                        position[1] = j*hy + ymin;
                        position[2] = k*hzed + zmin;

                        vec[IJK(i,j,k)] = Vacc_splineAcc(acc,position,parm,0);
                    }
                }   
            }
            break;

        case VDT_VDW:

            for (k=0; k<nz; k++) {
                for (j=0; j<ny; j++) {
                    for (i=0; i<nx; i++) {

                        position[0] = i*hx + xmin;
                        position[1] = j*hy + ymin;
                        position[2] = k*hzed + zmin;

                        vec[IJK(i,j,k)] = Vacc_vdwAcc(acc,position);
                    }
                }
            }
            break;

        case VDT_IVDW:

            for (k=0; k<nz; k++) {
                for (j=0; j<ny; j++) {
                    for (i=0; i<nx; i++) {

                        position[0] = i*hx + xmin;
                        position[1] = j*hy + ymin;
                        position[2] = k*hzed + zmin;

                        vec[IJK(i,j,k)] = Vacc_ivdwAcc(acc,position,parm);
                    }
                }
            }
            break;

        case VDT_LAP:

            grid = Vgrid_ctor(nx, ny, nz, hx, hy, hzed, xmin, ymin, zmin,
              thee->u);
            for (k=0; k<nz; k++) {
                for (j=0; j<ny; j++) {
                    for (i=0; i<nx; i++) {

                        if ((k==0) || (k==nz) ||
                            (j==0) || (j==ny) ||
                            (i==0) || (j==nz)) {

                            vec[IJK(i,j,k)] = 0;

                        } else { 

                                position[0] = i*hx + xmin;
                                position[1] = j*hy + ymin;
                                position[2] = k*hzed + zmin;
                                VASSERT(Vgrid_curvature(grid,position, 1,
                                  &(vec[IJK(i,j,k)])));
                        }
                    }
                }
            }
            Vgrid_dtor(&grid);
            break;

        case VDT_EDENS:

            grid = Vgrid_ctor(nx, ny, nz, hx, hy, hzed, xmin, ymin, zmin,
              thee->u);
            for (k=0; k<nz; k++) {
                for (j=0; j<ny; j++) {
                    for (i=0; i<nx; i++) {

                        position[0] = i*hx + xmin;
                        position[1] = j*hy + ymin;
                        position[2] = k*hzed + zmin;
                        VASSERT(Vgrid_gradient(grid, position, grad));
                        eps = epsp + (epss-epsp)*Vacc_molAcc(acc, position, 
                          pbe->solventRadius);
                        vec[IJK(i,j,k)] = 0.0;
                        for (l=0; l<3; l++) 
                          vec[IJK(i,j,k)] += eps*VSQR(grad[l]);
                    }
                }
            }
            Vgrid_dtor(&grid);
            break;

        case VDT_NDENS:

            for (k=0; k<nz; k++) {
                for (j=0; j<ny; j++) {
                    for (i=0; i<nx; i++) {

                        position[0] = i*hx + xmin;
                        position[1] = j*hy + ymin;
                        position[2] = k*hzed + zmin;
                        vec[IJK(i,j,k)] = 0.0;
                        if ( VABS(Vacc_ivdwAcc(acc, 
                               position, pbe->maxIonRadius) - 1.0) < VSMALL) {
                            for (l=0; l<pbe->numIon; l++) {
                                if (pbetype == PBE_NPBE) {
                                    vec[IJK(i,j,k)] += (pbe->ionConc[l]
                                        * Vcap_exp(-pbe->ionQ[l]*thee->u[IJK(i,j,k)], 
                                        &ichop));
                                } else if (pbetype == PBE_LPBE){
                                    vec[IJK(i,j,k)] += (pbe->ionConc[l]
                                        * (1 - pbe->ionQ[l]*thee->u[IJK(i,j,k)]));
                                }
                            }
                        } 
                    }
                }
            }
            break;

        case VDT_QDENS:

            for (k=0; k<nz; k++) {
                for (j=0; j<ny; j++) {
                    for (i=0; i<nx; i++) {

                        position[0] = i*hx + xmin;
                        position[1] = j*hy + ymin;
                        position[2] = k*hzed + zmin;
                        vec[IJK(i,j,k)] = 0.0;
                        if ( VABS(Vacc_ivdwAcc(acc, 
                               position, pbe->maxIonRadius) - 1.0) < VSMALL) {
                            for (l=0; l<pbe->numIon; l++) {
                                if (pbetype == PBE_NPBE) {
                                    vec[IJK(i,j,k)] += (pbe->ionConc[l] 
                                        * pbe->ionQ[l]
                                        * Vcap_exp(-pbe->ionQ[l]*thee->u[IJK(i,j,k)],
                                        &ichop));
                                } else if (pbetype == PBE_LPBE) {
                                    vec[IJK(i,j,k)] += (pbe->ionConc[l] 
                                        * pbe->ionQ[l]
                                        * (1 - pbe->ionQ[l]*thee->u[IJK(i,j,k)]));
                                }
                            }
                        }
                    }
                }
            }
            break;

        default:

            Vnm_print(2, "main:  Bogus data type (%d)!\n", type);
            break;

    }

}

VPUBLIC double Vpmg_energy(Vpmg *thee, int extFlag) {

    double totEnergy = 0.0;
    double dielEnergy = 0.0;
    double qmEnergy = 0.0;
    double qfEnergy = 0.0;
    double npEnergy = 0.0;

    VASSERT(thee != VNULL);
    /* VASSERT(thee->filled); */

    Vnm_print(0, "Vpmg_energy:  calculating apolar energy\n");
    npEnergy = Vpmg_npEnergy(thee, extFlag);
    Vnm_print(0, "Vpmg_energy:  npEnergy = %1.12E kT\n", npEnergy);

    if ((thee->pmgp->nonlin) && (Vpbe_getBulkIonicStrength(thee->pbe) > 0.)) {
        Vnm_print(0, "Vpmg_energy:  calculating full PBE energy\n");
        qmEnergy = Vpmg_qmEnergy(thee, extFlag);
        Vnm_print(0, "Vpmg_energy:  qmEnergy = %1.12E kT\n", qmEnergy);
        qfEnergy = Vpmg_qfEnergy(thee, extFlag);
        Vnm_print(0, "Vpmg_energy:  qfEnergy = %1.12E kT\n", qfEnergy);
        dielEnergy = Vpmg_dielEnergy(thee, extFlag);
        Vnm_print(0, "Vpmg_energy:  dielEnergy = %1.12E kT\n", dielEnergy);
        totEnergy = qfEnergy - dielEnergy - qmEnergy;
    } else {
        Vnm_print(0, "Vpmg_energy:  calculating only q-phi energy\n");
        qfEnergy = Vpmg_qfEnergy(thee, extFlag);
        Vnm_print(0, "Vpmg_energy:  qfEnergy = %1.12E kT\n", qfEnergy);
        totEnergy = 0.5*qfEnergy;
    }

    return totEnergy;

}

VPUBLIC double Vpmg_dielEnergy(Vpmg *thee, int extFlag) {

    double hx, hy, hzed, energy, nrgx, nrgy, nrgz, pvecx, pvecy, pvecz;
    int i, j, k, nx, ny, nz;
 
    VASSERT(thee != VNULL);

    /* Get the mesh information */
    nx = thee->pmgp->nx;
    ny = thee->pmgp->ny;
    nz = thee->pmgp->nz;
    hx = thee->pmgp->hx;
    hy = thee->pmgp->hy;
    hzed = thee->pmgp->hzed;

    energy = 0.0;

    /* Refill the dieletric coefficient arrays */
    if (!thee->filled) Vpmg_fillco(thee, 
      thee->surfMeth, thee->splineWin, thee->chargeMeth,
      thee->useDielXMap, thee->dielXMap,
      thee->useDielYMap, thee->dielYMap,
      thee->useDielZMap, thee->dielZMap,
      thee->useKappaMap, thee->kappaMap,
      thee->useChargeMap, thee->chargeMap);

    for (k=0; k<(nz-1); k++) {
        for (j=0; j<(ny-1); j++) {
            for (i=0; i<(nx-1); i++) {
                pvecx = 0.5*(thee->pvec[IJK(i,j,k)]+thee->pvec[IJK(i+1,j,k)]);
                pvecy = 0.5*(thee->pvec[IJK(i,j,k)]+thee->pvec[IJK(i,j+1,k)]);
                pvecz = 0.5*(thee->pvec[IJK(i,j,k)]+thee->pvec[IJK(i,j,k+1)]);
                nrgx = thee->a1cf[IJK(i,j,k)]*pvecx
                  * VSQR((thee->u[IJK(i,j,k)]-thee->u[IJK(i+1,j,k)])/hx);
                nrgy = thee->a2cf[IJK(i,j,k)]*pvecy
                  * VSQR((thee->u[IJK(i,j,k)]-thee->u[IJK(i,j+1,k)])/hy);
                nrgz = thee->a3cf[IJK(i,j,k)]*pvecz
                  * VSQR((thee->u[IJK(i,j,k)]-thee->u[IJK(i,j,k+1)])/hzed);
                energy += (nrgx + nrgy + nrgz);
            }
        }
    }

    energy = 0.5*energy*hx*hy*hzed;
    energy = energy/Vpbe_getZmagic(thee->pbe);

    if (extFlag == 1) energy += (thee->extDiEnergy);

    return energy;
}

VPUBLIC double Vpmg_dielGradNorm(Vpmg *thee) {

    double hx, hy, hzed, energy, nrgx, nrgy, nrgz, pvecx, pvecy, pvecz;
    int i, j, k, nx, ny, nz;
 
    VASSERT(thee != VNULL);

    /* Get the mesh information */
    nx = thee->pmgp->nx;
    ny = thee->pmgp->ny;
    nz = thee->pmgp->nz;
    hx = thee->pmgp->hx;
    hy = thee->pmgp->hy;
    hzed = thee->pmgp->hzed;

    energy = 0.0;

    /* Refill the dieletric coefficient arrays */
    if (!thee->filled) Vpmg_fillco(thee, 
      thee->surfMeth, thee->splineWin, thee->chargeMeth,
      thee->useDielXMap, thee->dielXMap,
      thee->useDielYMap, thee->dielYMap,
      thee->useDielZMap, thee->dielZMap,
      thee->useKappaMap, thee->kappaMap,
      thee->useChargeMap, thee->chargeMap);

    for (k=1; k<nz; k++) {
        for (j=1; j<ny; j++) {
            for (i=1; i<nx; i++) {
                pvecx = 0.5*(thee->pvec[IJK(i,j,k)]+thee->pvec[IJK(i-1,j,k)]);
                pvecy = 0.5*(thee->pvec[IJK(i,j,k)]+thee->pvec[IJK(i,j-1,k)]);
                pvecz = 0.5*(thee->pvec[IJK(i,j,k)]+thee->pvec[IJK(i,j,k-1)]);
                nrgx = pvecx
                 * VSQR((thee->a1cf[IJK(i,j,k)]-thee->a1cf[IJK(i-1,j,k)])/hx);
                nrgy = pvecy
                 * VSQR((thee->a2cf[IJK(i,j,k)]-thee->a2cf[IJK(i,j-1,k)])/hy);
                nrgz = pvecz
                 * VSQR((thee->a3cf[IJK(i,j,k)]-thee->a3cf[IJK(i,j,k-1)])/hzed);
                energy += VSQRT(nrgx + nrgy + nrgz);
            }
        }
    }

    energy = energy*hx*hy*hzed;

    return energy;
}

VPUBLIC double Vpmg_npEnergy(Vpmg *thee, int extFlag) {

    double area, energy, epsp, epss, gamma, temp;

    epsp = Vpbe_getSoluteDiel(thee->pbe);
    epss = Vpbe_getSolventDiel(thee->pbe);
    gamma = Vpbe_getGamma(thee->pbe);
    temp = Vpbe_getTemperature(thee->pbe);
    gamma = gamma/(1e-3*Vunit_Na*Vunit_kb*temp);

    if ((VABS(epsp-epss) < VSMALL) || (gamma < VSMALL)) {
        return 0.0;
    } 

    area = Vpmg_dielGradNorm(thee);
    energy = gamma*area/(epss-epsp);
   
    if (extFlag == 1) energy += (thee->extNpEnergy); 

    return energy;

}
    
VPUBLIC double Vpmg_qmEnergy(Vpmg *thee, int extFlag) {

    double hx, hy, hzed, energy, ionConc[MAXION], ionRadii[MAXION];
    double ionQ[MAXION], zkappa2, ionstr, zks2;
    int i, j, nx, ny, nz, nion, ichop, nchop;
 
    VASSERT(thee != VNULL);

    /* Get the mesh information */
    nx = thee->pmgp->nx;
    ny = thee->pmgp->ny;
    nz = thee->pmgp->nz;
    hx = thee->pmgp->hx;
    hy = thee->pmgp->hy;
    hzed = thee->pmgp->hzed;
    zkappa2 = Vpbe_getZkappa2(thee->pbe);
    ionstr = Vpbe_getBulkIonicStrength(thee->pbe);

    /* Bail if we're at zero ionic strength */
    if (zkappa2 < VSMALL) {

#ifndef VAPBSQUIET
        Vnm_print(0, "Vpmg_qmEnergy:  Zero energy for zero ionic strength!\n");
#endif

        return 0.0;
    }
    zks2 = 0.5*zkappa2/ionstr;

    /* Because PMG seems to overwrite some of the coefficient arrays... */
    if (!thee->filled) Vpmg_fillco(thee, 
      thee->surfMeth, thee->splineWin, thee->chargeMeth,
      thee->useDielXMap, thee->dielXMap,
      thee->useDielYMap, thee->dielYMap,
      thee->useDielZMap, thee->dielZMap,
      thee->useKappaMap, thee->kappaMap,
      thee->useChargeMap, thee->chargeMap);

    energy = 0.0;
    nchop = 0;
    Vpbe_getIons(thee->pbe, &nion, ionConc, ionRadii, ionQ);
    if (thee->pmgp->nonlin) {
        Vnm_print(0, "Vpmg_qmEnergy:  Calculating nonlinear energy\n");
        for (i=0; i<(nx*ny*nz); i++) {
            if (thee->pvec[i]*thee->ccf[i] > VSMALL) {
                for (j=0; j<nion; j++) {
                    energy += (thee->pvec[i]*thee->ccf[i]*zks2
                      * ionConc[j] * VSQR(ionQ[j]) 
                      * (Vcap_exp(-ionQ[j]*thee->u[i], &ichop)-1.0));
                    nchop += ichop;
                }
            }
        }
        if (nchop > 0) Vnm_print(2, "Vpmg_qmEnergy:  Chopped EXP %d times!\n",
          nchop);
    } else {
        /* Zkappa2 OK here b/c LPBE approx */
        Vnm_print(0, "Vpmg_qmEnergy:  Calculating linear energy\n");
        for (i=0; i<(nx*ny*nz); i++) {
            if (thee->pvec[i]*thee->ccf[i] > VSMALL) 
              energy += (thee->pvec[i]*zkappa2*thee->ccf[i]*VSQR(thee->u[i]));
        }
        energy = 0.5*energy;
    }
    energy = energy*hx*hy*hzed;
    energy = energy/Vpbe_getZmagic(thee->pbe);

    if (extFlag == 1) energy += thee->extQmEnergy;

    return energy;
}
    
VPUBLIC double Vpmg_qfEnergy(Vpmg *thee, int extFlag) {

    double energy = 0.0;

    VASSERT(thee != VNULL);

    if ((thee->useChargeMap) || (thee->chargeMeth == VCM_BSPL2)) { 
        energy = Vpmg_qfEnergyVolume(thee, extFlag); 
    } else { 
        energy = Vpmg_qfEnergyPoint(thee, extFlag); 
    } 
 
    return energy;
}

VPRIVATE double Vpmg_qfEnergyPoint(Vpmg *thee, int extFlag) {

    int iatom, nx, ny, nz, ihi, ilo, jhi, jlo, khi, klo;
    double xmax, ymax, zmax, xmin, ymin, zmin, hx, hy, hzed, ifloat, jfloat;
    double charge, kfloat, dx, dy, dz, energy, uval, *position;
    double *u;
    double *pvec;
    Valist *alist;
    Vatom *atom; 
    Vpbe *pbe;

    pbe = thee->pbe;
    alist = pbe->alist;
    VASSERT(alist != VNULL);

    /* Get the mesh information */
    nx = thee->pmgp->nx;
    ny = thee->pmgp->ny;
    nz = thee->pmgp->nz;
    hx = thee->pmgp->hx;
    hy = thee->pmgp->hy;
    hzed = thee->pmgp->hzed;
    xmax = thee->pmgp->xmax;
    ymax = thee->pmgp->ymax;
    zmax = thee->pmgp->zmax;
    xmin = thee->pmgp->xmin;
    ymin = thee->pmgp->ymin;
    zmin = thee->pmgp->zmin;

    u = thee->u;
    pvec = thee->pvec;
  
    energy = 0.0;

    for (iatom=0; iatom<Valist_getNumberAtoms(alist); iatom++) {

        /* Get atomic information */
        atom = Valist_getAtom(alist, iatom);

        position = Vatom_getPosition(atom);
        charge = Vatom_getCharge(atom);

        /* Figure out which vertices we're next to */
        ifloat = (position[0] - xmin)/hx;
        jfloat = (position[1] - ymin)/hy;
        kfloat = (position[2] - zmin)/hzed;
        ihi = (int)ceil(ifloat);
        ilo = (int)floor(ifloat);
        jhi = (int)ceil(jfloat);
        jlo = (int)floor(jfloat);
        khi = (int)ceil(kfloat);
        klo = (int)floor(kfloat);

        if (atom->partID > 0) {

            if ((ihi<nx) && (jhi<ny) && (khi<nz) &&
                (ilo>=0) && (jlo>=0) && (klo>=0)) {

                /* Now get trilinear interpolation constants */
                dx = ifloat - (double)(ilo);
                dy = jfloat - (double)(jlo);
                dz = kfloat - (double)(klo);
                uval =  
                  dx*dy*dz*u[IJK(ihi,jhi,khi)]
                + dx*(1.0-dy)*dz*u[IJK(ihi,jlo,khi)]
                + dx*dy*(1.0-dz)*u[IJK(ihi,jhi,klo)]
                + dx*(1.0-dy)*(1.0-dz)*u[IJK(ihi,jlo,klo)]
                + (1.0-dx)*dy*dz*u[IJK(ilo,jhi,khi)]
                + (1.0-dx)*(1.0-dy)*dz*u[IJK(ilo,jlo,khi)]
                + (1.0-dx)*dy*(1.0-dz)*u[IJK(ilo,jhi,klo)]
                + (1.0-dx)*(1.0-dy)*(1.0-dz)*u[IJK(ilo,jlo,klo)];
                energy += (uval*charge*atom->partID);
            } else if (thee->pmgp->bcfl != BCFL_FOCUS) {
                Vnm_print(2, "Vpmg_qfEnergy:  Atom #%d at (%4.3f, %4.3f, \
%4.3f) is off the mesh (ignoring)!\n",
                iatom, position[0], position[1], position[2]);
            }
        } 
    }

    if (extFlag) energy += thee->extQfEnergy;
 
    return energy;
}

VPUBLIC double Vpmg_qfAtomEnergy(Vpmg *thee, Vatom *atom) {

    int nx, ny, nz, ihi, ilo, jhi, jlo, khi, klo;
    double xmax, xmin, ymax, ymin, zmax, zmin, hx, hy, hzed, ifloat, jfloat;
    double charge, kfloat, dx, dy, dz, energy, uval, *position;
    double *u;


    /* Get the mesh information */
    nx = thee->pmgp->nx;
    ny = thee->pmgp->ny;
    nz = thee->pmgp->nz;
    hx = thee->pmgp->hx;
    hy = thee->pmgp->hy;
    hzed = thee->pmgp->hzed;
    xmax = thee->xf[nx-1];
    ymax = thee->yf[ny-1];
    zmax = thee->zf[nz-1];
    xmin = thee->xf[0];
    ymin = thee->yf[0];
    zmin = thee->zf[0];

    u = thee->u;

    energy = 0.0;


    position = Vatom_getPosition(atom);
    charge = Vatom_getCharge(atom);

    /* Figure out which vertices we're next to */
    ifloat = (position[0] - xmin)/hx;
    jfloat = (position[1] - ymin)/hy;
    kfloat = (position[2] - zmin)/hzed;
    ihi = (int)ceil(ifloat);
    ilo = (int)floor(ifloat);
    jhi = (int)ceil(jfloat);
    jlo = (int)floor(jfloat);
    khi = (int)ceil(kfloat);
    klo = (int)floor(kfloat);

    if (atom->partID > 0) {

        if ((ihi<nx) && (jhi<ny) && (khi<nz) &&
            (ilo>=0) && (jlo>=0) && (klo>=0)) {

            /* Now get trilinear interpolation constants */
            dx = ifloat - (double)(ilo);
            dy = jfloat - (double)(jlo);
            dz = kfloat - (double)(klo);
            uval =
              dx*dy*dz*u[IJK(ihi,jhi,khi)]
            + dx*(1.0-dy)*dz*u[IJK(ihi,jlo,khi)]
            + dx*dy*(1.0-dz)*u[IJK(ihi,jhi,klo)]
            + dx*(1.0-dy)*(1.0-dz)*u[IJK(ihi,jlo,klo)]
            + (1.0-dx)*dy*dz*u[IJK(ilo,jhi,khi)]
            + (1.0-dx)*(1.0-dy)*dz*u[IJK(ilo,jlo,khi)]
            + (1.0-dx)*dy*(1.0-dz)*u[IJK(ilo,jhi,klo)]
            + (1.0-dx)*(1.0-dy)*(1.0-dz)*u[IJK(ilo,jlo,klo)];
            energy += (uval*charge*atom->partID);
        } else if (thee->pmgp->bcfl != BCFL_FOCUS) {
            Vnm_print(2, "Vpmg_qfAtomEnergy:  Atom at (%4.3f, %4.3f, \
%4.3f) is off the mesh (ignoring)!\n",
            position[0], position[1], position[2]);
        }
    } 

    return energy; 
}
    
VPRIVATE double Vpmg_qfEnergyVolume(Vpmg *thee, int extFlag) {

    double hx, hy, hzed, energy;
    int i, nx, ny, nz;
 
    VASSERT(thee != VNULL);

    /* Get the mesh information */
    nx = thee->pmgp->nx;
    ny = thee->pmgp->ny;
    nz = thee->pmgp->nz;
    hx = thee->pmgp->hx;
    hy = thee->pmgp->hy;
    hzed = thee->pmgp->hzed;

    /* Because PMG seems to overwrite some of the coefficient arrays... 
     * NAB:  it would be useful to figure out which particular arrays need to
     * be recreated. */
    if (!thee->filled) Vpmg_fillco(thee, 
      thee->surfMeth, thee->splineWin, thee->chargeMeth,
      thee->useDielXMap, thee->dielXMap,
      thee->useDielYMap, thee->dielYMap,
      thee->useDielZMap, thee->dielZMap,
      thee->useKappaMap, thee->kappaMap,
      thee->useChargeMap, thee->chargeMap);

    energy = 0.0;
    Vnm_print(0, "Vpmg_qfEnergyVolume:  Calculating energy\n");
    for (i=0; i<(nx*ny*nz); i++) {
        energy += (thee->pvec[i]*thee->u[i]*thee->fcf[i]);
    }
    energy = energy*hx*hy*hzed/Vpbe_getZmagic(thee->pbe);

    if (extFlag == 1) energy += thee->extQfEnergy;

    return energy;
}

VPUBLIC int Vpmg_ctor2(Vpmg *thee, Vpmgp *pmgp, Vpbe *pbe, int focusFlag,
        Vpmg *pmgOLD, MGparm *mgparm, PBEparm_calcEnergy energyFlag) {

    int i, j, nion;
    double ionConc[MAXION], ionQ[MAXION], ionRadii[MAXION], zkappa2, zks2;
    double ionstr, partMin[3], partMax[3];

    /* Get the parameters */
    VASSERT(pmgp != VNULL);
    VASSERT(pbe != VNULL);
    thee->pmgp = pmgp;
    thee->pbe = pbe;

    /* Set up the memory */
    thee->vmem = Vmem_ctor("APBS:VPMG");

    /* Calculate storage requirements */
    F77MGSZ(&(thee->pmgp->mgcoar), &(thee->pmgp->mgdisc),
          &(thee->pmgp->mgsolv), &(thee->pmgp->nx), &(thee->pmgp->ny),
      &(thee->pmgp->nz),
      &(thee->pmgp->nlev), &(thee->pmgp->nxc), &(thee->pmgp->nyc),
      &(thee->pmgp->nzc), &(thee->pmgp->nf), &(thee->pmgp->nc),
      &(thee->pmgp->narr), &(thee->pmgp->narrc), &(thee->pmgp->n_rpc),
      &(thee->pmgp->n_iz), &(thee->pmgp->n_ipc), &(thee->pmgp->nrwk),
      &(thee->pmgp->niwk));

    /* We need some additional storage if: nonlinear & newton OR cgmg */
    if (((thee->pmgp->nonlin == 1) && (thee->pmgp->meth == 1))
        || (thee->pmgp->meth == 0)) { thee->pmgp->nrwk += (2*(thee->pmgp->nf));
    }

    Vnm_print(0, "Vpmg_ctor2: PMG chose nx = %d, ny = %d, nz = %d, nlev = %d\n",
       thee->pmgp->nx, thee->pmgp->ny, thee->pmgp->nz, thee->pmgp->nlev);

    /* Allocate boundary storage */
    thee->gxcf = (double *)Vmem_malloc(thee->vmem,
      10*(thee->pmgp->ny)*(thee->pmgp->nz), sizeof(double));
    thee->gycf = (double *)Vmem_malloc(thee->vmem,
      10*(thee->pmgp->nx)*(thee->pmgp->nz), sizeof(double));
    thee->gzcf = (double *)Vmem_malloc(thee->vmem,
      10*(thee->pmgp->nx)*(thee->pmgp->ny), sizeof(double));

    /* Allocate partition vector storage */
    thee->pvec = (double *)Vmem_malloc(thee->vmem,
      (thee->pmgp->nx)*(thee->pmgp->ny)*(thee->pmgp->nz), sizeof(double));

    /* Allocate remaining storage */
    thee->iparm = (int *)Vmem_malloc(thee->vmem, 100, sizeof(int));
    thee->rparm = (double *)Vmem_malloc(thee->vmem, 100, sizeof(double));
    thee->iwork = (int *)Vmem_malloc(thee->vmem, thee->pmgp->niwk,
      sizeof(int));
    thee->rwork = (double *)Vmem_malloc(thee->vmem, thee->pmgp->nrwk,
      sizeof(double));
    thee->a1cf = (double *)Vmem_malloc(thee->vmem, thee->pmgp->narr,
      sizeof(double));
    thee->a2cf = (double *)Vmem_malloc(thee->vmem, thee->pmgp->narr,
      sizeof(double));
    thee->a3cf = (double *)Vmem_malloc(thee->vmem, thee->pmgp->narr,
      sizeof(double));
    thee->ccf = (double *)Vmem_malloc(thee->vmem, thee->pmgp->narr,
      sizeof(double));
    thee->fcf = (double *)Vmem_malloc(thee->vmem, thee->pmgp->narr,
      sizeof(double));
    thee->tcf = (double *)Vmem_malloc(thee->vmem, thee->pmgp->narr,
      sizeof(double));
    thee->u = (double *)Vmem_malloc(thee->vmem, thee->pmgp->narr,
      sizeof(double));
    thee->xf = (double *)Vmem_malloc(thee->vmem, 5*(thee->pmgp->nx),
      sizeof(double));
    thee->yf = (double *)Vmem_malloc(thee->vmem, 5*(thee->pmgp->ny),
      sizeof(double));
    thee->zf = (double *)Vmem_malloc(thee->vmem, 5*(thee->pmgp->nz),
      sizeof(double));

    if (focusFlag) {
        /* Overwrite any default or user-specified boundary condition
         * arguments; we are now committed to a calculation via focusing */
        if (thee->pmgp->bcfl != BCFL_FOCUS) {
            Vnm_print(2, 
            "Vpmg_ctor2: reset boundary condition flag to BCFL_FOCUS!\n");
            thee->pmgp->bcfl = BCFL_FOCUS;
        }

        /* Fill boundaries */
        Vnm_print(0, "Vpmg_ctor2:  Filling boundary with old solution!\n");
        focusFillBound(thee, pmgOLD);

        /* Ignore old maps */
        if (pmgOLD->useDielXMap || pmgOLD->useDielYMap || pmgOLD->useDielZMap ||
            pmgOLD->useKappaMap || pmgOLD->useChargeMap)
           Vnm_print(2, "Vpmg_ctor2:  WARNING!  Ignoring coefficient and charge distribution maps during focusing!\n");

        /* Calculate energetic contributions from region outside focusing
         * domain */
        if (energyFlag != PCE_NO) {

            if (mgparm->type == MCT_PAR) {

                for (j=0; j<3; j++) {
                    partMin[j] = mgparm->center[j] 
                        + mgparm->partDisjCenterShift[j]
                        - 0.5*mgparm->partDisjLength[j];
                    partMax[j] = mgparm->center[j] 
                        + mgparm->partDisjCenterShift[j]
                        + 0.5*mgparm->partDisjLength[j];
                }

            } else {
                for (j=0; j<3; j++) {
                    partMin[j] = mgparm->center[j] - 0.5*mgparm->glen[j];
                    partMax[j] = mgparm->center[j] + 0.5*mgparm->glen[j];
                }
            }
            extEnergy(thee, pmgOLD, energyFlag, partMin, partMax, 
                    mgparm->partDisjOwnSide);
        }
        /* Destroy old Vpmg object */
        Vpmg_dtor(&pmgOLD);

    } else {

        /* Ignore external energy contributions */
        thee->extQmEnergy = 0;
        thee->extDiEnergy = 0;
        thee->extQfEnergy = 0;
        thee->extNpEnergy = 0;
    }

    /* Plop some of the parameters into the iparm and rparm arrays */
    F77PACKMG(thee->iparm, thee->rparm, &(thee->pmgp->nrwk),
&(thee->pmgp->niwk),
      &(thee->pmgp->nx), &(thee->pmgp->ny), &(thee->pmgp->nz),
      &(thee->pmgp->nlev), &(thee->pmgp->nu1), &(thee->pmgp->nu2),
      &(thee->pmgp->mgkey), &(thee->pmgp->itmax), &(thee->pmgp->istop),
      &(thee->pmgp->ipcon), &(thee->pmgp->nonlin), &(thee->pmgp->mgsmoo),
      &(thee->pmgp->mgprol), &(thee->pmgp->mgcoar), &(thee->pmgp->mgsolv),
      &(thee->pmgp->mgdisc), &(thee->pmgp->iinfo), &(thee->pmgp->errtol),
      &(thee->pmgp->ipkey), &(thee->pmgp->omegal), &(thee->pmgp->omegan),
      &(thee->pmgp->irite), &(thee->pmgp->iperf));


    /* Initialize ion concentrations and valencies in PMG routines */
    zkappa2 = Vpbe_getZkappa2(thee->pbe);
    ionstr = Vpbe_getBulkIonicStrength(thee->pbe);
    if (ionstr > 0.0) zks2 = 0.5/ionstr;
    else zks2 = 0.0;
    Vpbe_getIons(thee->pbe, &nion, ionConc, ionRadii, ionQ);
    for (i=0; i<nion; i++) {
        ionConc[i] = zks2 * ionConc[i];
    }
    F77MYPDEFINIT(&nion, ionQ, ionConc);

    /* Turn off restriction of observable calculations to a specific
     * partition */
    Vpmg_unsetPart(thee);

    /* The coefficient arrays have not been filled */
    thee->filled = 0;

    return 1;
}

VPRIVATE void focusFillBound(Vpmg *thee, Vpmg *pmgOLD) {

    Vpbe *pbe;
    double hxOLD, hyOLD, hzOLD, xminOLD, yminOLD, zminOLD, xmaxOLD, ymaxOLD;
    double zmaxOLD;
    int nxOLD, nyOLD, nzOLD;
    double hxNEW, hyNEW, hzNEW, xminNEW, yminNEW, zminNEW, xmaxNEW, ymaxNEW;
    double zmaxNEW;
    int nxNEW, nyNEW, nzNEW;
    int i, j, k, ihi, ilo, jhi, jlo, khi, klo, nx, ny, nz;
    double x, y, z, dx, dy, dz, ifloat, jfloat, kfloat, uval;
    double eps_w, T, pre1, xkappa, size, *apos, charge, pos[3];

    /* Calculate new problem dimensions */
    hxNEW = thee->pmgp->hx;
    hyNEW = thee->pmgp->hy;
    hzNEW = thee->pmgp->hzed;
    nx =  thee->pmgp->nx;
    ny =  thee->pmgp->ny;
    nz =  thee->pmgp->nz;
    nxNEW = thee->pmgp->nx;
    nyNEW = thee->pmgp->ny;
    nzNEW = thee->pmgp->nz;
    xminNEW = thee->pmgp->xcent - ((double)(nxNEW-1)*hxNEW)/2.0;
    xmaxNEW = thee->pmgp->xcent + ((double)(nxNEW-1)*hxNEW)/2.0;
    yminNEW = thee->pmgp->ycent - ((double)(nyNEW-1)*hyNEW)/2.0;
    ymaxNEW = thee->pmgp->ycent + ((double)(nyNEW-1)*hyNEW)/2.0;
    zminNEW = thee->pmgp->zcent - ((double)(nzNEW-1)*hzNEW)/2.0;
    zmaxNEW = thee->pmgp->zcent + ((double)(nzNEW-1)*hzNEW)/2.0;

    /* Relevant old problem parameters */
    hxOLD = pmgOLD->pmgp->hx;
    hyOLD = pmgOLD->pmgp->hy;
    hzOLD = pmgOLD->pmgp->hzed;
    nxOLD = pmgOLD->pmgp->nx;
    nyOLD = pmgOLD->pmgp->ny;
    nzOLD = pmgOLD->pmgp->nz;
    xminOLD = pmgOLD->pmgp->xcent - ((double)(nxOLD-1)*hxOLD)/2.0;
    xmaxOLD = pmgOLD->pmgp->xcent + ((double)(nxOLD-1)*hxOLD)/2.0;
    yminOLD = pmgOLD->pmgp->ycent - ((double)(nyOLD-1)*hyOLD)/2.0;
    ymaxOLD = pmgOLD->pmgp->ycent + ((double)(nyOLD-1)*hyOLD)/2.0;
    zminOLD = pmgOLD->pmgp->zcent - ((double)(nzOLD-1)*hzOLD)/2.0;
    zmaxOLD = pmgOLD->pmgp->zcent + ((double)(nzOLD-1)*hzOLD)/2.0;

    /* BOUNDARY CONDITION SETUP FOR POINTS OFF OLD MESH:
     * For each "atom" (only one for bcfl=1), we use the following formula to
     * calculate the boundary conditions:
     *    g(x) = \frac{q e_c}{4*\pi*\eps_0*\eps_w*k_b*T}
     *          * \frac{exp(-xkappa*(d - a))}{1+xkappa*a}
     *          * 1/d
     * where d = ||x - x_0|| (in m) and a is the size of the atom (in m).
     * We only need to evaluate some of these prefactors once:
     *    pre1 = \frac{e_c}{4*\pi*\eps_0*\eps_w*k_b*T}
     * which gives the potential as
     *    g(x) = pre1 * q/d * \frac{exp(-xkappa*(d - a))}{1+xkappa*a}
     */
    pbe = thee->pbe;
    eps_w = Vpbe_getSolventDiel(pbe);           /* Dimensionless */
    T = Vpbe_getTemperature(pbe);               /* K             */
    pre1 = (Vunit_ec)/(4*VPI*Vunit_eps0*eps_w*Vunit_kb*T);

    /* Finally, if we convert keep xkappa in A^{-1} and scale pre1 by
     * m/A, then we will only need to deal with distances and sizes in
     * Angstroms rather than meters.                                       */
    xkappa = Vpbe_getXkappa(pbe);              /* A^{-1}        */
    pre1 = pre1*(1.0e10);
    size = Vpbe_getSoluteRadius(pbe);
    apos = Vpbe_getSoluteCenter(pbe);
    charge = Vunit_ec*Vpbe_getSoluteCharge(pbe);

    /* Check for rounding error */
    if (VABS(xminOLD-xminNEW) < VSMALL) xminNEW = xminOLD;
    if (VABS(xmaxOLD-xmaxNEW) < VSMALL) xmaxNEW = xmaxOLD;
    if (VABS(yminOLD-yminNEW) < VSMALL) yminNEW = yminOLD;
    if (VABS(ymaxOLD-ymaxNEW) < VSMALL) ymaxNEW = ymaxOLD;
    if (VABS(zminOLD-zminNEW) < VSMALL) zminNEW = zminOLD;
    if (VABS(zmaxOLD-zmaxNEW) < VSMALL) zmaxNEW = zmaxOLD;
    

    /* Sanity check: make sure we're within the old mesh */
    Vnm_print(0, "VPMG::focusFillBound -- New mesh mins = %g, %g, %g\n",
      xminNEW, yminNEW, zminNEW);
    Vnm_print(0, "VPMG::focusFillBound -- New mesh maxs = %g, %g, %g\n",
      xmaxNEW, ymaxNEW, zmaxNEW);
    Vnm_print(0, "VPMG::focusFillBound -- Old mesh mins = %g, %g, %g\n",
      xminOLD, yminOLD, zminOLD);
    Vnm_print(0, "VPMG::focusFillBound -- Old mesh maxs = %g, %g, %g\n",
      xmaxOLD, ymaxOLD, zmaxOLD);

    /* The following is obsolete; we'll substitute analytical boundary
     * condition values when the new mesh falls outside the old */
    if ((xmaxNEW>xmaxOLD) || (ymaxNEW>ymaxOLD) || (zmaxNEW>zmaxOLD) ||
        (xminOLD>xminNEW) || (yminOLD>yminNEW) || (zminOLD>zminNEW)) {

        Vnm_print(2, "VPMG::focusFillBound -- new mesh not contained in old!\n");
        fflush(stderr);
        VASSERT(0);
    }

    
    /* Fill the "i" boundaries (dirichlet) */
    for (k=0; k<nzNEW; k++) {
        for (j=0; j<nyNEW; j++) {
            /* Low X face */
            x = xminNEW;
            y = yminNEW + j*hyNEW;
            z = zminNEW + k*hzNEW;
            if ((x >= (xminOLD-VSMALL)) && (y >= (yminOLD-VSMALL)) && (z >= (zminOLD-VSMALL)) &&
                (x <= (xmaxOLD+VSMALL)) && (y <= (ymaxOLD+VSMALL)) && (z <= (zmaxOLD+VSMALL))) {
                ifloat = (x - xminOLD)/hxOLD;
                jfloat = (y - yminOLD)/hyOLD;
                kfloat = (z - zminOLD)/hzOLD;
                ihi = (int)ceil(ifloat);
                if (ihi > (nxOLD-1)) ihi = nxOLD-1;
                ilo = (int)floor(ifloat);
                if (ilo < 0) ilo = 0;
                jhi = (int)ceil(jfloat);
                if (jhi > (nyOLD-1)) jhi = nyOLD-1;
                jlo = (int)floor(jfloat);
                if (jlo < 0) jlo = 0;
                khi = (int)ceil(kfloat);
                if (khi > (nzOLD-1)) khi = nzOLD-1;
                klo = (int)floor(kfloat);
                if (klo < 0) klo = 0;
                dx = ifloat - (double)(ilo);
                dy = jfloat - (double)(jlo);
                dz = kfloat - (double)(klo);
                nx = nxOLD; ny = nyOLD; nz = nzOLD;
                uval =  dx*dy*dz*(pmgOLD->u[IJK(ihi,jhi,khi)])
                  + dx*(1.0-dy)*dz*(pmgOLD->u[IJK(ihi,jlo,khi)])
                  + dx*dy*(1.0-dz)*(pmgOLD->u[IJK(ihi,jhi,klo)])
                  + dx*(1.0-dy)*(1.0-dz)*(pmgOLD->u[IJK(ihi,jlo,klo)])
                  + (1.0-dx)*dy*dz*(pmgOLD->u[IJK(ilo,jhi,khi)])
                  + (1.0-dx)*(1.0-dy)*dz*(pmgOLD->u[IJK(ilo,jlo,khi)])
                  + (1.0-dx)*dy*(1.0-dz)*(pmgOLD->u[IJK(ilo,jhi,klo)])
                  + (1.0-dx)*(1.0-dy)*(1.0-dz)*(pmgOLD->u[IJK(ilo,jlo,klo)]);
                nx = nxNEW; ny = nyNEW; nz = nzNEW;
            } else {
                pos[0] = x; pos[1] = y; pos[2] = z;
                Vnm_print(1, "focusFillBound -- DEBUG:  CALLING BCFL1 for %g, \
%g, %g!\n", x, y, z);
                uval = bcfl1sp(size, apos, charge, xkappa, pre1, pos);
            }
            nx = nxNEW; ny = nyNEW; nz = nzNEW;
            thee->gxcf[IJKx(j,k,0)] = uval;

            /* High X face */
            x = xmaxNEW;
            if ((x >= xminOLD) && (y >= yminOLD) && (z >= zminOLD) &&
                (x <= xmaxOLD) && (y <= ymaxOLD) && (z <= zmaxOLD)) {
                ifloat = (x - xminOLD)/hxOLD;
                jfloat = (y - yminOLD)/hyOLD;
                kfloat = (z - zminOLD)/hzOLD;
                ihi = (int)ceil(ifloat);
                if (ihi > (nxOLD-1)) ihi = nxOLD-1;
                ilo = (int)floor(ifloat);
                if (ilo < 0) ilo = 0;
                jhi = (int)ceil(jfloat);
                if (jhi > (nyOLD-1)) jhi = nyOLD-1;
                jlo = (int)floor(jfloat);
                if (jlo < 0) jlo = 0;
                khi = (int)ceil(kfloat);
                if (khi > (nzOLD-1)) khi = nzOLD-1;
                klo = (int)floor(kfloat);
                if (klo < 0) klo = 0;
                dx = ifloat - (double)(ilo);
                dy = jfloat - (double)(jlo);
                dz = kfloat - (double)(klo);
                nx = nxOLD; ny = nyOLD; nz = nzOLD;
                uval =  dx*dy*dz*(pmgOLD->u[IJK(ihi,jhi,khi)])
                  + dx*(1.0-dy)*dz*(pmgOLD->u[IJK(ihi,jlo,khi)])
                  + dx*dy*(1.0-dz)*(pmgOLD->u[IJK(ihi,jhi,klo)])
                  + dx*(1.0-dy)*(1.0-dz)*(pmgOLD->u[IJK(ihi,jlo,klo)])
                  + (1.0-dx)*dy*dz*(pmgOLD->u[IJK(ilo,jhi,khi)])
                  + (1.0-dx)*(1.0-dy)*dz*(pmgOLD->u[IJK(ilo,jlo,khi)])
                  + (1.0-dx)*dy*(1.0-dz)*(pmgOLD->u[IJK(ilo,jhi,klo)])
                  + (1.0-dx)*(1.0-dy)*(1.0-dz)*(pmgOLD->u[IJK(ilo,jlo,klo)]);
                nx = nxNEW; ny = nyNEW; nz = nzNEW;
            } else {
                Vnm_print(1, "focusFillBound -- DEBUG:  CALLING BCFL1 for %g, \
%g, %g!\n", x, y, z);
                pos[0] = x; pos[1] = y; pos[2] = z;
                uval = bcfl1sp(size, apos, charge, xkappa, pre1, pos);
            }
            nx = nxNEW; ny = nyNEW; nz = nzNEW;
            thee->gxcf[IJKx(j,k,1)] = uval;
            
            /* Zero Neumann conditions */             
            nx = nxNEW; ny = nyNEW; nz = nzNEW;
            thee->gxcf[IJKx(j,k,2)] = 0.0;
            nx = nxNEW; ny = nyNEW; nz = nzNEW;
            thee->gxcf[IJKx(j,k,3)] = 0.0;
        }
    }

    /* Fill the "j" boundaries (dirichlet) */
    for (k=0; k<nzNEW; k++) {
        for (i=0; i<nxNEW; i++) {
            /* Low Y face */
            x = xminNEW + i*hxNEW;
            y = yminNEW;
            z = zminNEW + k*hzNEW;
            if ((x >= xminOLD) && (y >= yminOLD) && (z >= zminOLD) &&
                (x <= xmaxOLD) && (y <= ymaxOLD) && (z <= zmaxOLD)) {
                ifloat = (x - xminOLD)/hxOLD;
                jfloat = (y - yminOLD)/hyOLD;
                kfloat = (z - zminOLD)/hzOLD;
                ihi = (int)ceil(ifloat);
                if (ihi > (nxOLD-1)) ihi = nxOLD-1;
                ilo = (int)floor(ifloat);
                if (ilo < 0) ilo = 0;
                jhi = (int)ceil(jfloat);
                if (jhi > (nyOLD-1)) jhi = nyOLD-1;
                jlo = (int)floor(jfloat);
                if (jlo < 0) jlo = 0;
                khi = (int)ceil(kfloat);
                if (khi > (nzOLD-1)) khi = nzOLD-1;
                klo = (int)floor(kfloat);
                if (klo < 0) klo = 0;
                dx = ifloat - (double)(ilo);
                dy = jfloat - (double)(jlo);
                dz = kfloat - (double)(klo);
                nx = nxOLD; ny = nyOLD; nz = nzOLD;
                uval =  dx*dy*dz*(pmgOLD->u[IJK(ihi,jhi,khi)])
                  + dx*(1.0-dy)*dz*(pmgOLD->u[IJK(ihi,jlo,khi)])
                  + dx*dy*(1.0-dz)*(pmgOLD->u[IJK(ihi,jhi,klo)])
                  + dx*(1.0-dy)*(1.0-dz)*(pmgOLD->u[IJK(ihi,jlo,klo)])
                  + (1.0-dx)*dy*dz*(pmgOLD->u[IJK(ilo,jhi,khi)])
                  + (1.0-dx)*(1.0-dy)*dz*(pmgOLD->u[IJK(ilo,jlo,khi)])
                  + (1.0-dx)*dy*(1.0-dz)*(pmgOLD->u[IJK(ilo,jhi,klo)])
                  + (1.0-dx)*(1.0-dy)*(1.0-dz)*(pmgOLD->u[IJK(ilo,jlo,klo)]);
                nx = nxNEW; ny = nyNEW; nz = nzNEW;
            } else {
                Vnm_print(1, "focusFillBound -- DEBUG:  CALLING BCFL1 for %g, \
%g, %g!\n", x, y, z);
                pos[0] = x; pos[1] = y; pos[2] = z;
                uval = bcfl1sp(size, apos, charge, xkappa, pre1, pos);
            }
            nx = nxNEW; ny = nyNEW; nz = nzNEW;
            thee->gycf[IJKy(i,k,0)] = uval;

            /* High Y face */
            y = ymaxNEW;
            if ((x >= xminOLD) && (y >= yminOLD) && (z >= zminOLD) &&
                (x <= xmaxOLD) && (y <= ymaxOLD) && (z <= zmaxOLD)) {
                ifloat = (x - xminOLD)/hxOLD;
                jfloat = (y - yminOLD)/hyOLD;
                kfloat = (z - zminOLD)/hzOLD;
                ihi = (int)ceil(ifloat);
                if (ihi > (nxOLD-1)) ihi = nxOLD-1;
                ilo = (int)floor(ifloat);
                if (ilo < 0) ilo = 0;
                jhi = (int)ceil(jfloat);
                if (jhi > (nyOLD-1)) jhi = nyOLD-1;
                jlo = (int)floor(jfloat);
                if (jlo < 0) jlo = 0;
                khi = (int)ceil(kfloat);
                if (khi > (nzOLD-1)) khi = nzOLD-1;
                klo = (int)floor(kfloat);
                if (klo < 0) klo = 0;
                dx = ifloat - (double)(ilo);
                dy = jfloat - (double)(jlo);
                dz = kfloat - (double)(klo);
                nx = nxOLD; ny = nyOLD; nz = nzOLD;
                uval =  dx*dy*dz*(pmgOLD->u[IJK(ihi,jhi,khi)])
                  + dx*(1.0-dy)*dz*(pmgOLD->u[IJK(ihi,jlo,khi)])
                  + dx*dy*(1.0-dz)*(pmgOLD->u[IJK(ihi,jhi,klo)])
                  + dx*(1.0-dy)*(1.0-dz)*(pmgOLD->u[IJK(ihi,jlo,klo)])
                  + (1.0-dx)*dy*dz*(pmgOLD->u[IJK(ilo,jhi,khi)])
                  + (1.0-dx)*(1.0-dy)*dz*(pmgOLD->u[IJK(ilo,jlo,khi)])
                  + (1.0-dx)*dy*(1.0-dz)*(pmgOLD->u[IJK(ilo,jhi,klo)])
                  + (1.0-dx)*(1.0-dy)*(1.0-dz)*(pmgOLD->u[IJK(ilo,jlo,klo)]);
                nx = nxNEW; ny = nyNEW; nz = nzNEW;
            } else {
                pos[0] = x; pos[1] = y; pos[2] = z;
                Vnm_print(1, "focusFillBound -- DEBUG:  CALLING BCFL1 for %g, \
%g, %g!\n", x, y, z);
                uval = bcfl1sp(size, apos, charge, xkappa, pre1, pos);
            }
            nx = nxNEW; ny = nyNEW; nz = nzNEW;
            thee->gycf[IJKy(i,k,1)] = uval;

            /* Zero Neumann conditions */
            nx = nxNEW; ny = nyNEW; nz = nzNEW;
            thee->gycf[IJKy(i,k,2)] = 0.0;
            nx = nxNEW; ny = nyNEW; nz = nzNEW;
            thee->gycf[IJKy(i,k,3)] = 0.0;
        }
    }

    /* Fill the "k" boundaries (dirichlet) */
    for (j=0; j<nyNEW; j++) {
        for (i=0; i<nxNEW; i++) {
            /* Low Z face */
            x = xminNEW + i*hxNEW;
            y = yminNEW + j*hyNEW;
            z = zminNEW;
            if ((x >= xminOLD) && (y >= yminOLD) && (z >= zminOLD) &&
                (x <= xmaxOLD) && (y <= ymaxOLD) && (z <= zmaxOLD)) {
                ifloat = (x - xminOLD)/hxOLD;
                jfloat = (y - yminOLD)/hyOLD;
                kfloat = (z - zminOLD)/hzOLD;
                ihi = (int)ceil(ifloat);
                if (ihi > (nxOLD-1)) ihi = nxOLD-1;
                ilo = (int)floor(ifloat);
                if (ilo < 0) ilo = 0;
                jhi = (int)ceil(jfloat);
                if (jhi > (nyOLD-1)) jhi = nyOLD-1;
                jlo = (int)floor(jfloat);
                if (jlo < 0) jlo = 0;
                khi = (int)ceil(kfloat);
                if (khi > (nzOLD-1)) khi = nzOLD-1;
                klo = (int)floor(kfloat);
                if (klo < 0) klo = 0;
                dx = ifloat - (double)(ilo);
                dy = jfloat - (double)(jlo);
                dz = kfloat - (double)(klo);
                nx = nxOLD; ny = nyOLD; nz = nzOLD;
                uval =  dx*dy*dz*(pmgOLD->u[IJK(ihi,jhi,khi)])
                  + dx*(1.0-dy)*dz*(pmgOLD->u[IJK(ihi,jlo,khi)])
                  + dx*dy*(1.0-dz)*(pmgOLD->u[IJK(ihi,jhi,klo)])
                  + dx*(1.0-dy)*(1.0-dz)*(pmgOLD->u[IJK(ihi,jlo,klo)])
                  + (1.0-dx)*dy*dz*(pmgOLD->u[IJK(ilo,jhi,khi)])
                  + (1.0-dx)*(1.0-dy)*dz*(pmgOLD->u[IJK(ilo,jlo,khi)])
                  + (1.0-dx)*dy*(1.0-dz)*(pmgOLD->u[IJK(ilo,jhi,klo)])
                  + (1.0-dx)*(1.0-dy)*(1.0-dz)*(pmgOLD->u[IJK(ilo,jlo,klo)]);
                nx = nxNEW; ny = nyNEW; nz = nzNEW;
            } else {
                pos[0] = x; pos[1] = y; pos[2] = z;
                Vnm_print(1, "focusFillBound -- DEBUG:  CALLING BCFL1 for %g, \
%g, %g!\n", x, y, z);
                uval = bcfl1sp(size, apos, charge, xkappa, pre1, pos);
            }
            nx = nxNEW; ny = nyNEW; nz = nzNEW;
            thee->gzcf[IJKz(i,j,0)] = uval;

            /* High Z face */
            z = zmaxNEW;
            if ((x >= xminOLD) && (y >= yminOLD) && (z >= zminOLD) &&
                (x <= xmaxOLD) && (y <= ymaxOLD) && (z <= zmaxOLD)) {
                ifloat = (x - xminOLD)/hxOLD;
                jfloat = (y - yminOLD)/hyOLD;
                kfloat = (z - zminOLD)/hzOLD;
                ihi = (int)ceil(ifloat);
                if (ihi > (nxOLD-1)) ihi = nxOLD-1;
                ilo = (int)floor(ifloat);
                if (ilo < 0) ilo = 0;
                jhi = (int)ceil(jfloat);
                if (jhi > (nyOLD-1)) jhi = nyOLD-1;
                jlo = (int)floor(jfloat);
                if (jlo < 0) jlo = 0;
                khi = (int)ceil(kfloat);
                if (khi > (nzOLD-1)) khi = nzOLD-1;
                klo = (int)floor(kfloat);
                if (klo < 0) klo = 0;
                dx = ifloat - (double)(ilo);
                dy = jfloat - (double)(jlo);
                dz = kfloat - (double)(klo);
                nx = nxOLD; ny = nyOLD; nz = nzOLD;
                uval =  dx*dy*dz*(pmgOLD->u[IJK(ihi,jhi,khi)])
                  + dx*(1.0-dy)*dz*(pmgOLD->u[IJK(ihi,jlo,khi)])
                  + dx*dy*(1.0-dz)*(pmgOLD->u[IJK(ihi,jhi,klo)])
                  + dx*(1.0-dy)*(1.0-dz)*(pmgOLD->u[IJK(ihi,jlo,klo)])
                  + (1.0-dx)*dy*dz*(pmgOLD->u[IJK(ilo,jhi,khi)])
                  + (1.0-dx)*(1.0-dy)*dz*(pmgOLD->u[IJK(ilo,jlo,khi)])
                  + (1.0-dx)*dy*(1.0-dz)*(pmgOLD->u[IJK(ilo,jhi,klo)])
                  + (1.0-dx)*(1.0-dy)*(1.0-dz)*(pmgOLD->u[IJK(ilo,jlo,klo)]);
                nx = nxNEW; ny = nyNEW; nz = nzNEW;
            } else {
                pos[0] = x; pos[1] = y; pos[2] = z;
                Vnm_print(1, "focusFillBound -- DEBUG:  CALLING BCFL1 for %g, \
%g, %g!\n", x, y, z);
                uval = bcfl1sp(size, apos, charge, xkappa, pre1, pos);
            }
            nx = nxNEW; ny = nyNEW; nz = nzNEW;
            thee->gzcf[IJKz(i,j,1)] = uval;

            /* Zero Neumann conditions */
            nx = nxNEW; ny = nyNEW; nz = nzNEW;
            thee->gzcf[IJKz(i,j,2)] = 0.0;
            nx = nxNEW; ny = nyNEW; nz = nzNEW;
            thee->gzcf[IJKz(i,j,3)] = 0.0;
        }
    }
}

VPRIVATE void extEnergy(Vpmg *thee, Vpmg *pmgOLD, PBEparm_calcEnergy extFlag, 
        double partMin[3], double partMax[3], int bflags[6]) {

    Vatom *atom;
    double hxNEW, hyNEW, hzNEW;
    double lowerCorner[3], upperCorner[3];
    int nxNEW, nyNEW, nzNEW;
    int nxOLD, nyOLD, nzOLD;
    int i,j,k;
    double xmin, xmax, ymin, ymax, zmin, zmax;
    double hxOLD, hyOLD, hzOLD;
    double xval, yval, zval;
    double x,y,z;
    int nx, ny, nz;
    
    /* Set the new external energy contribution to zero.  Any external
     * contributions from higher levels will be included in the appropriate
     * energy function call. */
    thee->extQmEnergy = 0;
    thee->extQfEnergy = 0;
    thee->extDiEnergy = 0;
    thee->extNpEnergy = 0;

    /* New problem dimensions */
    hxNEW = thee->pmgp->hx;
    hyNEW = thee->pmgp->hy;
    hzNEW = thee->pmgp->hzed;
    nxNEW = thee->pmgp->nx;
    nyNEW = thee->pmgp->ny;
    nzNEW = thee->pmgp->nz;
    lowerCorner[0] = thee->pmgp->xcent - ((double)(nxNEW-1)*hxNEW)/2.0;
    upperCorner[0] = thee->pmgp->xcent + ((double)(nxNEW-1)*hxNEW)/2.0;
    lowerCorner[1] = thee->pmgp->ycent - ((double)(nyNEW-1)*hyNEW)/2.0;
    upperCorner[1] = thee->pmgp->ycent + ((double)(nyNEW-1)*hyNEW)/2.0;
    lowerCorner[2] = thee->pmgp->zcent - ((double)(nzNEW-1)*hzNEW)/2.0;
    upperCorner[2] = thee->pmgp->zcent + ((double)(nzNEW-1)*hzNEW)/2.0;

    Vnm_print(0, "VPMG::extEnergy:  energy flag = %d\n", extFlag);

    /* Old problem dimensions */
    nxOLD = pmgOLD->pmgp->nx;
    nyOLD = pmgOLD->pmgp->ny;
    nzOLD = pmgOLD->pmgp->nz;

    /* Create a partition based on the new problem dimensions */
   
    Vpmg_setPart(pmgOLD, lowerCorner, upperCorner, bflags);

    
    Vnm_print(0,"VPMG::extEnergy:   Finding extEnergy dimensions...\n");
    Vnm_print(0,"VPMG::extEnergy    Disj part lower corner = (%g, %g, %g)\n",
               partMin[0], partMin[1], partMin[2]);
    Vnm_print(0,"VPMG::extEnergy    Disj part upper corner = (%g, %g, %g)\n",
               partMax[0], partMax[1], partMax[2]);
    
    /* Find the old dimensions */

    hxOLD = pmgOLD->pmgp->hx;
    hyOLD = pmgOLD->pmgp->hy;
    hzOLD = pmgOLD->pmgp->hzed;
    xmin =  pmgOLD->pmgp->xcent - 0.5*hxOLD*(nxOLD-1);
    ymin =  pmgOLD->pmgp->ycent - 0.5*hyOLD*(nyOLD-1);
    zmin =  pmgOLD->pmgp->zcent - 0.5*hzOLD*(nzOLD-1);
    xmax =  xmin+hxOLD*(nxOLD-1);
    ymax =  ymin+hyOLD*(nyOLD-1);
    zmax =  zmin+hzOLD*(nzOLD-1);
    
    Vnm_print(0,"VPMG::extEnergy    Old lower corner = (%g, %g, %g)\n",
               xmin, ymin, zmin);
    Vnm_print(0,"VPMG::extEnergy    Old upper corner = (%g, %g, %g)\n",
               xmax, ymax, zmax);

    /* Flip the partition, but do not include any points that will
       be included by another processor */
    
    nx = nxOLD;
    ny = nyOLD;
    nz = nzOLD;
    
    for(i=0; i<nx; i++) {
        xval = 1;
        x = i*hxOLD + xmin;
        if (x < partMin[0] && bflags[VAPBS_LEFT] == 1) xval = 0;
        else if (x > partMax[0] && bflags[VAPBS_RIGHT] == 1) xval = 0;
   
        for(j=0; j<ny; j++) {
            yval = 1;
            y = j*hyOLD + ymin;  
            if (y < partMin[1] && bflags[VAPBS_BACK] == 1) yval = 0;
            else if (y > partMax[1] && bflags[VAPBS_FRONT] == 1) yval = 0;

            for(k=0; k<nz; k++) {
                zval = 1;
                z = k*hzOLD + zmin;
                if (z < partMin[2] && bflags[VAPBS_DOWN] == 1) zval = 0;
                else if (z > partMax[2] && bflags[VAPBS_UP] == 1) zval = 0;
               
                if (pmgOLD->pvec[IJK(i,j,k)] > VSMALL) pmgOLD->pvec[IJK(i,j,k)] = 1.0;
                pmgOLD->pvec[IJK(i,j,k)] = (1 - (pmgOLD->pvec[IJK(i,j,k)])) * (xval*yval*zval);
            }
        }
    }

    for (i=0; i<Valist_getNumberAtoms(thee->pbe->alist); i++) {
        xval=1;
        yval=1;
        zval=1;
        atom = Valist_getAtom(thee->pbe->alist, i);
        x = atom->position[0];
        y = atom->position[1];
        z = atom->position[2];
        if (x < partMin[0] && bflags[VAPBS_LEFT] == 1) xval = 0;
        else if (x > partMax[0] && bflags[VAPBS_RIGHT] == 1) xval = 0;
        if (y < partMin[1] && bflags[VAPBS_BACK] == 1) yval = 0;
        else if (y > partMax[1] && bflags[VAPBS_FRONT] == 1) yval = 0;
        if (z < partMin[2] && bflags[VAPBS_DOWN] == 1) zval = 0;
        else if (z > partMax[2] && bflags[VAPBS_UP] == 1) zval = 0;
        if (atom->partID > VSMALL) atom->partID = 1.0;
        atom->partID = (1 - atom->partID) * (xval*yval*zval);
    }

    /* Now calculate the energy on inverted subset of the domain */
    thee->extQmEnergy = Vpmg_qmEnergy(pmgOLD, 1);
    Vnm_print(0, "VPMG::extEnergy: extQmEnergy = %g kT\n", thee->extQmEnergy);
    thee->extQfEnergy = Vpmg_qfEnergy(pmgOLD, 1);
    Vnm_print(0, "VPMG::extEnergy: extQfEnergy = %g kT\n", thee->extQfEnergy);
    thee->extNpEnergy = Vpmg_npEnergy(pmgOLD, 1);
    Vnm_print(0, "VPMG::extEnergy: extNpEnergy = %g kT\n", thee->extNpEnergy);
    thee->extDiEnergy = Vpmg_dielEnergy(pmgOLD, 1);
    Vnm_print(0, "VPMG::extEnergy: extDiEnergy = %g kT\n", thee->extDiEnergy);
    Vpmg_unsetPart(pmgOLD);
}

VPRIVATE double bcfl1sp(double size, double *apos, double charge, 
  double xkappa, double pre1, double *pos) {

    double dist, val;

    dist = VSQRT(VSQR(pos[0]-apos[0]) + VSQR(pos[1]-apos[1])
      + VSQR(pos[2]-apos[2]));
    if (xkappa > VSMALL) {
        val = pre1*(charge/dist)*VEXP(-xkappa*(dist-size))
          / (1+xkappa*size);
    } else {
        val = pre1*(charge/dist);
    } 

    return val;
}

VPRIVATE void bcfl1(double size, double *apos, double charge, 
  double xkappa, double pre1, double *gxcf, double *gycf, double *gzcf,
  double *xf, double *yf, double *zf, int nx, int ny, int nz) {

    int i, j, k;
    double dist, val;
    double gpos[3];

    /* the "i" boundaries (dirichlet) */
    for (k=0; k<nz; k++) {
        gpos[2] = zf[k];
        for (j=0; j<ny; j++) {
            gpos[1] = yf[j];
            gpos[0] = xf[0];
            dist = VSQRT(VSQR(gpos[0]-apos[0]) + VSQR(gpos[1]-apos[1])
              + VSQR(gpos[2]-apos[2]));
            if (xkappa > VSMALL) {
                val = pre1*(charge/dist)*VEXP(-xkappa*(dist-size))
                       / (1+xkappa*size);
            } else {
                val = pre1*(charge/dist);
            } 
            gxcf[IJKx(j,k,0)] += val;
            gpos[0] = xf[nx-1];
            dist = VSQRT(VSQR(gpos[0]-apos[0]) + VSQR(gpos[1]-apos[1])
              + VSQR(gpos[2]-apos[2]));
            if (xkappa > VSMALL) {
                val = pre1*(charge/dist)*VEXP(-xkappa*(dist-size))
                       / (1+xkappa*size);
            } else {
                val = pre1*(charge/dist);
            }
            gxcf[IJKx(j,k,1)] += val;
        }
    }

    /* the "j" boundaries (dirichlet) */
    for (k=0; k<nz; k++) {
        gpos[2] = zf[k];
        for (i=0; i<nx; i++) {
            gpos[0] = xf[i];
            gpos[1] = yf[0];
            dist = VSQRT(VSQR(gpos[0]-apos[0]) + VSQR(gpos[1]-apos[1])
              + VSQR(gpos[2]-apos[2]));
            if (xkappa > VSMALL) {
                val = pre1*(charge/dist)*VEXP(-xkappa*(dist-size))
                       / (1+xkappa*size);
            } else {
                val = pre1*(charge/dist);
            }
            gycf[IJKy(i,k,0)] += val;
            gpos[1] = yf[ny-1];
            dist = VSQRT(VSQR(gpos[0]-apos[0]) + VSQR(gpos[1]-apos[1])
              + VSQR(gpos[2]-apos[2]));
            if (xkappa > VSMALL) {
                val = pre1*(charge/dist)*VEXP(-xkappa*(dist-size))
                       / (1+xkappa*size);
            } else {
                val = pre1*(charge/dist);
            }
            gycf[IJKy(i,k,1)] += val;
        }
    }

    /* the "k" boundaries (dirichlet) */
    for (j=0; j<ny; j++) {
        gpos[1] = yf[j];
        for (i=0; i<nx; i++) {
            gpos[0] = xf[i];
            gpos[2] = zf[0];
            dist = VSQRT(VSQR(gpos[0]-apos[0]) + VSQR(gpos[1]-apos[1])
              + VSQR(gpos[2]-apos[2]));
            if (xkappa > VSMALL) {
                val = pre1*(charge/dist)*VEXP(-xkappa*(dist-size))
                       / (1+xkappa*size);
            } else {
                val = pre1*(charge/dist);
            }
            gzcf[IJKz(i,j,0)] += val;
            gpos[2] = zf[nz-1];
            dist = VSQRT(VSQR(gpos[0]-apos[0]) + VSQR(gpos[1]-apos[1])
              + VSQR(gpos[2]-apos[2]));
            if (xkappa > VSMALL) {
                val = pre1*(charge/dist)*VEXP(-xkappa*(dist-size))
                       / (1+xkappa*size);
            } else {
                val = pre1*(charge/dist);
            }
            gzcf[IJKz(i,j,1)] += val;
        }
    }
}


VPRIVATE void bcCalc(Vpmg *thee) {

    int nx, ny, nz;
    double size, *position, charge, xkappa, eps_w, T, pre1;
    int i, j, k, iatom;
    Vpbe *pbe;
    Vatom *atom;
    Valist *alist;
    
    pbe = thee->pbe;
    alist = thee->pbe->alist;
    nx = thee->pmgp->nx;
    ny = thee->pmgp->ny;
    nz = thee->pmgp->nz;

    /* Zero out the boundaries */
    /* the "i" boundaries (dirichlet) */
    for (k=0; k<nz; k++) {
        for (j=0; j<ny; j++) {
            thee->gxcf[IJKx(j,k,0)] = 0.0;
            thee->gxcf[IJKx(j,k,1)] = 0.0;
            thee->gxcf[IJKx(j,k,2)] = 0.0;
            thee->gxcf[IJKx(j,k,3)] = 0.0;
        }
    }

    /* the "j" boundaries (dirichlet) */
    for (k=0; k<nz; k++) {
        for (i=0; i<nx; i++) {
            thee->gycf[IJKy(i,k,0)] = 0.0;
            thee->gycf[IJKy(i,k,1)] = 0.0;
            thee->gycf[IJKy(i,k,2)] = 0.0;
            thee->gycf[IJKy(i,k,3)] = 0.0;
        }
    }

    /* the "k" boundaries (dirichlet) */
    for (j=0; j<ny; j++) {
        for (i=0; i<nx; i++) {
            thee->gzcf[IJKz(i,j,0)] = 0.0;
            thee->gzcf[IJKz(i,j,1)] = 0.0;
            thee->gzcf[IJKz(i,j,2)] = 0.0;
            thee->gzcf[IJKz(i,j,3)] = 0.0;
        }
    }

    /* For each "atom" (only one for bcfl=1), we use the following formula to
     * calculate the boundary conditions: 
     *    g(x) = \frac{q e_c}{4*\pi*\eps_0*\eps_w*k_b*T}
     *          * \frac{exp(-xkappa*(d - a))}{1+xkappa*a}
     *          * 1/d
     * where d = ||x - x_0|| (in m) and a is the size of the atom (in m).
     * We only need to evaluate some of these prefactors once:
     *    pre1 = \frac{e_c}{4*\pi*\eps_0*\eps_w*k_b*T}
     * which gives the potential as
     *    g(x) = pre1 * q/d * \frac{exp(-xkappa*(d - a))}{1+xkappa*a} 
     */
    eps_w = Vpbe_getSolventDiel(pbe);           /* Dimensionless */
    T = Vpbe_getTemperature(pbe);               /* K             */
    pre1 = (Vunit_ec)/(4*VPI*Vunit_eps0*eps_w*Vunit_kb*T);

    /* Finally, if we convert keep xkappa in A^{-1} and scale pre1 by
     * m/A, then we will only need to deal with distances and sizes in
     * Angstroms rather than meters.                                       */
    xkappa = Vpbe_getXkappa(pbe);              /* A^{-1}        */
    pre1 = pre1*(1.0e10);
   
    switch (thee->pmgp->bcfl) {
        /*  If we have zero boundary conditions, we're done */
        case BCFL_ZERO: 
            return;

        /*  For single DH sphere BC's, we only have one "atom" to deal with;
         *  get its information and */
        case BCFL_SDH:
            size = Vpbe_getSoluteRadius(pbe);
            position = Vpbe_getSoluteCenter(pbe);
            charge = Vunit_ec*Vpbe_getSoluteCharge(pbe);
    
            bcfl1(size, position, charge, xkappa, pre1,
              thee->gxcf, thee->gycf, thee->gzcf, 
              thee->xf, thee->yf, thee->zf, nx, ny, nz);
            break;

        case BCFL_MDH:
            for (iatom=0; iatom<Valist_getNumberAtoms(alist); iatom++) {
                atom = Valist_getAtom(alist, iatom);
                position = Vatom_getPosition(atom);
                charge = Vunit_ec*Vatom_getCharge(atom);
                size = Vatom_getRadius(atom);
                bcfl1(size, position, charge, xkappa, pre1,
                  thee->gxcf, thee->gycf, thee->gzcf, 
                  thee->xf, thee->yf, thee->zf, nx, ny, nz);
            }
            break;

        case BCFL_UNUSED:
            Vnm_print(2, "bcCalc:  Invalid bcfl (%d)!\n", thee->pmgp->bcfl);
            VASSERT(0);

        case BCFL_FOCUS:
            Vnm_print(2, "VPMG::bcCalc -- not appropriate for focusing!\n");
            VASSERT(0);

        default:
            Vnm_print(2, "VPMG::bcCalc -- invalid boundary condition \
flag (%d)!\n", thee->pmgp->bcfl);
            VASSERT(0);
    }
}

VPRIVATE void fillcoCoefMap(Vpmg *thee) {

    Vpbe *pbe;
    double ionstr, position[3], tkappa, eps, hx, hy, hzed;
    int i, j, k, nx, ny, nz;
    double kappamax;
    VASSERT(thee != VNULL);

    /* Get PBE info */
    pbe = thee->pbe;
    ionstr = Vpbe_getBulkIonicStrength(pbe);

    /* Mesh info */
    nx = thee->pmgp->nx;
    ny = thee->pmgp->ny;
    nz = thee->pmgp->nz;
    hx = thee->pmgp->hx;
    hy = thee->pmgp->hy;
    hzed = thee->pmgp->hzed;

    if ((!thee->useDielXMap) || (!thee->useDielYMap) || (!thee->useDielZMap) ||
      ((!thee->useKappaMap) && (ionstr>VPMGSMALL))) {

        Vnm_print(2, "fillcoCoefMap:  You need to use all coefficient maps!\n");
        VASSERT(0);

    }
    
    /* Scale the kappa map to values between 0 and 1 
       Thus get the maximum value in the map - this 
       is theoretically unnecessary, but a good check.*/
    kappamax = -1.00;
    for (k=0; k<nz; k++) {
        for (j=0; j<ny; j++) {
            for (i=0; i<nx; i++) {
                if (ionstr > VPMGSMALL) {
                     position[0] = thee->xf[i];
                     position[1] = thee->yf[j];
                     position[2] = thee->zf[k];
                     if (!Vgrid_value(thee->kappaMap, position, &tkappa)) {
                         Vnm_print(2, "Vpmg_fillco:  Off kappaMap at:\n");
                         Vnm_print(2, "Vpmg_fillco:  (x,y,z) = (%g,%g %g)\n",
                                   position[0], position[1], position[2]);
                         VASSERT(0);
                     }
                     if (tkappa > kappamax) {
                         kappamax = tkappa;
                     }
                     if (tkappa < 0.0){
                       Vnm_print(2, "Vpmg_fillcoCoefMap: Kappa map less than 0\n");
                       Vnm_print(2, "Vpmg_fillcoCoefMap: at (x,y,z) = (%g,%g %g)\n",
                                 position[0], position[1], position[2]);
                       VASSERT(0);
                     }
                }
            }
        }
    }
    
    if (kappamax > 1.0){
      Vnm_print(2, "Vpmg_fillcoCoefMap:  Maximum Kappa value\n");
      Vnm_print(2, "%g is greater than 1 - will scale appropriately!\n",
                kappamax);
    }
    else {
      kappamax = 1.0;
    }

    for (k=0; k<nz; k++) {
        for (j=0; j<ny; j++) {
            for (i=0; i<nx; i++) {

                if (ionstr > VPMGSMALL) {
                     position[0] = thee->xf[i];
                     position[1] = thee->yf[j];
                     position[2] = thee->zf[k];
                     if (!Vgrid_value(thee->kappaMap, position, &tkappa)) {
                         Vnm_print(2, "Vpmg_fillco:  Off kappaMap at:\n");
                         Vnm_print(2, "Vpmg_fillco:  (x,y,z) = (%g,%g %g)\n",
                           position[0], position[1], position[2]);
                         VASSERT(0);
                     }
                     if (tkappa < VPMGSMALL) tkappa = 0.0;
                     thee->ccf[IJK(i,j,k)] = (tkappa / kappamax);
                }

                position[0] = thee->xf[i] + 0.5*hx;
                position[1] = thee->yf[j];
                position[2] = thee->zf[k];
                if (!Vgrid_value(thee->dielXMap, position, &eps)) {
                    Vnm_print(2, "Vpmg_fillco:  Off dielXMap at:\n");
                    Vnm_print(2, "Vpmg_fillco:  (x,y,z) = (%g,%g %g)\n",
                      position[0], position[1], position[2]);
                    VASSERT(0);
                 }
                 thee->a1cf[IJK(i,j,k)] = eps;
        
                 position[0] = thee->xf[i];
                 position[1] = thee->yf[j] + 0.5*hy;
                 position[2] = thee->zf[k];
                 if (!Vgrid_value(thee->dielYMap, position, &eps)) {
                    Vnm_print(2, "Vpmg_fillco:  Off dielYMap at:\n");
                    Vnm_print(2, "Vpmg_fillco:  (x,y,z) = (%g,%g %g)\n",
                      position[0], position[1], position[2]);
                    VASSERT(0);
                 }
                 thee->a2cf[IJK(i,j,k)] = eps;
            
                 position[0] = thee->xf[i];
                 position[1] = thee->yf[j];
                 position[2] = thee->zf[k] + 0.5*hzed;
                 if (!Vgrid_value(thee->dielZMap, position, &eps)) {
                    Vnm_print(2, "Vpmg_fillco:  Off dielZMap at:\n");
                    Vnm_print(2, "Vpmg_fillco:  (x,y,z) = (%g,%g %g)\n",
                      position[0], position[1], position[2]);
                    VASSERT(0);
                 }
                 thee->a3cf[IJK(i,j,k)] = eps;
            }
        }
    }
}

VPRIVATE void fillcoCoefMol(Vpmg *thee) {

    if (thee->useDielXMap || thee->useDielYMap || thee->useDielZMap ||
      thee->useKappaMap)  {

        fillcoCoefMap(thee);

    } else { 

        fillcoCoefMolDiel(thee); 
        fillcoCoefMolIon(thee);

    }

}

VPRIVATE void fillcoCoefMolIon(Vpmg *thee) {

    Vacc *acc;
    Valist *alist;
    Vpbe *pbe;
    Vatom *atom;
    double xmin, xmax, ymin, ymax, zmin, zmax, ionmask, ionstr;
    double xlen, ylen, zlen, irad;
    double hx, hy, hzed, *apos, arad;
    int i, nx, ny, nz, iatom;
    Vsurf_Meth surfMeth;

    VASSERT(thee != VNULL);
    surfMeth = thee->surfMeth;

    /* Get PBE info */
    pbe = thee->pbe;
    acc = pbe->acc;
    alist = pbe->alist;
    irad = Vpbe_getMaxIonRadius(pbe);
    ionstr = Vpbe_getBulkIonicStrength(pbe);

    /* Mesh info */
    nx = thee->pmgp->nx;
    ny = thee->pmgp->ny;
    nz = thee->pmgp->nz;
    hx = thee->pmgp->hx;
    hy = thee->pmgp->hy;
    hzed = thee->pmgp->hzed;

    /* Define the total domain size */
    xlen = thee->pmgp->xlen;
    ylen = thee->pmgp->ylen;
    zlen = thee->pmgp->zlen;

    /* Define the min/max dimensions */
    xmin = thee->pmgp->xcent - (xlen/2.0);
    ymin = thee->pmgp->ycent - (ylen/2.0);
    zmin = thee->pmgp->zcent - (zlen/2.0);
    xmax = thee->pmgp->xcent + (xlen/2.0);
    ymax = thee->pmgp->ycent + (ylen/2.0);
    zmax = thee->pmgp->zcent + (zlen/2.0);

    /* This is a floating point parameter related to the non-zero nature of the
     * bulk ionic strength.  If the ionic strength is greater than zero; this
     * parameter is set to 1.0 and later scaled by the appropriate pre-factors.
     * Otherwise, this parameter is set to 0.0 */
    if (ionstr > VPMGSMALL) ionmask = 1.0;
    else ionmask = 0.0;

    /* Reset the ccf array, marking everything accessible */
    for (i=0; i<(nx*ny*nz); i++) thee->ccf[i] = ionmask;

    /* Loop through the atoms and set ccf = 0.0 (inaccessible) if a point
     * is inside the ion-inflated van der Waals radii */
    for (iatom=0; iatom<Valist_getNumberAtoms(alist); iatom++) {

        atom = Valist_getAtom(alist, iatom);
        apos = Vatom_getPosition(atom);
        arad = Vatom_getRadius(atom);

        /* Make sure we're on the grid */
        if ((apos[0]<=xmin) || (apos[0]>=xmax)  || \
            (apos[1]<=ymin) || (apos[1]>=ymax)  || \
            (apos[2]<=zmin) || (apos[2]>=zmax)) {
            if (thee->pmgp->bcfl != BCFL_FOCUS) {
                Vnm_print(2, 
"Vpmg_fillco:  Atom #%d at (%4.3f, %4.3f, %4.3f) is off the mesh (ignoring):\n",
                  iatom, apos[0], apos[1], apos[2]);
                Vnm_print(2, "Vpmg_fillco:  xmin = %g, xmax = %g\n", 
                  xmin, xmax);
                Vnm_print(2, "Vpmg_fillco:  ymin = %g, ymax = %g\n", 
                  ymin, ymax);
                Vnm_print(2, "Vpmg_fillco:  zmin = %g, zmax = %g\n", 
                  zmin, zmax);
            }
            fflush(stderr);

        } else { /* if we're on the mesh */

            /* Mark ions */
            markSphere((irad+arad), apos, 
                    nx, ny, nz,
                    hx, hy, hzed,
                    xmin, ymin, zmin,
                    thee->ccf, 0.0);

        } /* endif (on the mesh) */
    } /* endfor (over all atoms) */

}

VPRIVATE void fillcoCoefMolDiel(Vpmg *thee) {

    /* Now figure out what to do with the solvent accessibility */
    switch (thee->surfMeth) {
        case VSM_MOL:
            fillcoCoefMolDielNoSmooth(thee);
            break;
        case VSM_MOLSMOOTH:
            fillcoCoefMolDielSmooth(thee);
            break;
        default:
            Vnm_print(2, "Error in surfMeth!\n");
            VASSERT(0);
    }
}

VPRIVATE void fillcoCoefMolDielNoSmooth(Vpmg *thee) {

    Vacc *acc;
    VaccSurf *asurf;
    Valist *alist;
    Vpbe *pbe;
    Vatom *atom;
    double xmin, xmax, ymin, ymax, zmin, zmax;
    double xlen, ylen, zlen, position[3];
    double srad, epsw, epsp, deps;
    double hx, hy, hzed, *apos, arad;
    int i, nx, ny, nz, iatom, ipt;

    /* Get PBE info */
    pbe = thee->pbe;
    acc = pbe->acc;
    alist = pbe->alist;
    srad = Vpbe_getSolventRadius(pbe);
    epsw = Vpbe_getSolventDiel(pbe);
    epsp = Vpbe_getSoluteDiel(pbe);

    /* Mesh info */
    nx = thee->pmgp->nx;
    ny = thee->pmgp->ny;
    nz = thee->pmgp->nz;
    hx = thee->pmgp->hx;
    hy = thee->pmgp->hy;
    hzed = thee->pmgp->hzed;

    /* Define the total domain size */
    xlen = thee->pmgp->xlen;
    ylen = thee->pmgp->ylen;
    zlen = thee->pmgp->zlen;

    /* Define the min/max dimensions */
    xmin = thee->pmgp->xcent - (xlen/2.0);
    ymin = thee->pmgp->ycent - (ylen/2.0);
    zmin = thee->pmgp->zcent - (zlen/2.0);
    xmax = thee->pmgp->xcent + (xlen/2.0);
    ymax = thee->pmgp->ycent + (ylen/2.0);
    zmax = thee->pmgp->zcent + (zlen/2.0);

    /* Reset the ccf, a1cf, a2cf, and a3cf arrays */
    for (i=0; i<(nx*ny*nz); i++) {
        thee->a1cf[i] = 1.0;
        thee->a2cf[i] = 1.0;
        thee->a3cf[i] = 1.0;
    }

    /* Loop through the atoms and set a{123}cf = 0.0 (inaccessible)
     * if a point is inside the solvent-inflated van der Waals radii */
    for (iatom=0; iatom<Valist_getNumberAtoms(alist); iatom++) {

        atom = Valist_getAtom(alist, iatom);
        apos = Vatom_getPosition(atom);
        arad = Vatom_getRadius(atom);

        /* Make sure we're on the grid */
        if ((apos[0]<=xmin) || (apos[0]>=xmax)  || \
            (apos[1]<=ymin) || (apos[1]>=ymax)  || \
            (apos[2]<=zmin) || (apos[2]>=zmax)) {
            if (thee->pmgp->bcfl != BCFL_FOCUS) {
                Vnm_print(2, "Vpmg_fillco:  Atom #%d at (%4.3f, %4.3f,\
 %4.3f) is off the mesh (ignoring):\n",
                  iatom, apos[0], apos[1], apos[2]);
                Vnm_print(2, "Vpmg_fillco:  xmin = %g, xmax = %g\n", 
                  xmin, xmax);
                Vnm_print(2, "Vpmg_fillco:  ymin = %g, ymax = %g\n", 
                  ymin, ymax);
                Vnm_print(2, "Vpmg_fillco:  zmin = %g, zmax = %g\n", 
                  zmin, zmax);
            }
            fflush(stderr);

        } else { /* if we're on the mesh */

            if (arad > VSMALL) {
                /* Mark x-shifted dielectric */
                markSphere((arad+srad), apos, 
                        nx, ny, nz,
                        hx, hy, hzed,
                        (xmin+0.5*hx), ymin, zmin,
                        thee->a1cf, 0.0);

                /* Mark y-shifted dielectric */
                markSphere((arad+srad), apos, 
                        nx, ny, nz,
                        hx, hy, hzed,
                        xmin, (ymin+0.5*hy), zmin,
                        thee->a2cf, 0.0);

                /* Mark z-shifted dielectric */
                markSphere((arad+srad), apos, 
                        nx, ny, nz,
                        hx, hy, hzed,
                        xmin, ymin, (zmin+0.5*hzed),
                        thee->a3cf, 0.0);
            }

        } /* endif (on the mesh) */
    } /* endfor (over all atoms) */

    /* We only need to do the next step for non-zero solvent radii */
    if (srad > VSMALL) {

        /* Now loop over the solvent accessible surface points */
        for (iatom=0; iatom<Valist_getNumberAtoms(alist); iatom++) {
            atom = Valist_getAtom(alist, iatom);
            asurf = Vacc_atomSASPoints(acc, srad, atom);
    
            /* Use each point on the SAS to reset the solvent accessibility */
            for (ipt=0; ipt<(asurf->npts); ipt++) {
    
                position[0] = asurf->xpts[ipt];
                position[1] = asurf->ypts[ipt];
                position[2] = asurf->zpts[ipt];
    
                /* Mark x-shifted dielectric */
                markSphere(srad, position, 
                        nx, ny, nz,
                        hx, hy, hzed,
                        (xmin+0.5*hx), ymin, zmin,
                        thee->a1cf, 1.0);
    
                /* Mark y-shifted dielectric */
                markSphere(srad, position, 
                        nx, ny, nz,
                        hx, hy, hzed,
                        xmin, (ymin+0.5*hy), zmin,
                        thee->a2cf, 1.0);
    
                /* Mark z-shifted dielectric */
                markSphere(srad, position, 
                        nx, ny, nz,
                        hx, hy, hzed,
                        xmin, ymin, (zmin+0.5*hzed),
                        thee->a3cf, 1.0);
    
            }

        }
    }

    /* Set the dielectric values from the accessibility */
    deps = epsw - epsp;
    for (i=0; i<(nx*ny*nz); i++) {
        thee->a1cf[i] = deps*(thee->a1cf[i]) + epsp;
        thee->a2cf[i] = deps*(thee->a2cf[i]) + epsp;
        thee->a3cf[i] = deps*(thee->a3cf[i]) + epsp;
    }

}

VPRIVATE void fillcoCoefMolDielSmooth(Vpmg *thee) {

    Vacc *acc;
    VaccSurf *asurf;
    Valist *alist;
    Vpbe *pbe;
    Vatom *atom;
    double xmin, xmax, ymin, ymax, zmin, zmax;
    double xlen, ylen, zlen, position[3], frac;
    double srad, epsw, epsp, ap, am, a;
    double hx, hy, hzed, *apos, arad;
    int i, j, k, nx, ny, nz, iatom, ipt;

    /* Get PBE info */
    pbe = thee->pbe;
    acc = pbe->acc;
    alist = pbe->alist;
    srad = Vpbe_getSolventRadius(pbe);
    epsw = Vpbe_getSolventDiel(pbe);
    epsp = Vpbe_getSoluteDiel(pbe);

    /* Mesh info */
    nx = thee->pmgp->nx;
    ny = thee->pmgp->ny;
    nz = thee->pmgp->nz;
    hx = thee->pmgp->hx;
    hy = thee->pmgp->hy;
    hzed = thee->pmgp->hzed;

    /* Define the total domain size */
    xlen = thee->pmgp->xlen;
    ylen = thee->pmgp->ylen;
    zlen = thee->pmgp->zlen;

    /* Define the min/max dimensions */
    xmin = thee->pmgp->xcent - (xlen/2.0);
    ymin = thee->pmgp->ycent - (ylen/2.0);
    zmin = thee->pmgp->zcent - (zlen/2.0);
    xmax = thee->pmgp->xcent + (xlen/2.0);
    ymax = thee->pmgp->ycent + (ylen/2.0);
    zmax = thee->pmgp->zcent + (zlen/2.0);

    /* Reset the ccf, a1cf, a2cf, and a3cf arrays */
    for (i=0; i<(nx*ny*nz); i++) {
        thee->ccf[i]  = 1.0;
        thee->a1cf[i] = 1.0;
        thee->a2cf[i] = 1.0;
        thee->a3cf[i] = 1.0;
    }

    /* Loop through the atoms and set ccf to 0.0 (inaccessible) if a point is
     * inside the solvent-inflated van der Waals radii */
    for (iatom=0; iatom<Valist_getNumberAtoms(alist); iatom++) {

        atom = Valist_getAtom(alist, iatom);
        apos = Vatom_getPosition(atom);
        arad = Vatom_getRadius(atom);

        /* Make sure we're on the grid */
        if ((apos[0]<=xmin) || (apos[0]>=xmax)  || \
            (apos[1]<=ymin) || (apos[1]>=ymax)  || \
            (apos[2]<=zmin) || (apos[2]>=zmax)) {
            if (thee->pmgp->bcfl != BCFL_FOCUS) {
                Vnm_print(2, "Vpmg_fillco:  Atom #%d at (%4.3f, %4.3f,\
 %4.3f) is off the mesh (ignoring):\n",
                  iatom, apos[0], apos[1], apos[2]);
                Vnm_print(2, "Vpmg_fillco:  xmin = %g, xmax = %g\n", 
                  xmin, xmax);
                Vnm_print(2, "Vpmg_fillco:  ymin = %g, ymax = %g\n", 
                  ymin, ymax);
                Vnm_print(2, "Vpmg_fillco:  zmin = %g, zmax = %g\n", 
                  zmin, zmax);
            }
            fflush(stderr);

        } else { /* if we're on the mesh */

            if (arad > VSMALL) {
                /* Mark unshifted dielectric */
                markSphere((arad+srad), apos, 
                        nx, ny, nz,
                        hx, hy, hzed,
                        xmin, ymin, zmin,
                        thee->ccf, 0.0);

                /* Mark shifted dielectrics */
                markSphere((arad+srad), apos, 
                        nx, ny, nz,
                        hx, hy, hzed,
                        xmin+0.5*hx, ymin, zmin,
                        thee->a1cf, 0.0);
                markSphere((arad+srad), apos, 
                        nx, ny, nz,
                        hx, hy, hzed,
                        xmin, ymin+0.5*hy, zmin,
                        thee->a2cf, 0.0);
                markSphere((arad+srad), apos, 
                        nx, ny, nz,
                        hx, hy, hzed,
                        xmin, ymin, zmin+0.5*hzed,
                        thee->a3cf, 0.0);
            }

        } /* endif (on the mesh) */

    } /* endfor (over all atoms) */

    /* We only need to do the next step for non-zero solvent radii */
    if (srad > VSMALL) {
        /* Now use the SAS points to reset the dielectric */
        for (iatom=0; iatom<Valist_getNumberAtoms(alist); iatom++) {
    
            atom = Valist_getAtom(alist, iatom);
            asurf = Vacc_atomSASPoints(acc, srad, atom);
    
            /* Use each point on the SAS to reset the solvent
             * accessibility and set grid length fractions */
            for (ipt=0; ipt<(asurf->npts); ipt++) {
                position[0] = asurf->xpts[ipt];
                position[1] = asurf->ypts[ipt];
                position[2] = asurf->zpts[ipt];
    
                /* Mark unshifted dielectric */
                markSphere(srad, position, 
                        nx, ny, nz,
                        hx, hy, hzed,
                        xmin, ymin, zmin,
                        thee->ccf, 1.0);
    
                /* Mark shifted dielectrics */
                markSphere(srad, position, 
                        nx, ny, nz,
                        hx, hy, hzed,
                        xmin+0.5*hx, ymin, zmin,
                        thee->a1cf, 1.0);
                markSphere(srad, position, 
                        nx, ny, nz,
                        hx, hy, hzed,
                        xmin, ymin+0.5*hy, zmin,
                        thee->a2cf, 1.0);
                markSphere(srad, position, 
                        nx, ny, nz,
                        hx, hy, hzed,
                        xmin, ymin, zmin+0.5*hzed,
                        thee->a3cf, 1.0);
            }
        }
    }

    /* Convert to dielectric values */
    for (i=0; i<nx; i++) {
        for (j=0; j<ny; j++) {
            for (k=0; k<nz; k++) {

                /* X-shifted */
                a = thee->a1cf[IJK(i,j,k)];
                am = thee->ccf[IJK(i,j,k)];
                if (i < (nx-1)) ap = thee->ccf[IJK(i+1,j,k)];
                else ap = a;
                frac = (a + ap + am)/3.0;
                thee->a1cf[IJK(i,j,k)] = molSmoothHarm(epsw, epsp, frac);

                /* Y-shifted */
                a = thee->a2cf[IJK(i,j,k)];
                am = thee->ccf[IJK(i,j,k)];
                if (j < (ny-1)) ap = thee->ccf[IJK(i,j+1,k)];
                else ap = a;
                frac = (a + ap + am)/3.0;
                thee->a2cf[IJK(i,j,k)] = molSmoothHarm(epsw, epsp, frac);

                /* Y-shifted */
                a = thee->a3cf[IJK(i,j,k)];
                am = thee->ccf[IJK(i,j,k)];
                if (k < (nz-1)) ap = thee->ccf[IJK(i,j,k+1)];
                else ap = a;
                frac = (a + ap + am)/3.0;
                thee->a3cf[IJK(i,j,k)] = molSmoothHarm(epsw, epsp, frac);

            }
        }
    }

}

VPRIVATE double molSmoothHarm(double epsw, double epsp, double frac) {

    double num, den, eps;

    num = epsw*epsp;
    den = epsp*frac + epsw*(1.0-frac);
    eps = num/den;

    return eps;
}


VPRIVATE void fillcoCoefSpline(Vpmg *thee) {

    Valist *alist;
    Vpbe *pbe;
    Vatom *atom;
    double xmin, xmax, ymin, ymax, zmin, zmax, ionmask, ionstr, dist2;
    double xlen, ylen, zlen, position[3], itot, stot, ictot, ictot2, sctot;
    double irad, dx, dy, dz, epsw, epsp, w2i;
    double hx, hy, hzed, *apos, arad, sctot2;
    double dx2, dy2, dz2, stot2, itot2, rtot, rtot2, splineWin, w3i;
    double dist, value, sm, sm2;
    int i, j, k, nx, ny, nz, iatom;
    int imin, imax, jmin, jmax, kmin, kmax;

    VASSERT(thee != VNULL);
    splineWin = thee->splineWin;
    w2i = 1.0/(splineWin*splineWin);
    w3i = 1.0/(splineWin*splineWin*splineWin);

    /* Get PBE info */
    pbe = thee->pbe;
    alist = pbe->alist;
    irad = Vpbe_getMaxIonRadius(pbe);
    ionstr = Vpbe_getBulkIonicStrength(pbe);
    epsw = Vpbe_getSolventDiel(pbe);
    epsp = Vpbe_getSoluteDiel(pbe);

    /* Mesh info */
    nx = thee->pmgp->nx;
    ny = thee->pmgp->ny;
    nz = thee->pmgp->nz;
    hx = thee->pmgp->hx;
    hy = thee->pmgp->hy;
    hzed = thee->pmgp->hzed;

    /* Define the total domain size */
    xlen = thee->pmgp->xlen;
    ylen = thee->pmgp->ylen;
    zlen = thee->pmgp->zlen;

    /* Define the min/max dimensions */
    xmin = thee->pmgp->xcent - (xlen/2.0);
    ymin = thee->pmgp->ycent - (ylen/2.0);
    zmin = thee->pmgp->zcent - (zlen/2.0);
    xmax = thee->pmgp->xcent + (xlen/2.0);
    ymax = thee->pmgp->ycent + (ylen/2.0);
    zmax = thee->pmgp->zcent + (zlen/2.0);

    /* This is a floating point parameter related to the non-zero nature of the
     * bulk ionic strength.  If the ionic strength is greater than zero; this
     * parameter is set to 1.0 and later scaled by the appropriate pre-factors.
     * Otherwise, this parameter is set to 0.0 */
    if (ionstr > VPMGSMALL) ionmask = 1.0;
    else ionmask = 0.0;

    /* Reset the fcf, tcf, ccf, a1cf, a2cf, and a3cf arrays */
    for (i=0; i<(nx*ny*nz); i++) {
        thee->ccf[i] = 1.0;
        thee->a1cf[i] = 1.0;
        thee->a2cf[i] = 1.0;
        thee->a3cf[i] = 1.0;
    }

    /* Loop through the atoms and do assign the dielectric */
    for (iatom=0; iatom<Valist_getNumberAtoms(alist); iatom++) {

        atom = Valist_getAtom(alist, iatom);
        apos = Vatom_getPosition(atom);
        arad = Vatom_getRadius(atom);

        /* Make sure we're on the grid */
        if ((apos[0]<=xmin) || (apos[0]>=xmax)  || \
            (apos[1]<=ymin) || (apos[1]>=ymax)  || \
            (apos[2]<=zmin) || (apos[2]>=zmax)) {
            if (thee->pmgp->bcfl != BCFL_FOCUS) {
                Vnm_print(2, "Vpmg_fillco:  Atom #%d at (%4.3f, %4.3f,\
 %4.3f) is off the mesh (ignoring):\n",
                  iatom, apos[0], apos[1], apos[2]);
                Vnm_print(2, "Vpmg_fillco:    xmin = %g, xmax = %g\n", 
                  xmin, xmax);
                Vnm_print(2, "Vpmg_fillco:    ymin = %g, ymax = %g\n", 
                  ymin, ymax);
                Vnm_print(2, "Vpmg_fillco:    zmin = %g, zmax = %g\n", 
                  zmin, zmax);
            }
            fflush(stderr);

        } else if (arad > VPMGSMALL ) { /* if we're on the mesh */

            /* Convert the atom position to grid reference frame */
            position[0] = apos[0] - xmin;
            position[1] = apos[1] - ymin;
            position[2] = apos[2] - zmin;

            /* MARK ION ACCESSIBILITY AND DIELECTRIC VALUES FOR LATER
             * ASSIGNMENT (Steps #1-3) */
            itot = irad + arad + splineWin;
            itot2 = VSQR(itot2);     
            ictot = VMAX2(0, (irad + arad - splineWin));
            ictot2 = VSQR(ictot2);
            stot = arad + splineWin;
            stot2 = VSQR(stot);
            sctot = VMAX2(0, (arad - splineWin));
            sctot2 = VSQR(sctot);

           /* We'll search over grid points which are in the greater of
             * these two radii */
            rtot = VMAX2(itot, stot);
            rtot2 = VMAX2(itot2, stot2);
            dx = rtot + 0.5*hx;
            dy = rtot + 0.5*hy;
            dz = rtot + 0.5*hzed;
            imin = VMAX2(0,(int)floor((position[0] - dx)/hx));
            imax = VMIN2(nx-1,(int)ceil((position[0] + dx)/hx));
            jmin = VMAX2(0,(int)floor((position[1] - dy)/hy));
            jmax = VMIN2(ny-1,(int)ceil((position[1] + dy)/hy));
            kmin = VMAX2(0,(int)floor((position[2] - dz)/hzed));
            kmax = VMIN2(nz-1,(int)ceil((position[2] + dz)/hzed));
            for (i=imin; i<=imax; i++) {
                dx2 = VSQR(position[0] - hx*i);
                for (j=jmin; j<=jmax; j++) {
                    dy2 = VSQR(position[1] - hy*j);
                    for (k=kmin; k<=kmax; k++) {
                        dz2 = VSQR(position[2] - k*hzed);

                        /* ASSIGN CCF */
                        if (thee->ccf[IJK(i,j,k)] > VPMGSMALL) {
                            dist2 = dz2 + dy2 + dx2;
                            if (dist2 >= itot2) {
                                thee->ccf[IJK(i,j,k)] *= 1.0;
                            } 
                            if (dist2 <= ictot2) {
                                thee->ccf[IJK(i,j,k)] = 0.0;
                            }
                            if ((dist2 < itot2) && (dist2 > ictot2)) {
                                dist = VSQRT(dist2);
                                sm = dist - (arad + irad) + splineWin;
                                sm2 = VSQR(sm);
                                value = 0.75*sm2*w2i - 0.25*sm*sm2*w3i;
                                thee->ccf[IJK(i,j,k)] *= value;
                            }
                        }

                        /* ASSIGN A1CF */
                        if (thee->a1cf[IJK(i,j,k)] > VPMGSMALL) {
                            dist2 = dz2+dy2+VSQR(position[0]-(i+0.5)*hx);
                            if (dist2 >= stot2) {
                                thee->a1cf[IJK(i,j,k)] *= 1.0;
                            } 
                            if (dist2 <= sctot2) {
                                thee->a1cf[IJK(i,j,k)] = 0.0;
                            } 
                            if ((dist2 > sctot2) && (dist2 < stot2)) {
                                dist = VSQRT(dist2);
                                sm = dist - arad + splineWin;
                                sm2 = VSQR(sm);
                                value = 0.75*sm2*w2i - 0.25*sm*sm2*w3i;
                                thee->a1cf[IJK(i,j,k)] *= value;
                            } 
                        }

                        /* ASSIGN A2CF */
                        if (thee->a2cf[IJK(i,j,k)] > VPMGSMALL) {
                            dist2 = dz2+dx2+VSQR(position[1]-(j+0.5)*hy);
                            if (dist2 >= stot2) {
                                thee->a2cf[IJK(i,j,k)] *= 1.0;
                            } 
                            if (dist2 <= sctot2) {
                                thee->a2cf[IJK(i,j,k)] = 0.0;
                            }
                            if ((dist2 > sctot2) && (dist2 < stot2)) {
                                dist = VSQRT(dist2);
                                sm = dist - arad + splineWin;
                                sm2 = VSQR(sm);
                                value = 0.75*sm2*w2i - 0.25*sm*sm2*w3i;
                                thee->a2cf[IJK(i,j,k)] *= value;
                            }
                        }

                        /* ASSIGN A3CF */
                        if (thee->a3cf[IJK(i,j,k)] > VPMGSMALL) {
                            dist2 = dy2+dx2+VSQR(position[2]-(k+0.5)*hzed);
                            if (dist2 >= stot2) {
                                thee->a3cf[IJK(i,j,k)] *= 1.0;
                            } 
                            if (dist2 <= sctot2) {
                                thee->a3cf[IJK(i,j,k)] = 0.0;
                            } 
                            if ((dist2 > sctot2) && (dist2 < stot2)) {
                                dist = VSQRT(dist2);
                                sm = dist - arad + splineWin;
                                sm2 = VSQR(sm);
                                value = 0.75*sm2*w2i - 0.25*sm*sm2*w3i;
                                thee->a3cf[IJK(i,j,k)] *= value;
                            }
                        }


                    } /* k loop */
                } /* j loop */
            } /* i loop */
        } /* endif (on the mesh) */
    } /* endfor (over all atoms) */

    Vnm_print(0, "Vpmg_fillco:  filling coefficient arrays\n");
    /* Interpret markings and fill the coefficient arrays */
    for (k=0; k<nz; k++) {
        for (j=0; j<ny; j++) {
            for (i=0; i<nx; i++) {

                thee->ccf[IJK(i,j,k)] = ionmask*thee->ccf[IJK(i,j,k)];
                thee->a1cf[IJK(i,j,k)] = (epsw-epsp)*thee->a1cf[IJK(i,j,k)] 
                  + epsp;
                thee->a2cf[IJK(i,j,k)] = (epsw-epsp)*thee->a2cf[IJK(i,j,k)] 
                  + epsp;
                thee->a3cf[IJK(i,j,k)] = (epsw-epsp)*thee->a3cf[IJK(i,j,k)] 
                  + epsp;

#if 0
                Vnm_print(2, "ccf = %g\n", thee->ccf[IJK(i,j,k)]);
                Vnm_print(2, "a2cf = %g\n", thee->a2cf[IJK(i,j,k)]);
                Vnm_print(2, "a1cf = %g\n", thee->a1cf[IJK(i,j,k)]);
                Vnm_print(2, "a3cf = %g\n", thee->a3cf[IJK(i,j,k)]);
#endif

            } /* i loop */
        } /* j loop */
    } /* k loop */

}

VPRIVATE void fillcoCoef(Vpmg *thee) {

    VASSERT(thee != VNULL);

    if (thee->useDielXMap || thee->useDielYMap || thee->useDielZMap ||
      thee->useKappaMap) {
        fillcoCoefMap(thee);
        return;
    }

    switch(thee->surfMeth) {
        case VSM_MOL:
            Vnm_print(0, "fillcoCoef:  Calling fillcoCoefMol...\n");
            fillcoCoefMol(thee);
            break;
        case VSM_MOLSMOOTH:
            Vnm_print(0, "fillcoCoef:  Calling fillcoCoefMol...\n");
            fillcoCoefMol(thee);
            break;
        case VSM_SPLINE:
            Vnm_print(0, "fillcoCoef:  Calling fillcoCoefSpline...\n");
            fillcoCoefSpline(thee);
            break;
        default:
            Vnm_print(2, "fillcoCoef:  Invalid surfMeth (%d)!\n",
              thee->surfMeth);
            VASSERT(0);
            break;
    }
}


VPRIVATE void fillcoCharge(Vpmg *thee) {

    VASSERT(thee != VNULL);

    if (thee->useChargeMap) {
        fillcoChargeMap(thee);
        return;
    }

    switch(thee->chargeMeth) {
        case VCM_TRIL:
            Vnm_print(0, "fillcoCharge:  Calling fillcoChargeSpline1...\n");
            fillcoChargeSpline1(thee);
            break;
        case VCM_BSPL2:
            Vnm_print(0, "fillcoCharge:  Calling fillcoChargeSpline2...\n");
            fillcoChargeSpline2(thee);
            break;
        default:
            Vnm_print(2, "fillcoCharge:  Invalid chargeMeth (%d)!\n",
              thee->chargeMeth);
            VASSERT(0);
            break;
    }
}

VPRIVATE void fillcoChargeMap(Vpmg *thee) {

    Vpbe *pbe;
    double position[3], charge, zmagic, hx, hy, hzed;
    int i, j, k, nx, ny, nz;


    VASSERT(thee != VNULL);

    /* Get PBE info */
    pbe = thee->pbe;
    zmagic = Vpbe_getZmagic(pbe);

    /* Mesh info */
    nx = thee->pmgp->nx;
    ny = thee->pmgp->ny;
    nz = thee->pmgp->nz;
    hx = thee->pmgp->hx;
    hy = thee->pmgp->hy;
    hzed = thee->pmgp->hzed;
   
    /* Reset the fcf, tcf, ccf, a1cf, a2cf, and a3cf arrays */
    for (i=0; i<(nx*ny*nz); i++) thee->fcf[i] = 0.0;

    /* Fill in the source term (atomic charges) */
    Vnm_print(0, "Vpmg_fillco:  filling in source term.\n");
    for (k=0; k<nz; k++) {
        for (j=0; j<ny; j++) {
            for (i=0; i<nx; i++) {
                position[0] = thee->xf[i];
                position[1] = thee->yf[j];
                position[2] = thee->zf[k];
                VASSERT(Vgrid_value(thee->chargeMap, position, &charge));
                /* Scale the charge to internal units */
                charge = charge*zmagic;
                thee->fcf[IJK(i,j,k)] = charge;
            }
        }
    }
}

VPRIVATE void fillcoChargeSpline1(Vpmg *thee) {

    Valist *alist;
    Vpbe *pbe;
    Vatom *atom;
    double xmin, xmax, ymin, ymax, zmin, zmax;
    double xlen, ylen, zlen, position[3], ifloat, jfloat, kfloat;
    double charge, dx, dy, dz, zmagic, hx, hy, hzed, *apos;
    int i, nx, ny, nz, iatom, ihi, ilo, jhi, jlo, khi, klo;


    VASSERT(thee != VNULL);

    /* Get PBE info */
    pbe = thee->pbe;
    alist = pbe->alist;
    zmagic = Vpbe_getZmagic(pbe);

    /* Mesh info */
    nx = thee->pmgp->nx;
    ny = thee->pmgp->ny;
    nz = thee->pmgp->nz;
    hx = thee->pmgp->hx;
    hy = thee->pmgp->hy;
    hzed = thee->pmgp->hzed;
   
    /* Define the total domain size */
    xlen = thee->pmgp->xlen;
    ylen = thee->pmgp->ylen;
    zlen = thee->pmgp->zlen;

    /* Define the min/max dimensions */
    xmin = thee->pmgp->xcent - (xlen/2.0);
    ymin = thee->pmgp->ycent - (ylen/2.0);
    zmin = thee->pmgp->zcent - (zlen/2.0);
    xmax = thee->pmgp->xcent + (xlen/2.0);
    ymax = thee->pmgp->ycent + (ylen/2.0);
    zmax = thee->pmgp->zcent + (zlen/2.0);

    /* Reset the fcf, tcf, ccf, a1cf, a2cf, and a3cf arrays */
    for (i=0; i<(nx*ny*nz); i++) thee->fcf[i] = 0.0;

    /* Fill in the source term (atomic charges) */
    Vnm_print(0, "Vpmg_fillco:  filling in source term.\n");
    for (iatom=0; iatom<Valist_getNumberAtoms(alist); iatom++) {

        atom = Valist_getAtom(alist, iatom);
        apos = Vatom_getPosition(atom);
        charge = Vatom_getCharge(atom);

        /* Make sure we're on the grid */
        if ((apos[0]<=xmin) || (apos[0]>=xmax)  || \
            (apos[1]<=ymin) || (apos[1]>=ymax)  || \
            (apos[2]<=zmin) || (apos[2]>=zmax)) {
            if (thee->pmgp->bcfl != BCFL_FOCUS) {
                Vnm_print(2, "Vpmg_fillco:  Atom #%d at (%4.3f, %4.3f, \
%4.3f) is off the mesh (ignoring):\n",
                  iatom, apos[0], apos[1], apos[2]);
                Vnm_print(2, "Vpmg_fillco:    xmin = %g, xmax = %g\n", 
                  xmin, xmax);
                Vnm_print(2, "Vpmg_fillco:    ymin = %g, ymax = %g\n", 
                  ymin, ymax);
                Vnm_print(2, "Vpmg_fillco:    zmin = %g, zmax = %g\n", 
                  zmin, zmax);
            }
            fflush(stderr);
        } else {

            /* Convert the atom position to grid reference frame */
            position[0] = apos[0] - xmin;
            position[1] = apos[1] - ymin;
            position[2] = apos[2] - zmin;

            /* Scale the charge to be a delta function */
            charge = charge*zmagic/(hx*hy*hzed);

            /* Figure out which vertices we're next to */
            ifloat = position[0]/hx;
            jfloat = position[1]/hy;
            kfloat = position[2]/hzed;

            ihi = (int)ceil(ifloat);
            ilo = (int)floor(ifloat);
            jhi = (int)ceil(jfloat);
            jlo = (int)floor(jfloat);
            khi = (int)ceil(kfloat);
            klo = (int)floor(kfloat);

            /* Now assign fractions of the charge to the nearby verts */
            dx = ifloat - (double)(ilo);
            dy = jfloat - (double)(jlo);
            dz = kfloat - (double)(klo);
            thee->fcf[IJK(ihi,jhi,khi)] += (dx      *     dy *     dz *charge);
            thee->fcf[IJK(ihi,jlo,khi)] += (dx      *(1.0-dy)*     dz *charge);
            thee->fcf[IJK(ihi,jhi,klo)] += (dx      *     dy *(1.0-dz)*charge);
            thee->fcf[IJK(ihi,jlo,klo)] += (dx      *(1.0-dy)*(1.0-dz)*charge);
            thee->fcf[IJK(ilo,jhi,khi)] += ((1.0-dx)*     dy *     dz *charge);
            thee->fcf[IJK(ilo,jlo,khi)] += ((1.0-dx)*(1.0-dy)*     dz *charge);
            thee->fcf[IJK(ilo,jhi,klo)] += ((1.0-dx)*     dy *(1.0-dz)*charge);
            thee->fcf[IJK(ilo,jlo,klo)] += ((1.0-dx)*(1.0-dy)*(1.0-dz)*charge);
        } /* endif (on the mesh) */
    } /* endfor (each atom) */
}

VPRIVATE double bspline2(double x) {

    double m2m, m2, m3; 

    if ((x >= 0.0) && (x <= 2.0)) m2m = 1.0 - VABS(x - 1.0);
    else m2m = 0.0;
    if ((x >= 1.0) && (x <= 3.0)) m2 = 1.0 - VABS(x - 2.0);
    else m2 = 0.0;

    if ((x >= 0.0) && (x <= 3.0)) m3 = 0.5*x*m2m + 0.5*(3.0-x)*m2;
    else m3 = 0.0;

    return m3;

}

VPRIVATE double dbspline2(double x) {

    double m2m, m2, dm3;

    if ((x >= 0.0) && (x <= 2.0)) m2m = 1.0 - VABS(x - 1.0);
    else m2m = 0.0;
    if ((x >= 1.0) && (x <= 3.0)) m2 = 1.0 - VABS(x - 2.0);
    else m2 = 0.0;

    dm3 = m2m - m2;

    return dm3;

}


VPRIVATE void fillcoChargeSpline2(Vpmg *thee) {

    Valist *alist;
    Vpbe *pbe;
    Vatom *atom;
    double xmin, xmax, ymin, ymax, zmin, zmax, zmagic;
    double xlen, ylen, zlen, position[3], ifloat, jfloat, kfloat;
    double charge, hx, hy, hzed, *apos, mx, my, mz;
    int i, ii, jj, kk, nx, ny, nz, iatom; 
    int im2, im1, ip1, ip2, jm2, jm1, jp1, jp2, km2, km1, kp1, kp2;


    VASSERT(thee != VNULL);

    /* Get PBE info */
    pbe = thee->pbe;
    alist = pbe->alist;
    zmagic = Vpbe_getZmagic(pbe);

    /* Mesh info */
    nx = thee->pmgp->nx;
    ny = thee->pmgp->ny;
    nz = thee->pmgp->nz;
    hx = thee->pmgp->hx;
    hy = thee->pmgp->hy;
    hzed = thee->pmgp->hzed;
   
    /* Define the total domain size */
    xlen = thee->pmgp->xlen;
    ylen = thee->pmgp->ylen;
    zlen = thee->pmgp->zlen;

    /* Define the min/max dimensions */
    xmin = thee->pmgp->xcent - (xlen/2.0);
    ymin = thee->pmgp->ycent - (ylen/2.0);
    zmin = thee->pmgp->zcent - (zlen/2.0);
    xmax = thee->pmgp->xcent + (xlen/2.0);
    ymax = thee->pmgp->ycent + (ylen/2.0);
    zmax = thee->pmgp->zcent + (zlen/2.0);

    /* Reset the fcf, tcf, ccf, a1cf, a2cf, and a3cf arrays */
    for (i=0; i<(nx*ny*nz); i++) thee->fcf[i] = 0.0;

    /* Fill in the source term (atomic charges) */
    Vnm_print(0, "Vpmg_fillco:  filling in source term.\n");
    for (iatom=0; iatom<Valist_getNumberAtoms(alist); iatom++) {

        atom = Valist_getAtom(alist, iatom);
        apos = Vatom_getPosition(atom);
        charge = Vatom_getCharge(atom);

        /* Make sure we're on the grid */
        if ((apos[0]<=(xmin-hx)) || (apos[0]>=(xmax+hx))  || \
            (apos[1]<=(ymin-hy)) || (apos[1]>=(ymax+hy))  || \
            (apos[2]<=(zmin-hzed)) || (apos[2]>=(zmax+hzed))) {
            if (thee->pmgp->bcfl != BCFL_FOCUS) {
                Vnm_print(2, "Vpmg_fillco:  Atom #%d at (%4.3f, %4.3f, \
%4.3f) is off the mesh (for cubic splines!!) (ignoring this atom):\n",
                  iatom, apos[0], apos[1], apos[2]);
                Vnm_print(2, "Vpmg_fillco:    xmin = %g, xmax = %g\n", 
                  xmin, xmax);
                Vnm_print(2, "Vpmg_fillco:    ymin = %g, ymax = %g\n", 
                  ymin, ymax);
                Vnm_print(2, "Vpmg_fillco:    zmin = %g, zmax = %g\n", 
                  zmin, zmax);
            }
            fflush(stderr);
        } else {

            /* Convert the atom position to grid reference frame */
            position[0] = apos[0] - xmin;
            position[1] = apos[1] - ymin;
            position[2] = apos[2] - zmin;

            /* Scale the charge to be a delta function */
            charge = charge*zmagic/(hx*hy*hzed);

            /* Figure out which vertices we're next to */
            ifloat = position[0]/hx;
            jfloat = position[1]/hy;
            kfloat = position[2]/hzed;

            ip1   = (int)ceil(ifloat);
            ip2   = ip1 + 1;
            im1   = (int)floor(ifloat);
            im2   = im1 - 1;
            jp1   = (int)ceil(jfloat);
            jp2   = jp1 + 1;
            jm1   = (int)floor(jfloat);
            jm2   = jm1 - 1;
            kp1   = (int)ceil(kfloat);
            kp2   = kp1 + 1;
            km1   = (int)floor(kfloat);
            km2   = km1 - 1;

            /* This step shouldn't be necessary, but it saves nasty debugging
             * later on if something goes wrong */
            ip2 = VMIN2(ip2,nx-1);
            ip1 = VMIN2(ip1,nx-1);
            im1 = VMAX2(im1,0);
            im2 = VMAX2(im2,0);
            jp2 = VMIN2(jp2,ny-1);
            jp1 = VMIN2(jp1,ny-1);
            jm1 = VMAX2(jm1,0);
            jm2 = VMAX2(jm2,0);
            kp2 = VMIN2(kp2,nz-1);
            kp1 = VMIN2(kp1,nz-1);
            km1 = VMAX2(km1,0);
            km2 = VMAX2(km2,0);

            /* Now assign fractions of the charge to the nearby verts */
            for (ii=im2; ii<=ip2; ii++) {
                mx = bspline2(VFCHI(ii,ifloat));
                for (jj=jm2; jj<=jp2; jj++) {
                    my = bspline2(VFCHI(jj,jfloat));
                    for (kk=km2; kk<=kp2; kk++) {
                        mz = bspline2(VFCHI(kk,kfloat));
                        thee->fcf[IJK(ii,jj,kk)] += (charge*mx*my*mz);
                    }
                }
            }

        } /* endif (on the mesh) */
    } /* endfor (each atom) */
}

VPUBLIC void Vpmg_fillco(Vpmg *thee, 
  Vsurf_Meth surfMeth, double splineWin, Vchrg_Meth chargeMeth,
  int useDielXMap,   Vgrid *dielXMap, 
  int useDielYMap,   Vgrid *dielYMap, 
  int useDielZMap,   Vgrid *dielZMap, 
  int useKappaMap,   Vgrid *kappaMap, 
  int useChargeMap,  Vgrid *chargeMap) {

    Vpbe *pbe;
    double xmin, xmax, ymin, ymax, zmin, zmax;
    double xlen, ylen, zlen, hx, hy, hzed;
    double epsw, epsp, ionstr;
    int i, nx, ny, nz, islap;

    VASSERT(thee != VNULL);
    thee->surfMeth = surfMeth;
    thee->splineWin = splineWin;
    thee->chargeMeth = chargeMeth;
    thee->useDielXMap = useDielXMap;
    if (thee->useDielXMap) thee->dielXMap = dielXMap;
    thee->useDielYMap = useDielYMap;
    if (thee->useDielYMap) thee->dielYMap = dielYMap;
    thee->useDielZMap = useDielZMap;
    if (thee->useDielZMap) thee->dielZMap = dielZMap;
    thee->useKappaMap = useKappaMap;
    if (thee->useKappaMap) thee->kappaMap = kappaMap;
    thee->useChargeMap = useChargeMap;
    if (thee->useChargeMap) thee->chargeMap = chargeMap;

    /* Get PBE info */
    pbe = thee->pbe;
    ionstr = Vpbe_getBulkIonicStrength(pbe);
    epsw = Vpbe_getSolventDiel(pbe);
    epsp = Vpbe_getSoluteDiel(pbe);

    /* Mesh info */
    nx = thee->pmgp->nx;
    ny = thee->pmgp->ny;
    nz = thee->pmgp->nz;
    hx = thee->pmgp->hx;
    hy = thee->pmgp->hy;
    hzed = thee->pmgp->hzed;
   
    /* Define the total domain size */
    xlen = thee->pmgp->xlen;
    ylen = thee->pmgp->ylen;
    zlen = thee->pmgp->zlen;

    /* Define the min/max dimensions */
    xmin = thee->pmgp->xcent - (xlen/2.0);
    thee->pmgp->xmin = xmin;
    ymin = thee->pmgp->ycent - (ylen/2.0);
    thee->pmgp->ymin = ymin;
    zmin = thee->pmgp->zcent - (zlen/2.0);
    thee->pmgp->zmin = zmin;
    xmax = thee->pmgp->xcent + (xlen/2.0);
    thee->pmgp->xmax = xmax;
    ymax = thee->pmgp->ycent + (ylen/2.0);
    thee->pmgp->ymax = ymax;
    zmax = thee->pmgp->zcent + (zlen/2.0);
    thee->pmgp->zmax = zmax;
    thee->rparm[2] = xmin;
    thee->rparm[3] = xmax;
    thee->rparm[4] = ymin;
    thee->rparm[5] = ymax;
    thee->rparm[6] = zmin;
    thee->rparm[7] = zmax;

    /* This is a flag that gets set if the operator is a simple Laplacian;
     * i.e., in the case of a homogenous dielectric and zero ionic strength */
    if ((ionstr < VPMGSMALL) && (VABS(epsp-epsw) < VPMGSMALL)) islap = 1;
    else islap = 0;

    /* Fill the mesh point coordinate arrays */
    for (i=0; i<nx; i++) thee->xf[i] = xmin + i*hx;
    for (i=0; i<ny; i++) thee->yf[i] = ymin + i*hy;
    for (i=0; i<nz; i++) thee->zf[i] = zmin + i*hzed;

    /* Reset the fcf, tcf, ccf, a1cf, a2cf, and a3cf arrays */
    for (i=0; i<(nx*ny*nz); i++) thee->tcf[i] = 0.0;

    /* Fill in the source term (atomic charges) */
    Vnm_print(0, "Vpmg_fillco:  filling in source term.\n");
    fillcoCharge(thee);

    /* THE FOLLOWING NEEDS TO BE DONE IF WE'RE NOT USING A SIMPLE LAPLACIAN
     * OPERATOR */
    if (!islap) {

        Vnm_print(0, "Vpmg_fillco:  marking ion and solvent accessibility.\n");
        fillcoCoef(thee);
        Vnm_print(0, "Vpmg_fillco:  done filling coefficient arrays\n");

    } else { /* else (!islap) ==> It's a Laplacian operator! */

        for (i=0; i<(nx*ny*nz); i++) {
            thee->ccf[i] = 0.0;
            thee->a1cf[i] = epsp;
            thee->a2cf[i] = epsp;
            thee->a3cf[i] = epsp;
        }

    } /* endif (!islap) */

    /* Fill the boundary arrays (except when focusing, bcfl = 4) */
    if (thee->pmgp->bcfl != BCFL_FOCUS) {
        Vnm_print(0, "Vpmg_fillco:  filling boundary arrays\n");
        bcCalc(thee);
        Vnm_print(0, "Vpmg_fillco:  done filling boundary arrays\n");
    }

    thee->filled = 1;
}


VPUBLIC void Vpmg_force(Vpmg *thee, double *force, int atomID, 
  Vsurf_Meth srfm, Vchrg_Meth chgm) {

    double qfF[3];                  /* Charge-field force */  
    double dbF[3];                  /* Dielectric boundary force */
    double ibF[3];                  /* Ion boundary force */
    double npF[3];                  /* Non-polar boundary force */

    VASSERT(thee != VNULL);
 
    Vpmg_dbnpForce(thee, qfF, npF, atomID, srfm);
    Vpmg_ibForce(thee, dbF, atomID, srfm); 
    Vpmg_qfForce(thee, ibF, atomID, chgm); 

    force[0] = qfF[0] + dbF[0] + npF[0] + ibF[0];
    force[1] = qfF[1] + dbF[1] + npF[1] + ibF[1];
    force[2] = qfF[2] + dbF[2] + npF[2] + ibF[2];

}

VPUBLIC void Vpmg_ibForce(Vpmg *thee, double *force, int atomID, 
  Vsurf_Meth srfm) {

    Valist *alist;
    Vacc *acc;
    Vpbe *pbe;
    Vatom *atom;

    double *apos, position[3], arad, irad, zkappa2, hx, hy, hzed;
    double xlen, ylen, zlen, xmin, ymin, zmin, xmax, ymax, zmax, rtot2;
    double rtot, dx, dx2, dy, dy2, dz, dz2, gpos[3], tgrad[3], fmag;
    double izmagic;
    int i, j, k, nx, ny, nz, imin, imax, jmin, jmax, kmin, kmax;
   
    VASSERT(thee != VNULL);
    /* VASSERT(thee->filled); */

   
    acc = thee->pbe->acc;
    atom = Valist_getAtom(thee->pbe->alist, atomID);
    apos = Vatom_getPosition(atom);
    arad = Vatom_getRadius(atom);

    /* Reset force */
    force[0] = 0.0;
    force[1] = 0.0;
    force[2] = 0.0;

    /* Check surface definition */
    if (srfm != VSM_SPLINE) {
        Vnm_print(2, "Vpmg_ibForce:  Forces *must* be calculated with \
spline-based surfaces!\n");
        Vnm_print(2, "Vpmg_ibForce:  Skipping ionic boundary force \
calculation!\n");
        return;
    }

    /* If we aren't in the current position, then we're done */
    if (atom->partID == 0) return;

    /* Get PBE info */
    pbe = thee->pbe;
    acc = pbe->acc;
    alist = pbe->alist;
    irad = Vpbe_getMaxIonRadius(pbe);
    zkappa2 = Vpbe_getZkappa2(pbe);
    izmagic = 1.0/Vpbe_getZmagic(pbe);

    /* Mesh info */
    nx = thee->pmgp->nx;
    ny = thee->pmgp->ny;
    nz = thee->pmgp->nz;
    hx = thee->pmgp->hx;
    hy = thee->pmgp->hy;
    hzed = thee->pmgp->hzed;
    xlen = thee->pmgp->xlen;
    ylen = thee->pmgp->ylen;
    zlen = thee->pmgp->zlen;
    xmin = thee->pmgp->xmin;
    ymin = thee->pmgp->ymin;
    zmin = thee->pmgp->zmin;
    xmax = thee->pmgp->xmax;
    ymax = thee->pmgp->ymax;
    zmax = thee->pmgp->zmax;

    /* Sanity check: there is no force if there is zero ionic strength */
    if (zkappa2 < VPMGSMALL) {
        Vnm_print(2, "Vpmg_ibForce:  No force for zero ionic strength!\n");
        return;
    }

    /* Make sure we're on the grid */
    if ((apos[0]<=xmin) || (apos[0]>=xmax)  || \
      (apos[1]<=ymin) || (apos[1]>=ymax)  || \
      (apos[2]<=zmin) || (apos[2]>=zmax)) {
        if (thee->pmgp->bcfl != BCFL_FOCUS) {
            Vnm_print(2, "Vpmg_ibForce:  Atom #%d at (%4.3f, %4.3f, %4.3f) is off the mesh (ignoring):\n",
                  atom, position[0], position[1], position[2]);
            Vnm_print(2, "Vpmg_ibForce:    xmin = %g, xmax = %g\n",
              xmin, xmax);
            Vnm_print(2, "Vpmg_ibForce:    ymin = %g, ymax = %g\n",
              ymin, ymax);
            Vnm_print(2, "Vpmg_ibForce:    zmin = %g, zmax = %g\n",
              zmin, zmax);
        }
        fflush(stderr);
    } else {

        /* Convert the atom position to grid reference frame */
        position[0] = apos[0] - xmin;
        position[1] = apos[1] - ymin;
        position[2] = apos[2] - zmin;

        /* Integrate over points within this atom's (inflated) radius */
        rtot = (irad + arad + thee->splineWin);
        rtot2 = VSQR(rtot);
        dx = rtot + 0.5*hx;
        imin = VMAX2(0,(int)ceil((position[0] - dx)/hx));
        imax = VMIN2(nx-1,(int)floor((position[0] + dx)/hx));
        for (i=imin; i<=imax; i++) { 
            dx2 = VSQR(position[0] - hx*i);
            if (rtot2 > dx2) dy = VSQRT(rtot2 - dx2) + 0.5*hy;
            else dy = 0.5*hy;
            jmin = VMAX2(0,(int)ceil((position[1] - dy)/hy));
            jmax = VMIN2(ny-1,(int)floor((position[1] + dy)/hy));
            for (j=jmin; j<=jmax; j++) { 
                dy2 = VSQR(position[1] - hy*j);
                if (rtot2 > (dx2+dy2)) dz = VSQRT(rtot2-dx2-dy2)+0.5*hzed;
                else dz = 0.5*hzed;
                kmin = VMAX2(0,(int)ceil((position[2] - dz)/hzed));
                kmax = VMIN2(nz-1,(int)floor((position[2] + dz)/hzed));
                for (k=kmin; k<=kmax; k++) {
                    dz2 = VSQR(k*hzed - position[2]);
                    /* See if grid point is inside ivdw radius and set ccf
                     * accordingly (do spline assignment here) */
                    if ((dz2 + dy2 + dx2) <= rtot2) {
                        gpos[0] = i*hx + xmin;
                        gpos[1] = j*hy + ymin;
                        gpos[2] = k*hzed + zmin;
                        Vacc_splineAccGradAtom(acc, gpos, thee->splineWin, irad,
                          atom, tgrad);
                        if (thee->pmgp->nonlin) {
                            /* Nonlinear forces not done */
                            Vnm_print(2, "Vpmg_ibForce:  No NPBE forces yet!\n");
                            force[0] = 0.0;
                            force[1] = 0.0;
                            force[2] = 0.0;
                            return;
                        } else {
                            /* Use of bulk factor (zkappa2) OK here becuase
                             * LPBE force approximation */
                            fmag = VSQR(thee->u[IJK(i,j,k)]);
                            force[0] += (zkappa2*fmag*tgrad[0]);
                            force[1] += (zkappa2*fmag*tgrad[1]);
                            force[2] += (zkappa2*fmag*tgrad[2]);
                        }
                    }
                } /* k loop */
            } /* j loop */
        } /* i loop */
    } 
    force[0] = force[0] * 0.5 * hx * hy * hzed * izmagic;
    force[1] = force[1] * 0.5 * hx * hy * hzed * izmagic;
    force[2] = force[2] * 0.5 * hx * hy * hzed * izmagic;
}

VPUBLIC void Vpmg_dbnpForce(Vpmg *thee, double *dbForce, double *npForce, 
  int atomID, Vsurf_Meth srfm) {

    Vacc *acc;
    Vpbe *pbe;
    Vatom *atom;

    double *apos, position[3], arad, hx, hy, hzed, izmagic, deps, depsi;
    double xlen, ylen, zlen, xmin, ymin, zmin, xmax, ymax, zmax, rtot2, epsp;
    double rtot, dx, gpos[3], tgrad[3], dbFmag, epsw, gamma, kT;
    double npFmag, *u, Hxijk, Hyijk, Hzijk, Hxim1jk, Hyijm1k, Hzijkm1;
    double dHxijk[3], dHyijk[3], dHzijk[3], dHxim1jk[3], dHyijm1k[3]; 
    double dHzijkm1[3];
    int i, j, k, l, nx, ny, nz, imin, imax, jmin, jmax, kmin, kmax;

    VASSERT(thee != VNULL);
    if (!thee->filled) Vpmg_fillco(thee,
      thee->surfMeth, thee->splineWin, thee->chargeMeth,
      thee->useDielXMap, thee->dielXMap,
      thee->useDielYMap, thee->dielYMap,
      thee->useDielZMap, thee->dielZMap,
      thee->useKappaMap, thee->kappaMap,
      thee->useChargeMap, thee->chargeMap);

    acc = thee->pbe->acc;
    atom = Valist_getAtom(thee->pbe->alist, atomID);
    apos = Vatom_getPosition(atom);
    arad = Vatom_getRadius(atom);

    /* Reset force */
    dbForce[0] = 0.0;
    dbForce[1] = 0.0;
    dbForce[2] = 0.0;
    npForce[0] = 0.0;
    npForce[1] = 0.0;
    npForce[2] = 0.0;

    /* Check surface definition */
    if (srfm != VSM_SPLINE) {
        Vnm_print(2, "Vpmg_dbnpForce:  Forces *must* be calculated with \
spline-based surfaces!\n");
        Vnm_print(2, "Vpmg_dbnpForce:  Skipping dielectric/apolar boundary \
force calculation!\n");
        return;
    }


    /* If we aren't in the current position, then we're done */
    if (atom->partID == 0) return;

    /* Get PBE info */
    pbe = thee->pbe;
    acc = pbe->acc;
    epsp = Vpbe_getSoluteDiel(pbe);
    epsw = Vpbe_getSolventDiel(pbe);
    kT = Vpbe_getTemperature(pbe)*(1e-3)*Vunit_Na*Vunit_kb;
    gamma = Vpbe_getGamma(pbe)/kT;
    izmagic = 1.0/Vpbe_getZmagic(pbe);

    /* Mesh info */
    nx = thee->pmgp->nx;
    ny = thee->pmgp->ny;
    nz = thee->pmgp->nz;
    hx = thee->pmgp->hx;
    hy = thee->pmgp->hy;
    hzed = thee->pmgp->hzed;
    xlen = thee->pmgp->xlen;
    ylen = thee->pmgp->ylen;
    zlen = thee->pmgp->zlen;
    xmin = thee->pmgp->xmin;
    ymin = thee->pmgp->ymin;
    zmin = thee->pmgp->zmin;
    xmax = thee->pmgp->xmax;
    ymax = thee->pmgp->ymax;
    zmax = thee->pmgp->zmax;
    u = thee->u;

    /* Sanity check: there is no force if there is zero ionic strength */
    if (VABS(epsp-epsw) < VPMGSMALL) {
       Vnm_print(0, "Vpmg_dbnpForce: No force for uniform dielectric!\n");
       return;
    }
    deps = (epsw - epsp);
    depsi = 1.0/deps;

    /* Make sure we're on the grid */
    if ((apos[0]<=xmin) || (apos[0]>=xmax)  || \
      (apos[1]<=ymin) || (apos[1]>=ymax)  || \
      (apos[2]<=zmin) || (apos[2]>=zmax)) {
        if (thee->pmgp->bcfl != BCFL_FOCUS) {
            Vnm_print(2, "Vpmg_dbnpForce:  Atom #%d at (%4.3f, %4.3f, %4.3f) is off the mesh (ignoring):\n",
                  atomID, position[0], position[1], position[2]);
            Vnm_print(2, "Vpmg_dbnpForce:    xmin = %g, xmax = %g\n",
              xmin, xmax);
            Vnm_print(2, "Vpmg_dbnpForce:    ymin = %g, ymax = %g\n",
              ymin, ymax);
            Vnm_print(2, "Vpmg_dbnpForce:    zmin = %g, zmax = %g\n",
              zmin, zmax);
        }
        fflush(stderr);
    } else {

        /* Convert the atom position to grid reference frame */
        position[0] = apos[0] - xmin;
        position[1] = apos[1] - ymin;
        position[2] = apos[2] - zmin;

        /* Integrate over points within this atom's (inflated) radius */
        rtot = (arad + thee->splineWin);
        rtot2 = VSQR(rtot);
        dx = rtot/hx;
        imin = (int)floor((position[0]-rtot)/hx);
        if (imin < 1) {
            Vnm_print(2, "Vpmg_dbnpForce:  Atom %d off grid!\n", atomID); 
            return;
        }
        imax = (int)ceil((position[0]+rtot)/hx);
        if (imax > (nx-2)) {
            Vnm_print(2, "Vpmg_dbnpForce:  Atom %d off grid!\n", atomID); 
            return;
        }
        jmin = (int)floor((position[1]-rtot)/hy);
        if (jmin < 1) {
            Vnm_print(2, "Vpmg_dbnpForce:  Atom %d off grid!\n", atomID); 
            return;
        }
        jmax = (int)ceil((position[1]+rtot)/hy);
        if (jmax > (ny-2)) {
            Vnm_print(2, "Vpmg_dbnpForce:  Atom %d off grid!\n", atomID); 
            return;
        }
        kmin = (int)floor((position[2]-rtot)/hzed);
        if (kmin < 1) {
            Vnm_print(2, "Vpmg_dbnpForce:  Atom %d off grid!\n", atomID); 
            return;
        }
        kmax = (int)ceil((position[2]+rtot)/hzed);
        if (kmax > (nz-2)) {
            Vnm_print(2, "Vpmg_dbnpForce:  Atom %d off grid!\n", atomID); 
            return;
        }
        for (i=imin; i<=imax; i++) {
            for (j=jmin; j<=jmax; j++) {
                for (k=kmin; k<=kmax; k++) {
                    /* i,j,k */
                    gpos[0] = (i+0.5)*hx + xmin;
                    gpos[1] = j*hy + ymin;
                    gpos[2] = k*hzed + zmin;
                    Hxijk = (thee->a1cf[IJK(i,j,k)] - epsp)*depsi;
                    Vacc_splineAccGradAtom(acc, gpos, thee->splineWin, 0., 
                            atom, dHxijk);
                    for (l=0; l<3; l++) dHxijk[l] *= Hxijk;
                    gpos[0] = i*hx + xmin;
                    gpos[1] = (j+0.5)*hy + ymin;
                    gpos[2] = k*hzed + zmin;
                    Hyijk = (thee->a2cf[IJK(i,j,k)] - epsp)*depsi;
                    Vacc_splineAccGradAtom(acc, gpos, thee->splineWin, 0., 
                            atom, dHyijk);
                    for (l=0; l<3; l++) dHyijk[l] *= Hyijk;
                    gpos[0] = i*hx + xmin;
                    gpos[1] = j*hy + ymin;
                    gpos[2] = (k+0.5)*hzed + zmin;
                    Hzijk = (thee->a3cf[IJK(i,j,k)] - epsp)*depsi;
                    Vacc_splineAccGradAtom(acc, gpos, thee->splineWin, 0., 
                            atom, dHzijk);
                    for (l=0; l<3; l++) dHzijk[l] *= Hzijk;
                    /* i-1,j,k */
                    gpos[0] = (i-0.5)*hx + xmin;
                    gpos[1] = j*hy + ymin;
                    gpos[2] = k*hzed + zmin;
                    Hxim1jk = (thee->a1cf[IJK(i-1,j,k)] - epsp)*depsi;
                    Vacc_splineAccGradAtom(acc, gpos, thee->splineWin, 0.,
                            atom, dHxim1jk);
                    for (l=0; l<3; l++) dHxim1jk[l] *= Hxim1jk;
                    /* i,j-1,k */
                    gpos[0] = i*hx + xmin;
                    gpos[1] = (j-0.5)*hy + ymin;
                    gpos[2] = k*hzed + zmin;
                    Hyijm1k = (thee->a2cf[IJK(i,j-1,k)] - epsp)*depsi;
                    Vacc_splineAccGradAtom(acc, gpos, thee->splineWin, 0.,
                            atom, dHyijm1k);
                    for (l=0; l<3; l++) dHyijm1k[l] *= Hyijm1k;
                    /* i,j,k-1 */
                    gpos[0] = i*hx + xmin;
                    gpos[1] = j*hy + ymin;
                    gpos[2] = (k-0.5)*hzed + zmin;
                    Hzijkm1 = (thee->a3cf[IJK(i,j,k-1)] - epsp)*depsi;
                    Vacc_splineAccGradAtom(acc, gpos, thee->splineWin, 0.,
                            atom, dHzijkm1);
                    for (l=0; l<3; l++) dHzijkm1[l] *= Hzijkm1;
                    /* *** CALCULATE DIELECTRIC BOUNDARY FORCES *** */
                    dbFmag = u[IJK(i,j,k)];
                    tgrad[0] = 
                       (dHxijk[0]  *(u[IJK(i+1,j,k)]-u[IJK(i,j,k)])
                     +  dHxim1jk[0]*(u[IJK(i-1,j,k)]-u[IJK(i,j,k)]))/VSQR(hx)
                     + (dHyijk[0]  *(u[IJK(i,j+1,k)]-u[IJK(i,j,k)])
                     +  dHyijm1k[0]*(u[IJK(i,j-1,k)]-u[IJK(i,j,k)]))/VSQR(hy)
                     + (dHzijk[0]  *(u[IJK(i,j,k+1)]-u[IJK(i,j,k)])
                     + dHzijkm1[0]*(u[IJK(i,j,k-1)]-u[IJK(i,j,k)]))/VSQR(hzed);
                    tgrad[1] = 
                       (dHxijk[1]  *(u[IJK(i+1,j,k)]-u[IJK(i,j,k)])
                     +  dHxim1jk[1]*(u[IJK(i-1,j,k)]-u[IJK(i,j,k)]))/VSQR(hx)
                     + (dHyijk[1]  *(u[IJK(i,j+1,k)]-u[IJK(i,j,k)])
                     +  dHyijm1k[1]*(u[IJK(i,j-1,k)]-u[IJK(i,j,k)]))/VSQR(hy)
                     + (dHzijk[1]  *(u[IJK(i,j,k+1)]-u[IJK(i,j,k)])
                     + dHzijkm1[1]*(u[IJK(i,j,k-1)]-u[IJK(i,j,k)]))/VSQR(hzed);
                    tgrad[2] = 
                       (dHxijk[2]  *(u[IJK(i+1,j,k)]-u[IJK(i,j,k)])
                     +  dHxim1jk[2]*(u[IJK(i-1,j,k)]-u[IJK(i,j,k)]))/VSQR(hx)
                     + (dHyijk[2]  *(u[IJK(i,j+1,k)]-u[IJK(i,j,k)])
                     +  dHyijm1k[2]*(u[IJK(i,j-1,k)]-u[IJK(i,j,k)]))/VSQR(hy)
                     + (dHzijk[2]  *(u[IJK(i,j,k+1)]-u[IJK(i,j,k)])
                     + dHzijkm1[2]*(u[IJK(i,j,k-1)]-u[IJK(i,j,k)]))/VSQR(hzed);
                     dbForce[0] += (dbFmag*tgrad[0]);
                     dbForce[1] += (dbFmag*tgrad[1]);
                     dbForce[2] += (dbFmag*tgrad[2]);
                    /* *** CALCULATE NONPOLAR FORCES *** */
                    /* First we calculate the local H1-seminorm of the
                     * characteristic function */
                    npFmag =  VSQR((Hxijk - Hxim1jk)/hx)
                            + VSQR((Hyijk - Hyijm1k)/hy)
                            + VSQR((Hzijk - Hzijkm1)/hzed);
                    npFmag = VSQRT(npFmag);
                    if (npFmag > VPMGSMALL) {
                        tgrad[0] = 
                          (Hxijk-Hxim1jk)*(dHxijk[0]-dHxim1jk[0])/VSQR(hx)
                        + (Hyijk-Hyijm1k)*(dHyijk[0]-dHyijm1k[0])/VSQR(hy)
                        + (Hzijk-Hzijkm1)*(dHzijk[0]-dHzijkm1[0])/VSQR(hzed);
                        tgrad[1] = 
                          (Hxijk-Hxim1jk)*(dHxijk[1]-dHxim1jk[1])/VSQR(hx)
                        + (Hyijk-Hyijm1k)*(dHyijk[1]-dHyijm1k[1])/VSQR(hy)
                        + (Hzijk-Hzijkm1)*(dHzijk[1]-dHzijkm1[1])/VSQR(hzed);
                        tgrad[2] = 
                          (Hxijk-Hxim1jk)*(dHxijk[2]-dHxim1jk[2])/VSQR(hx)
                        + (Hyijk-Hyijm1k)*(dHyijk[2]-dHyijm1k[2])/VSQR(hy)
                        + (Hzijk-Hzijkm1)*(dHzijk[2]-dHzijkm1[2])/VSQR(hzed);
                        npForce[0] += (tgrad[0]/npFmag);
                        npForce[1] += (tgrad[1]/npFmag);
                        npForce[2] += (tgrad[2]/npFmag);
                    } 
                } /* k loop */
            } /* j loop */
        } /* i loop */
        
        dbForce[0] = -dbForce[0]*hx*hy*hzed*deps*0.5*izmagic;
        dbForce[1] = -dbForce[1]*hx*hy*hzed*deps*0.5*izmagic;
        dbForce[2] = -dbForce[2]*hx*hy*hzed*deps*0.5*izmagic;
        npForce[0] = -npForce[0]*hx*hy*hzed*gamma;
        npForce[1] = -npForce[1]*hx*hy*hzed*gamma;
        npForce[2] = -npForce[2]*hx*hy*hzed*gamma;
    }
}

VPUBLIC void Vpmg_qfForce(Vpmg *thee, double *force, int atomID, 
  Vchrg_Meth chgm) {

    double tforce[3];

    /* Reset force */
    force[0] = 0.0;
    force[1] = 0.0;
    force[2] = 0.0;

    /* Check surface definition */
    if (chgm != VCM_BSPL2) {
        Vnm_print(2, "Vpmg_qfForce:  It is recommended that forces be \
calculated with the\n");
        Vnm_print(2, "Vpmg_qfForce:  cubic spline charge discretization \
scheme\n");
	}

    switch (chgm) {
        case VCM_TRIL:
            qfForceSpline1(thee, tforce, atomID);
            break;
        case VCM_BSPL2:
            qfForceSpline2(thee, tforce, atomID);
            break;
        default:
            Vnm_print(2, "Vpmg_qfForce:  Undefined charge discretization \
method (%d)!\n", chgm);
            Vnm_print(2, "Vpmg_qfForce:  Forces not calculated!\n");
            return;
    }

    /* Assign forces */
    force[0] = tforce[0];
    force[1] = tforce[1];
    force[2] = tforce[2];

    return;
}


VPRIVATE void qfForceSpline1(Vpmg *thee, double *force, int atomID) {

    Vatom *atom;
    
    double *apos, position[3], hx, hy, hzed;
    double xmin, ymin, zmin, xmax, ymax, zmax;
    double dx, dy, dz;
    double *u, charge, ifloat, jfloat, kfloat;
    int nx, ny, nz, ihi, ilo, jhi, jlo, khi, klo;

    VASSERT(thee != VNULL);
    /* VASSERT(thee->filled); */

    atom = Valist_getAtom(thee->pbe->alist, atomID);
    apos = Vatom_getPosition(atom);
    charge = Vatom_getCharge(atom);

    /* Reset force */
    force[0] = 0.0;
    force[1] = 0.0;
    force[2] = 0.0;

    /* If we aren't in the current position, then we're done */
    if (atom->partID == 0) return;

    /* Mesh info */
    nx = thee->pmgp->nx;
    ny = thee->pmgp->ny;
    nz = thee->pmgp->nz;
    hx = thee->pmgp->hx;
    hy = thee->pmgp->hy;
    hzed = thee->pmgp->hzed;
    xmin = thee->pmgp->xmin;
    ymin = thee->pmgp->ymin;
    zmin = thee->pmgp->zmin;
    xmax = thee->pmgp->xmax;
    ymax = thee->pmgp->ymax;
    zmax = thee->pmgp->zmax;
    u = thee->u;

    /* Make sure we're on the grid */
    if ((apos[0]<=xmin) || (apos[0]>=xmax) || (apos[1]<=ymin) || \
        (apos[1]>=ymax) || (apos[2]<=zmin) || (apos[2]>=zmax)) {
        if (thee->pmgp->bcfl != BCFL_FOCUS) {
            Vnm_print(2, "Vpmg_qfForce:  Atom #%d at (%4.3f, %4.3f, %4.3f) is off the mesh (ignoring):\n", atomID, position[0], position[1], position[2]);
            Vnm_print(2, "Vpmg_qfForce:    xmin = %g, xmax = %g\n", xmin, xmax);
            Vnm_print(2, "Vpmg_qfForce:    ymin = %g, ymax = %g\n", ymin, ymax);
            Vnm_print(2, "Vpmg_qfForce:    zmin = %g, zmax = %g\n", zmin, zmax);
        }
        fflush(stderr);
    } else {
    
        /* Convert the atom position to grid coordinates */
        position[0] = apos[0] - xmin;
        position[1] = apos[1] - ymin;
        position[2] = apos[2] - zmin;
        ifloat = position[0]/hx;
        jfloat = position[1]/hy;
        kfloat = position[2]/hzed;
        ihi = (int)ceil(ifloat);
        ilo = (int)floor(ifloat);
        jhi = (int)ceil(jfloat);
        jlo = (int)floor(jfloat);
        khi = (int)ceil(kfloat);
        klo = (int)floor(kfloat);
        VASSERT((ihi < nx) && (ihi >=0));
        VASSERT((ilo < nx) && (ilo >=0));
        VASSERT((jhi < ny) && (jhi >=0));
        VASSERT((jlo < ny) && (jlo >=0));
        VASSERT((khi < nz) && (khi >=0));
        VASSERT((klo < nz) && (klo >=0));
        dx = ifloat - (double)(ilo);
        dy = jfloat - (double)(jlo);
        dz = kfloat - (double)(klo);


#if 0
        Vnm_print(1, "Vpmg_qfForce: (DEBUG) u ~ %g\n", 
          dx    *dy    *dz    *u[IJK(ihi,jhi,khi)]
         +dx    *dy    *(1-dz)*u[IJK(ihi,jhi,klo)]
         +dx    *(1-dy)*dz    *u[IJK(ihi,jlo,khi)]
         +dx    *(1-dy)*(1-dz)*u[IJK(ihi,jlo,klo)]
         +(1-dx)*dy    *dz    *u[IJK(ilo,jhi,khi)]
         +(1-dx)*dy    *(1-dz)*u[IJK(ilo,jhi,klo)]
         +(1-dx)*(1-dy)*dz    *u[IJK(ilo,jlo,khi)]
         +(1-dx)*(1-dy)*(1-dz)*u[IJK(ilo,jlo,klo)]);
#endif


        if ((dx > VPMGSMALL) && (VABS(1.0-dx) > VPMGSMALL)) {
            force[0] = 
              -charge*(dy    *dz    *u[IJK(ihi,jhi,khi)]
                     + dy    *(1-dz)*u[IJK(ihi,jhi,klo)]
                     + (1-dy)*dz    *u[IJK(ihi,jlo,khi)]
                     + (1-dy)*(1-dz)*u[IJK(ihi,jlo,klo)]
                     - dy    *dz    *u[IJK(ilo,jhi,khi)]
                     - dy    *(1-dz)*u[IJK(ilo,jhi,klo)]
                     - (1-dy)*dz    *u[IJK(ilo,jlo,khi)]
                     - (1-dy)*(1-dz)*u[IJK(ilo,jlo,klo)])/hx;
        } else {
            force[0] = 0;
            Vnm_print(0, 
              "Vpmg_qfForce:  Atom %d on x gridline; zero x-force\n", atomID);
        }
        if ((dy > VPMGSMALL) && (VABS(1.0-dy) > VPMGSMALL)) {
            force[1] = 
              -charge*(dx    *dz    *u[IJK(ihi,jhi,khi)]
                     + dx    *(1-dz)*u[IJK(ihi,jhi,klo)]
                     - dx    *dz    *u[IJK(ihi,jlo,khi)]
                     - dx    *(1-dz)*u[IJK(ihi,jlo,klo)]
                     + (1-dx)*dz    *u[IJK(ilo,jhi,khi)]
                     + (1-dx)*(1-dz)*u[IJK(ilo,jhi,klo)]
                     - (1-dx)*dz    *u[IJK(ilo,jlo,khi)]
                     - (1-dx)*(1-dz)*u[IJK(ilo,jlo,klo)])/hy;
        } else {
            force[1] = 0;
            Vnm_print(0, 
              "Vpmg_qfForce:  Atom %d on y gridline; zero y-force\n", atomID);
        }
        if ((dz > VPMGSMALL) && (VABS(1.0-dz) > VPMGSMALL)) {
            force[2] = 
              -charge*(dy    *dx    *u[IJK(ihi,jhi,khi)]
                     - dy    *dx    *u[IJK(ihi,jhi,klo)]
                     + (1-dy)*dx    *u[IJK(ihi,jlo,khi)]
                     - (1-dy)*dx    *u[IJK(ihi,jlo,klo)]
                     + dy    *(1-dx)*u[IJK(ilo,jhi,khi)]
                     - dy    *(1-dx)*u[IJK(ilo,jhi,klo)]
                     + (1-dy)*(1-dx)*u[IJK(ilo,jlo,khi)]
                     - (1-dy)*(1-dx)*u[IJK(ilo,jlo,klo)])/hzed;
        } else {
            force[2] = 0;
            Vnm_print(0, 
              "Vpmg_qfForce:  Atom %d on z gridline; zero z-force\n", atomID);
        }
    }
}

VPRIVATE void qfForceSpline2(Vpmg *thee, double *force, int atomID) {

    Vatom *atom;
    
    double *apos, position[3], hx, hy, hzed;
    double xlen, ylen, zlen, xmin, ymin, zmin, xmax, ymax, zmax;
    double mx, my, mz, dmx, dmy, dmz;
    double *u, charge, ifloat, jfloat, kfloat;
    int nx, ny, nz, im2, im1, ip1, ip2, jm2, jm1, jp1, jp2, km2, km1; 
    int kp1, kp2, ii, jj, kk;

    VASSERT(thee != VNULL);
    /* VASSERT(thee->filled); */

    atom = Valist_getAtom(thee->pbe->alist, atomID);
    apos = Vatom_getPosition(atom);
    charge = Vatom_getCharge(atom);

    /* Reset force */
    force[0] = 0.0;
    force[1] = 0.0;
    force[2] = 0.0;

    /* If we aren't in the current position, then we're done */
    if (atom->partID == 0) return;

    /* Mesh info */
    nx = thee->pmgp->nx;
    ny = thee->pmgp->ny;
    nz = thee->pmgp->nz;
    hx = thee->pmgp->hx;
    hy = thee->pmgp->hy;
    hzed = thee->pmgp->hzed;
    xlen = thee->pmgp->xlen;
    ylen = thee->pmgp->ylen;
    zlen = thee->pmgp->zlen;
    xmin = thee->pmgp->xmin;
    ymin = thee->pmgp->ymin;
    zmin = thee->pmgp->zmin;
    xmax = thee->pmgp->xmax;
    ymax = thee->pmgp->ymax;
    zmax = thee->pmgp->zmax;
    u = thee->u;

    /* Make sure we're on the grid */
    if ((apos[0]<=(xmin+hx))   || (apos[0]>=(xmax-hx)) \
     || (apos[1]<=(ymin+hy))   || (apos[1]>=(ymax-hy)) \
     || (apos[2]<=(zmin+hzed)) || (apos[2]>=(zmax-hzed))) {
        if (thee->pmgp->bcfl != BCFL_FOCUS) {
            Vnm_print(2, "qfForceSpline2:  Atom #%d off the mesh \
(ignoring)\n", atomID);
        }
        fflush(stderr);

    } else {
    
        /* Convert the atom position to grid coordinates */
        position[0] = apos[0] - xmin;
        position[1] = apos[1] - ymin;
        position[2] = apos[2] - zmin;
        ifloat = position[0]/hx;
        jfloat = position[1]/hy;
        kfloat = position[2]/hzed;
        ip1 = (int)ceil(ifloat);
        ip2 = ip1 + 1;
        im1 = (int)floor(ifloat);
        im2 = im1 - 1;
        jp1 = (int)ceil(jfloat);
        jp2 = jp1 + 1;
        jm1 = (int)floor(jfloat);
        jm2 = jm1 - 1;
        kp1 = (int)ceil(kfloat);
        kp2 = kp1 + 1;
        km1 = (int)floor(kfloat);
        km2 = km1 - 1;

        /* This step shouldn't be necessary, but it saves nasty debugging
         * later on if something goes wrong */
        ip2 = VMIN2(ip2,nx-1);
        ip1 = VMIN2(ip1,nx-1);
        im1 = VMAX2(im1,0);
        im2 = VMAX2(im2,0);
        jp2 = VMIN2(jp2,ny-1);
        jp1 = VMIN2(jp1,ny-1);
        jm1 = VMAX2(jm1,0);
        jm2 = VMAX2(jm2,0);
        kp2 = VMIN2(kp2,nz-1);
        kp1 = VMIN2(kp1,nz-1);
        km1 = VMAX2(km1,0);
        km2 = VMAX2(km2,0);


        for (ii=im2; ii<=ip2; ii++) {
            mx = bspline2(VFCHI(ii,ifloat));
            dmx = dbspline2(VFCHI(ii,ifloat));
            for (jj=jm2; jj<=jp2; jj++) {
                my = bspline2(VFCHI(jj,jfloat));
                dmy = dbspline2(VFCHI(jj,jfloat));
                for (kk=km2; kk<=kp2; kk++) {
                    mz = bspline2(VFCHI(kk,kfloat));
                    dmz = dbspline2(VFCHI(kk,kfloat));

                    force[0] += (charge*dmx*my*mz*u[IJK(ii,jj,kk)])/hx;
                    force[1] += (charge*mx*dmy*mz*u[IJK(ii,jj,kk)])/hy;
                    force[2] += (charge*mx*my*dmz*u[IJK(ii,jj,kk)])/hzed;

                }
            }
        }

    }
}

VPRIVATE void markFrac(
        double rtot, double *tpos,
        int nx, int ny, int nz,
        double hx, double hy, double hzed,
        double xmin, double ymin, double zmin,
        double *xarray, double *yarray, double *zarray) {

    int i, j, k, imin, imax, jmin, jmax, kmin, kmax;
    double dx, dx2, dy, dy2, dz, dz2, a000, a001, a010, a100, r2;
    double x, xp, xm, y, yp, ym, zp, z, zm, xspan, yspan, zspan;
    double rtot2, pos[3];

    /* Convert to grid reference frame */
    pos[0] = tpos[0] - xmin;
    pos[1] = tpos[1] - ymin;
    pos[2] = tpos[2] - zmin;

    rtot2 = VSQR(rtot);

    xspan = rtot + 2*hx;
    imin = VMAX2(0, (int)ceil((pos[0] - xspan)/hx));
    imax = VMIN2(nx-1, (int)floor((pos[0] + xspan)/hx));
    for (i=imin; i<=imax; i++) {
        x = hx*i;
        dx2 = VSQR(pos[0] - x);
        if (rtot2 > dx2) {
            yspan = VSQRT(rtot2 - dx2) + 2*hy;
        } else {
            yspan = 2*hy;
        }
        jmin = VMAX2(0,(int)ceil((pos[1] - yspan)/hy));
        jmax = VMIN2(ny-1,(int)floor((pos[1] + yspan)/hy));
        for (j=jmin; j<=jmax; j++) {
            y = hy*j;
            dy2 = VSQR(pos[1] - y);
            if (rtot2 > (dx2+dy2)) { 
                zspan = VSQRT(rtot2-dx2-dy2) + 2*hzed; 
            } else {
                zspan = 2*hzed;
            }
            kmin = VMAX2(0,(int)ceil((pos[2] - zspan)/hzed));
            kmax = VMIN2(nz-1,(int)floor((pos[2] + zspan)/hzed));
            for (k=kmin; k<=kmax; k++) {
                z = hzed*k;
                dz2 = VSQR(pos[2] - z);

                r2 = dx2 + dy2 + dz2;

                /* We need to determine the inclusion value a000 at (i,j,k) */
                if (r2 < rtot2) a000 = 1.0;
                else a000 = 0.0;

                /* We need to evaluate the values of x which intersect the
                 * sphere and determine if these are in the interval 
                 * [(i,j,k), (i+1,j,k)] */
                if (r2 < (rtot2 - hx*hx)) a100 = 1.0;
                else if (r2 > (rtot2 + hx*hx)) a100 = 0.0;
                else if (rtot2 > (dy2 + dz2)) {
                    dx = VSQRT(rtot2 - dy2 - dz2);
                    xm = pos[0] - dx;
                    xp = pos[0] + dx;
                    if ((xm < x+hx) && (xm > x)) {
                        a100 = (xm - x)/hx;
                    } else if ((xp < x+hx) && (xp > x)) {
                        a100 = (xp - x)/hx;
                    }
                } else a100 = 0.0;

                /* We need to evaluate the values of y which intersect the
                 * sphere and determine if these are in the interval 
                 * [(i,j,k), (i,j+1,k)] */
                if (r2 < (rtot2 - hy*hy)) a010 = 1.0;
                else if (r2 > (rtot2 + hy*hy)) a010 = 0.0;
                else if (rtot2 > (dx2 + dz2)) {
                    dy = VSQRT(rtot2 - dx2 - dz2);
                    ym = pos[1] - dy;
                    yp = pos[1] + dy;
                    if ((ym < y+hy) && (ym > y)) {
                        a010 = (ym - y)/hy;
                    } else if ((yp < y+hy) && (yp > y)) {
                        a010 = (yp - y)/hy;
                    }
                } else a010 = 0.0;

                /* We need to evaluate the values of y which intersect the
                 * sphere and determine if these are in the interval 
                 * [(i,j,k), (i,j,k+1)] */
                if (r2 < (rtot2 - hzed*hzed)) a001 = 1.0;
                else if (r2 > (rtot2 + hzed*hzed)) a001 = 0.0;
                else if (rtot2 > (dx2 + dy2)) {
                    dz = VSQRT(rtot2 - dx2 - dy2);
                    zm = pos[2] - dz;
                    zp = pos[2] + dz;
                    if ((zm < z+hzed) && (zm > z)) {
                        a001 = (zm - z)/hzed;
                    } else if ((zp < z+hzed) && (zp > z)) {
                        a001 = (zp - z)/hzed;
                    }
                } else a001 = 0.0;

                if (a100 < xarray[IJK(i,j,k)]) xarray[IJK(i,j,k)] = a100;
                if (a010 < yarray[IJK(i,j,k)]) yarray[IJK(i,j,k)] = a010;
                if (a001 < zarray[IJK(i,j,k)]) zarray[IJK(i,j,k)] = a001;

            } /* k loop */
        } /* j loop */
    } /* i loop */
}

VPRIVATE void markSphere(
        double rtot, double *tpos,
        int nx, int ny, int nz,
        double hx, double hy, double hzed,
        double xmin, double ymin, double zmin,
        double *array, double markVal) {

    int i, j, k, imin, imax, jmin, jmax, kmin, kmax;
    double dx, dx2, dy, dy2, dz, dz2;
    double rtot2, pos[3];

    /* Convert to grid reference frame */
    pos[0] = tpos[0] - xmin;
    pos[1] = tpos[1] - ymin;
    pos[2] = tpos[2] - zmin;

    rtot2 = VSQR(rtot);

    dx = rtot + 0.5*hx;
    imin = VMAX2(0,(int)ceil((pos[0] - dx)/hx));
    imax = VMIN2(nx-1,(int)floor((pos[0] + dx)/hx));
    for (i=imin; i<=imax; i++) {
        dx2 = VSQR(pos[0] - hx*i);
        if (rtot2 > dx2) {
            dy = VSQRT(rtot2 - dx2) + 0.5*hy;
        } else {
            dy = 0.5*hy;
        }
        jmin = VMAX2(0,(int)ceil((pos[1] - dy)/hy));
        jmax = VMIN2(ny-1,(int)floor((pos[1] + dy)/hy));
        for (j=jmin; j<=jmax; j++) {
            dy2 = VSQR(pos[1] - hy*j);
            if (rtot2 > (dx2+dy2)) { 
                dz = VSQRT(rtot2-dx2-dy2)+0.5*hzed; 
            } else {
                dz = 0.5*hzed;
            }
            kmin = VMAX2(0,(int)ceil((pos[2] - dz)/hzed));
            kmax = VMIN2(nz-1,(int)floor((pos[2] + dz)/hzed));
            for (k=kmin; k<=kmax; k++) {
                dz2 = VSQR(k*hzed - pos[2]);
                if ((dz2 + dy2 + dx2) <= rtot2) {
                    array[IJK(i,j,k)] = markVal;
                }
            } /* k loop */
        } /* j loop */
    } /* i loop */
}
