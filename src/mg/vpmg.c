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
 * Washington University in St. Louis
 *
 * Additional contributing authors listed in the code documentation.
 *
 * Copyright (c) 2002.  Washington University in St. Louis.
 * All Rights Reserved.
 *
 * Portions Copyright (c) 1999-2002.  The Regents of the University of
 * California.  
 * Portions Copyright (c) 1995.  Michael Holst.
 *
 * Permission to use, copy, modify, and distribute this software and its
 * documentation for educational, research, and not-for-profit purposes,
 * without fee and without a signed licensing agreement, is hereby granted,
 * provided that the above copyright notice, this paragraph and the
 * following two paragraphs appear in all copies, modifications, and
 * distributions.
 *
 * IN NO EVENT SHALL THE AUTHORS BE LIABLE TO ANY PARTY FOR DIRECT, INDIRECT,
 * SPECIAL, INCIDENTAL, OR CONSEQUENTIAL DAMAGES, INCLUDING LOST PROFITS,
 * ARISING OUT OF THE USE OF THIS SOFTWARE AND ITS DOCUMENTATION, EVEN IF THE
 * AUTHORS HAVE BEEN ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 *
 * THE AUTHORS SPECIFICALLY DISCLAIM ANY WARRANTIES, INCLUDING, BUT NOT LIMITED
 * TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
 * PURPOSE.  THE SOFTWARE AND ACCOMPANYING DOCUMENTATION, IF ANY, PROVIDED
 * HEREUNDER IS PROVIDED "AS IS".  THE AUTHORS HAVE NO OBLIGATION TO PROVIDE
 * MAINTENANCE, SUPPORT, UPDATES, ENHANCEMENTS, OR MODIFICATIONS.

 * @endverbatim
 */

#include "apbscfg.h"
#include "vpmg-private.h"
#include "apbs/vpmg.h"

VEMBED(rcsid="$Id$")

/* ///////////////////////////////////////////////////////////////////////////
// Class Vpmg: Inlineable methods
/////////////////////////////////////////////////////////////////////////// */
#if !defined(VINLINE_VPMG)

/* ///////////////////////////////////////////////////////////////////////////
// Routine:  Vpmg_memChk
//
// Author:   Nathan Baker
/////////////////////////////////////////////////////////////////////////// */
VPUBLIC int Vpmg_memChk(Vpmg *thee) {
    if (thee == VNULL) return 0;
    return Vmem_bytes(thee->vmem);
}

#endif /* if !defined(VINLINE_VPMG) */


/* ///////////////////////////////////////////////////////////////////////////
// Routine:  Vpmg_printColComp
//
// Author:   Nathan Baker
/////////////////////////////////////////////////////////////////////////// */
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

/* ///////////////////////////////////////////////////////////////////////////
// Routine:  Vpmg_ctor
// Author:   Nathan Baker
/////////////////////////////////////////////////////////////////////////// */
VPUBLIC Vpmg* Vpmg_ctor(Vpmgp *pmgp, Vpbe *pbe) {

    Vpmg *thee = VNULL;

    /* Set up the structure */
    thee = Vmem_malloc(VNULL, 1, sizeof(Vpmg) );
    VASSERT( thee != VNULL);
    VASSERT(Vpmg_ctor2(thee, pmgp, pbe));

    return thee;
}

/* ///////////////////////////////////////////////////////////////////////////
// Routine:  Vpmg_ctorFocus
// Author:   Nathan Baker
/////////////////////////////////////////////////////////////////////////// */
VPUBLIC Vpmg* Vpmg_ctorFocus(Vpmgp *pmgp, Vpbe *pbe, Vpmg *pmgOLD, 
  int energyFlag) {

    Vpmg *thee = VNULL;

    /* Set up the structure */
    thee = Vmem_malloc(VNULL, 1, sizeof(Vpmg) );
    VASSERT( thee != VNULL);
    VASSERT(Vpmg_ctor2Focus(thee, pmgp, pbe, pmgOLD, energyFlag));

    return thee;
}


/* ///////////////////////////////////////////////////////////////////////////
// Routine:  Vpmg_solve
// Author:   Nathan Baker
/////////////////////////////////////////////////////////////////////////// */
VPUBLIC void Vpmg_solve(Vpmg *thee) {

    thee->filled = 0;

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

}

    
/* ///////////////////////////////////////////////////////////////////////////
// Routine:  Vpmg_dtor
// Author:   Nathan Baker
/////////////////////////////////////////////////////////////////////////// */
VPUBLIC void Vpmg_dtor(Vpmg **thee) {
    
    if ((*thee) != VNULL) {
        Vpmg_dtor2(*thee);
        Vmem_free(VNULL, 1, sizeof(Vpmg), (void **)thee);
        (*thee) = VNULL;
    }

}

/* ///////////////////////////////////////////////////////////////////////////
// Routine:  Vpmg_dtor2
// Author:   Nathan Baker
/////////////////////////////////////////////////////////////////////////// */
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
      sizeof(int), (void **)&(thee->pvec));

    Vmem_dtor(&(thee->vmem));
}

/* ///////////////////////////////////////////////////////////////////////////
// Routine:  Vpmg_setPart
//
// Author:   Nathan Baker
/////////////////////////////////////////////////////////////////////////// */
VPUBLIC void Vpmg_setPart(Vpmg *thee, double lowerCorner[3],
  double upperCorner[3], int bflags[6]) {

    Valist *alist;
    Vatom *atom;
    int i, j, k, nx, ny, nz, xok, yok, zok;
    double xmin, ymin, zmin, x, y, z, hx, hy, hzed;

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

    /* Identify atoms as inside or outside */
    for (i=0; i<Valist_getNumberAtoms(alist); i++) {
        atom = Valist_getAtom(alist, i);
        if ((atom->position[0] < upperCorner[0]) &&
            (atom->position[0] > lowerCorner[0])) xok = 1;
        else {
            if ((atom->position[0] == lowerCorner[0]) &&
                (bflags[VAPBS_LEFT] == 1)) xok = 1;
            else if ((atom->position[0] == upperCorner[0]) &&
                     (bflags[VAPBS_RIGHT] == 1)) xok = 1;
            else xok = 0; 
        }
        if ((atom->position[1] < upperCorner[1]) &&
            (atom->position[1] > lowerCorner[1])) yok = 1;
        else {
            if ((atom->position[1] == lowerCorner[1]) &&
                (bflags[VAPBS_BACK] == 1)) yok = 1;
            else if ((atom->position[1] == upperCorner[1]) &&
                     (bflags[VAPBS_FRONT] == 1)) yok = 1;
            else yok = 0;
        }
        if ((atom->position[2] < upperCorner[2]) &&
            (atom->position[2] > lowerCorner[2])) zok = 1;
        else {
            if ((atom->position[2] == lowerCorner[2]) &&
                (bflags[VAPBS_DOWN] == 1)) zok = 1;
            else if ((atom->position[2] == upperCorner[2]) &&
                     (bflags[VAPBS_UP] == 1)) zok = 1;
            else zok = 0; 
        }

        if ((xok && yok) && zok) atom->partID = 1;
        else atom->partID = 0;
    }


    /* Load up pvec */
    for (i=0; i<(nx*ny*nz); i++) thee->pvec[i] = 0;
    for (i=0; i<nx; i++) {
        xok = 0;
        x = i*hx + xmin;
        if ((x < upperCorner[0]) && (x > lowerCorner[0])) xok = 1;
        else { 
            if ((x == lowerCorner[0]) && 
                (bflags[VAPBS_LEFT] == 1)) xok = 1;
            else if ((x == upperCorner[0]) &&
                (bflags[VAPBS_RIGHT] == 1)) xok = 1;
            else xok = 0;
        }
        if (xok) {
            for (j=0; j<ny; j++) {
                yok = 0;
                y = j*hy + ymin;
                if ((y < upperCorner[1]) && 
                    (y > lowerCorner[1])) yok = 1;
                else {
                    if ((y == lowerCorner[1]) &&
                        (bflags[VAPBS_BACK] == 1)) yok = 1;
                    else if ((y == upperCorner[1]) &&
                        (bflags[VAPBS_FRONT] == 1)) yok = 1;
                    else yok = 0;
                }
                if (yok) {
                    for (k=0; k<nz; k++) {
                        zok = 0; 
                        z = k*hzed + zmin;
                        if ((z < upperCorner[2]) && 
                            (z > lowerCorner[2])) zok = 1;
                        else {
                            if ((z == lowerCorner[2]) &&
                                (bflags[VAPBS_DOWN] == 1)) zok = 1;
                            else if ((z == upperCorner[2]) &&
                                (bflags[VAPBS_UP] == 1)) zok = 1;
                            else zok = 0;

                        }
                        if (zok) thee->pvec[IJK(i,j,k)] = 1;
                        else thee->pvec[IJK(i,j,k)] = 0;
                    }
                }
            }
        }
    }
}

/* ///////////////////////////////////////////////////////////////////////////
// Routine:  Vpmg_unsetPart
//
// Author:   Nathan Baker
/////////////////////////////////////////////////////////////////////////// */
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

/* ///////////////////////////////////////////////////////////////////////////
// Routine:  Vpmg_fillArray
//
// Author:   Nathan Baker
/////////////////////////////////////////////////////////////////////////// */
VPUBLIC void Vpmg_fillArray(Vpmg *thee, double *vec, Vdata_Type type, 
  double parm) {

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
                        if (Vacc_molAcc(acc, position, pbe->maxIonRadius)) {
                            for (l=0; l<pbe->numIon; l++) {
                              vec[IJK(i,j,k)] += (pbe->ionConc[l]
                                * Vcap_exp(-pbe->ionQ[l]*thee->u[IJK(i,j,k)], 
                                &ichop));
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
                        if (Vacc_molAcc(acc, position, pbe->maxIonRadius)) {
                            for (l=0; l<pbe->numIon; l++) {
                              vec[IJK(i,j,k)] += (pbe->ionConc[l] 
                                * pbe->ionQ[l]
                                * Vcap_exp(-pbe->ionQ[l]*thee->u[IJK(i,j,k)],
                                &ichop));
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
