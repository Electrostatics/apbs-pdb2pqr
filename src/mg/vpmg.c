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
 * Nathan A. Baker (nbaker@wasabi.ucsd.edu)
 * Dept. of Chemistry and Biochemistry
 * University of California, San Diego 
 *
 * Additional contributing authors listed in the code documentation.
 *
 * Copyright (c) 1999-2002. The Regents of the University of California
 *                          (Regents).  All Rights Reserved.
 *
 * Permission to use, copy, modify, and distribute this software and its
 * documentation for educational, research, and not-for-profit purposes,
 * without fee and without a signed licensing agreement, is hereby granted,
 * provided that the above copyright notice, this paragraph and the
 * following two paragraphs appear in all copies, modifications, and
 * distributions.
 *
 * IN NO EVENT SHALL REGENTS BE LIABLE TO ANY PARTY FOR DIRECT, INDIRECT,
 * SPECIAL, INCIDENTAL, OR CONSEQUENTIAL DAMAGES, INCLUDING LOST PROFITS,
 * ARISING OUT OF THE USE OF THIS SOFTWARE AND ITS DOCUMENTATION, EVEN IF
 * REGENTS HAS BEEN ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 *
 * REGENTS SPECIFICALLY DISCLAIMS ANY WARRANTIES, INCLUDING, BUT NOT
 * LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A
 * PARTICULAR PURPOSE.  THE SOFTWARE AND ACCOMPANYING DOCUMENTATION, IF
 * ANY, PROVIDED HEREUNDER IS PROVIDED "AS IS".  REGENTS HAS NO OBLIGATION
 * TO PROVIDE MAINTENANCE, SUPPORT, UPDATES, ENHANCEMENTS, OR
 * MODIFICATIONS.
 *
 * @endverbatim
 */

#include "apbscfg.h"
#include "vpmg-private.h"
#include "apbs/vpmg.h"

#define VPMGSMALL 1e-14

VEMBED(rcsid="$Id$")

/* ///////////////////////////////////////////////////////////////////////////
// Class Vpmg: Inlineable methods
/////////////////////////////////////////////////////////////////////////// */
#if !defined(VINLINE_VACC)
#endif /* if !defined(VINLINE_VACC) */


/* ///////////////////////////////////////////////////////////////////////////
// Routine:  Vpmg_printColComp
//
// Author:   Nathan Baker
/////////////////////////////////////////////////////////////////////////// */
VPUBLIC void Vpmg_printColComp(Vpmg *thee, char path[72], char title[72], 
  char mxtype[3], int flag) {

    int i, nn, nxm2, nym2, nzm2, ncol, nrow, nonz; 
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

    Vnm_print(1, "Vpmg_printColComp:  Allocated space for %d nonzeros\n",
      nonz);

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
// Routine:  focusFillBound
//
// Purpose:  Fill boundaries with old values before destroying
//
// Author:   Nathan Baker
/////////////////////////////////////////////////////////////////////////// */
VPRIVATE void focusFillBound(Vpmg *thee, Vpmg *pmgOLD) {

    Vpbe *pbe;
    Valist *alist;
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
    alist = thee->pbe->alist;
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

/* ///////////////////////////////////////////////////////////////////////////
// Routine:  extEnergy
//
// Purpose:  Calculate energy from region outside of current (focused) domain
//
// Arguments:  extFlag (1 => calculate total energy only, 2 => calculate energy
//             components)
//
// Author:   Nathan Baker
/////////////////////////////////////////////////////////////////////////// */
VPRIVATE void extEnergy(Vpmg *thee, Vpmg *pmgOLD, int extFlag) {

    Vatom *atom;
    double hxNEW, hyNEW, hzNEW;
    double lowerCorner[3], upperCorner[3];
    int nxNEW, nyNEW, nzNEW;
    int nxOLD, nyOLD, nzOLD, bflags[6];
    int i;

    /* Set the new external energy contribution to zero.  Any external
     * contributions from higher levels will be included in the appropriate
     * energy function call. */
    thee->extQmEnergy = 0;
    thee->extQfEnergy = 0;
    thee->extDiEnergy = 0;

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

    /* Old problem dimensions */
    nxOLD = pmgOLD->pmgp->nx;
    nyOLD = pmgOLD->pmgp->ny;
    nzOLD = pmgOLD->pmgp->nz;

    /* Create a partition based on the new problem dimensions */
    for (i=0; i<6; i++) bflags[i] = 1;
    Vpmg_setPart(pmgOLD, lowerCorner, upperCorner, bflags);
    /* Invert partition mask */
    for (i=0; i<(nxOLD*nyOLD*nzOLD); i++) {
        pmgOLD->pvec[i] = (!(pmgOLD->pvec[i]));
    }
    for (i=0; i<Valist_getNumberAtoms(thee->pbe->alist); i++) {
        atom = Valist_getAtom(thee->pbe->alist, i);
        atom->partID = !(atom->partID);
    }
    /* Now calculate the energy on inverted subset of the domain */
    thee->extQmEnergy = Vpmg_qmEnergy(pmgOLD, 1);
    Vnm_print(0, "VPMG::extEnergy: extQmEnergy = %g kT\n", thee->extQmEnergy);
    thee->extQfEnergy = Vpmg_qfEnergy(pmgOLD, 1);
    Vnm_print(0, "VPMG::extEnergy: extQfEnergy = %g kT\n", thee->extQfEnergy);
    thee->extDiEnergy = Vpmg_dielEnergy(pmgOLD, 1);
    Vnm_print(0, "VPMG::extEnergy: extDiEnergy = %g kT\n", thee->extDiEnergy);
    Vpmg_unsetPart(pmgOLD);
}

/* ///////////////////////////////////////////////////////////////////////////
// Routine:  bcfl1sp
//
// Purpose:  Return the value of
//              pre1*(charge/d)*(exp(-xkappa*(d-size))/(1+xkappa*size)
//
// Args:     apos and pos are 3-vectors
//
// Author:   Nathan Baker and Michael Holst
/////////////////////////////////////////////////////////////////////////// */
VPRIVATE double bcfl1sp(double size, double *apos, double charge, 
  double xkappa, double pre1, double *pos) {

    double dist, val;

    dist = VSQRT(VSQR(pos[0]-apos[0]) + VSQR(pos[1]-apos[1])
      + VSQR(pos[2]-apos[2]));
    if (xkappa != 0.0) {
        val = pre1*(charge/dist)*VEXP(-xkappa*(dist-size))
          / (1+xkappa*size);
    } else {
        val = pre1*(charge/dist);
    } 

    return val;
}

/* ///////////////////////////////////////////////////////////////////////////
// Routine:  bcfl1
//
// Purpose:  Increment all the boundary points by 
//              pre1*(charge/d)*(exp(-xkappa*(d-size))/(1+xkappa*size)
//
// Args:     apos is a 3-vector
//
// Author:   Nathan Baker and Michael Holst
/////////////////////////////////////////////////////////////////////////// */
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
            if (xkappa != 0.0) {
                val = pre1*(charge/dist)*VEXP(-xkappa*(dist-size))
                       / (1+xkappa*size);
            } else {
                val = pre1*(charge/dist);
            } 
            gxcf[IJKx(j,k,0)] += val;
            gpos[0] = xf[nx-1];
            dist = VSQRT(VSQR(gpos[0]-apos[0]) + VSQR(gpos[1]-apos[1])
              + VSQR(gpos[2]-apos[2]));
            if (xkappa != 0.0) {
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
            if (xkappa != 0.0) {
                val = pre1*(charge/dist)*VEXP(-xkappa*(dist-size))
                       / (1+xkappa*size);
            } else {
                val = pre1*(charge/dist);
            }
            gycf[IJKy(i,k,0)] += val;
            gpos[1] = yf[ny-1];
            dist = VSQRT(VSQR(gpos[0]-apos[0]) + VSQR(gpos[1]-apos[1])
              + VSQR(gpos[2]-apos[2]));
            if (xkappa != 0.0) {
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
            if (xkappa != 0.0) {
                val = pre1*(charge/dist)*VEXP(-xkappa*(dist-size))
                       / (1+xkappa*size);
            } else {
                val = pre1*(charge/dist);
            }
            gzcf[IJKz(i,j,0)] += val;
            gpos[2] = zf[nz-1];
            gzcf[IJKz(i,j,1)] += val;
        }
    }
}


/* ///////////////////////////////////////////////////////////////////////////
// Routine:  bcCalc
//
// Purpose:  Dirichlet boundary function and initial approximation function.
//
// Args:     x    = position vector
//           flag = evaluation flag 
//                    0 => zero B.C.
//                    1 => single Debye-Huckel sphere
//                    2 => Debye-Huckel sphere for each atom
//                    4 => focusing
//
// Author:   Nathan Baker and Michael Holst
/////////////////////////////////////////////////////////////////////////// */
VPRIVATE void bcCalc(Vpmg *thee) {

    int flag, nx, ny, nz;
    double size, *position, charge, xkappa, eps_w, T, pre1;
    int i, j, k, iatom;
    Vpbe *pbe;
    Vatom *atom;
    Valist *alist;
    
    pbe = thee->pbe;
    flag = thee->pmgp->bcfl;
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
   
    /*  If we have zero boundary conditions, we're done */
    if (flag == 0) return;

    /*  For single DH sphere BC's, we only have one "atom" to deal with; get
     *  its information and */
    else if (flag == 1) {

        size = Vpbe_getSoluteRadius(pbe);
        position = Vpbe_getSoluteCenter(pbe);
        charge = Vunit_ec*Vpbe_getSoluteCharge(pbe);

        bcfl1(size, position, charge, xkappa, pre1,
          thee->gxcf, thee->gycf, thee->gzcf, 
          thee->xf, thee->yf, thee->zf, nx, ny, nz);

    } else if (flag == 2) {
        for (iatom=0; iatom<Valist_getNumberAtoms(alist); iatom++) {
            atom = Valist_getAtom(alist, iatom);
            position = Vatom_getPosition(atom);
            charge = Vunit_ec*Vatom_getCharge(atom);
            size = Vatom_getRadius(atom);
            bcfl1(size, position, charge, xkappa, pre1,
              thee->gxcf, thee->gycf, thee->gzcf, 
              thee->xf, thee->yf, thee->zf, nx, ny, nz);
        }
    } else if (flag == 4) {
        Vnm_print(2, "VPMG::bcCalc -- not appropriate for focusing!\n");
        VASSERT(0);
    } else {
        Vnm_print(2, "VPMG::bcCalc -- invalid boundary condition flag (%d)!\n",
          flag);
        VASSERT(0);
    }
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
// Routine:  Vpmg_ctor2
// Author:   Nathan Baker
/////////////////////////////////////////////////////////////////////////// */
VPUBLIC int Vpmg_ctor2(Vpmg *thee, Vpmgp *pmgp, Vpbe *pbe) {

    int i, nion;
    double ionConc[MAXION], ionQ[MAXION], ionRadii[MAXION], zkappa2, zks2;
    double ionstr;

    /* Get the parameters */    
    VASSERT(pmgp != VNULL); 
    VASSERT(pbe != VNULL); 
    thee->pmgp = pmgp;
    thee->pbe = pbe;

    /* Set up the memory */
    thee->vmem = Vmem_ctor("APBS:VPMG");

    /* Calculate storage requirements */
    F77MGSZ(&(thee->pmgp->mgcoar), &(thee->pmgp->mgdisc), &(thee->pmgp->mgsolv),
      &(thee->pmgp->nx), &(thee->pmgp->ny), &(thee->pmgp->nz),
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

    /* Allocate storage */
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
    thee->gxcf = (double *)Vmem_malloc(thee->vmem, 
      10*(thee->pmgp->ny)*(thee->pmgp->nz), sizeof(double));
    thee->gycf = (double *)Vmem_malloc(thee->vmem, 
      10*(thee->pmgp->nx)*(thee->pmgp->nz), sizeof(double));
    thee->gzcf = (double *)Vmem_malloc(thee->vmem, 
      10*(thee->pmgp->nx)*(thee->pmgp->ny), sizeof(double));
    thee->pvec = (int *)Vmem_malloc(thee->vmem, 
      (thee->pmgp->nx)*(thee->pmgp->ny)*(thee->pmgp->nz), sizeof(int));

    /* Plop some of the parameters into the iparm and rparm arrays */
    F77PACKMG(thee->iparm, thee->rparm, &(thee->pmgp->nrwk), &(thee->pmgp->niwk),
      &(thee->pmgp->nx), &(thee->pmgp->ny), &(thee->pmgp->nz),
      &(thee->pmgp->nlev), &(thee->pmgp->nu1), &(thee->pmgp->nu2),
      &(thee->pmgp->mgkey), &(thee->pmgp->itmax), &(thee->pmgp->istop),
      &(thee->pmgp->ipcon), &(thee->pmgp->nonlin), &(thee->pmgp->mgsmoo),
      &(thee->pmgp->mgprol), &(thee->pmgp->mgcoar), &(thee->pmgp->mgsolv),
      &(thee->pmgp->mgdisc), &(thee->pmgp->iinfo), &(thee->pmgp->errtol),
      &(thee->pmgp->ipkey), &(thee->pmgp->omegal), &(thee->pmgp->omegan),
      &(thee->pmgp->irite), &(thee->pmgp->iperf));

    /* Turn off restriction of observable calculations to a specific 
     * partition */
    Vpmg_unsetPart(thee);

    /* Initialize ion concentrations and valencies in PMG routines */
    zkappa2 = Vpbe_getZkappa2(thee->pbe);
    ionstr = Vpbe_getBulkIonicStrength(thee->pbe);
    if (ionstr > 0.0) zks2 = 0.5*zkappa2/ionstr;
    else zks2 = 0.0;
    Vpbe_getIons(thee->pbe, &nion, ionConc, ionRadii, ionQ);
    for (i=0; i<nion; i++) {
        ionConc[i] = zks2 * ionConc[i] * ionQ[i];
    }
    F77MYPDEFINIT(&nion, ionQ, ionConc);

    /* Ignore external energy contributions */
    thee->extQmEnergy = 0;
    thee->extDiEnergy = 0;
    thee->extQfEnergy = 0;

    thee->filled = 0;

    return 1;
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
// Routine:  Vpmg_ctor2Focus
// Author:   Nathan Baker
/////////////////////////////////////////////////////////////////////////// */
VPUBLIC int Vpmg_ctor2Focus(Vpmg *thee, Vpmgp *pmgp, Vpbe *pbe, Vpmg *pmgOLD,
  int energyFlag) {

    int i, nion;
    double ionstr, zkappa2, zks2, ionQ[MAXION], ionConc[MAXION];
    double ionRadii[MAXION];

    /* Get the parameters */    
    VASSERT(pmgp != VNULL); 
    VASSERT(pbe != VNULL); 
    VASSERT(pmgOLD != VNULL); 
    thee->pmgp = pmgp;
    thee->pbe = pbe;

    /* Set up the memory */
    thee->vmem = Vmem_ctor("APBS:VPMG");

    /* Calculate storage requirements */
    F77MGSZ(&(thee->pmgp->mgcoar), &(thee->pmgp->mgdisc),
      &(thee->pmgp->mgsolv), &(thee->pmgp->nx), &(thee->pmgp->ny),
      &(thee->pmgp->nz), &(thee->pmgp->nlev), &(thee->pmgp->nxc),
      &(thee->pmgp->nyc), &(thee->pmgp->nzc), &(thee->pmgp->nf),
      &(thee->pmgp->nc), &(thee->pmgp->narr), &(thee->pmgp->narrc),
      &(thee->pmgp->n_rpc), &(thee->pmgp->n_iz), &(thee->pmgp->n_ipc),
      &(thee->pmgp->nrwk), &(thee->pmgp->niwk));


    /* We need some additional storage if: nonlinear & newton OR cgmg */
    if (((thee->pmgp->nonlin == 1) && (thee->pmgp->meth == 1))
        || (thee->pmgp->meth == 0)) { thee->pmgp->nrwk += (2*(thee->pmgp->nf));
    }


    /* Overwrite any default or user-specified boundary condition arguments; we
     * are now committed to a calculation via focusing */
    if (thee->pmgp->bcfl != 4) {
        Vnm_print(2, "Vpmg_ctor2Focus: reset boundary condition flag to 4!\n");
        thee->pmgp->bcfl = 4;
    }

    /* Allocate storage for boundaries */
    thee->gxcf = (double *)Vmem_malloc(thee->vmem, 
      10*(thee->pmgp->ny)*(thee->pmgp->nz), sizeof(double));
    thee->gycf = (double *)Vmem_malloc(thee->vmem, 
      10*(thee->pmgp->nx)*(thee->pmgp->nz), sizeof(double));
    thee->gzcf = (double *)Vmem_malloc(thee->vmem, 
      10*(thee->pmgp->nx)*(thee->pmgp->ny), sizeof(double));

    /* Fill boundaries */
    Vnm_print(0, "Vpmg_ctor2Focus:  Filling boundary with old solution!\n");
    focusFillBound(thee, pmgOLD);

    /* Calculate energetic contributions from region outside focusing domain */
    if (energyFlag != 0) extEnergy(thee, pmgOLD, energyFlag);

    /* Destroy old Vpmg object */
    Vpmg_dtor(&pmgOLD);

    /* Allocate storage for everything else */
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
    thee->pvec = (int *)Vmem_malloc(thee->vmem, 
      (thee->pmgp->nz)*(thee->pmgp->nx)*(thee->pmgp->ny), sizeof(int));

    /* Plop some of the parameters into the iparm and rparm arrays */
    F77PACKMG(thee->iparm, thee->rparm, &(thee->pmgp->nrwk), &(thee->pmgp->niwk),
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
    if (ionstr > 0.0) zks2 = 0.5*zkappa2/ionstr;
    else zks2 = 0.0;
    Vpbe_getIons(thee->pbe, &nion, ionConc, ionRadii, ionQ);
    for (i=0; i<nion; i++) {
        ionConc[i] = zks2 * ionConc[i] * ionQ[i];
    }
    F77MYPDEFINIT(&nion, ionQ, ionConc);

    /* Turn off restriction of observable calculations to a specific 
     * partition */
    Vpmg_unsetPart(thee);

    thee->filled = 0;

    return 1;
}

/* ///////////////////////////////////////////////////////////////////////////
// Routine:  Vpmg_solve
// Author:   Nathan Baker
/////////////////////////////////////////////////////////////////////////// */
VPUBLIC void Vpmg_solve(Vpmg *thee) {

    switch(thee->pmgp->meth) {
        /* CGMG (linear) */
        case 0:
            F77CGMGDRIV(thee->iparm, thee->rparm, thee->iwork, thee->rwork,
              thee->u, thee->xf, thee->yf, thee->zf, thee->gxcf, thee->gycf,
              thee->gzcf, thee->a1cf, thee->a2cf, thee->a3cf, thee->ccf,
              thee->fcf, thee->tcf);
            break;
        /* Newton (nonlinear) */
        case 1:
            F77NEWDRIV(thee->iparm, thee->rparm, thee->iwork, thee->rwork, 
              thee->u, thee->xf, thee->yf, thee->zf, thee->gxcf, thee->gycf,
              thee->gzcf, thee->a1cf, thee->a2cf, thee->a3cf, thee->ccf, 
              thee->fcf, thee->tcf);
            break;
        /* MG (linear/nonlinear) */
        case 2:
	    F77MGDRIV(thee->iparm, thee->rparm, thee->iwork, thee->rwork,
	      thee->u, thee->xf, thee->yf, thee->zf, thee->gxcf, thee->gycf,
	      thee->gzcf, thee->a1cf, thee->a2cf, thee->a3cf, thee->ccf,
              thee->fcf, thee->tcf);
            break;
        /* CGHS (linear/nonlinear) */
        case 3: 
	    F77NCGHSDRIV(thee->iparm, thee->rparm, thee->iwork, thee->rwork,
	      thee->u, thee->xf, thee->yf, thee->zf, thee->gxcf, thee->gycf,
	      thee->gzcf, thee->a1cf, thee->a2cf, thee->a3cf, thee->ccf,
              thee->fcf, thee->tcf);
            break;
        /* SOR (linear/nonlinear) */
        case 4:
	    F77NSORDRIV(thee->iparm, thee->rparm, thee->iwork, thee->rwork,
	      thee->u, thee->xf, thee->yf, thee->zf, thee->gxcf, thee->gycf,
	      thee->gzcf, thee->a1cf, thee->a2cf, thee->a3cf, thee->ccf,
              thee->fcf, thee->tcf);
            break;
        /* GSRB (linear/nonlinear) */
        case 5:
	    F77NGSRBDRIV(thee->iparm, thee->rparm, thee->iwork, thee->rwork,
	      thee->u, thee->xf, thee->yf, thee->zf, thee->gxcf, thee->gycf,
	      thee->gzcf, thee->a1cf, thee->a2cf, thee->a3cf, thee->ccf,
              thee->fcf, thee->tcf); 
            break;
        /* WJAC (linear/nonlinear) */
        case 6:
	    F77NWJACDRIV(thee->iparm, thee->rparm, thee->iwork, thee->rwork,
	      thee->u, thee->xf, thee->yf, thee->zf, thee->gxcf, thee->gycf,
	      thee->gzcf, thee->a1cf, thee->a2cf, thee->a3cf, thee->ccf,
              thee->fcf, thee->tcf);
            break;
        /* RICH (linear/nonlinear) */
        case 7:
	    F77NRICHDRIV(thee->iparm, thee->rparm, thee->iwork, thee->rwork,
	      thee->u, thee->xf, thee->yf, thee->zf, thee->gxcf, thee->gycf,
	      thee->gzcf, thee->a1cf, thee->a2cf, thee->a3cf, thee->ccf,
              thee->fcf, thee->tcf);
            break;
        /* Error handling */
        default: 
            Vnm_print(2, "Vpgm_solve: invalid solver method key (%d)\n",
              thee->pmgp->key);
            break;
    }

}

/* ///////////////////////////////////////////////////////////////////////////
// Routine:  Vpmg_fillco
// Author:   Nathan Baker
/////////////////////////////////////////////////////////////////////////// */
VPUBLIC void Vpmg_fillco(Vpmg *thee, int surfMeth, double splineWin) {

    Vacc *acc;
    Valist *alist;
    Vpbe *pbe;
    Vatom *atom;
    double xmin, xmax, ymin, ymax, zmin, zmax, chi, ionmask, ionstr;
    double xlen, ylen, zlen, position[3], ifloat, jfloat, kfloat, accf;
    double zmagic, irad, srad, charge, dx, dy, dz, zkappa2, epsw, epsp;
    double hx, hy, hzed, *apos, arad, gpos[3];
    int i, j, k, nx, ny, nz, iatom, ihi, ilo, jhi, jlo, khi, klo;
    int imin, imax, jmin, jmax, kmin, kmax;
    double dx2, dy2, dz2, arad2, stot2, itot2, rtot, rtot2;
    int acclo, accmid, acchi, a000, islap;

    VASSERT(thee != VNULL);
    thee->surfMeth = surfMeth;
    thee->splineWin = splineWin;

    /* Get PBE info */
    pbe = thee->pbe;
    acc = pbe->acc;
    alist = pbe->alist;
    irad = Vpbe_getMaxIonRadius(pbe);
    srad = Vpbe_getSolventRadius(pbe);
    zmagic = Vpbe_getZmagic(pbe);
    zkappa2 = Vpbe_getZkappa2(pbe);
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

    /* This is a floating point parameter related to the non-zero nature of the
     * bulk ionic strength.  If the ionic strength is greater than zero; this
     * parameter is set to 1.0 and later scaled by the appropriate pre-factors.
     * Otherwise, this parameter is set to 0.0 */
    if (ionstr > VPMGSMALL) ionmask = 1.0;
    else ionmask = 0.0;

    /* This is a flag that gets set if the operator is a simple Laplacian;
     * i.e., in the case of a homogenous dielectric and zero ionic strength */
    if ((ionmask == 0.0) && (VABS(epsp-epsw) < VPMGSMALL)) islap = 1;
    else islap = 0;

    /* Fill the mesh point coordinate arrays */
    for (i=0; i<nx; i++) thee->xf[i] = xmin + i*hx;
    for (i=0; i<ny; i++) thee->yf[i] = ymin + i*hy;
    for (i=0; i<nz; i++) thee->zf[i] = zmin + i*hzed;

    /* Reset the fcf, tcf, ccf, a1cf, a2cf, and a3cf arrays */
    for (i=0; i<(nx*ny*nz); i++) {
        thee->tcf[i] = 0.0;
        thee->fcf[i] = 0.0;
        thee->ccf[i] = 0.0;
        thee->a1cf[i] = 0.0;
        thee->a2cf[i] = 0.0;
        thee->a3cf[i] = 0.0;
    }

    /* Fill in the source term (atomic charges) */
    Vnm_print(0, "Vpmg_fillco:  filling in source term.\n");
    for (iatom=0; iatom<Valist_getNumberAtoms(alist); iatom++) {

        atom = Valist_getAtom(alist, iatom);
        apos = Vatom_getPosition(atom);
        arad = Vatom_getRadius(atom);
        charge = Vatom_getCharge(atom);

        /* Make sure we're on the grid */
        if ((apos[0]<=xmin) || (apos[0]>=xmax)  || \
            (apos[1]<=ymin) || (apos[1]>=ymax)  || \
            (apos[2]<=zmin) || (apos[2]>=zmax)) {
            if (thee->pmgp->bcfl != 4) {
                Vnm_print(2, "Vpmg_fillco:  Atom #%d at (%4.3f, %4.3f, %4.3f)\
 is off the mesh (ignoring):\n",
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
            thee->fcf[IJK(ihi,jhi,khi)] += (dx*dy*dz*charge);
            thee->fcf[IJK(ihi,jlo,khi)] += (dx*(1.0-dy)*dz*charge);
            thee->fcf[IJK(ihi,jhi,klo)] += (dx*dy*(1.0-dz)*charge);
            thee->fcf[IJK(ihi,jlo,klo)] += (dx*(1.0-dy)*(1.0-dz)*charge);
            thee->fcf[IJK(ilo,jhi,khi)] += ((1.0-dx)*dy*dz*charge);
            thee->fcf[IJK(ilo,jlo,khi)] += ((1.0-dx)*(1.0-dy)*dz*charge);
            thee->fcf[IJK(ilo,jhi,klo)] += ((1.0-dx)*dy*(1.0-dz)*charge);
            thee->fcf[IJK(ilo,jlo,klo)] += ((1.0-dx)*(1.0-dy)*(1.0-dz)*charge);
        } /* endif (on the mesh) */
    } /* endfor (each atom) */

    /* THE FOLLOWING NEEDS TO BE DONE IF WE'RE NOT USING A SIMPLE LAPLACIAN
     * OPERATOR */
    if (!islap) {

        Vnm_print(0, "Vpmg_fillco:  marking ion and solvent accessibility.\n");

        /* Loop through the atoms and do the following:
         * 1.  Set ccf = -1.0, for all points inside the
         *     (possibly spline-based) inflated van der Waals surface
	 * 2.  Set a{123}cf = -1.0 if a point is inside the inflated van der
	 *     Waals radii
	 * 3.  Set a{123}cf = -2.0 if a point is inside the van der Waals radii
	 * 4.  Fill in the source term array
         */
        for (iatom=0; iatom<Valist_getNumberAtoms(alist); iatom++) {

            atom = Valist_getAtom(alist, iatom);
            apos = Vatom_getPosition(atom);
            arad = Vatom_getRadius(atom);
            charge = Vatom_getCharge(atom);

            /* Make sure we're on the grid */
            if ((apos[0]<=xmin) || (apos[0]>=xmax)  || \
                (apos[1]<=ymin) || (apos[1]>=ymax)  || \
                (apos[2]<=zmin) || (apos[2]>=zmax)) {
                if (thee->pmgp->bcfl != 4) {
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
            } else { /* if we're on the mesh */

                /* Convert the atom position to grid reference frame */
                position[0] = apos[0] - xmin;
                position[1] = apos[1] - ymin;
                position[2] = apos[2] - zmin;
   
                /* MARK ION ACCESSIBILITY AND DIELECTRIC VALUES FOR LATER
                 * ASSIGNMENT (Steps #1-3) */
                if (surfMeth == 2) itot2 = VSQR(irad + arad + splineWin);     
                else itot2 = VSQR(irad + arad);
                if (surfMeth == 2) stot2 = VSQR(arad + splineWin);
                else stot2 = VSQR(srad + arad);
                arad2 = VSQR(arad);
		/* We'll search over grid points which are in the greater of
                 * these two radii */
                rtot = VMAX2((irad + arad), (srad + arad));
                rtot2 = VMAX2(itot2, stot2);
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
                        if (rtot2 > (dx2+dy2)) 
                          dz = VSQRT(rtot2-dx2-dy2)+0.5*hzed;
                        else dz = 0.5*hzed;
                        kmin = VMAX2(0,(int)ceil((position[2] - dz)/hzed));
                        kmax = VMIN2(nz-1,(int)floor((position[2] + dz)/hzed));
                        for (k=kmin; k<=kmax; k++) {
                            dz2 = VSQR(k*hzed - position[2]);
			    /* See if grid point is inside ivdw radius and set
                             * ccf accordingly (do spline assignment here) */
                            if ((dz2 + dy2 + dx2) <= itot2) 
                              thee->ccf[IJK(i,j,k)] = -1.0;
			    /* See if x-shifted grid point is inside ivdw rad */
                            if (thee->a1cf[IJK(i,j,k)] != -2.0) {
                                if ((dz2+dy2+VSQR((i+0.5)*hx-position[0]))
                                     <=stot2) {
                                    /* See if inside vdw rad */
                                    if ((dz2+dy2+VSQR((i+0.5)*hx-position[0]))
                                         <=arad2) 
                                      thee->a1cf[IJK(i,j,k)] = -2.0;
                                     else thee->a1cf[IJK(i,j,k)] = -1.0;
                                } 
                            }
                            /* See if y-shifted grid point is inside ivdw rad */
                            if (thee->a2cf[IJK(i,j,k)] != -2.0) {
                                if ((dz2+VSQR((j+0.5)*hy-position[1])+dx2) 
                                     <= stot2) {
                                    /* See if inside vdw rad */
                                    if ((dz2+VSQR((j+0.5)*hy-position[1])+dx2)
                                         <=arad2) 
                                      thee->a2cf[IJK(i,j,k)] = -2.0;
                                    else thee->a2cf[IJK(i,j,k)] = -1.0;
                                }        
                            }
                            /* See if z-shifted grid point is inside ivdw rad */
                            if (thee->a3cf[IJK(i,j,k)] != -2.0) {
                                if ((VSQR((k+0.5)*hzed-position[2])+dy2+dx2)
                                     <=stot2) {
                                    /* See if inside vdw rad */
                                    if ((VSQR((k+0.5)*hzed-position[2])+dy2+dx2)
                                         <=arad2)
                                      thee->a3cf[IJK(i,j,k)] = -2.0;
                                    else thee->a3cf[IJK(i,j,k)] = -1.0;
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
                    position[0] = thee->xf[i];
                    position[1] = thee->yf[j];
                    position[2] = thee->zf[k];

		    /* the scalar (0th derivative) entry.  This is simply a
		     * number between 0 and 1; the actual coefficent value is
		     * calculated in mypde.f */
                    if (surfMeth == 2) {
                        if (thee->ccf[IJK(i,j,k)] == -1.0) {
                           gpos[0] = i*hx + xmin;
                           gpos[1] = j*hy + ymin;
                           gpos[2] = k*hzed + zmin;
                           thee->ccf[IJK(i,j,k)] = 
                             ionmask * Vacc_splineAcc(acc, gpos, splineWin, 
                             irad);
                        } else thee->ccf[IJK(i,j,k)] = ionmask;
                    } else {
                        if (thee->ccf[IJK(i,j,k)] == -1.0) 
                          thee->ccf[IJK(i,j,k)] = 0.0;
                        else thee->ccf[IJK(i,j,k)] = ionmask;
                    }

		    /* The diagonal tensor (2nd derivative) entries.  Each of
                     * these entries is evaluated ad the grid edges midpoints */
                    switch (surfMeth) {

                      /* No dielectric smoothing */
                      case 0: 
                        /* x-direction */
                        if (thee->a1cf[IJK(i,j,k)] == -1.0) {
                            position[0] = thee->xf[i] + 0.5*hx;
                            position[1] = thee->yf[j];
                            position[2] = thee->zf[k];
                            if (Vacc_fastMolAcc(acc, position, srad) == 0) 
                              thee->a1cf[IJK(i,j,k)] = epsp; 
                            else thee->a1cf[IJK(i,j,k)] = epsw; 
                        }
                        /* y-direction */
                        if (thee->a2cf[IJK(i,j,k)] == -1.0) {
                            position[0] = thee->xf[i];
                            position[1] = thee->yf[j] + 0.5*hy;
                            position[2] = thee->zf[k];
                            if (Vacc_fastMolAcc(acc, position, srad) == 0) 
                              thee->a2cf[IJK(i,j,k)] = epsp; 
                            else thee->a2cf[IJK(i,j,k)] = epsw; 
                        }
                        /* z-direction */
                        if (thee->a3cf[IJK(i,j,k)] == -1.0) {
                            position[0] = thee->xf[i];
                            position[1] = thee->yf[j];
                            position[2] = thee->zf[k] + 0.5*hzed;
                            if (Vacc_fastMolAcc(acc, position, srad) == 0) 
                              thee->a3cf[IJK(i,j,k)] = epsp; 
                            else thee->a3cf[IJK(i,j,k)] = epsw; 
                        }
                        break; 
    
                      /* A very rudimentary form of dielectric smoothing.
		       * Specifically, the dielectric will be evaluated at the
		       * mid point and the two flanking mesh points.  The
		       * fraction of the grid edge in the solvent will then be
		       * calculated from these three values (i.e., either 0,
		       * 1/3, 2/3, or 1).  The dielectric value at the midpoint
		       * will then be assigned based on the usual dielectric
		       * smoothing formula: \epsilon_s\epsilon_i/(a\epsilon_s +
		       * (1-a)\epsilon_i)  */
                      case 1:
                        a000 = 1.0;
                        if ((thee->a1cf[IJK(i,j,k)] == -1.0) ||
                            (thee->a2cf[IJK(i,j,k)] == -1.0) ||
                            (thee->a3cf[IJK(i,j,k)] == -1.0)) {
                            position[0] = thee->xf[i];
                            position[1] = thee->yf[j];
                            position[2] = thee->zf[k];
                            a000 = Vacc_fastMolAcc(acc, position, srad);
                        }
                        /* x-direction */
                        if (thee->a1cf[IJK(i,j,k)] == -1.0) {
                            position[0] = thee->xf[i] + 0.5*hx;
                            position[1] = thee->yf[j];
                            position[2] = thee->zf[k];
                            acclo = a000;
                            accmid = Vacc_molAcc(acc, position, srad);
                            position[0] = thee->xf[i] + hx;
                            acchi = Vacc_molAcc(acc, position, srad);
                            accf = ((double)acchi+(double)accmid
                              +(double)acclo)/3.0;
                            thee->a1cf[IJK(i,j,k)] = 
                              epsw*epsp/((1-accf)*epsw + accf*epsp);
                        }
                        /* y-direction */
                        if (thee->a2cf[IJK(i,j,k)] == -1.0) {
                            position[0] = thee->xf[i];
                            position[1] = thee->yf[j] + 0.5*hy;
                            position[2] = thee->zf[k];
                            accmid = Vacc_molAcc(acc, position, srad);
                            acclo = a000;
                            position[1] = thee->yf[j] + hy;
                            acchi = Vacc_molAcc(acc, position, srad);
                            accf = ((double)acchi+(double)accmid
                              +(double)acclo)/3.0;
                            thee->a2cf[IJK(i,j,k)] = 
                              epsw*epsp/((1-accf)*epsw + accf*epsp);
                        }
                        /* z-direction */
                        if (thee->a3cf[IJK(i,j,k)] == -1.0) {
                            position[0] = thee->xf[i];
                            position[1] = thee->yf[j];
                            position[2] = thee->zf[k] + 0.5*hzed;
                            accmid = Vacc_molAcc(acc, position, srad);
                            acclo = a000;
                            position[2] = thee->zf[k] + hzed;
                            acchi = Vacc_molAcc(acc, position, srad);
                            accf = ((double)acchi+(double)accmid
                              +(double)acclo)/3.0;
                            thee->a3cf[IJK(i,j,k)] = 
                              epsw*epsp/((1-accf)*epsw + accf*epsp);
                        }
                        break;

                      /* Spline-based accessibility function for force
		       * calculations.  See Im et al, Comp Phys Comm 111,
		       * (1998) 59--75. */
                      case 2:
                        /* x-direction */
                        if (thee->a1cf[IJK(i,j,k)] < 0.0) {
                            position[0] = thee->xf[i] + 0.5*hx;
                            position[1] = thee->yf[j];
                            position[2] = thee->zf[k];
                            chi = Vacc_splineAcc(acc, position, splineWin, 0.0);
                            thee->a1cf[IJK(i,j,k)] = epsp + (epsw - epsp)*chi;
                        }
                        /* y-direction */
                        if (thee->a2cf[IJK(i,j,k)] < 0.0) {
                            position[0] = thee->xf[i];
                            position[1] = thee->yf[j] + 0.5*hy;
                            position[2] = thee->zf[k];
                            chi = Vacc_splineAcc(acc, position, splineWin, 0.0);
                            thee->a2cf[IJK(i,j,k)] = epsp + (epsw - epsp)*chi;
                        }
                        /* z-direction */
                        if (thee->a3cf[IJK(i,j,k)] < 0.0) {
                            position[0] = thee->xf[i];
                            position[1] = thee->yf[j];
                            position[2] = thee->zf[k] + 0.5*hzed;
                            chi = Vacc_splineAcc(acc, position, splineWin, 0.0);
                            thee->a3cf[IJK(i,j,k)] = epsp + (epsw - epsp)*chi;
                        }
                        break;


                      /* Oops, invalid surfMeth */
                      default:
                        Vnm_print(2, "Vpmg_fillco:  Bad surfMeth (%d)!\n", 
                          surfMeth);
                        VASSERT(0);
                    } /* end switch(surfMeth) */

                    /* Fill in the remaining dielectric values */
                    if (thee->a1cf[IJK(i,j,k)] == -2.0)
                      thee->a1cf[IJK(i,j,k)] = epsp;
                    if (thee->a2cf[IJK(i,j,k)] == -2.0)
                      thee->a2cf[IJK(i,j,k)] = epsp;
                    if (thee->a3cf[IJK(i,j,k)] == -2.0)
                      thee->a3cf[IJK(i,j,k)] = epsp;
                    if (thee->a1cf[IJK(i,j,k)] == 0.0)
                      thee->a1cf[IJK(i,j,k)] = epsw;
                    if (thee->a2cf[IJK(i,j,k)] == 0.0)
                      thee->a2cf[IJK(i,j,k)] = epsw;
                    if (thee->a3cf[IJK(i,j,k)] == 0.0)
                      thee->a3cf[IJK(i,j,k)] = epsw;
                } /* i loop */
            } /* j loop */
        } /* k loop */

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
    if (thee->pmgp->bcfl != 4) {
        Vnm_print(0, "Vpmg_fillco:  filling boundary arrays\n");
        bcCalc(thee);
        Vnm_print(0, "Vpmg_fillco:  done filling boundary arrays\n");
    }

    thee->filled = 1;
}

/* ///////////////////////////////////////////////////////////////////////////
// Routine:  Vpmg_force
// Author:   Nathan Baker
/////////////////////////////////////////////////////////////////////////// */
VPUBLIC void Vpmg_force(Vpmg *thee, double *force, double gamma, 
  int atomID) {

    double qfF[3];                  /* Charge-field force */  
    double dbF[3];                  /* Dielectric boundary force */
    double ibF[3];                  /* Ion boundary force */
    double npF[3];                  /* Non-polar boundary force */

    VASSERT(thee != VNULL);
    VASSERT(thee->filled);
 
    Vpmg_dbnpForce(thee, qfF, npF, gamma, atomID);
    Vpmg_ibForce(thee, dbF, atomID); 
    Vpmg_qfForce(thee, ibF, atomID); 

    force[0] = qfF[0] + dbF[0] + npF[0] + ibF[0];
    force[1] = qfF[1] + dbF[1] + npF[1] + ibF[1];
    force[2] = qfF[2] + dbF[2] + npF[2] + ibF[2];

}

/* ///////////////////////////////////////////////////////////////////////////
// Routine:  Vpmg_ibForce
// Author:   Nathan Baker
/////////////////////////////////////////////////////////////////////////// */
VPUBLIC void Vpmg_ibForce(Vpmg *thee, double *force, int atomID) {

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
    VASSERT(thee->filled);
   
    acc = thee->pbe->acc;
    atom = Valist_getAtom(thee->pbe->alist, atomID);
    apos = Vatom_getPosition(atom);
    arad = Vatom_getRadius(atom);

    /* Reset force */
    force[0] = 0.0;
    force[1] = 0.0;
    force[2] = 0.0;

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

    /* Bail on focusing cases right now; not sure what to do */
    if (thee->pmgp->bcfl == 4) {
        Vnm_print(2, "Vpmg_ibForce:  Sorry, but force evaluation doesn't work with focusing (yet).\n");
        VASSERT(0);
    }


    /* Make sure we're on the grid */
    if ((apos[0]<=xmin) || (apos[0]>=xmax)  || \
      (apos[1]<=ymin) || (apos[1]>=ymax)  || \
      (apos[2]<=zmin) || (apos[2]>=zmax)) {
        if (thee->pmgp->bcfl != 4) {
            Vnm_print(2, "Vpmg_fillco:  Atom #%d at (%4.3f, %4.3f, %4.3f) is off the mesh (ignoring):\n",
                  atomID, position[0], position[1], position[2]);
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
                        Vacc_splineAccGrad(acc, gpos, thee->splineWin, irad,
                          atomID, tgrad);
                        if (thee->pmgp->nonlin) {
                            /* Nonlinear forces not done */
                            VASSERT(0);
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

/* ///////////////////////////////////////////////////////////////////////////
// Routine:  Vpmg_dbnpForce
// Author:   Nathan Baker
/////////////////////////////////////////////////////////////////////////// */
VPUBLIC void Vpmg_dbnpForce(Vpmg *thee, double *dbForce, double *npForce, 
  double gamma, int atomID) {

    Vacc *acc;
    Vpbe *pbe;
    Vatom *atom;

    double *apos, position[3], arad, hx, hy, hzed, izmagic;
    double xlen, ylen, zlen, xmin, ymin, zmin, xmax, ymax, zmax, rtot2, epsp;
    double rtot, dx, gpos[3], tgrad[3], dbFmag, epsw;
    double npFmag, *u, Hxijk, Hyijk, Hzijk, Hxim1jk, Hyijm1k, Hzijkm1;
    double dHxijk[3], dHyijk[3], dHzijk[3], dHxim1jk[3], dHyijm1k[3]; 
    double dHzijkm1[3];
    int i, j, k, nx, ny, nz, imin, imax, jmin, jmax, kmin, kmax;

    VASSERT(thee != VNULL);
    VASSERT(thee->filled);

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

    /* Get PBE info */
    pbe = thee->pbe;
    acc = pbe->acc;
    epsp = Vpbe_getSoluteDiel(pbe);
    epsw = Vpbe_getSolventDiel(pbe);
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

    /* Bail on focusing cases right now; not sure what to do */
    if (thee->pmgp->bcfl == 4) {
        Vnm_print(2, "Vpmg_ibForce:  Sorry, but force evaluation doesn't work with focusing (yet).\n");
        VASSERT(0);
    }

    /* Make sure we're on the grid */
    if ((apos[0]<=xmin) || (apos[0]>=xmax)  || \
      (apos[1]<=ymin) || (apos[1]>=ymax)  || \
      (apos[2]<=zmin) || (apos[2]>=zmax)) {
        if (thee->pmgp->bcfl != 4) {
            Vnm_print(2, "Vpmg_fillco:  Atom #%d at (%4.3f, %4.3f, %4.3f) is off the mesh (ignoring):\n",
                  atomID, position[0], position[1], position[2]);
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

        /* Integrate over points within this atom's (inflated) radius */
        rtot = (arad + thee->splineWin);
        rtot2 = VSQR(rtot);
        dx = rtot/hx;
        imin = VMAX2(1,(int)floor((position[0]-rtot)/hx));
        imax = VMIN2(nx-2,(int)ceil((position[0]+rtot)/hx));
        jmin = VMAX2(1,(int)floor((position[1]-rtot)/hy));
        jmax = VMIN2(ny-2,(int)ceil((position[1]+rtot)/hy));
        kmin = VMAX2(1,(int)floor((position[2]-rtot)/hzed));
        kmax = VMIN2(nz-2,(int)ceil((position[2]+rtot)/hzed));
        for (i=imin; i<=imax; i++) {
            for (j=jmin; j<=jmax; j++) {
                for (k=kmin; k<=kmax; k++) {
                    /* i,j,k */
                    gpos[0] = (i+0.5)*hx + xmin;
                    gpos[1] = j*hy + ymin;
                    gpos[2] = k*hzed + zmin;
                    Hxijk = Vacc_splineAcc(acc, gpos, thee->splineWin, 0.);
                    Vacc_splineAccGrad(acc, gpos, thee->splineWin, 0., atomID,
                      dHxijk);
                    gpos[0] = i*hx + xmin;
                    gpos[1] = (j+0.5)*hy + ymin;
                    gpos[2] = k*hzed + zmin;
                    Hyijk = Vacc_splineAcc(acc, gpos, thee->splineWin, 0.);
                    Vacc_splineAccGrad(acc, gpos, thee->splineWin, 0., atomID,
                      dHyijk);
                    gpos[0] = i*hx + xmin;
                    gpos[1] = j*hy + ymin;
                    gpos[2] = (k+0.5)*hzed + zmin;
                    Hzijk = Vacc_splineAcc(acc, gpos, thee->splineWin, 0.);
                    Vacc_splineAccGrad(acc, gpos, thee->splineWin, 0., atomID,
                      dHzijk);
                    /* i-1,j,k */
                    gpos[0] = (i-0.5)*hx + xmin;
                    gpos[1] = j*hy + ymin;
                    gpos[2] = k*hzed + zmin;
                    Hxim1jk = Vacc_splineAcc(acc, gpos, thee->splineWin, 0.);
                    Vacc_splineAccGrad(acc,gpos,thee->splineWin,0.,atomID,
                      dHxim1jk);
                    /* i,j-1,k */
                    gpos[0] = i*hx + xmin;
                    gpos[1] = (j-0.5)*hy + ymin;
                    gpos[2] = k*hzed + zmin;
                    Hyijm1k = Vacc_splineAcc(acc, gpos, thee->splineWin, 0.);
                    Vacc_splineAccGrad(acc, gpos,thee->splineWin,0.,atomID,
                      dHyijm1k);
                    /* i,j,k-1 */
                    gpos[0] = i*hx + xmin;
                    gpos[1] = j*hy + ymin;
                    gpos[2] = (k-0.5)*hzed + zmin;
                    Hzijkm1 = Vacc_splineAcc(acc, gpos, thee->splineWin, 0.);
                    Vacc_splineAccGrad(acc, gpos, thee->splineWin,0.,atomID,
                      dHzijkm1);
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
        
        dbForce[0] = -dbForce[0]*hx*hy*hzed*(epsw-epsp)*0.5*izmagic;
        dbForce[1] = -dbForce[1]*hx*hy*hzed*(epsw-epsp)*0.5*izmagic;
        dbForce[2] = -dbForce[2]*hx*hy*hzed*(epsw-epsp)*0.5*izmagic;
        npForce[0] = -npForce[0]*hx*hy*hzed*gamma;
        npForce[1] = -npForce[1]*hx*hy*hzed*gamma;
        npForce[2] = -npForce[2]*hx*hy*hzed*gamma;
    }
}

/* ///////////////////////////////////////////////////////////////////////////
// Routine:  Vpmg_qfForce
// Author:   Nathan Baker
/////////////////////////////////////////////////////////////////////////// */
VPUBLIC void Vpmg_qfForce(Vpmg *thee, double *force, int atomID) {

    Vatom *atom;
    
    double *apos, position[3], hx, hy, hzed;
    double xlen, ylen, zlen, xmin, ymin, zmin, xmax, ymax, zmax;
    double dx, dy, dz;
    double *u, charge, ifloat, jfloat, kfloat;
    int nx, ny, nz, ihi, ilo, jhi, jlo, khi, klo;

    VASSERT(thee != VNULL);
    VASSERT(thee->filled);

    atom = Valist_getAtom(thee->pbe->alist, atomID);
    apos = Vatom_getPosition(atom);
    charge = Vatom_getCharge(atom);

    /* Reset force */
    force[0] = 0.0;
    force[1] = 0.0;
    force[2] = 0.0;

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

    /* Bail on focusing cases right now; not sure what to do */
    if (thee->pmgp->bcfl == 4) {
        Vnm_print(2, "Vpmg_ibForce:  Sorry, but force evaluation doesn't work with focusing (yet).\n");
        VASSERT(0);
    }

    /* Make sure we're on the grid */
    if ((apos[0]<=xmin) || (apos[0]>=xmax) || (apos[1]<=ymin) || \
        (apos[1]>=ymax) || (apos[2]<=zmin) || (apos[2]>=zmax)) {
        if (thee->pmgp->bcfl != 4) {
            Vnm_print(2, "Vpmg_fillco:  Atom #%d at (%4.3f, %4.3f, %4.3f) is off the mesh (ignoring):\n", atomID, position[0], position[1], position[2]);
            Vnm_print(2, "Vpmg_fillco:    xmin = %g, xmax = %g\n", xmin, xmax);
            Vnm_print(2, "Vpmg_fillco:    ymin = %g, ymax = %g\n", ymin, ymax);
            Vnm_print(2, "Vpmg_fillco:    zmin = %g, zmax = %g\n", zmin, zmax);
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


        if (dx > VPMGSMALL) {
            force[0] = 
              -charge*(dy    *dz    *u[IJK(ihi,jhi,khi)]
                     + dy    *(1-dz)*u[IJK(ihi,jhi,klo)]
                     + (1-dy)*dz    *u[IJK(ihi,jlo,khi)]
                     + (1-dy)*(1-dz)*u[IJK(ihi,jlo,klo)]
                     - dy    *dz    *u[IJK(ilo,jhi,khi)]
                     - dy    *(1-dz)*u[IJK(ilo,jhi,klo)]
                     - (1-dy)*dz    *u[IJK(ilo,jlo,khi)]
                     - (1-dy)*(1-dz)*u[IJK(ilo,jlo,klo)])/hx;
        } else force[0] = 0;
        if (dy > VPMGSMALL) {
            force[1] = 
              -charge*(dx    *dz    *u[IJK(ihi,jhi,khi)]
                     + dx    *(1-dz)*u[IJK(ihi,jhi,klo)]
                     - dx    *dz    *u[IJK(ihi,jlo,khi)]
                     - dx    *(1-dz)*u[IJK(ihi,jlo,klo)]
                     + (1-dx)*dz    *u[IJK(ilo,jhi,khi)]
                     + (1-dx)*(1-dz)*u[IJK(ilo,jhi,klo)]
                     - (1-dx)*dz    *u[IJK(ilo,jlo,khi)]
                     - (1-dx)*(1-dz)*u[IJK(ilo,jlo,klo)])/hy;
        } else force[1] = 0;
        if (dz > VPMGSMALL) {
            force[2] = 
              -charge*(dy    *dx    *u[IJK(ihi,jhi,khi)]
                     - dy    *dx    *u[IJK(ihi,jhi,klo)]
                     + (1-dy)*dx    *u[IJK(ihi,jlo,khi)]
                     - (1-dy)*dx    *u[IJK(ihi,jlo,klo)]
                     + dy    *(1-dx)*u[IJK(ilo,jhi,khi)]
                     - dy    *(1-dx)*u[IJK(ilo,jhi,klo)]
                     + (1-dy)*(1-dx)*u[IJK(ilo,jlo,khi)]
                     - (1-dy)*(1-dx)*u[IJK(ilo,jlo,klo)])/hzed;
        } else force[2] = 0;
    }
}

/* ///////////////////////////////////////////////////////////////////////////
// Routine:  Vpmg_energy
// Author:   Nathan Baker
/////////////////////////////////////////////////////////////////////////// */
VPUBLIC double Vpmg_energy(Vpmg *thee, int extFlag) {

    double totEnergy = 0.0;
    double dielEnergy = 0.0;
    double qmEnergy = 0.0;
    double qfEnergy = 0.0;

    VASSERT(thee != VNULL);
    VASSERT(thee->filled);

    if ((thee->pmgp->nonlin) && (Vpbe_getBulkIonicStrength(thee->pbe) > 0.)) {
        Vnm_print(0, "Vpmg_energy:  calculating full PBE energy\n");
        qmEnergy = Vpmg_qmEnergy(thee, extFlag);
        Vnm_print(0, "Vpmg_energy:  qmEnergy = %g kT\n", qmEnergy);
        qfEnergy = Vpmg_qfEnergy(thee, extFlag);
        Vnm_print(0, "Vpmg_energy:  qfEnergy = %g kT\n", qfEnergy);
        dielEnergy = Vpmg_dielEnergy(thee, extFlag);
        Vnm_print(0, "Vpmg_energy:  dielEnergy = %g kT\n", dielEnergy);
        totEnergy = qfEnergy - dielEnergy - qmEnergy;
    } else {
        Vnm_print(0, "Vpmg_energy:  calculating only q-phi energy\n");
        qfEnergy = Vpmg_qfEnergy(thee, extFlag);
        Vnm_print(0, "Vpmg_energy:  qfEnergy = %g kT\n", qfEnergy);
        totEnergy = 0.5*qfEnergy;
    }

    return totEnergy;

}



/* ///////////////////////////////////////////////////////////////////////////
// Routine:  Vpmg_dielEnergy
// Author:   Nathan Baker
/////////////////////////////////////////////////////////////////////////// */
VPUBLIC double Vpmg_dielEnergy(Vpmg *thee, int extFlag) {

    double hx, hy, hzed, energy, nrgx, nrgy, nrgz, pvecx, pvecy, pvecz;
    int i, j, k, nx, ny, nz;
 
    VASSERT(thee != VNULL);
    VASSERT(thee->filled);

    /* Get the mesh information */
    nx = thee->pmgp->nx;
    ny = thee->pmgp->ny;
    nz = thee->pmgp->nz;
    hx = thee->pmgp->hx;
    hy = thee->pmgp->hy;
    hzed = thee->pmgp->hzed;

    energy = 0.0;

    /* Refill the dieletric coefficient arrays */
    Vpmg_fillco(thee, thee->surfMeth, thee->splineWin);

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

/* ///////////////////////////////////////////////////////////////////////////
// Routine:  Vpmg_qmEnergy
// Author:   Nathan Baker
/////////////////////////////////////////////////////////////////////////// */
VPUBLIC double Vpmg_qmEnergy(Vpmg *thee, int extFlag) {

    double hx, hy, hzed, energy, ionConc[MAXION], ionRadii[MAXION];
    double ionQ[MAXION], zkappa2, ionstr, zks2;
    int i, j, nx, ny, nz, nion, ichop, nchop;
 
    VASSERT(thee != VNULL);
    VASSERT(thee->filled);

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
    if (zkappa2 == 0.0) {
        Vnm_print(0, "Vpmg_qmEnergy:  Zero energy for zero ionic strength!\n");
        return 0.0;
    }
    zks2 = 0.5*zkappa2/ionstr;

    /* Because PMG seems to overwrite some of the coefficient arrays... */
    Vpmg_fillco(thee, thee->surfMeth, thee->splineWin);

    energy = 0.0;
    nchop = 0;
    Vpbe_getIons(thee->pbe, &nion, ionConc, ionRadii, ionQ);
    if (thee->pmgp->nonlin) {
        Vnm_print(0, "Vpmg_qmEnergy:  Calculating nonlinear energy\n");
        for (i=0; i<(nx*ny*nz); i++) {
            if (thee->pvec[i]*thee->ccf[i] > 0) {
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
            if (thee->pvec[i]*thee->ccf[i] > 0) 
              energy += (thee->pvec[i]*zkappa2*thee->ccf[i]*VSQR(thee->u[i]));
        }
        energy = 0.5*energy;
    }
    energy = energy*hx*hy*hzed;
    energy = energy/Vpbe_getZmagic(thee->pbe);

    if (extFlag == 1) energy += thee->extQmEnergy;

    return energy;
}
    
/* ///////////////////////////////////////////////////////////////////////////
// Routine:  Vpmg_qfEnergy
// Author:   Nathan Baker
/////////////////////////////////////////////////////////////////////////// */
VPUBLIC double Vpmg_qfEnergy(Vpmg *thee, int extFlag) {

    int iatom, nx, ny, nz, ihi, ilo, jhi, jlo, khi, klo;
    double xmax, xmin, ymax, ymin, zmax, zmin, hx, hy, hzed, ifloat, jfloat;
    double charge, kfloat, dx, dy, dz, energy, uval, *position;
    double *u;
    int *pvec;
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
    xmax = thee->xf[nx-1];
    ymax = thee->yf[ny-1];
    zmax = thee->zf[nz-1];
    xmin = thee->xf[0];
    ymin = thee->yf[0];
    zmin = thee->zf[0];

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

        if (atom->partID) {

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
                energy += (uval*charge);
            } else if (thee->pmgp->bcfl != 4) {
                Vnm_print(2, "Vpmg_qfEnergy:  Atom #%d at (%4.3f, %4.3f, \
%4.3f) is off the mesh (ignoring)!\n",
                iatom, position[0], position[1], position[2]);
            }
        } 
    }

    if (extFlag) energy += thee->extQfEnergy;
 
    return energy;
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
// Routine:  Vpmg_writeUHBD
// Author:   Nathan Baker
/////////////////////////////////////////////////////////////////////////// */
VPUBLIC void Vpmg_writeUHBD(Vpmg *thee, const char *iodev, const char *iofmt, 
  const char *thost, const char *fname, char *title, double *data) {

    int icol, i, j, k, u, nx, ny, nz, gotit;
    double xmin, ymin, zmin, hzed, hy, hx;
    Vio *sock;

    VASSERT(thee != VNULL);
    if ((thee->pmgp->hx!=thee->pmgp->hy) || (thee->pmgp->hy!=thee->pmgp->hzed) 
      || (thee->pmgp->hx!=thee->pmgp->hzed)) {
        Vnm_print(2, "Vpmg_writeUHBD: can't write UHBD mesh with non-uniform spacing\n");
        return;
    }

    /* Set up the virtual socket */
    sock = Vio_ctor(iodev,iofmt,thost,fname,"w");
    if (sock == VNULL) {
        Vnm_print(2, "Vpmg_writeUHBD: Problem opening virtual socket %s\n",
          fname);
        return;
    }
    if (Vio_connect(sock, 0) < 0) {
        Vnm_print(2, "Vpmg_writeUHBD: Problem connecting virtual socket %s\n",
          fname);
        return;
    }

    /* Get the lower corner and number of grid points for the local 
     * partition */
    hx = thee->pmgp->hx;
    hy = thee->pmgp->hy;
    hzed = thee->pmgp->hzed;
    nx = thee->pmgp->nx;
    ny = thee->pmgp->ny;
    nz = thee->pmgp->nz;
    xmin = thee->pmgp->xcent - 0.5*hx*(nx-1);
    ymin = thee->pmgp->ycent - 0.5*hy*(ny-1);
    zmin = thee->pmgp->zcent - 0.5*hzed*(nz-1);

    /* Let interested folks know that partition information is ignored */
    gotit = 0;
    for (i=0; i<(nx*ny*nz); i++) {
        if (thee->pvec[i] == 0) {
            gotit = 1;
            break;
        }
    }
    if (gotit) { 
        Vnm_print(2, "Vpmg_writeUHBD:  IGNORING PARTITION INFORMATION!\n");
        Vnm_print(2, "Vpmg_writeUHBD:  This means I/O from parallel runs will\
 have significant overlap.\n");
    }
 
    /* Write out the header */
    Vio_printf(sock, "%72s\n", title);
    Vio_printf(sock, "%12.5E%12.5E%7d%7d%7d%7d%7d\n", 1.0, 0.0, -1, 0, 
      nz, 1, nz);
    Vio_printf(sock, "%7d%7d%7d%12.5E%12.5E%12.5E%12.5E\n", nx, ny, nz, 
      hx, xmin, ymin, zmin);
    Vio_printf(sock, "%12.5E%12.5E%12.5E%12.5E\n", 0.0, 0.0, 0.0, 0.0);
    Vio_printf(sock, "%12.5E%12.5E%7d%7d", 0.0, 0.0, 0, 0);

    /* Write out the entries */
    icol = 0;
    for (k=0; k<nz; k++) {
        Vio_printf(sock, "\n%7d%7d%7d\n", k+1, thee->pmgp->nx, thee->pmgp->ny);
        icol = 0;
        for (j=0; j<ny; j++) {
            for (i=0; i<nx; i++) {
                u = k*(nx)*(ny)+j*(nx)+i;
                icol++;
                Vio_printf(sock, " %12.5E", data[u]);
                if (icol == 6) {
                    icol = 0;
                    Vio_printf(sock, "\n");
                }
            }
        }
    } 
    if (icol != 0) Vio_printf(sock, "\n");

    /* Close off the socket */
    Vio_connectFree(sock);
    Vio_dtor(&sock);
}

/* ///////////////////////////////////////////////////////////////////////////
// Routine:  Vpmg_readDX
//
// Author:   Nathan Baker
/////////////////////////////////////////////////////////////////////////// */
VPUBLIC void Vpmg_readDX(const char *iodev, const char *iofmt,
  const char *thost, const char *fname, 
  int *nx, int *ny, int *nz, 
  double *hx, double *hy, double *hzed, 
  double *xmin, double *ymin, double *zmin,
  double **data) {

    int i, j, k, itmp, u;
    double dtmp;
    char tok[VMAX_BUFSIZE];
    Vio *sock;
    char *MCwhiteChars = " =,;\t\n";
    char *MCcommChars  = "#%";


    /* Set up the virtual socket */
    sock = Vio_ctor(iodev,iofmt,thost,fname,"r");
    if (sock == VNULL) {
        Vnm_print(2, "Vpmg_readDX: Problem opening virtual socket %s\n",
          fname);
        return;
    }
    if (Vio_accept(sock, 0) < 0) {
        Vnm_print(2, "Vpmg_readDX: Problem accepting virtual socket %s\n",
          fname);
        return;
    }

    Vio_setWhiteChars(sock, MCwhiteChars);
    Vio_setCommChars(sock, MCcommChars);
    

    /* Read in the DX regular positions */
    /* Get "object" */
    VJMPERR2(1 == Vio_scanf(sock, "%s", tok));
    VJMPERR1(!strcmp(tok, "object"));
    /* Get "1" */
    VJMPERR2(1 == Vio_scanf(sock, "%s", tok));
    /* Get "class" */
    VJMPERR2(1 == Vio_scanf(sock, "%s", tok));
    VJMPERR1(!strcmp(tok, "class"));
    /* Get "gridpositions" */
    VJMPERR2(1 == Vio_scanf(sock, "%s", tok));
    VJMPERR1(!strcmp(tok, "gridpositions"));
    /* Get "counts" */
    VJMPERR2(1 == Vio_scanf(sock, "%s", tok));
    VJMPERR1(!strcmp(tok, "counts"));
    /* Get nz */
    VJMPERR2(1 == Vio_scanf(sock, "%s", tok));
    VJMPERR1(1 == sscanf(tok, "%d", nz));
    /* Get ny */
    VJMPERR2(1 == Vio_scanf(sock, "%s", tok));
    VJMPERR1(1 == sscanf(tok, "%d", ny));
    /* Get nx */
    VJMPERR2(1 == Vio_scanf(sock, "%s", tok));
    VJMPERR1(1 == sscanf(tok, "%d", nx));
    Vnm_print(0, "Vpmg_readDX:  Grid dimensions %d x %d x %d grid\n", *nx, *ny, *nz);
    /* Get "origin" */ 
    VJMPERR2(1 == Vio_scanf(sock, "%s", tok));
    VJMPERR1(!strcmp(tok, "origin"));
    /* Get zmin */
    VJMPERR2(1 == Vio_scanf(sock, "%s", tok));
    VJMPERR1(1 == sscanf(tok, "%lf", zmin));
    /* Get ymin */
    VJMPERR2(1 == Vio_scanf(sock, "%s", tok));
    VJMPERR1(1 == sscanf(tok, "%lf", ymin));
    /* Get xmin */
    VJMPERR2(1 == Vio_scanf(sock, "%s", tok));
    VJMPERR1(1 == sscanf(tok, "%lf", xmin));
    Vnm_print(0, "Vpmg_readDX:  Grid origin = (%g, %g, %g)\n", 
      *xmin, *ymin, *zmin);
    /* Get "delta" */ 
    VJMPERR2(1 == Vio_scanf(sock, "%s", tok));
    VJMPERR1(!strcmp(tok, "delta"));
    /* Get 0.0 */
    VJMPERR2(1 == Vio_scanf(sock, "%s", tok));
    VJMPERR1(1 == sscanf(tok, "%lf", &dtmp));
    VJMPERR1(dtmp == 0.0);
    /* Get 0.0 */
    VJMPERR2(1 == Vio_scanf(sock, "%s", tok));
    VJMPERR1(1 == sscanf(tok, "%lf", &dtmp));
    VJMPERR1(dtmp == 0.0);
    /* Get hz */
    VJMPERR2(1 == Vio_scanf(sock, "%s", tok));
    VJMPERR1(1 == sscanf(tok, "%lf", hzed));
    /* Get "delta" */ 
    VJMPERR2(1 == Vio_scanf(sock, "%s", tok));
    VJMPERR1(!strcmp(tok, "delta"));
    /* Get 0.0 */
    VJMPERR2(1 == Vio_scanf(sock, "%s", tok));
    VJMPERR1(1 == sscanf(tok, "%lf", &dtmp));
    VJMPERR1(dtmp == 0.0);
    /* Get hy */
    VJMPERR2(1 == Vio_scanf(sock, "%s", tok));
    VJMPERR1(1 == sscanf(tok, "%lf", hy));
    /* Get 0.0 */
    VJMPERR2(1 == Vio_scanf(sock, "%s", tok));
    VJMPERR1(1 == sscanf(tok, "%lf", &dtmp));
    VJMPERR1(dtmp == 0.0);
    /* Get "delta" */ 
    VJMPERR2(1 == Vio_scanf(sock, "%s", tok));
    VJMPERR1(!strcmp(tok, "delta"));
    /* Get hx */
    VJMPERR2(1 == Vio_scanf(sock, "%s", tok));
    VJMPERR1(1 == sscanf(tok, "%lf", hx));
    Vnm_print(0, "Vpmg_readDX:  Grid spacings = (%g, %g, %g)\n", 
      *hx, *hy, *hzed);
    /* Get 0.0 */
    VJMPERR2(1 == Vio_scanf(sock, "%s", tok));
    VJMPERR1(1 == sscanf(tok, "%lf", &dtmp));
    VJMPERR1(dtmp == 0.0);
    /* Get 0.0 */
    VJMPERR2(1 == Vio_scanf(sock, "%s", tok));
    VJMPERR1(1 == sscanf(tok, "%lf", &dtmp));
    VJMPERR1(dtmp == 0.0);
    /* Get "object" */ 
    VJMPERR2(1 == Vio_scanf(sock, "%s", tok));
    VJMPERR1(!strcmp(tok, "object"));
    /* Get "2" */ 
    VJMPERR2(1 == Vio_scanf(sock, "%s", tok));
    /* Get "class" */ 
    VJMPERR2(1 == Vio_scanf(sock, "%s", tok));
    VJMPERR1(!strcmp(tok, "class"));
    /* Get "gridconnections" */ 
    VJMPERR2(1 == Vio_scanf(sock, "%s", tok));
    VJMPERR1(!strcmp(tok, "gridconnections"));
    /* Get "counts" */ 
    VJMPERR2(1 == Vio_scanf(sock, "%s", tok));
    VJMPERR1(!strcmp(tok, "counts"));
    /* Get the dimensions again */ 
    VJMPERR2(1 == Vio_scanf(sock, "%s", tok));
    VJMPERR2(1 == Vio_scanf(sock, "%s", tok));
    VJMPERR2(1 == Vio_scanf(sock, "%s", tok));
    /* Get "object" */ 
    VJMPERR2(1 == Vio_scanf(sock, "%s", tok));
    VJMPERR1(!strcmp(tok, "object"));
    /* Get # */ 
    VJMPERR2(1 == Vio_scanf(sock, "%s", tok));
    /* Get "class" */ 
    VJMPERR2(1 == Vio_scanf(sock, "%s", tok));
    VJMPERR1(!strcmp(tok, "class"));
    /* Get "array" */ 
    VJMPERR2(1 == Vio_scanf(sock, "%s", tok));
    VJMPERR1(!strcmp(tok, "array"));
    /* Get "type" */ 
    VJMPERR2(1 == Vio_scanf(sock, "%s", tok));
    VJMPERR1(!strcmp(tok, "type"));
    /* Get "double" */ 
    VJMPERR2(1 == Vio_scanf(sock, "%s", tok));
    VJMPERR1(!strcmp(tok, "double"));
    /* Get "rank" */ 
    VJMPERR2(1 == Vio_scanf(sock, "%s", tok));
    VJMPERR1(!strcmp(tok, "rank"));
    /* Get # */ 
    VJMPERR2(1 == Vio_scanf(sock, "%s", tok));
    /* Get "items" */ 
    VJMPERR2(1 == Vio_scanf(sock, "%s", tok));
    VJMPERR1(!strcmp(tok, "items"));
    /* Get # */ 
    VJMPERR2(1 == Vio_scanf(sock, "%s", tok));
    VJMPERR1(1 == sscanf(tok, "%d", &itmp));
    VJMPERR1(((*nx)*(*ny)*(*nz)) == itmp);
    /* Get "data" */ 
    VJMPERR2(1 == Vio_scanf(sock, "%s", tok));
    VJMPERR1(!strcmp(tok, "data"));
    /* Get "follows" */ 
    VJMPERR2(1 == Vio_scanf(sock, "%s", tok));
    VJMPERR1(!strcmp(tok, "follows"));

    /* Allocate space for the data */
    Vnm_print(0, "Vpmg_readDX:  allocating %d x %d x %d doubles for storage\n",
      *nx, *ny, *nz);
    *data = VNULL;
    *data = Vmem_malloc(VNULL, (*nx)*(*ny)*(*nz), sizeof(double));
    if (*data == VNULL) {
        Vnm_print(2, "Vpmg_readDX:  Unable to allocate space for data!\n");
        return;
    }

    for (k=0; k<*nz; k++) {
        for (j=0; j<*ny; j++) {
            for (i=0; i<*nx; i++) {
                u = k*(*nx)*(*ny)+j*(*nx)+i;
                VJMPERR2(1 == Vio_scanf(sock, "%s", tok));
                VJMPERR1(1 == sscanf(tok, "%lf", &dtmp));
                (*data)[u] = dtmp;
            }
        }
    }

    /* Close off the socket */
    Vio_acceptFree(sock);
    Vio_dtor(&sock);

    return;

  VERROR1:
    Vio_dtor(&sock);
    Vnm_print(2, "Vpmg_readDX:  Format problem with input file <%s>\n", 
      fname);
    return;

  VERROR2:
    Vio_dtor(&sock);
    Vnm_print(2, "Vpmg_readDX:  I/O problem with input file <%s>\n", 
      fname);
    return;
   


}

/* ///////////////////////////////////////////////////////////////////////////
// Routine:  Vpmg_writeDX
//
// Author:   Nathan Baker
/////////////////////////////////////////////////////////////////////////// */
VPUBLIC void Vpmg_writeDX(Vpmg *thee, const char *iodev, const char *iofmt,
  const char *thost, const char *fname, char *title, double *data) {

    double xmin, ymin, zmin, hx, hy, hzed;
    int nx, ny, nz;

    hx = thee->pmgp->hx;
    hy = thee->pmgp->hy;
    hzed = thee->pmgp->hzed;
    nx = thee->pmgp->nx;
    ny = thee->pmgp->ny;
    nz = thee->pmgp->nz;
    xmin = thee->pmgp->xcent - 0.5*hx*(nx-1);
    ymin = thee->pmgp->ycent - 0.5*hy*(ny-1);
    zmin = thee->pmgp->zcent - 0.5*hzed*(nz-1);

    Vpmg_writeDX2(iodev, iofmt, thost, fname, title, data, thee->pvec,
      hx, hy, hzed, nx, ny, nz, xmin, ymin, zmin);
}


/* ///////////////////////////////////////////////////////////////////////////
// Routine:  Vpmg_writeDX2
//
// Author:   Nathan Baker
/////////////////////////////////////////////////////////////////////////// */
VPUBLIC void Vpmg_writeDX2(const char *iodev, const char *iofmt,
  const char *thost, const char *fname, char *title, double *data, int *pvec,
  double hx, double hy, double hzed, int nx, int ny, int nz, 
  double xmin, double ymin, double zmin) {

    int icol, i, j, k, u, usepart, nxPART, nyPART, nzPART, gotit;
    double x, y, z, xminPART, yminPART, zminPART;
    Vio *sock;

    if (pvec == VNULL) usepart = 0;
    else usepart = 1;

    /* Set up the virtual socket */
    sock = Vio_ctor(iodev,iofmt,thost,fname,"w");
    if (sock == VNULL) {
        Vnm_print(2, "Vpmg_writeDX: Problem opening virtual socket %s\n",
          fname);
        return;
    }
    if (Vio_connect(sock, 0) < 0) {
        Vnm_print(2, "Vpmg_writeDX: Problem connecting virtual socket %s\n",
          fname);
        return;
    }

    if (usepart) {
        /* Get the lower corner and number of grid points for the local
         * partition */
        xminPART = VLARGE;
        yminPART = VLARGE;
        zminPART = VLARGE;
        nxPART = 0;
        nyPART = 0;
        nzPART = 0;
        /* First, search for the lower corner */
        for (k=0; k<nz; k++) {
            z = k*hzed + zmin;
            for (j=0; j<ny; j++) {
                y = j*hy + ymin;
                for (i=0; i<nx; i++) {
                    x = i*hx + xmin;
                    if (pvec[IJK(i,j,k)] != 0) {
                        if (x < xminPART) xminPART = x;
                        if (y < yminPART) yminPART = y;
                        if (z < zminPART) zminPART = z;
                    }
                }
            }
        }
        /* Now search for the number of grid points in the z direction */
        for (k=0; k<nz; k++) {
            gotit = 0;
            for (j=0; j<ny; j++) {
                for (i=0; i<nx; i++) {
                    if (pvec[IJK(i,j,k)] != 0) {
                        gotit = 1;
                        break;
                    }
                }
                if (gotit) break;
            }
            if (gotit) nzPART++;
        }
        /* Now search for the number of grid points in the y direction */
        for (j=0; j<ny; j++) {
            gotit = 0;
            for (k=0; k<nz; k++) {
                for (i=0; i<nx; i++) {
                    if (pvec[IJK(i,j,k)] != 0) {
                        gotit = 1;
                        break;
                    }
                }
                if (gotit) break;
            }
            if (gotit) nyPART++;
        }
        /* Now search for the number of grid points in the x direction */
        for (i=0; i<nx; i++) {
            gotit = 0;
            for (k=0; k<nz; k++) {
                for (j=0; j<ny; j++) {
                    if (pvec[IJK(i,j,k)] != 0) {
                        gotit = 1;
                        break;
                    }
                }
                if (gotit) break;
            }
            if (gotit) nxPART++;
        }

        if ((nxPART != nx) || (nyPART != ny) || (nzPART != nz)) {
            Vnm_print(0, "Vpmg_writeUHBD:  printing only subset of domain\n");
        }
            

        /* Write off the title */
        Vio_printf(sock, "# Electrostatic potential data from APBS/PMG\n");
        Vio_printf(sock, "# \n");
        Vio_printf(sock, "# %s\n", title);
        Vio_printf(sock, "# \n");
    
        /* Write off the DX regular positions */
        Vio_printf(sock, "object 1 class gridpositions counts %d %d %d\n",
          nzPART, nyPART, nxPART);
        Vio_printf(sock, "origin %12.6E %12.6E %12.6E\n", xminPART, yminPART,
          zminPART);
        Vio_printf(sock, "delta %12.6E %12.6E %12.6E\n", 0.0, 0.0, hzed);
        Vio_printf(sock, "delta %12.6E %12.6E %12.6E\n", 0.0, hy, 0.0);
        Vio_printf(sock, "delta %12.6E %12.6E %12.6E\n", hx, 0.0, 0.0);
    
        /* Write off the DX regular connections */
        Vio_printf(sock, "object 2 class gridconnections counts %d %d %d\n",
          nzPART, nyPART, nxPART);
    
        /* Write off the DX data */
        Vio_printf(sock, "object 3 class array type double rank 0 items %d data follows\n",
          (nxPART*nyPART*nzPART));
        icol = 0;
        for (k=0; k<nz; k++) {
            for (j=0; j<ny; j++) {
                for (i=0; i<nx; i++) {
                    u = k*(nx)*(ny)+j*(nx)+i;
                    if (pvec[u] != 0) {
                        Vio_printf(sock, "%12.6E ", data[u]);
                        icol++;
                        if (icol == 3) {
                            icol = 0;
                            Vio_printf(sock, "\n");
                        }
                    }
                }
            }
        }

        if (icol != 0) Vio_printf(sock, "\n");
                    
        /* Create the field */
        Vio_printf(sock, "attribute \"dep\" string \"positions\"\n");
        Vio_printf(sock, "object \"regular positions regular connections\" class field\n");
        Vio_printf(sock, "component \"positions\" value 1\n");
        Vio_printf(sock, "component \"connections\" value 2\n");
        Vio_printf(sock, "component \"data\" value 3\n");

    } else { 
    
        /* Write off the title */
        Vio_printf(sock, "# Electrostatic potential data from APBS/PMG\n");
        Vio_printf(sock, "# \n");
        Vio_printf(sock, "# %s\n", title);
        Vio_printf(sock, "# \n");
    
        /* Write off the DX regular positions */
        Vio_printf(sock, "object 1 class gridpositions counts %d %d %d\n",
          nz, ny, nx);
        Vio_printf(sock, "origin %12.6E %12.6E %12.6E\n", xmin, ymin, zmin);
        Vio_printf(sock, "delta %12.6E %12.6E %12.6E\n", 0.0, 0.0, hzed);
        Vio_printf(sock, "delta %12.6E %12.6E %12.6E\n", 0.0, hy, 0.0);
        Vio_printf(sock, "delta %12.6E %12.6E %12.6E\n", hx, 0.0, 0.0);
    
        /* Write off the DX regular connections */
        Vio_printf(sock, "object 2 class gridconnections counts %d %d %d\n",
          nz, ny, nx);
    
        /* Write off the DX data */
        Vio_printf(sock, "object 3 class array type double rank 0 items %d data follows\n",
          (nx*ny*nz));
        icol = 0;
        for (k=0; k<nz; k++) {
            for (j=0; j<ny; j++) {
                for (i=0; i<nx; i++) {
                    u = k*(nx)*(ny)+j*(nx)+i;
                    Vio_printf(sock, "%12.6E ", data[u]);
                    icol++;
                    if (icol == 3) {
                        icol = 0;
                        Vio_printf(sock, "\n");
                    }
                }
            }
        }
        if (icol != 0) Vio_printf(sock, "\n");
                    
        /* Create the field */
        Vio_printf(sock, "attribute \"dep\" string \"positions\"\n");
        Vio_printf(sock, "object \"regular positions regular connections\" class field\n");
        Vio_printf(sock, "component \"positions\" value 1\n");
        Vio_printf(sock, "component \"connections\" value 2\n");
        Vio_printf(sock, "component \"data\" value 3\n");
    }
    
    /* Close off the socket */
    Vio_connectFree(sock);
    Vio_dtor(&sock);

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
    VASSERT(thee->filled == 1);

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
// Routine:  Vpmg_fillAcc
//
// Author:   Nathan Baker
/////////////////////////////////////////////////////////////////////////// */
VPUBLIC void Vpmg_fillAcc(Vpmg *thee, double *vec, int meth, double parm) {

    Vacc *acc = VNULL;
    Vpbe *pbe = VNULL;
    double position[3], hx, hy, hzed, xmin, ymin, zmin;
    int i, j, k, nx, ny, nz;

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

    if (meth == 0) {
        Vnm_print(0, "Vpmg_fillAcc: using molecular surface with %g A probe\n",
          parm);
    } else if (meth == 1) {
        Vnm_print(0, "Vpmg_fillAcc: using van der Waals surface\n");
    } else if (meth == 2) {
        Vnm_print(0, "Vpmg_fillAcc: using inflated van der Waals surface with %g A probe\n",
          parm);
    } else if (meth == 3) {
        Vnm_print(0, "Vpmg_fillAcc: using spline surface with %g window\n",
          parm);
    } else {
        Vnm_print(2, "Vpmg_fillAcc: invalid surface method (%d)!\n", meth);
        VASSERT(0);
    }


    for (k=0; k<nz; k++) {
        for (j=0; j<ny; j++) {
            for (i=0; i<nx; i++) {

                position[0] = i*hx + xmin;
                position[1] = j*hy + ymin;
                position[2] = k*hzed + zmin;

                /* the scalar (0th derivative) entry */
                if (meth == 0) {
                    vec[IJK(i,j,k)] = (Vacc_molAcc(acc,position,parm));
                } else if (meth == 1) {
                    vec[IJK(i,j,k)] = (Vacc_vdwAcc(acc,position));
                } else if (meth == 2) {
                    vec[IJK(i,j,k)] = (Vacc_ivdwAcc(acc,position,parm));
                } else if (meth == 3) {
                    vec[IJK(i,j,k)] = (Vacc_splineAcc(acc,position,parm,0.0));
                }
            }
        }
    }
}

#undef VPMGSMALL
