/**
 *  @file    vopot.c
 *  @author  Nathan Baker
 *  @brief   Class Vopot methods
 *  @ingroup Vopot
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
 * Copyright (c) 1999-2002.  Nathan A. Baker.  All Rights Reserved.
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
#include "apbs/vopot.h"

VEMBED(rcsid="$Id$")

/* ///////////////////////////////////////////////////////////////////////////
// Routine:  Vopot_ctor
// Author:   Nathan Baker
/////////////////////////////////////////////////////////////////////////// */
VPUBLIC Vopot* Vopot_ctor(int nx, int ny, int nz,
                  double hx, double hy, double hzed,
                  double xmin, double ymin, double zmin,
                  double *data, Vpbe *pbe) {

    Vopot *thee = VNULL;

    thee = Vmem_malloc(VNULL, 1, sizeof(Vopot));
    VASSERT(thee != VNULL);
    VASSERT(Vopot_ctor2(thee, nx, ny, nz, hx, hy, hzed,
                  xmin, ymin, zmin, data, pbe));

    return thee;
}

/* ///////////////////////////////////////////////////////////////////////////
// Routine:  Vopot_ctor2
// Author:   Nathan Baker
/////////////////////////////////////////////////////////////////////////// */
VPUBLIC int Vopot_ctor2(Vopot *thee, int nx, int ny, int nz,
                  double hx, double hy, double hzed,
                  double xmin, double ymin, double zmin,
                  double *data, Vpbe *pbe) {

    if (thee == VNULL) return 0;
    thee->nx = nx;
    thee->ny = ny;
    thee->nz = nz;
    thee->hx = hx;
    thee->hy = hy;
    thee->hzed = hzed;
    thee->xmin = xmin;
    thee->ymin = ymin;
    thee->zmin = zmin;
    thee->data = data;
    thee->pbe = pbe;

    return 1;
}

/* ///////////////////////////////////////////////////////////////////////////
// Routine:  Vopot_dtor
// Author:   Nathan Baker
/////////////////////////////////////////////////////////////////////////// */
VPUBLIC void Vopot_dtor(Vopot **thee) {

    if ((*thee) != VNULL) {
        Vopot_dtor2(*thee);
        Vmem_free(VNULL, 1, sizeof(Vopot), (void **)thee);
        (*thee) = VNULL;
    }
}

/* ///////////////////////////////////////////////////////////////////////////
// Routine:  Vopot_dtor2
// Author:   Nathan Baker
/////////////////////////////////////////////////////////////////////////// */
VPUBLIC void Vopot_dtor2(Vopot *thee) { return; }

/* ///////////////////////////////////////////////////////////////////////////
// Routine:  Vopot_pot
// Author:   Nathan Baker
/////////////////////////////////////////////////////////////////////////// */
#define IJK(i,j,k)  (((k)*(nx)*(ny))+((j)*(nx))+(i))
VPUBLIC double Vopot_pot(Vopot *thee, double pt[3]) {

    int nx, ny, nz, ihi, jhi, khi, ilo, jlo, klo;
    double hx, hy, hzed, xmin, ymin, zmin, ifloat, jfloat, kfloat;
    double u, dx, dy, dz, T, charge, eps_w, xkappa, dist;
    Valist *alist;

    VASSERT(thee != VNULL);

    eps_w = Vpbe_getSolventDiel(thee->pbe);
    xkappa = (1.0e10)*Vpbe_getXkappa(thee->pbe);
    T = Vpbe_getTemperature(thee->pbe);
    alist = Vpbe_getValist(thee->pbe);

    nx = thee->nx;
    ny = thee->ny;
    nz = thee->nz;
    hx = thee->hx;
    hy = thee->hy;
    hzed = thee->hzed;
    xmin = thee->xmin;
    ymin = thee->ymin;
    zmin = thee->zmin;

    u = 0;

    ifloat = (pt[0] - xmin)/hx;
    jfloat = (pt[1] - ymin)/hy;
    kfloat = (pt[2] - zmin)/hzed;
    ihi = (int)ceil(ifloat);
    jhi = (int)ceil(jfloat);
    khi = (int)ceil(kfloat);
    ilo = (int)floor(ifloat);
    jlo = (int)floor(jfloat);
    klo = (int)floor(kfloat);

    /* See if we're on the mesh */
    if ((ihi<nx) && (jhi<ny) && (khi<nz) &&
        (ilo>=0) && (jlo>=0) && (klo>=0)) {

        dx = ifloat - (double)(ilo);
        dy = jfloat - (double)(jlo);
        dz = kfloat - (double)(klo);
        u = dx      *dy      *dz      *(thee->data[IJK(ihi,jhi,khi)])
          + dx      *(1.0-dy)*dz      *(thee->data[IJK(ihi,jlo,khi)])
          + dx      *dy      *(1.0-dz)*(thee->data[IJK(ihi,jhi,klo)])
          + dx      *(1.0-dy)*(1.0-dz)*(thee->data[IJK(ihi,jlo,klo)])
          + (1.0-dx)*dy      *dz      *(thee->data[IJK(ilo,jhi,khi)])
          + (1.0-dx)*(1.0-dy)*dz      *(thee->data[IJK(ilo,jlo,khi)])
          + (1.0-dx)*dy      *(1.0-dz)*(thee->data[IJK(ilo,jhi,klo)])
          + (1.0-dx)*(1.0-dy)*(1.0-dz)*(thee->data[IJK(ilo,jlo,klo)]);

    } else {

#if (VOPOT_BCFL == 0)

        u = 0;

#elseif (VOPOT_BCFL == 1)

        size = (1.0e-10)*Vpbe_getSoluteRadius(thee->pbe);
        position = Vpbe_getSoluteCenter(thee->pbe);
        charge = Vunit_ec*Vpbe_getSoluteCharge(thee->pbe);
        dist = 0;
        for (i=0; i<d; i++)
          dist += ((position[i] - pt[i])*(position[i] - pt[i]));
        dist = (1.0e-10)*VSQRT(dist);
        val = (charge)/(4*VPI*Vunit_eps0*eps_w*dist);
        if (xkappa != 0.0) val = val*(exp(-xkappa*(dist-size))/(1+xkappa*size));
        val = val*Vunit_ec/(Vunit_kb*T);

        u = val;


#elseif (VOPOT_BCFL == 2)

        u = 0;
        for (iatom=0; iatom<Valist_getNumberAtoms(alist); iatom++) {
            atom = Valist_getAtom(alist, iatom);
            position = Vatom_getPosition(atom);
            charge = Vunit_ec*Vatom_getCharge(atom);
            size = (1e-10)*Vatom_getRadius(atom);
            dist = 0;
            for (i=0; i<d; i++)
              dist += ((position[i] - pt[i])*(position[i] - pt[i]));
            dist = (1.0e-10)*VSQRT(dist);
            val = (charge)/(4*VPI*Vunit_eps0*eps_w*dist);
            if (xkappa != 0.0)
              val = val*(exp(-xkappa*(dist-size))/(1+xkappa*size));
            val = val*Vunit_ec/(Vunit_kb*T);
            u = u + val;
        }

#endif

    }

    return u;

}

/* ///////////////////////////////////////////////////////////////////////////
// Routine:  Vopot_curvature
//
//   Notes:  cflag=0 ==> Reduced Maximal Curvature
//           cflag=1 ==> Mean Curvature (Laplace)
//           cflag=2 ==> Gauss Curvature
//           cflag=3 ==> True Maximal Curvature
//
// Authors:  Stephen Bond and Nathan Baker
/////////////////////////////////////////////////////////////////////////// */
VPUBLIC double Vopot_curvature(Vopot *thee, double pt[3], int cflag)
{
    double hx, hy, hzed, curv;
    double dxx, dyy, dzz;
    double uleft, umid, uright, testpt[3];

    hx = thee->hx;
    hy = thee->hy;
    hzed = thee->hzed;

    curv = 0.0;

    testpt[0] = pt[0];
    testpt[1] = pt[1];
    testpt[2] = pt[2];

    /* Compute 2nd derivative in the x-direction */
    umid = Vopot_pot( thee, testpt );
    testpt[0] = pt[0] - hx;
    uleft = Vopot_pot( thee, testpt );
    testpt[0] = pt[0] + hx;
    uright = Vopot_pot( thee, testpt );
    testpt[0] = pt[0];

    dxx = (uright - 2*umid + uleft)/(hx*hx);

    /* Compute 2nd derivative in the y-direction */
    umid = Vopot_pot( thee, testpt );
    testpt[1] = pt[1] - hy;
    uleft = Vopot_pot( thee, testpt );
    testpt[1] = pt[1] + hy;
    uright = Vopot_pot( thee, testpt );
    testpt[1] = pt[1];

    dyy = (uright - 2*umid + uleft)/(hy*hy);

    /* Compute 2nd derivative in the z-direction */
    umid = Vopot_pot( thee, testpt );
    testpt[2] = pt[2] - hzed;
    uleft = Vopot_pot( thee, testpt );
    testpt[2] = pt[2] + hzed;
    uright = Vopot_pot( thee, testpt );

    dzz = (uright - 2*umid + uleft)/(hzed*hzed);


    if ( cflag == 0 ) {
        curv = fabs(dxx);
        curv = ( curv > fabs(dyy) ) ? curv : fabs(dyy);
        curv = ( curv > fabs(dzz) ) ? curv : fabs(dzz);
    } else if ( cflag == 1 ) {
        curv = (dxx + dyy + dzz)/3.0;
    } else {
        VASSERT( 0 ); /* Feature Not Coded Yet! */
    }

    return curv;

}

/* ///////////////////////////////////////////////////////////////////////////
// Routine:  Vopot_gradient
//
// Authors:  Nathan Baker and Stephen Bond
/////////////////////////////////////////////////////////////////////////// */
VPUBLIC void Vopot_gradient(Vopot *thee, double pt[3], double grad[3]) {

    double hx, hy, hzed;
    double uleft, umid, uright, testpt[3];

    hx = thee->hx;
    hy = thee->hy;
    hzed = thee->hzed;

    testpt[0] = pt[0];
    testpt[1] = pt[1];
    testpt[2] = pt[2];

    /* Compute derivative in the x-direction */
    umid = Vopot_pot( thee, testpt );
    testpt[0] = pt[0] - hx;
    uleft = Vopot_pot( thee, testpt );
    testpt[0] = pt[0] + hx;
    uright = Vopot_pot( thee, testpt );
    testpt[0] = pt[0];

    grad[0] = (uright - uleft)/(2*hx);

    /* Compute derivative in the y-direction */
    umid = Vopot_pot( thee, testpt );
    testpt[1] = pt[1] - hy;
    uleft = Vopot_pot( thee, testpt );
    testpt[1] = pt[1] + hy;
    uright = Vopot_pot( thee, testpt );
    testpt[1] = pt[1];

    grad[1] = (uright - uleft)/(2*hy);

    /* Compute derivative in the z-direction */
    umid = Vopot_pot( thee, testpt );
    testpt[2] = pt[2] - hzed;
    uleft = Vopot_pot( thee, testpt );
    testpt[2] = pt[2] + hzed;
    uright = Vopot_pot( thee, testpt );

    grad[2] = (uright - uleft)/(2*hzed);

}

