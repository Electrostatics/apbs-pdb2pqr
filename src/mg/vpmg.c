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
// File:     vpmg.c
//
// Purpose:  Class Vpmg: methods.
//
// Author:   Nathan Baker
/////////////////////////////////////////////////////////////////////////// */

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
// Class Vpmg: Non-inlineable methods
/////////////////////////////////////////////////////////////////////////// */

/* ///////////////////////////////////////////////////////////////////////////
// Routine:  focusFillBound
//
// Purpose:  Fill boundaries with old values before destroying
//
// Author:   Nathan Baker
/////////////////////////////////////////////////////////////////////////// */
VPRIVATE void focusFillBound(Vpmg *thee, Vpmg *pmgOLD) {

    double hxOLD, hyOLD, hzOLD, xminOLD, yminOLD, zminOLD, xmaxOLD, ymaxOLD;
    double zmaxOLD;
    int nxOLD, nyOLD, nzOLD;
    double hxNEW, hyNEW, hzNEW, xminNEW, yminNEW, zminNEW, xmaxNEW, ymaxNEW;
    double zmaxNEW;
    int nxNEW, nyNEW, nzNEW;
    int i, j, k, ihi, ilo, jhi, jlo, khi, klo, nx, ny, nz;
    double x, y, z, dx, dy, dz, ifloat, jfloat, kfloat, uval;


    /* Calculate new problem dimensions */
    hxNEW = thee->pmgp->hx;
    hyNEW = thee->pmgp->hy;
    hzNEW = thee->pmgp->hzed;
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

    /* Sanity check: make sure we're within the old mesh */
    if (((xmaxNEW-xmaxOLD)>VPMGSMALL) || 
        ((ymaxNEW-ymaxOLD)>VPMGSMALL) || 
        ((zmaxNEW-zmaxOLD)>VPMGSMALL) ||
        ((xminOLD-xminNEW)>VPMGSMALL) || 
        ((yminOLD-yminNEW)>VPMGSMALL) || 
        ((zminOLD-zminNEW)>VPMGSMALL)) {
        Vnm_print(2, "VPMG::focusFillBound -- new mesh not contained in old!\n");
        Vnm_print(2, "VPMG::focusFillBound -- New mesh mins = %g, %g, %g\n",
          xminNEW, yminNEW, zminNEW);
        Vnm_print(2, "VPMG::focusFillBound -- New mesh maxs = %g, %g, %g\n",
          xmaxNEW, ymaxNEW, zmaxNEW);
        Vnm_print(2, "VPMG::focusFillBound -- Old mesh mins = %g, %g, %g\n",
          xminOLD, yminOLD, zminOLD);
        Vnm_print(2, "VPMG::focusFillBound -- Old mesh maxs = %g, %g, %g\n",
          xmaxOLD, ymaxOLD, zmaxOLD);
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
            thee->gxcf[IJKx(j,k,0)] = uval;

            /* High X face */
            x = xmaxNEW;
            ifloat = (x - xminOLD)/hxOLD;
            ihi = (int)ceil(ifloat);
            if (ihi > (nxOLD-1)) ihi = nxOLD-1;
            ilo = (int)floor(ifloat);
            if (ilo < 0) ilo = 0;
            dx = ifloat - (double)(ilo);
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
            thee->gycf[IJKy(i,k,0)] = uval;

            /* High Y face */
            y = ymaxNEW;
            jfloat = (y - yminOLD)/hyOLD;
            jhi = (int)ceil(jfloat);
            if (jhi > (nyOLD-1)) jhi = nyOLD-1;
            jlo = (int)floor(jfloat);
            if (jlo < 0) jlo = 0;
            dy = jfloat - (double)(jlo);
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
            thee->gzcf[IJKz(i,j,0)] = uval;

            /* High Z face */
            z = zmaxNEW;
            kfloat = (z - zminOLD)/hzOLD;
            khi = (int)ceil(kfloat);
            if (khi > (nzOLD-1)) khi = nzOLD-1;
            klo = (int)floor(kfloat);
            if (klo < 0) klo = 0;
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
// Author:   Nathan Baker
/////////////////////////////////////////////////////////////////////////// */
VPRIVATE void extEnergy(Vpmg *thee, Vpmg *pmgOLD) {

    double hxNEW, hyNEW, hzNEW, xminNEW, yminNEW, zminNEW, xmaxNEW, ymaxNEW;
    double zmaxNEW;
    int nxNEW, nyNEW, nzNEW;
    int ihi, ilo, jhi, jlo, khi, klo, iatom;
    double x, y, z, ifloat, jfloat, kfloat, *position;
    Valist *alist; 
    Vatom *atom;

    /* Set the new external energy contribution to zero.  Any external
     * contributions from higher levels will be included in the appropriate
     * energy function call. */
    thee->extEnergy = 0;

    /* Calculate new problem dimensions */
    hxNEW = thee->pmgp->hx;
    hyNEW = thee->pmgp->hy;
    hzNEW = thee->pmgp->hzed;
    nxNEW = thee->pmgp->nx;
    nyNEW = thee->pmgp->ny;
    nzNEW = thee->pmgp->nz;
    xminNEW = thee->pmgp->xcent - ((double)(nxNEW-1)*hxNEW)/2.0;
    xmaxNEW = thee->pmgp->xcent + ((double)(nxNEW-1)*hxNEW)/2.0;
    yminNEW = thee->pmgp->ycent - ((double)(nyNEW-1)*hyNEW)/2.0;
    ymaxNEW = thee->pmgp->ycent + ((double)(nyNEW-1)*hyNEW)/2.0;
    zminNEW = thee->pmgp->zcent - ((double)(nzNEW-1)*hzNEW)/2.0;
    zmaxNEW = thee->pmgp->zcent + ((double)(nzNEW-1)*hzNEW)/2.0;

    /* Loop through the atoms, marking those outside the current domain */
    alist = Vpbe_getValist(pmgOLD->pbe);
    for (iatom=0; iatom<Valist_getNumberAtoms(alist); iatom++) {
        atom = Valist_getAtom(alist, iatom);
        position = Vatom_getPosition(atom);
        x = position[0];
        y = position[1];
        z = position[2];
        ifloat = (x - xminNEW)/hxNEW;
        jfloat = (y - yminNEW)/hyNEW;
        kfloat = (z - zminNEW)/hzNEW;
        ihi = (int)ceil(ifloat);
        ilo = (int)floor(ifloat);
        jhi = (int)ceil(jfloat);
        jlo = (int)floor(jfloat);
        khi = (int)ceil(kfloat);
        klo = (int)floor(kfloat);

        /* See if this atom is outside the new problem domain and mark it if it
         * is */
        if ((ihi<nxNEW) && (jhi<nyNEW) && (khi<nzNEW) &&
            (ilo>=0) && (jlo>=0) && (klo>=0)) Vatom_setPartID(atom, 0);
        else Vatom_setPartID(atom, 1);
    }

    /* Now calculate the energy on that subset of the domain */
    pmgOLD->partFlag = 1;
    if (pmgOLD->pmgp->nonlin == 0) {
        /* For linear calculations, we can just use the subset of atoms we
         * just marked. */
        thee->extEnergy = Vpmg_getLinearEnergy1(pmgOLD, 1);
    } else {
        /* For nonlinear calculations, we need to do a volume integral */
        Vnm_print(1, "extEnergy:  Focusing does not work with NPBE yet!\n");
        VASSERT(0);
    } 
    pmgOLD->partFlag = 0;
    for (iatom=0; iatom<Valist_getNumberAtoms(alist); iatom++) {
        atom = Valist_getAtom(alist, iatom);
        Vatom_setPartID(atom, 0);
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
VPRIVATE double bcCalc(Vpbe *pbe, double x[], int flag) {

    double size, *position, charge, xkappa, eps_w, dist, T, val, pot;
    int i, iatom;
    Vatom *atom;
    Valist *alist;

    if (flag == 0) {
        return 0.0;
    } else if (flag == 1) {
        /* Get the solute radius in meters and position in angstroms */
        size = (1.0e-10)*Vpbe_getSoluteRadius(pbe);
        position = Vpbe_getSoluteCenter(pbe);

        /* We keep the charge relative to units of ec that are factored out;
         * this term should be dimensionless. The dielectric is unitless. */
        charge = Vunit_ec*Vpbe_getSoluteCharge(pbe);
        eps_w = Vpbe_getSolventDiel(pbe);

        /* Get xkappa in units of inverse meters */
        xkappa = (1.0e10)*Vpbe_getXkappa(pbe);

        /* The temperature is in units of K */
        T = Vpbe_getTemperature(pbe);

        /* Compute the distance (in units of m) */
        dist = 0;
        for (i=0; i<3; i++) 
            dist += ((position[i] - x[i])*(position[i] - x[i])); 
        dist = (1.0e-10)*VSQRT(dist);

        /* Compute the potential in J/electron */
        val = (charge)/(4*VPI*Vunit_eps0*eps_w*dist);
        if (xkappa != 0.0) val = val*(exp(-xkappa*(dist-size))/(1+xkappa*size));
        /* Scale the potential to be dimensionless */
        val = val*Vunit_ec/(Vunit_kb*T);
        return val;
    } else if (flag == 2) {
        pot = 0.0;
        eps_w = Vpbe_getSolventDiel(pbe);
        xkappa = (1.0e10)*Vpbe_getXkappa(pbe);
        T = Vpbe_getTemperature(pbe);
        alist = Vpbe_getValist(pbe);
        for (iatom=0; iatom<Valist_getNumberAtoms(alist); iatom++) {
            atom = Valist_getAtom(alist, iatom);
            position = Vatom_getPosition(atom);
            charge = Vunit_ec*Vatom_getCharge(atom);
            size = (1e-10)*Vatom_getRadius(atom);
            dist = 0;
            for (i=0; i<3; i++)
               dist += ((position[i] - x[i])*(position[i] - x[i]));
            dist = (1.0e-10)*VSQRT(dist);
            val = (charge)/(4*VPI*Vunit_eps0*eps_w*dist);
            if (xkappa != 0.0)
              val = val*(exp(-xkappa*(dist-size))/(1+xkappa*size));
            val = val*Vunit_ec/(Vunit_kb*T);
            pot = pot + val;
        }

        return pot;
    } else if (flag == 4) {
        Vnm_print(2, "VPMG::bcCalc -- not appropriate for focusing!\n");
        VASSERT(0);
    } else {
        Vnm_print(2, "VPMG::bcCalc -- invalid boundary condition flag (%d)!\n",
          flag);
        VASSERT(0);
    }

    return 0;
}


/* ///////////////////////////////////////////////////////////////////////////
// Routine:  Vpmg_ctor
//
// Purpose:  Construct the PMG object from scratch using the solver parameters
//           specifed by the passed Vpmgp object and the equation data from
//           the Vpbe object
//
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
//
// Purpose:  Construct the PMG object
//
// Notes:    See header files for default parameter values
//
// Author:   Nathan Baker
/////////////////////////////////////////////////////////////////////////// */
VPUBLIC int Vpmg_ctor2(Vpmg *thee, Vpmgp *pmgp, Vpbe *pbe) {

    int nxc, nyc, nzc, nf, nc, narrc, n_rpc, n_iz, n_ipc;

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
      &(thee->pmgp->nlev), &nxc, &nyc, &nzc, &nf, &nc, &(thee->pmgp->narr),
      &narrc, &n_rpc, &n_iz, &n_ipc, &(thee->pmgp->nrwk), &(thee->pmgp->niwk));

    Vnm_print(2, "Vpmg_ctor2: PMG chose nx = %d, ny = %d, nz = %d, nlev = %d\n",
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
    thee->partFlag = 0;
    thee->partLower[0] = 0;
    thee->partLower[1] = 0;
    thee->partLower[2] = 0;
    thee->partUpper[0] = 0;
    thee->partUpper[1] = 0;
    thee->partUpper[2] = 0;

    /* Ignore external energy contributions */
    thee->extEnergy = 0;

    return 1;
}

/* ///////////////////////////////////////////////////////////////////////////
// Routine:  Vpmg_ctorFocus
//
// Purpose:  Construct the PMG object by focusing.  In other words, use the
//           solution from the passed Vpmg object to set the boundary
//           conditions for the new Vpmg object.  IN THE PROCESS, THE OLD VPMG
//           OBJECT IS DESTROYED.  The solver parameters specifed by the passed
//           Vpmgp object and the equation data from the Vpbe object are also
//           used.
//
// Author:   Nathan Baker
/////////////////////////////////////////////////////////////////////////// */
VPUBLIC Vpmg* Vpmg_ctorFocus(Vpmgp *pmgp, Vpbe *pbe, Vpmg *pmgOLD) {

    Vpmg *thee = VNULL;

    /* Set up the structure */
    thee = Vmem_malloc(VNULL, 1, sizeof(Vpmg) );
    VASSERT( thee != VNULL);
    VASSERT(Vpmg_ctor2Focus(thee, pmgp, pbe, pmgOLD));

    return thee;
}

/* ///////////////////////////////////////////////////////////////////////////
// Routine:  Vpmg_ctor2Focus
//
// Purpose:  Construct the PMG object
//
// Notes:    See Vpmg_ctor2Focus description
//
// Author:   Nathan Baker
/////////////////////////////////////////////////////////////////////////// */
VPUBLIC int Vpmg_ctor2Focus(Vpmg *thee, Vpmgp *pmgp, Vpbe *pbe, Vpmg *pmgOLD) {

    int nxc, nyc, nzc, nf, nc, narrc, n_rpc, n_iz, n_ipc;

    /* Get the parameters */    
    VASSERT(pmgp != VNULL); 
    VASSERT(pbe != VNULL); 
    VASSERT(pmgOLD != VNULL); 
    thee->pmgp = pmgp;
    thee->pbe = pbe;

    /* Set up the memory */
    thee->vmem = Vmem_ctor("APBS:VPMG");

    /* Calculate storage requirements */
    F77MGSZ(&(thee->pmgp->mgcoar), &(thee->pmgp->mgdisc), &(thee->pmgp->mgsolv),
      &(thee->pmgp->nx), &(thee->pmgp->ny), &(thee->pmgp->nz),
      &(thee->pmgp->nlev), &nxc, &nyc, &nzc, &nf, &nc, &(thee->pmgp->narr),
      &narrc, &n_rpc, &n_iz, &n_ipc, &(thee->pmgp->nrwk), &(thee->pmgp->niwk));

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
    focusFillBound(thee, pmgOLD);

    /* Calculate energetic contributions from region outside focusing domain */
    extEnergy(thee, pmgOLD);

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
    thee->partFlag = 0;
    thee->partLower[0] = 0;
    thee->partLower[1] = 0;
    thee->partLower[2] = 0;
    thee->partUpper[0] = 0;
    thee->partUpper[1] = 0;
    thee->partUpper[2] = 0;


    return 1;
}

/* ///////////////////////////////////////////////////////////////////////////
// Routine:  Vpmg_solve
//
// Purpose:  Solve the equation
//
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
//
// Purpose:  Fill the coefficient arrays prior to solving the equation
//
// Author:   Nathan Baker
/////////////////////////////////////////////////////////////////////////// */
VPUBLIC void Vpmg_fillco(Vpmg *thee) {

    Vacc *acc;
    Valist *alist;
    Vpbe *pbe;
    Vatom *atom;
    double xmin, xmax, ymin, ymax, zmin, zmax;
    double xlen, ylen, zlen, position[3], ifloat, jfloat, kfloat, accf;
    double zmagic, irad, srad, charge, dx, dy, dz, zkappa2, epsw, epsp;
    double hx, hy, hzed;
    int i, j, k, nx, ny, nz, iatom, ihi, ilo, jhi, jlo, khi, klo;
    int acclo, accmid, acchi;

    /* Get PBE info */
    pbe = thee->pbe;
    acc = pbe->acc;
    alist = pbe->alist;
    irad = Vpbe_getIonRadius(pbe);
    srad = Vpbe_getSolventRadius(pbe);
    zmagic = Vpbe_getZmagic(pbe);
    zkappa2 = Vpbe_getZkappa2(pbe);
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
    xlen = hx*(nx - 1);
    ylen = hy*(ny - 1);
    zlen = hzed*(nz - 1);

    /* Define the min/max dimensions */
    xmin = thee->pmgp->xcent - (xlen/2.0);
    ymin = thee->pmgp->ycent - (ylen/2.0);
    zmin = thee->pmgp->zcent - (zlen/2.0);
    xmax = thee->pmgp->xcent + (xlen/2.0);
    ymax = thee->pmgp->ycent + (ylen/2.0);
    zmax = thee->pmgp->zcent + (zlen/2.0);
    thee->rparm[2] = xmin;
    thee->rparm[3] = xmax;
    thee->rparm[4] = ymin;
    thee->rparm[5] = ymax;
    thee->rparm[6] = zmin;
    thee->rparm[7] = zmax;

    /* Fill the mesh point coordinate arrays */
    for (i=0; i<nx; i++) thee->xf[i] = xmin + i*hx;
    for (i=0; i<ny; i++) thee->yf[i] = ymin + i*hy;
    for (i=0; i<nz; i++) thee->zf[i] = zmin + i*hzed;

    /* Fill the coefficient arrays */
    for (k=0; k<nz; k++) {
        for (j=0; j<ny; j++) {
            for (i=0; i<nx; i++) {

                position[0] = thee->xf[i];
                position[1] = thee->yf[j];
                position[2] = thee->zf[k];

                /* the scalar (0th derivative) entry */
                if (Vacc_ivdwAcc(acc, position, irad) == 1) 
                  thee->ccf[IJK(i,j,k)] = zkappa2;
                else thee->ccf[IJK(i,j,k)] = 0.0;

                /* the chosen true solution */
                thee->tcf[IJK(i,j,k)] = 0.0;

                /* Clear out the load vector */
                thee->fcf[IJK(i,j,k)] = 0.0;

                /* The diagonal tensor (2nd derivative) entries.  Each of these
                 * entries is evaluated ad the grid edges midpoints.  We will
                 * implement a very rudimentary form of dielectric smoothing.
                 * Specifically, the dielectric will be evaluated at the mid
                 * point and the two flanking mesh points.  The fraction of the
                 * grid edge in the solvent will then be calculated from these
                 * three values (i.e., either 0, 1/3, 2/3, or 1).  The
                 * dielectric value at the midpoint will then be assigned based
                 * on the usual dielectric smoothing formula:
                 * \epsilon_s\epsilon_i/(a\epsilon_s + (1-a)\epsilon_i)  */
                /* x-direction */
                position[0] = thee->xf[i] + hx/2.0;
                position[1] = thee->yf[j];
                position[2] = thee->zf[k];
                accmid = Vacc_molAcc(acc, position, srad);
                position[0] = thee->xf[i];
                acclo = Vacc_molAcc(acc, position, srad);
                position[0] = thee->xf[i] + hx;
                acchi = Vacc_molAcc(acc, position, srad);
                accf = ((double)acchi + (double)accmid + (double)acclo)/3.0;
                thee->a1cf[IJK(i,j,k)] = 
                  epsw*epsp/((1-accf)*epsw + accf*epsp);
                /* y-direction */
                position[0] = thee->xf[i];
                position[1] = thee->yf[j] + hy/2.0;
                position[2] = thee->zf[k];
                accmid = Vacc_molAcc(acc, position, srad);
                position[1] = thee->yf[j];
                acclo = Vacc_molAcc(acc, position, srad);
                position[1] = thee->yf[j] + hy;
                acchi = Vacc_molAcc(acc, position, srad);
                accf = ((double)acchi + (double)accmid + (double)acclo)/3.0;
                thee->a2cf[IJK(i,j,k)] = 
                  epsw*epsp/((1-accf)*epsw + accf*epsp);
                /* y-direction */
                position[0] = thee->xf[i];
                position[1] = thee->yf[j];
                position[2] = thee->zf[k] + hzed/2.0;
                accmid = Vacc_molAcc(acc, position, srad);
                position[2] = thee->zf[k];
                acclo = Vacc_molAcc(acc, position, srad);
                position[2] = thee->zf[k] + hzed;
                acchi = Vacc_molAcc(acc, position, srad);
                accf = ((double)acchi + (double)accmid + (double)acclo)/3.0;
                thee->a3cf[IJK(i,j,k)] = 
                  epsw*epsp/((1-accf)*epsw + accf*epsp);
            }
        }
    }

    /* Fill source term array */
    for (iatom=0; iatom<Valist_getNumberAtoms(alist); iatom++) {
        atom = Valist_getAtom(alist, iatom);

        position[0] = Vatom_getPosition(atom)[0];
        position[1] = Vatom_getPosition(atom)[1];
        position[2] = Vatom_getPosition(atom)[2];

        /* Make sure we're on the grid */
        if ((position[0]<=xmin) || (position[0]>=xmax)  || \
            (position[1]<=ymin) || (position[1]>=ymax)  || \
            (position[2]<=zmin) || (position[2]>=zmax)) {
            if (thee->pmgp->bcfl != 4) {
                Vnm_print(2, "Vpmg_fillco:  Atom #%d at (%4.3f, %4.3f, %4.3f) is off the mesh (ignoring):\n",
                  iatom, position[0], position[1], position[2]);
                Vnm_print(2, "Vpmg_fillco:    xmin = %g, xmax = %g\n", 
                  xmin, xmax);
                Vnm_print(2, "Vpmg_fillco:    ymin = %g, ymax = %g\n", 
                  ymin, ymax);
                Vnm_print(2, "Vpmg_fillco:    zmin = %g, zmax = %g\n", 
                  zmin, zmax);
            }
            fflush(stderr);
        } else {

            /* Figure out which vertices we're next to */
            ifloat = (position[0] - xmin)/(hx);
            jfloat = (position[1] - ymin)/(hy);
            kfloat = (position[2] - zmin)/(hzed);

            ihi = (int)ceil(ifloat);
            ilo = (int)floor(ifloat);
            jhi = (int)ceil(jfloat);
            jlo = (int)floor(jfloat);
            khi = (int)ceil(kfloat);
            klo = (int)floor(kfloat);

            /* Now assign fractions of the charge to the nearby verts */
            charge = zmagic*Vatom_getCharge(atom)/hx/hy/hzed;
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

        }
    }

    /* Fill the boundary arrays (except when focusing, bcfl = 4) */
    if (thee->pmgp->bcfl != 4) {
        /* the "i" boundaries (dirichlet) */
        for (k=0; k<nz; k++) {
            for (j=0; j<ny; j++) {
                position[0] = thee->xf[0];
                position[1] = thee->yf[j];
                position[2] = thee->zf[k];
                thee->gxcf[IJKx(j,k,0)] = bcCalc(pbe, position, 
                  thee->pmgp->bcfl);
                position[0] = thee->xf[nx-1];
                thee->gxcf[IJKx(j,k,1)] = bcCalc(pbe, position, 
                  thee->pmgp->bcfl);
                thee->gxcf[IJKx(j,k,2)] = 0.0;
                thee->gxcf[IJKx(j,k,3)] = 0.0;
            }
        }

        /* the "j" boundaries (dirichlet) */
        for (k=0; k<nz; k++) {
            for (i=0; i<nx; i++) {
                position[0] = thee->xf[i];
                position[1] = thee->yf[0];
                position[2] = thee->zf[k];
                thee->gycf[IJKy(i,k,0)] = bcCalc(pbe, position, 
                  thee->pmgp->bcfl);
                position[1] = thee->yf[ny-1];
                thee->gycf[IJKy(i,k,1)] = bcCalc(pbe, position, 
                  thee->pmgp->bcfl);
                thee->gycf[IJKy(i,k,2)] = 0.0;
                thee->gycf[IJKy(i,k,3)] = 0.0;
            }
        }

        /* the "k" boundaries (dirichlet) */
        for (j=0; j<ny; j++) {
            for (i=0; i<nx; i++) {
                position[0] = thee->xf[i];
                position[1] = thee->yf[j];
                position[2] = thee->zf[0];
                thee->gzcf[IJKz(i,j,0)] = bcCalc(pbe, position, 
                  thee->pmgp->bcfl);
                position[2] = thee->zf[nz-1];
                thee->gzcf[IJKz(i,j,1)] = bcCalc(pbe, position, 
                  thee->pmgp->bcfl);
                thee->gzcf[IJKz(i,j,2)] = 0.0;
                thee->gzcf[IJKz(i,j,3)] = 0.0;
            }
        }
    }
}
    
/* ///////////////////////////////////////////////////////////////////////////
// Routine:  Vpmg_getLinearEnergy1
//
// Purpose:  using the solution at the finest mesh level, get the
//           electrostatic energy using the free energy functional for the
//           linearized Poisson-Boltzmann equation without removing any
//           self-interaction terms (i.e., removing the reference state of
//           isolated charges present in an infinite dielectric continuum with
//           the same relative permittivity as the interior of the protein).
//           In other words, we calculate
//             \[ G = \frac{1}{2} \sum_i q_i u(r_i) \]
//           and return the result in units of $k_B T$.  
//
// Args:     extFlag => If this was a focused calculation, then it is possible
//                      to include the energy contributions from the outside
//                      the focused domain.  This should be on (=1) for
//                      sequential focusing calculations and off (=0) for
//                      parallel calculations.
//     
// Notes:    The value of this observable may be modified by setting
//           restrictions on the subdomain over which it is calculated.  Such
//           limits can be set via Vpmg_setPart and are generally useful for
//           parallel runs.  In such cases, the values are calculated on the
//           interval [min, max).
//
// Author:   Nathan Baker
/////////////////////////////////////////////////////////////////////////// */
VPUBLIC double Vpmg_getLinearEnergy1(Vpmg *thee, int extFlag) {

    int iatom, nx, ny, nz, ihi, ilo, jhi, jlo, khi, klo;
    double xmax, xmin, ymax, ymin, zmax, zmin, hx, hy, hzed, ifloat, jfloat;
    double charge, kfloat, dx, dy, dz, energy, uval, *position;
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
  
    energy = 0.0;

    for (iatom=0; iatom<Valist_getNumberAtoms(alist); iatom++) {
        /* Get atomic information */
        atom = Valist_getAtom(alist, iatom);

        if ((thee->partFlag == 0) || (Vatom_getPartID(atom) == 1)) {

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
    
            if ((ihi<nx) && (jhi<ny) && (khi<nz) &&
                (ilo>=0) && (jlo>=0) && (klo>=0)) {
    
                /* Now get trilinear interpolation constants */
                dx = ifloat - (double)(ilo);
                dy = jfloat - (double)(jlo);
                dz = kfloat - (double)(klo);
                uval =  dx*dy*dz*(thee->u[IJK(ihi,jhi,khi)])
                      + dx*(1.0-dy)*dz*(thee->u[IJK(ihi,jlo,khi)])
                      + dx*dy*(1.0-dz)*(thee->u[IJK(ihi,jhi,klo)])
                      + dx*(1.0-dy)*(1.0-dz)*(thee->u[IJK(ihi,jlo,klo)])
                      + (1.0-dx)*dy*dz*(thee->u[IJK(ilo,jhi,khi)])
                      + (1.0-dx)*(1.0-dy)*dz*(thee->u[IJK(ilo,jlo,khi)])
                      + (1.0-dx)*dy*(1.0-dz)*(thee->u[IJK(ilo,jhi,klo)])
                      + (1.0-dx)*(1.0-dy)*(1.0-dz)*(thee->u[IJK(ilo,jlo,klo)]);
                energy += (uval*charge);
            } else if (thee->pmgp->bcfl != 4) {
                Vnm_print(2, "Vpmg_getLE1:  Atom #%d at (%4.3f, %4.3f, %4.3f) is off the mesh (ignoring)!\n",
                    iatom, position[0], position[1], position[2]);
            }
        }
    }

    energy = 0.5*energy;
    if (extFlag == 1) energy += (thee->extEnergy);
    return energy;
}
    
/* ///////////////////////////////////////////////////////////////////////////
// Routine:  Vpmg_dtor
//
// Purpose:  Clean up
//
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
//
// Purpose:  Clean up
//
// Author:   Nathan Baker
/////////////////////////////////////////////////////////////////////////// */
VPUBLIC void Vpmg_dtor2(Vpmg *thee) { 

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

    Vmem_dtor(&(thee->vmem));
}

/* ///////////////////////////////////////////////////////////////////////////
// Routine:  Vpmg_writeUHBD
//
// Purpose:  Write out a PMG array in UHBD grid format (ASCII)
//
// Args:     path => file name
//           title => title to be inserted in grid
//           data => nx*ny*nz length array of data
//
// Notes:    The mesh spacing should be uniform
//
// Author:   Nathan Baker
/////////////////////////////////////////////////////////////////////////// */
VPUBLIC void Vpmg_writeUHBD(Vpmg *thee, char *path, char *title, 
  double *data) {

    FILE *fp = VNULL;
    int icol, i, j, k, u;
    double xmin, ymin, zmin;

    VASSERT(thee != VNULL);
    if ((thee->pmgp->hx!=thee->pmgp->hy) || (thee->pmgp->hy!=thee->pmgp->hzed) 
      || (thee->pmgp->hx!=thee->pmgp->hzed)) {
        Vnm_print(2, "Vpmg_writeUHBD: can't write UHBD mesh with non-uniform spacing\n");
        return;
    }

    /* Open the file */
    fp = fopen(path, "w");
    if (fp == VNULL) {
        Vnm_print(2, "Vgrid_writeUHBD: Error opening %s for writing!\n", path);
        return;
    }

    /* Write out the header */
    xmin = thee->pmgp->xcent - (thee->pmgp->hx)*(thee->pmgp->nx-1);
    ymin = thee->pmgp->ycent - (thee->pmgp->hy)*(thee->pmgp->ny-1);
    zmin = thee->pmgp->zcent - (thee->pmgp->hzed)*(thee->pmgp->nz-1);
    fprintf(fp, "%72s\n", title);
    fprintf(fp, "%12.6E%12.6E%7d%7d%7d%7d%7d\n", 1.0, 0.0, -1, 0, 
      thee->pmgp->nz, 1, thee->pmgp->nz);
    fprintf(fp, "%7d%7d%7d%12.6E%12.6E%12.6E%12.6E\n", thee->pmgp->nx,
      thee->pmgp->ny, thee->pmgp->nz, thee->pmgp->hx, xmin, ymin, zmin);
    fprintf(fp, "%12.6E%12.6E%12.6E%12.6E\n", 0.0, 0.0, 0.0, 0.0);
    fprintf(fp, "%12.6E%12.6E%7d%7d", 0.0, 0.0, 0, 0);

    /* Write out the entries */
    for (k=0; k<thee->pmgp->nz; k++) {
        fprintf(fp, "\n%7d%7d%7d\n", k+1, thee->pmgp->nx, thee->pmgp->ny);
        icol = 0;
        for (j=0; j<thee->pmgp->ny; j++) {
            for (i=0; i<thee->pmgp->nx; i++) {
                u = k*(thee->pmgp->nx)*(thee->pmgp->ny)+j*(thee->pmgp->nx)+i;
                icol++;
                fprintf(fp, " %12.6E", data[u]);
                if (icol == 6) {
                    icol = 0;
                    fprintf(fp, "\n");
                }
            }
        }
    } 
    if (icol != 0) fprintf(fp, "\n");
    fclose(fp);
}

/* ///////////////////////////////////////////////////////////////////////////
// Routine:  Vpmg_setPart
//
// Purpose:  Set partition information which restricts the calculation of
//           observables to a (rectangular) subset of the problem domain
//
// Args:     [xyz]min => lower corner
//           [xyz]max => upper corner
//
// Author:   Nathan Baker
/////////////////////////////////////////////////////////////////////////// */
VPUBLIC void Vpmg_setPart(Vpmg *thee, double xmin, double ymin, double zmin,
           double xmax, double ymax, double zmax) {

    Valist *alist;
    Vatom *atom;
    int iatom;
    double *position;

    /* Store the partition information */
    thee->partFlag = 1;
    thee->partLower[0] = xmin;
    thee->partLower[1] = ymin;
    thee->partLower[2] = zmin;
    thee->partUpper[0] = xmax;
    thee->partUpper[1] = ymax;
    thee->partUpper[2] = zmax;

    /* Flag the atoms as one of the following:
     *   - Belonging to this partition (partID = 1)
     *   - Not belonging to this partition (partID = 0) */
    alist = Vpbe_getValist(thee->pbe);
    for (iatom=0; iatom<Valist_getNumberAtoms(alist); iatom++) {
        atom = Valist_getAtom(alist,iatom);
        position = Vatom_getPosition(atom);
        if ((position[0]>=xmin)&&(position[0]<xmax)
          &&(position[1]>=ymin)&&(position[1]<ymax)
          &&(position[2]>=zmin)&&(position[2]<zmax)) Vatom_setPartID(atom,1);
        else Vatom_setPartID(atom,0);
    }
}

#undef VPMGSMALL
