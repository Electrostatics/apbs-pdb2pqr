/**
 *  @file    vacc.c
 *  @ingroup Vacc
 *  @author  Nathan Baker
 *  @brief   Class Vacc methods
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
#include "apbs/vacc.h"

#if defined(HAVE_MC_H)
#include "mc/mc.h"
#endif

VEMBED(rcsid="$Id$")

#if !defined(VINLINE_VACC)

VPUBLIC unsigned long int Vacc_memChk(Vacc *thee) {
    if (thee == VNULL) return 0;
    return Vmem_bytes(thee->vmem);
}

#endif /* if !defined(VINLINE_VACC) */

/**
 * @brief  Determines if a point is within the union of the spheres centered
 *         at the atomic centers with radii equal to the sum of their van der
 *         Waals radii and the probe radius.  Does not include contributions
 *         from the specified atom.
 * @returns 1 if accessible (outside the inflated van der Waals radius), 0
 *          otherwise
 * @author  Nathan Baker
 */
VPRIVATE int ivdwAccExclus(
        Vacc *thee,  /** Accessibility object */
        double center[3],  /** Position to test */
        double radius,  /** Radius of probe */ 
        int atomID  /** ID of atom to ignore */
        ) {

    int iatom;
    double dist2, *apos;
    Vatom *atom;
    VclistCell *cell;

    VASSERT(thee != VNULL);

    /* We can only test probes with radii less than the max specified */
    if (radius > Vclist_maxRadius(thee->clist)) {
        Vnm_print(2, 
            "Vacc_ivdwAcc: got radius (%g) bigger than max radius (%g)\n", 
            radius, thee->max_radius);
         VASSERT(0);
    }

    /* Get the relevant cell from the cell list */
    cell = Vclist_getCell(thee->clist, center);

    /* If we have no cell, then no atoms are nearby and we're definitely
     * accessible */
    if (cell == VNULL) {
        return 1.0;
    }

    /* Otherwise, check for overlap with the atoms in the cell */
    for (iatom=0; iatom<cell->natoms; iatom++) {
        atom = cell->atoms[iatom];
        if (atom->id != atomID) {
            apos = Vatom_getPosition(atom);
            dist2 = VSQR(center[0]-apos[0]) + VSQR(center[1]-apos[1]) 
                + VSQR(center[2]-apos[2]); 
            if (dist2 < VSQR(Vatom_getRadius(atom)+radius)) return 0.0;
        }
    }

    /* If we're still here, then the point is accessible */
    return 1.0;

}


VPUBLIC Vacc* Vacc_ctor(Valist *alist, Vclist *clist, int nsphereSurf) {


    Vacc *thee = VNULL;

    /* Set up the structure */
    thee = Vmem_malloc(VNULL, 1, sizeof(Vacc) );
    VASSERT( thee != VNULL);
    VASSERT( Vacc_ctor2(thee, alist, clist, nsphereSurf));
    return thee;
}

/** Check and store parameters passed to constructor */
VPRIVATE int Vacc_storeParms(Vacc *thee, Valist *alist, Vclist *clist,
        int nsphereSurf) {

    if (alist == VNULL) {
        Vnm_print(2, "Vacc_storeParms:  Got NULL Valist!\n");
        return 0;
    } else thee->alist = alist;
    if (clist == VNULL) {
        Vnm_print(2, "Vacc_storeParms:  Got NULL Vclist!\n");
        return 0;
    } else thee->clist = clist;

    thee->nsphereSurf = nsphereSurf;
    Vnm_print(0, "Vacc_ctor2:  Using %d-point probe sphere\n", nsphereSurf);

    return 1;
}

/** Allocate (and clear) space for storage */
VPRIVATE int Vacc_allocate(Vacc *thee) {

    int i, natoms;

    natoms = Valist_getNumberAtoms(thee->alist);

    thee->atomFlags = Vmem_malloc(thee->vmem, natoms, sizeof(int));
    if (thee->atomFlags == VNULL) {
        Vnm_print(2, 
                "Vacc_ctor2:  Failed to allocate %d (int)s for atomFlags!\n", 
                natoms);
        return 0;
    }
    for (i=0; i<natoms; i++) (thee->atomFlags)[i] = 0;

    thee->area = Vmem_malloc(thee->vmem, natoms, sizeof(double));
    if (thee->area == VNULL) {
        Vnm_print(2, 
                "Vacc_ctor2:  Failed to allocate %d (double)s for area!\n", 
                natoms);
        return 0;
    }
    for (i=0; i<natoms; i++) thee->area[i] = 0;

    return 1;
}


VPUBLIC int Vacc_ctor2(Vacc *thee, Valist *alist, Vclist *clist,
    int nsphereSurf) {

    /* Check and store parameters */
    if (!Vacc_storeParms(thee, alist, clist, nsphereSurf)) {
        Vnm_print(2, "Vacc_ctor2:  parameter check failed!\n");
        return 0;
    }

    /* Set up memory management object */
    thee->vmem = Vmem_ctor("APBS::VACC");
    if (thee->vmem == VNULL) {
        Vnm_print(2, "Vacc_ctor2:  memory object setup failed!\n");
        return 0;
    }

    /* Setup and check probe */
    thee->sphereSurf = Vacc_sphereSurf(thee, &(thee->nsphereSurf));
    if (thee->sphereSurf == VNULL) {
        Vnm_print(2, "Vacc_ctor2:  probe sphere setup failed!\n");
        return 0;
    }
 
    /* Allocate space */
    if (!Vacc_allocate(thee)) {
        Vnm_print(2, "Vacc_ctor2:  memory allocation failed!\n");
        return 0;
    }

    return 1;
}


VPUBLIC void Vacc_dtor(Vacc **thee) {
    
    if ((*thee) != VNULL) {
        Vacc_dtor2(*thee);
        Vmem_free(VNULL, 1, sizeof(Vacc), (void **)thee);
        (*thee) = VNULL;
    }

}

VPUBLIC void Vacc_dtor2(Vacc *thee) {

    int i, natoms;

    natoms = Valist_getNumberAtoms(thee->alist);
    Vmem_free(thee->vmem, natoms, sizeof(double), (void **)&(thee->area));
    Vmem_free(thee->vmem, natoms, sizeof(int), (void **)&(thee->atomFlags));

    for (i=0; i<thee->nsphereSurf; i++) {
        Vmem_free(thee->vmem, 3, sizeof(double), 
                (void **)&(thee->sphereSurf[i]));
    }
    Vmem_free(thee->vmem, thee->nsphereSurf, sizeof(double *), 
      (void **)&(thee->sphereSurf));

    Vmem_dtor(&(thee->vmem));
}

VPUBLIC double Vacc_vdwAcc(Vacc *thee, double center[3]) {

    VclistCell *cell;
    Vatom *atom;
    int iatom;
    double *apos;
    double dist2;

    /* Get the relevant cell from the cell list */
    cell = Vclist_getCell(thee->clist, center);

    /* If we have no cell, then no atoms are nearby and we're definitely
     * accessible */
    if (cell == VNULL) return 1.0;

    /* Otherwise, check for overlap with the atoms in the cell */
    for (iatom=0; iatom<cell->natoms; iatom++) {
        atom = cell->atoms[iatom];
        apos = Vatom_getPosition(atom);
        dist2 = VSQR(center[0]-apos[0]) + VSQR(center[1]-apos[1])
               + VSQR(center[2]-apos[2]);
        if (dist2 < VSQR(Vatom_getRadius(atom))) return 0.0;
    }

    /* If we're still here, then the point is accessible */
    return 1.0;
}

VPUBLIC double Vacc_ivdwAcc(Vacc *thee, double center[3], double radius) {

    return (double)ivdwAccExclus(thee, center, radius, -1);

}

VPUBLIC void Vacc_splineAccGradAtom(Vacc *thee, double center[VAPBS_DIM], 
        double win, double infrad, Vatom *atom, double *grad) {

    int i;
    double dist, *apos, arad, sm, sm2, w2i, w3i, mygrad;
    double mychi = 1.0;           /* Char. func. value for given atom */

    VASSERT(thee != NULL);

    /* Inverse squared window parameter */
    w2i = 1.0/(win*win);
    w3i = 1.0/(win*win*win);

    /* The grad is zero by default */
    for (i=0; i<VAPBS_DIM; i++) grad[i] = 0.0;

    /* *** CALCULATE THE CHARACTERISTIC FUNCTION VALUE FOR THIS ATOM AND THE
     * *** MAGNITUDE OF THE FORCE *** */
    apos = Vatom_getPosition(atom);
    /* Zero-radius atoms don't contribute */
    if (Vatom_getRadius(atom) > 0.0) {
        arad = Vatom_getRadius(atom) + infrad;
        dist = VSQRT(VSQR(apos[0]-center[0]) + VSQR(apos[1]-center[1])
          + VSQR(apos[2]-center[2]));
        /* If we're inside an atom, the entire characteristic function
         * will be zero and the grad will be zero, so we can stop */
        if (dist < (arad - win)) return;
        /* Likewise, if we're outside the smoothing window, the characteristic
         * function is unity and the grad will be zero, so we can stop */
        else if (dist > (arad + win)) return;
        /* Account for floating point error at the border 
         * NAB:  COULDN'T THESE TESTS BE COMBINED AS BELOW
         * (Vacc_splineAccAtom)? */
        else if ((VABS(dist - (arad - win)) < VSMALL) || 
                 (VABS(dist - (arad + win)) < VSMALL)) return;
        /* If we're inside the smoothing window */
        else {
            sm = dist - arad + win;
            sm2 = VSQR(sm);
            mychi = 0.75*sm2*w2i -0.25*sm*sm2*w3i;
            mygrad = 1.5*sm*w2i - 0.75*sm2*w3i;
        }
        /* Now assemble the grad vector */
        VASSERT(mychi > 0.0);
        for (i=0; i<VAPBS_DIM; i++) 
            grad[i] = -(mygrad/mychi)*((center[i] - apos[i])/dist);
    }    
}

VPUBLIC double Vacc_splineAccAtom(Vacc *thee, double center[VAPBS_DIM], 
        double win, double infrad, Vatom *atom) {

    double dist, *apos, arad, sm, sm2, w2i, w3i, value, stot, sctot;

    VASSERT(thee != NULL);

    /* Inverse squared window parameter */
    w2i = 1.0/(win*win);
    w3i = 1.0/(win*win*win);

    apos = Vatom_getPosition(atom);
    /* Zero-radius atoms don't contribute */
    if (Vatom_getRadius(atom) > 0.0) {
        arad = Vatom_getRadius(atom) + infrad;
        stot = arad + win;
        sctot = VMAX2(0, (arad - win));
        dist = VSQRT(VSQR(apos[0]-center[0]) + VSQR(apos[1]-center[1])
          + VSQR(apos[2]-center[2]));
        /* If we're inside an atom, the entire characteristic function
         * will be zero */
        if ((dist < sctot) || (VABS(dist - sctot) < VSMALL)){
            value = 0.0;
        /* We're outside the smoothing window */
        } else if ((dist > stot) || (VABS(dist - stot) < VSMALL)) {
            value = 1.0;
        /* We're inside the smoothing window */
        } else {
            sm = dist - arad + win;
            sm2 = VSQR(sm);
            value = 0.75*sm2*w2i - 0.25*sm*sm2*w3i;
        }
    } else value = 1.0;
 
    return value;
}

/** 
 * @brief  Fast spline-based surface computation subroutine
 * @returns  Spline value
 * @author  Todd Dolinsky and Nathan Baker
 */
VPRIVATE double splineAcc(
        Vacc *thee,  /** Accessibility object */
        double center[VAPBS_DIM],  /** Point at which the acc is to be
                                    * evaluated */
        double win,  /** Spline window */
        double infrad,  /** Radius to inflate atomic radius */
        VclistCell *cell  /** Cell of atom objects */
        ) {

    int atomID, iatom;      
    Vatom *atom;
    double value = 1.0;          

    VASSERT(thee != NULL);

    /* Now loop through the atoms assembling the characteristic function */
    for (iatom=0; iatom<cell->natoms; iatom++) {

        atom = cell->atoms[iatom];
        atomID = atom->id;

        /* Check to see if we've counted this atom already */
        if ( !(thee->atomFlags[atomID]) ) {

            thee->atomFlags[atomID] = 1;
            value *= Vacc_splineAccAtom(thee, center, win, infrad, atom);
            
            if (value < VSMALL) return value;
        } 
    }
 
    return value;
}


VPUBLIC double Vacc_splineAcc(Vacc *thee, double center[VAPBS_DIM], double win, 
  double infrad) {

    VclistCell *cell;
    Vatom *atom;
    int iatom, atomID;      
    double value = 1.0;            


    VASSERT(thee != NULL);

    if (Vclist_maxRadius(thee->clist) < (win + infrad)) {
        Vnm_print(2, "Vacc_splineAcc:  Vclist has max_radius=%g;\n", 
                Vclist_maxRadius(thee->clist));
        Vnm_print(2, "Vacc_splineAcc:  Insufficient for win=%g, infrad=%g\n", 
                win, infrad);
        VASSERT(0);
    }

    /* Get a cell or VNULL; in the latter case return 1.0 */
    cell = Vclist_getCell(thee->clist, center);
    if (cell == VNULL) return 1.0;

    /* First, reset the list of atom flags 
     * NAB:  THIS SEEMS VERY INEFFICIENT */
    for (iatom=0; iatom<cell->natoms; iatom++) {
        atom = cell->atoms[iatom];
        atomID = atom->id; 
        thee->atomFlags[atomID] = 0;
    }

    return splineAcc(thee, center, win, infrad, cell);
}

VPUBLIC void Vacc_splineAccGrad(Vacc *thee, double center[VAPBS_DIM], 
        double win, double infrad, Vatom *atom, double *grad) {

    int iatom, i, atom2ID, atomID;
    double value = 1.0;            
    VclistCell *cell;
    Vatom *atom2;

    VASSERT(thee != NULL);

    if (thee->max_radius < (win + infrad)) {
        Vnm_print(2, "Vacc_splineAccGrad: Vclist max_radius=%g;\n", 
                thee->max_radius);
        Vnm_print(2, "Vacc_splineAccGrad: Insufficient for win=%g, infrad=%g\n", 
                win, infrad);
        VASSERT(0);
    }

    /* Reset the gradient */
    for (i=0; i<VAPBS_DIM; i++) grad[i] = 0.0;

    /* Get the cell; check for nullity */
    cell = Vclist_getCell(thee->clist, center);
    if (cell == VNULL) return;

    /* First, reset the list of atom flags for all atoms except the one of
     * interest.
     * NAB:  THIS SEEMS VERY INEFFICIENT */
    atomID = atom->id;
    for (iatom=0; iatom<cell->natoms; iatom++) {
        atom2 = cell->atoms[iatom];
        atom2ID = atom2->id;
        thee->atomFlags[atom2ID] = 0;
    }
    thee->atomFlags[atomID] = 1;

    value = splineAcc(thee, center, win, infrad, cell);

    if (value < VSMALL) {
        for (i=0; i<VAPBS_DIM; i++) grad[i] = 0.0;
    } else {
        Vacc_splineAccGradAtom(thee, center, win, infrad, atom, grad);
        for (i=0; i<VAPBS_DIM; i++) grad[i] *= value;
    }
}

VPUBLIC double Vacc_molAcc(Vacc *thee, double center[VAPBS_DIM], 
        double radius) {

    int ipt, i;
    double vec[VAPBS_DIM];

    /* ******* CHECK IF OUTSIDE ATOM+PROBE RADIUS SURFACE ***** */
    if (Vacc_ivdwAcc(thee, center, radius) == 1.0) return 1;

    /* ******* CHECK IF INSIDE ATOM RADIUS SURFACE ***** */
    if (Vacc_vdwAcc(thee, center) == 0.0) return 0;

    /* ******* CHECK IF OUTSIDE MOLECULAR SURFACE ***** */
    /* Let S be the sphere of radius radius centered at the point we are
     * testing.  We are outside the molecular surface if there is a point on
     * the surface of S that is outside the atom+probe radius surface */
    /* THIS IS INCORRECT:  the correct behavior should check all points within
     * S; not just one the surface.  Alternatively, we could check to see if
     * the point of interest is within a sphere radius of SAS.  THIS IS
     * AN OUTSTANDING BUG IN THE CODE!! 
     * (Thanks to John Mongan and Jessica Swanson for finding this) */
    VASSERT(thee->sphereSurf != VNULL);
    for (ipt=0; ipt<thee->nsphereSurf; ipt++) {
        for (i=0; i<VAPBS_DIM; i++) 
            vec[i] = radius*thee->sphereSurf[ipt][i] + center[i];
        if (Vacc_ivdwAcc(thee, vec, radius)) return 1.0;
    }

    /* If all else failed, we are not inside the molecular surface */
    return 0.0;
}

VPUBLIC double Vacc_fastMolAcc(Vacc *thee, double center[VAPBS_DIM], 
        double radius) {

    int ipt, i;
    double vec[VAPBS_DIM];

    /* ******* CHECK IF OUTSIDE MOLECULAR SURFACE ***** */
    /* Let S be the sphere of radius radius centered at the point we are
     * testing.  We are outside the molecular surface if there is a point on
     * the surface of S that is outside the atom+probe radius surface */
    /* THIS IS INCORRECT:  the correct behavior should check all points within
     * S; not just one the surface.  Alternatively, we could check to see if
     * the point of interest is within a sphere radius of SAS.  THIS IS
     * AN OUTSTANDING BUG IN THE CODE!! 
     * (Thanks to John Mongan and Jessica Swanson for finding this) */
    VASSERT(thee->sphereSurf != VNULL);
    for (ipt=0; ipt<thee->nsphereSurf; ipt++) {
        for (i=0; i<VAPBS_DIM; i++) 
            vec[i] = radius*thee->sphereSurf[ipt][i] + center[i];
        if (Vacc_ivdwAcc(thee, vec, radius)) return 1.0;
    }

    /* If all else failed, we are not inside the molecular surface */
    return 0.0;
}


#if defined(HAVE_MC_H)
VPUBLIC void Vacc_writeGMV(Vacc *thee, double radius, int meth, Gem *gm, 
  char *iodev, char *iofmt, char *iohost, char *iofile) {

    double *accVals[MAXV], coord[3];
    Vio *sock;
    int ivert, icoord;

    for (ivert=0; ivert<MAXV; ivert++) accVals[ivert] = VNULL;
    accVals[0] = (void *)Vmem_malloc(thee->vmem, Gem_numVV(gm), sizeof(double));
    accVals[1] = (void *)Vmem_malloc(thee->vmem, Gem_numVV(gm), sizeof(double));
    for (ivert=0; ivert<Gem_numVV(gm); ivert++) {
        for (icoord=0;icoord<3;icoord++) 
          coord[icoord] = VV_coord(Gem_VV(gm, ivert), icoord);
        if (meth == 0) {
            accVals[0][ivert] = Vacc_molAcc(thee, coord, radius);
            accVals[1][ivert] = Vacc_molAcc(thee, coord, radius);
        } else if (meth == 1) {
            accVals[0][ivert] = Vacc_ivdwAcc(thee, coord, radius);
            accVals[1][ivert] = Vacc_ivdwAcc(thee, coord, radius);
        } else if (meth == 2) {
            accVals[0][ivert] = Vacc_vdwAcc(thee, coord);
            accVals[1][ivert] = Vacc_vdwAcc(thee, coord);
        } else VASSERT(0);
    }
    sock = Vio_ctor(iodev, iofmt, iohost, iofile, "w");
    Gem_writeGMV(gm, sock, 1, accVals);
    Vio_dtor(&sock);
    Vmem_free(thee->vmem, Gem_numVV(gm), sizeof(double), 
      (void **)&(accVals[0]));
    Vmem_free(thee->vmem, Gem_numVV(gm), sizeof(double), 
      (void **)&(accVals[1]));
}
#endif /* defined(HAVE_MC_H) */

VPUBLIC double** Vacc_sphereSurf(Vacc *thee, int *npts) {

    double **points = VNULL;
    int nactual, i, itheta, ntheta, iphi, nphimax, nphi;
    double frac;
    double sintheta, costheta, theta, dtheta;
    double sinphi, cosphi, phi, dphi;

    frac = ((double)(*npts))/4.0;
    ntheta = VRINT(VSQRT(Vunit_pi*frac));
    dtheta = Vunit_pi/((double)(ntheta));
    nphimax = 2*ntheta;

    /* COUNT THE ACTUAL NUMBER OF POINTS TO BE USED */
    nactual = 0;
    for (itheta=0; itheta<ntheta; itheta++) {
        theta = dtheta*((double)(itheta));
        sintheta = VSIN(theta);
        costheta = VCOS(theta);
        nphi = VRINT(sintheta*nphimax);
        nactual += nphi;
    }

    /* ALLOCATE THE SPACE FOR THE POINTS */
    points = Vmem_malloc(thee->vmem, nactual, sizeof(double *));
    VASSERT(points != VNULL);
    for (i=0; i<nactual; i++) {
        points[i] = Vmem_malloc(thee->vmem, 3, sizeof(double));
        VASSERT(points[i] != VNULL);
    }

    /* ASSIGN THE POINTS */
    nactual = 0;
    for (itheta=0; itheta<ntheta; itheta++) {
        theta = dtheta*((double)(itheta));
        sintheta = VSIN(theta);
        costheta = VCOS(theta);
        nphi = VRINT(sintheta*nphimax);
        if (nphi != 0) {
            dphi = 2*Vunit_pi/((double)(nphi));
            for (iphi=0; iphi<nphi; iphi++) {
                phi = dphi*((double)(iphi));
                sinphi = VSIN(phi);
                cosphi = VCOS(phi);
                points[nactual][0] = cosphi * sintheta;
                points[nactual][1] = sinphi * sintheta;
                points[nactual][2] = costheta;
                nactual++;
            }
        }
    }

    *npts = nactual;
    return points;
}

VPUBLIC double Vacc_totalSASA(Vacc *thee, double radius) { 

    int i;
    double area = 0.0;
    Vatom *atom;

    for (i=0; i<Valist_getNumberAtoms(thee->alist); i++) {
        atom = Valist_getAtom(thee->alist, i);
        thee->area[i] = Vacc_atomSASA(thee, radius, atom);
        area += (thee->area[i]);
    }

    return area;

}

VPUBLIC double Vacc_atomSASA(Vacc *thee, double srad, 
        Vatom *thisAtom) { 

    int ipt, i;
    double area = 0.0;
    double *tPos, tRad, vec[VAPBS_DIM];

    /* Get the atom information */
    tPos = Vatom_getPosition(thisAtom);
    tRad = Vatom_getRadius(thisAtom);

    for (ipt=0; ipt<thee->nsphereSurf; ipt++) {
        for (i=0; i<VAPBS_DIM; i++) 
            vec[i] = (tRad+srad)*thee->sphereSurf[ipt][i] + tPos[i];
        if (ivdwAccExclus(thee, vec, srad, thisAtom->id)) area += 1.0;
    }

    /* We will return UHBD's asas2: probe-centered solvent-accessible surface
     * area */
    area = area/((double)(thee->nsphereSurf))*4.0*VPI*(tRad+srad)*(tRad+srad);

    return area;

}
