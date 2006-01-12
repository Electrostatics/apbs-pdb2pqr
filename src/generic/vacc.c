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
 * Copyright (c) 2002-2006.  Washington University in St. Louis.
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
 * code included in releases of ISIM, PMV, PyMOL, SMOL, VMD, and Vision.
 * Such combined software may be linked with APBS and redistributed together 
 * in original or modified form as mere aggregation without requirement that 
 * the entire work be under the scope of the GNU General Public License.
 * This special exception permission is also extended to any software listed
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
#include "apbs/vacc.h"

#if defined(HAVE_MC_H)
#include "mc/mc.h"
#endif

VEMBED(rcsid="$Id$")

#if !defined(VINLINE_VACC)

VPUBLIC unsigned long int Vacc_memChk(Vacc *thee) {
    if (thee == VNULL) return 0;
    return Vmem_bytes(thee->mem);
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
            radius, Vclist_maxRadius(thee->clist));
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


VPUBLIC Vacc* Vacc_ctor(Valist *alist, Vclist *clist, double surf_density) {


    Vacc *thee = VNULL;

    /* Set up the structure */
    thee = Vmem_malloc(VNULL, 1, sizeof(Vacc) );
    VASSERT( thee != VNULL);
    VASSERT( Vacc_ctor2(thee, alist, clist, surf_density));
    return thee;
}

/** Check and store parameters passed to constructor */
VPRIVATE int Vacc_storeParms(Vacc *thee, Valist *alist, Vclist *clist,
        double surf_density) {

    int nsphere, iatom;
    double maxrad, maxarea, rad;
    Vatom *atom;

    if (alist == VNULL) {
        Vnm_print(2, "Vacc_storeParms:  Got NULL Valist!\n");
        return 0;
    } else thee->alist = alist;
    if (clist == VNULL) {
        Vnm_print(2, "Vacc_storeParms:  Got NULL Vclist!\n");
        return 0;
    } else thee->clist = clist;
    thee->surf_density = surf_density;

    /* Loop through the atoms to determine the maximum radius */
    maxrad = 0.0;
    for (iatom=0; iatom<Valist_getNumberAtoms(alist); iatom++) {
        atom = Valist_getAtom(alist, iatom);
        rad = Vatom_getRadius(atom);
        if (rad > maxrad) maxrad = rad;
    }
    maxrad = maxrad + Vclist_maxRadius(thee->clist);

    maxarea = 4.0*VPI*maxrad*maxrad;
    nsphere = (int)ceil(maxarea*surf_density);

    Vnm_print(0, "Vacc_storeParms:  Surf. density = %g\n", surf_density);
    Vnm_print(0, "Vacc_storeParms:  Max area = %g\n", maxarea);
    thee->refSphere = VaccSurf_refSphere(thee->mem, nsphere);
    Vnm_print(0, "Vacc_storeParms:  Using %d-point reference sphere\n", 
            thee->refSphere->npts);

    return 1;
}

/** Allocate (and clear) space for storage */
VPRIVATE int Vacc_allocate(Vacc *thee) {

    int i, natoms;

    natoms = Valist_getNumberAtoms(thee->alist);

    thee->atomFlags = Vmem_malloc(thee->mem, natoms, sizeof(int));
    if (thee->atomFlags == VNULL) {
        Vnm_print(2, 
               "Vacc_allocate:  Failed to allocate %d (int)s for atomFlags!\n", 
                natoms);
        return 0;
    }
    for (i=0; i<natoms; i++) (thee->atomFlags)[i] = 0;

    return 1;
}


VPUBLIC int Vacc_ctor2(Vacc *thee, Valist *alist, Vclist *clist,
    double surf_density) {

    /* Check and store parameters */
    if (!Vacc_storeParms(thee, alist, clist, surf_density)) {
        Vnm_print(2, "Vacc_ctor2:  parameter check failed!\n");
        return 0;
    }

    /* Set up memory management object */
    thee->mem = Vmem_ctor("APBS::VACC");
    if (thee->mem == VNULL) {
        Vnm_print(2, "Vacc_ctor2:  memory object setup failed!\n");
        return 0;
    }

    /* Setup and check probe */
    thee->surf = VNULL;
 
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
    Vmem_free(thee->mem, natoms, sizeof(int), (void **)&(thee->atomFlags));

    if (thee->refSphere != VNULL) {
        VaccSurf_dtor(&(thee->refSphere));
        thee->refSphere = VNULL;
    }
    if (thee->surf != VNULL) {
        for (i=0; i<natoms; i++) VaccSurf_dtor(&(thee->surf[i]));
        Vmem_free(thee->mem, natoms, sizeof(VaccSurf *), 
                (void **)&(thee->surf));
        thee->surf = VNULL;
    }

    Vmem_dtor(&(thee->mem));
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

VPUBLIC void Vacc_splineAccGradAtomNorm(Vacc *thee, double center[VAPBS_DIM], 
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

VPUBLIC void Vacc_splineAccGradAtomUnnorm(Vacc *thee, double center[VAPBS_DIM], 
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
            grad[i] = -(mygrad)*((center[i] - apos[i])/dist);
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
        double win, double infrad, double *grad) {

    int iatom, i, atomID;
    double acc = 1.0;            
    double tgrad[VAPBS_DIM];
    VclistCell *cell;
    Vatom *atom = VNULL;

    VASSERT(thee != NULL);

    if (Vclist_maxRadius(thee->clist) < (win + infrad)) {
        Vnm_print(2, "Vacc_splineAccGrad: Vclist max_radius=%g;\n", 
                Vclist_maxRadius(thee->clist));
        Vnm_print(2, "Vacc_splineAccGrad: Insufficient for win=%g, infrad=%g\n", 
                win, infrad);
        VASSERT(0);
    }

    /* Reset the gradient */
    for (i=0; i<VAPBS_DIM; i++) grad[i] = 0.0;

    /* Get the cell; check for nullity */
    cell = Vclist_getCell(thee->clist, center);
    if (cell == VNULL) return;

    /* Reset the list of atom flags */
    for (iatom=0; iatom<cell->natoms; iatom++) {
        atom = cell->atoms[iatom];
        atomID = atom->id;
        thee->atomFlags[atomID] = 0;
    }

    /* Get the local accessibility */
    acc = splineAcc(thee, center, win, infrad, cell);

    /* Accumulate the gradient of all local atoms */
    if (acc > VSMALL) {
        for (iatom=0; iatom<cell->natoms; iatom++) {
            atom = cell->atoms[iatom];
            Vacc_splineAccGradAtomNorm(thee, center, win, infrad, atom, tgrad);
        }
        for (i=0; i<VAPBS_DIM; i++) grad[i] += tgrad[i];
    }
    for (i=0; i<VAPBS_DIM; i++) grad[i] *= -acc;
}

VPUBLIC double Vacc_molAcc(Vacc *thee, double center[VAPBS_DIM], 
        double radius) {
    
    double rc;

    /* ******* CHECK IF OUTSIDE ATOM+PROBE RADIUS SURFACE ***** */
    if (Vacc_ivdwAcc(thee, center, radius) == 1.0) {
       
        /* Vnm_print(2, "DEBUG:  ivdwAcc = 1.0\n"); */
        rc = 1.0;

    /* ******* CHECK IF INSIDE ATOM RADIUS SURFACE ***** */
    } else if (Vacc_vdwAcc(thee, center) == 0.0) {
       
        /* Vnm_print(2, "DEBUG:  vdwAcc = 0.0\n"); */
        rc = 0.0;

    /* ******* CHECK IF OUTSIDE MOLECULAR SURFACE ***** */
    } else {

        /* Vnm_print(2, "DEBUG:  calling fastMolAcc...\n"); */
        rc = Vacc_fastMolAcc(thee, center, radius);

    }

    return rc;

}

VPUBLIC double Vacc_fastMolAcc(Vacc *thee, double center[VAPBS_DIM], 
        double radius) {

    Vatom *atom;
    VaccSurf *surf;
    VclistCell *cell;
    int ipt, iatom, atomID;
    double *apos, dist2, rad2;

    rad2 = radius*radius;

    /* Check to see if the SAS has been defined */
    if (thee->surf == VNULL) Vacc_SASA(thee, radius);

    /* Get the cell associated with this point */
    cell = Vclist_getCell(thee->clist, center);
    if (cell == VNULL) {
        Vnm_print(2, "Vacc_fastMolAcc:  unexpected VNULL VclistCell!\n");
        return 1.0;
    }

    /* Loop through all the atoms in the cell */
    for (iatom=0; iatom<cell->natoms; iatom++) {
        atom = cell->atoms[iatom];
        atomID = Vatom_getAtomID(atom);
        surf = thee->surf[atomID];
        /* Loop through all SAS points associated with this atom */
        for (ipt=0; ipt<surf->npts; ipt++) {
            /* See if we're within a probe radius of the point */
            dist2 = VSQR(center[0]-(surf->xpts[ipt])) 
                + VSQR(center[1]-(surf->ypts[ipt])) 
                + VSQR(center[2]-(surf->zpts[ipt]));
            if (dist2 < rad2) return 1.0;
        }
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
    accVals[0] = (void *)Vmem_malloc(thee->mem, Gem_numVV(gm), sizeof(double));
    accVals[1] = (void *)Vmem_malloc(thee->mem, Gem_numVV(gm), sizeof(double));
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
    Vmem_free(thee->mem, Gem_numVV(gm), sizeof(double), 
      (void **)&(accVals[0]));
    Vmem_free(thee->mem, Gem_numVV(gm), sizeof(double), 
      (void **)&(accVals[1]));
}
#endif /* defined(HAVE_MC_H) */

VPUBLIC double Vacc_SASA(Vacc *thee, double radius) { 

    int i, natom;
    double area, *apos;
    Vatom *atom;
    VaccSurf *asurf;

    natom = Valist_getNumberAtoms(thee->alist);

    /* Check to see if we need to build the surface */
    if (thee->surf == VNULL) {
        thee->surf = Vmem_malloc(thee->mem, natom, sizeof(VaccSurf *));
        for (i=0; i<natom; i++) {
            atom = Valist_getAtom(thee->alist, i);
            apos = Vatom_getPosition(atom);
            /* NOTE:  RIGHT NOW WE DO THIS FOR THE ENTIRE MOLECULE WHICH IS
             * INCREDIBLY INEFFICIENT, PARTICULARLY DURING FOCUSING!!! */
            thee->surf[i] = Vacc_atomSurf(thee, atom, thee->refSphere, 
                    radius);
        }
    }

    /* Calculate the area */
    area = 0.0;
    for (i=0; i<natom; i++) {
        atom = Valist_getAtom(thee->alist, i);
        asurf = thee->surf[i];
        /* See if this surface needs to be rebuilt */
        if (asurf->probe_radius != radius) {
            Vnm_print(2, "Vacc_SASA:  Warning -- probe radius changed from %g to %g!\n", 
                    asurf->probe_radius, radius);
            VaccSurf_dtor2(asurf);
            thee->surf[i] = Vacc_atomSurf(thee, atom, thee->refSphere, radius);
            asurf = thee->surf[i];
        }
        area += (asurf->area);
    }

    return area;

}

VPUBLIC double Vacc_totalSASA(Vacc *thee, double radius) {

    return Vacc_SASA(thee, radius);

}

VPUBLIC double Vacc_atomSASA(Vacc *thee, double radius, Vatom *atom) {

    VaccSurf *asurf;
    int id;

    if (thee->surf == VNULL) Vacc_SASA(thee, radius);

    id = Vatom_getAtomID(atom);
    asurf = thee->surf[id];

    /* See if this surface needs to be rebuilt */
    if (asurf->probe_radius != radius) {
        Vnm_print(2, "Vacc_SASA:  Warning -- probe radius changed from %g to %g!\n", 
                asurf->probe_radius, radius);
        VaccSurf_dtor2(asurf);
        thee->surf[id] = Vacc_atomSurf(thee, atom, thee->refSphere, radius);
        asurf = thee->surf[id];
    }

    return asurf->area;

}

VPUBLIC VaccSurf* VaccSurf_ctor(Vmem *mem, double probe_radius, int nsphere) {
    VaccSurf *thee;

    thee = Vmem_malloc(mem, 1, sizeof(Vacc) );
    VASSERT( VaccSurf_ctor2(thee, mem, probe_radius, nsphere) );

    return thee;
}

VPUBLIC int VaccSurf_ctor2(VaccSurf *thee, Vmem *mem, double probe_radius,
        int nsphere) {

    if (thee == VNULL) return 0;

    thee->mem = mem;
    thee->npts = nsphere;
    thee->probe_radius = probe_radius;
    thee->area = 0.0;

    if (thee->npts > 0) {
        thee->xpts = Vmem_malloc(thee->mem, thee->npts, sizeof(double));
        thee->ypts = Vmem_malloc(thee->mem, thee->npts, sizeof(double));
        thee->zpts = Vmem_malloc(thee->mem, thee->npts, sizeof(double));
        thee->bpts = Vmem_malloc(thee->mem, thee->npts, sizeof(int));
    } else {
        thee->xpts = VNULL;
        thee->ypts = VNULL;
        thee->zpts = VNULL;
        thee->bpts = VNULL;
    } 

    return 1;
}

VPUBLIC void VaccSurf_dtor(VaccSurf **thee) {

    Vmem *mem;

    if ((*thee) != VNULL) {
        mem = (*thee)->mem;
        VaccSurf_dtor2(*thee);
        Vmem_free(mem, 1, sizeof(VaccSurf), (void **)thee);
        (*thee) = VNULL;
    }

}

VPUBLIC void VaccSurf_dtor2(VaccSurf *thee) {

    if (thee->npts > 0) {
        Vmem_free(thee->mem, thee->npts, sizeof(double), 
                (void **)&(thee->xpts));
        Vmem_free(thee->mem, thee->npts, sizeof(double), 
                (void **)&(thee->ypts));
        Vmem_free(thee->mem, thee->npts, sizeof(double), 
                (void **)&(thee->zpts));
        Vmem_free(thee->mem, thee->npts, sizeof(int), 
                (void **)&(thee->bpts));
    }
}

VPUBLIC VaccSurf* Vacc_atomSurf(Vacc *thee, Vatom *atom, 
        VaccSurf *ref, double prad) {

    VaccSurf *surf;
    int i, j, npts, atomID;
    double arad, rad, pos[3], *apos;

    /* Get atom information */
    arad = Vatom_getRadius(atom);
    apos = Vatom_getPosition(atom);
    atomID = Vatom_getAtomID(atom);

    if (arad < VSMALL) {
        return VaccSurf_ctor(thee->mem, prad, 0);
    }

    rad = arad + prad;

    /* Determine which points will contribute */
    npts = 0;
    for (i=0; i<ref->npts; i++) {
        /* Reset point flag: zero-radius atoms do not contribute */
        pos[0] = rad*(ref->xpts[i]) + apos[0];
        pos[1] = rad*(ref->ypts[i]) + apos[1];
        pos[2] = rad*(ref->zpts[i]) + apos[2];
        if (ivdwAccExclus(thee, pos, prad, atomID)) {
            npts++;
            ref->bpts[i] = 1;
        } else {
            ref->bpts[i] = 0;
        }
    }

    /* Allocate space for the points */
    surf = VaccSurf_ctor(thee->mem, prad, npts);

    /* Assign the points */
    j = 0;
    for (i=0; i<ref->npts; i++) {
        if (ref->bpts[i]) {
            surf->bpts[j] = 1;
            surf->xpts[j] = rad*(ref->xpts[i]) + apos[0];
            surf->ypts[j] = rad*(ref->ypts[i]) + apos[1];
            surf->zpts[j] = rad*(ref->zpts[i]) + apos[2];
            j++;
        }
    }

    /* Assign the area */
    surf->area = 4.0*VPI*rad*rad*((double)(surf->npts))/((double)(ref->npts));

    return surf;

}

VPUBLIC VaccSurf* VaccSurf_refSphere(Vmem *mem, int npts) {

    VaccSurf *surf;
    int nactual, i, itheta, ntheta, iphi, nphimax, nphi;
    double frac;
    double sintheta, costheta, theta, dtheta;
    double sinphi, cosphi, phi, dphi;

    /* Setup "constants" */
    frac = ((double)(npts))/4.0;
    ntheta = VRINT(VSQRT(Vunit_pi*frac));
    dtheta = Vunit_pi/((double)(ntheta));
    nphimax = 2*ntheta;

    /* Count the actual number of points to be used */
    nactual = 0;
    for (itheta=0; itheta<ntheta; itheta++) {
        theta = dtheta*((double)(itheta));
        sintheta = VSIN(theta);
        costheta = VCOS(theta);
        nphi = VRINT(sintheta*nphimax);
        nactual += nphi;
    }

    /* Allocate space for the points */
    surf = VaccSurf_ctor(mem, 1.0, nactual);

    /* Clear out the boolean array */
    for (i=0; i<nactual; i++) surf->bpts[i] = 1;

    /* Assign the points */
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
                surf->xpts[nactual] = cosphi * sintheta;
                surf->ypts[nactual] = sinphi * sintheta;
                surf->zpts[nactual] = costheta;
                nactual++;
            }
        }
    }

    surf->npts = nactual;

    return surf;
}

VPUBLIC VaccSurf* Vacc_atomSASPoints(Vacc *thee, double radius, 
        Vatom *atom) {

    VaccSurf *asurf = VNULL; 
    int id;

    if (thee->surf == VNULL) Vacc_SASA(thee, radius);
    id = Vatom_getAtomID(atom);

    asurf = thee->surf[id];

    /* See if this surface needs to be rebuilt */
    if (asurf->probe_radius != radius) {
        Vnm_print(2, "Vacc_SASA:  Warning -- probe radius changed from %g to %g!\n", 
                asurf->probe_radius, radius);
        VaccSurf_dtor2(asurf);
        thee->surf[id] = Vacc_atomSurf(thee, atom, thee->refSphere, radius);
        asurf = thee->surf[id];
    }

    return asurf;

}
