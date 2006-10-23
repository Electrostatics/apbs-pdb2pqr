/**
 *  @file    vclist.c
 *  @ingroup Vclist
 *  @author  Nathan Baker
 *  @brief   Class Vclist methods
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
#include "apbs/vclist.h"

#if defined(HAVE_MC_H)
#include "mc/mc.h"
#endif

VEMBED(rcsid="$Id$")

#if !defined(VINLINE_VCLIST)

VPUBLIC unsigned long int Vclist_memChk(Vclist *thee) {
    if (thee == VNULL) return 0;
    return Vmem_bytes(thee->vmem);
}

VPUBLIC double Vclist_maxRadius(Vclist *thee) {
    VASSERT(thee != VNULL);
    return thee->max_radius;
}

#endif /* if !defined(VINLINE_VCLIST) */

VPUBLIC Vclist* Vclist_ctor(Valist *alist, double max_radius, 
        int npts[VAPBS_DIM], Vclist_DomainMode mode, 
        double lower_corner[VAPBS_DIM],  double upper_corner[VAPBS_DIM]) {

    Vclist *thee = VNULL;

    /* Set up the structure */
    thee = Vmem_malloc(VNULL, 1, sizeof(Vclist) );
    VASSERT( thee != VNULL);
    VASSERT( Vclist_ctor2(thee, alist, max_radius, npts, mode, lower_corner,
                upper_corner) );
    return thee;
}

/** Get the dimensions of the molecule stored in thee->alist */
VPRIVATE void Vclist_getMolDims(
        Vclist *thee, 
        double lower_corner[VAPBS_DIM], /* Set to lower corner of molecule */
        double upper_corner[VAPBS_DIM], /* Set to lower corner of molecule */
        double *r_max /** Set to max atom radius */
        ) {

    int i, j;
    double pos;
    Valist *alist;
    Vatom *atom;

    alist = thee->alist;

    /* Initialize */
    for (i=0; i<VAPBS_DIM; i++) {
        lower_corner[i] = VLARGE;
        upper_corner[i] = -VLARGE;
    }
    *r_max = -1.0;

    /* Check each atom */
    for (i=0; i<Valist_getNumberAtoms(alist); i++) {
        atom = Valist_getAtom(alist, i);
		for (j=0; j<VAPBS_DIM; j++) {
            pos = (Vatom_getPosition(atom))[j];
            if ( pos < lower_corner[j] ) lower_corner[j] = pos;
            if ( pos > upper_corner[j] ) upper_corner[j] = pos;
        }
        if (Vatom_getRadius(atom) > *r_max) *r_max = Vatom_getRadius(atom);
    }

}

/** Setup lookup grid */
VPRIVATE int Vclist_setupGrid(Vclist *thee) {

    /* Inflation factor ~ sqrt(2)*/
    #define VCLIST_INFLATE 1.42

    int i;
    double length[VAPBS_DIM], r_max;

    /* Set up the grid corners */
    switch (thee->mode) {
        case CLIST_AUTO_DOMAIN:
            /* Get molecule dimensions */
            Vclist_getMolDims(thee, thee->lower_corner, thee->upper_corner, 
                    &r_max);
            /* Set up grid spacings */
            for (i=0; i<VAPBS_DIM; i++) {
                thee->upper_corner[i] = thee->upper_corner[i]
                    + VCLIST_INFLATE*(r_max+thee->max_radius);
                thee->lower_corner[i] = thee->lower_corner[i]
                    - VCLIST_INFLATE*(r_max+thee->max_radius);
            }
            break;
        case CLIST_MANUAL_DOMAIN:
            /* Grid corners established in constructor */
            break;
        default:
            Vnm_print(2, "Vclist_setupGrid:  invalid setup mode (%d)!\n",
                    thee->mode);
            return 0;
    }
   
    /* Set up the grid lengths and spacings */
    for (i=0; i<VAPBS_DIM; i++) { 
        length[i] = thee->upper_corner[i] - thee->lower_corner[i];
        thee->spacs[i] = length[i]/((double)(thee->npts[i] - 1));
    }
    Vnm_print(0, "Vclist_setupGrid:  Grid lengths = (%g, %g, %g)\n", 
            length[0], length[1], length[2]);

    Vnm_print(0, "Vclist_setupGrid:  Grid lower corner = (%g, %g, %g)\n",
      (thee->lower_corner)[0], (thee->lower_corner)[1],
      (thee->lower_corner)[2]);

    return 1;

    #undef VCLIST_INFLATE
}

/** Check and store parameters passed to constructor */
VPRIVATE int Vclist_storeParms(Vclist *thee, Valist *alist,
        double max_radius, int npts[VAPBS_DIM], Vclist_DomainMode mode,
        double lower_corner[VAPBS_DIM], double upper_corner[VAPBS_DIM] ) {

    int i = 0;

    if (alist == VNULL) {
        Vnm_print(2, "Vclist_ctor2:  Got NULL Valist!\n");
        return 0;
    } else thee->alist = alist;

    thee->n = 1;
    for (i=0; i<VAPBS_DIM; i++) {
        if (npts[i] < 3) {
            Vnm_print(2, 
                    "Vclist_ctor2:  n[%d] (%d) must be greater than 2!\n", 
                    i, npts[i]);
            return 0; 
        } 
        thee->npts[i] = npts[i];
        thee->n *= npts[i];
    }
    Vnm_print(0, "Vclist_ctor2:  Using %d x %d x %d hash table\n", 
            npts[0], npts[1], npts[2]);

    thee->mode = mode;
    switch (thee->mode) {
        case CLIST_AUTO_DOMAIN:
            Vnm_print(0, "Vclist_ctor2:  automatic domain setup.\n");
            break;
        case CLIST_MANUAL_DOMAIN:
            Vnm_print(0, "Vclist_ctor2:  manual domain setup.\n");
            Vnm_print(0, "Vclist_ctor2:  lower corner = [ \n");
            for (i=0; i<VAPBS_DIM; i++) {
                thee->lower_corner[i] = lower_corner[i];
                Vnm_print(0, "%g ", lower_corner[i]);
            }
            Vnm_print(0, "]\n");
            Vnm_print(0, "Vclist_ctor2:  upper corner = [ \n");
            for (i=0; i<VAPBS_DIM; i++) {
                thee->upper_corner[i] = upper_corner[i];
                Vnm_print(0, "%g ", upper_corner[i]);
            }
            Vnm_print(0, "]\n");
            break;
        default:
            Vnm_print(2, "Vclist_ctor2:  invalid setup mode (%d)!\n", mode);
            return 0;
    }

    thee->max_radius = max_radius;
    Vnm_print(0, "Vclist_ctor2:  Using %g max radius\n", max_radius);

    return 1;
}

/** Calculate the gridpoints an atom spans */
VPRIVATE void Vclist_gridSpan(Vclist *thee, 
        Vatom *atom, /** Atom */
        int imin[VAPBS_DIM], /** Set to min grid indices */
        int imax[VAPBS_DIM]  /** Set to max grid indices */
        ) {

    int i;
    double *coord, dc, idc, rtot;

    /* Get the position in the grid's frame of reference */
    coord = Vatom_getPosition(atom);

    /* Get the range the atom radius + probe radius spans */
    rtot = Vatom_getRadius(atom) + thee->max_radius;
    
    /* Calculate the range of grid points the inflated atom spans in the x 
     * direction. */
    for (i=0; i<VAPBS_DIM; i++) {
        dc = coord[i] - (thee->lower_corner)[i];
        idc = (dc + rtot)/(thee->spacs[i]);
        imax[i] = (int)(ceil(idc));
        imax[i] = VMIN2(imax[i], thee->npts[i]-1);
        idc = (dc - rtot)/(thee->spacs[i]);
        imin[i] = (int)(floor(idc));
        imin[i] = VMAX2(imin[i], 0);
    }

}

/** Get the array index for a particular cell based on its i,j,k 
 * coordinates */
VPRIVATE int Vclist_arrayIndex(Vclist *thee, int i, int j, int k) {

    return (thee->npts[2])*(thee->npts[1])*i + (thee->npts[2])*j + k;

}


/** Assign atoms to cells */
VPRIVATE int Vclist_assignAtoms(Vclist *thee) {

    int iatom, i, j, k, ui, inext;
    int imax[VAPBS_DIM], imin[VAPBS_DIM];
    int totatoms;
    Vatom *atom;
    VclistCell *cell;


    /* Find out how many atoms are associated with each grid point */
    totatoms = 0;
    for (iatom=0; iatom<Valist_getNumberAtoms(thee->alist); iatom++) { 

        /* Get grid span for atom */
        atom = Valist_getAtom(thee->alist, iatom);
        Vclist_gridSpan(thee, atom, imin, imax);

        /* Now find and assign the grid points */
        VASSERT(VAPBS_DIM == 3);
        for ( i = imin[0]; i <= imax[0]; i++) {
            for ( j = imin[1]; j <= imax[1]; j++) {
                for ( k = imin[2]; k <= imax[2]; k++) {
                    /* Get index to array */
                    ui = Vclist_arrayIndex(thee, i, j, k);
                    /* Increment number of atoms for this grid point */
                    cell = &(thee->cells[ui]);
                    (cell->natoms)++;
                    totatoms++;
                } 
            } 
        } 
    } 
    Vnm_print(0, "Vclist_assignAtoms:  Have %d atom entries\n", totatoms);

    /* Allocate the space to store the pointers to the atoms */
    for (ui=0; ui<thee->n; ui++) {
        cell = &(thee->cells[ui]);
        if ( !VclistCell_ctor2(cell, cell->natoms) ) {
            Vnm_print(2, "Vclist_assignAtoms:  cell error!\n");
            return 0;
        }
        /* Clear the counter for later use */
        cell->natoms = 0;
    }
 
    /* Assign the atoms to grid points */
    for (iatom=0; iatom<Valist_getNumberAtoms(thee->alist); iatom++) {

        /* Get grid span for atom */
        atom = Valist_getAtom(thee->alist, iatom);
        Vclist_gridSpan(thee, atom, imin, imax);

        /* Now find and assign the grid points */
        for (i = imin[0]; i <= imax[0]; i++) {
            for (j = imin[1]; j <= imax[1]; j++) {
                for (k = imin[2]; k <= imax[2]; k++) {
                    /* Get index to array */
                    ui = Vclist_arrayIndex(thee, i, j, k);
                    cell = &(thee->cells[ui]);
                    /* Index of next available array location */
                    inext = cell->natoms;
                    cell->atoms[inext] = atom;
                    /* Increment number of atoms */
                    (cell->natoms)++;
                }
            }
        }
    } 

    return 1;
}

/** Main (FORTRAN stub) constructor */
VPUBLIC int Vclist_ctor2(Vclist *thee, Valist *alist, double max_radius,
        int npts[VAPBS_DIM], Vclist_DomainMode mode, 
        double lower_corner[VAPBS_DIM], double upper_corner[VAPBS_DIM]) {

    int i;
    VclistCell *cell;

    /* Check and store parameters */
    if ( !Vclist_storeParms(thee, alist, max_radius, npts, mode, lower_corner,
                upper_corner) ) {
        Vnm_print(2, "Vclist_ctor2:  parameter check failed!\n");
        return 0;
    }

    /* Set up memory */
    thee->vmem = Vmem_ctor("APBS::VCLIST");
    if (thee->vmem == VNULL) {
        Vnm_print(2, "Vclist_ctor2:  memory object setup failed!\n");
        return 0;
    }

    /* Set up cells */
    thee->cells = Vmem_malloc( thee->vmem, thee->n, sizeof(VclistCell) );
    if (thee->cells == VNULL) {
        Vnm_print(2, 
                "Vclist_ctor2:  Failed allocating %d VclistCell objects!\n",
                thee->n);
        return 0;
    }
    for (i=0; i<thee->n; i++) {
        cell = &(thee->cells[i]);
        cell->natoms = 0;
    }

    /* Set up the grid */
    if (!Vclist_setupGrid(thee)) {
        Vnm_print(2, "Vclist_ctor2:  grid setup failed!\n");
        return 0;
    }

    /* Assign atoms to grid cells */
    if (!Vclist_assignAtoms(thee)) {
        Vnm_print(2, "Vclist_ctor2:  atom assignment failed!\n");
        return 0;
    }
	

	


    return 1;
}

/** Destructor */
VPUBLIC void Vclist_dtor(Vclist **thee) {
    
    if ((*thee) != VNULL) {
        Vclist_dtor2(*thee);
        Vmem_free(VNULL, 1, sizeof(Vclist), (void **)thee);
        (*thee) = VNULL;
    }

}

/** Main (stub) destructor */
VPUBLIC void Vclist_dtor2(Vclist *thee) {

    VclistCell *cell;
    int i;

    for (i=0; i<thee->n; i++) {
        cell = &(thee->cells[i]);
        VclistCell_dtor2(cell);
    }
    Vmem_free(thee->vmem, thee->n, sizeof(VclistCell), 
            (void **)&(thee->cells));
    Vmem_dtor(&(thee->vmem));

}

VPUBLIC VclistCell* Vclist_getCell(Vclist *thee, double pos[VAPBS_DIM]) {

    int i, ic[VAPBS_DIM], ui;
    double c[VAPBS_DIM];

    /* Convert to grid based coordinates */
    for (i=0; i<VAPBS_DIM; i++) {
        c[i] = pos[i] - (thee->lower_corner)[i];
        ic[i] = (int)(c[i]/thee->spacs[i]);
        if (ic[i] < 0) {
            /* printf("OFF LOWER CORNER!\n"); */
            return VNULL;
        } else if (ic[i] >= thee->npts[i]) {
            /* printf("OFF UPPER CORNER!\n"); */
            return VNULL;
        }
    }

    /* Get the array index */
    VASSERT(VAPBS_DIM == 3);
    ui = Vclist_arrayIndex(thee, ic[0], ic[1], ic[2]);

    return &(thee->cells[ui]);

}

VPUBLIC VclistCell* VclistCell_ctor(int natoms) {

    VclistCell *thee = VNULL;

    /* Set up the structure */
    thee = Vmem_malloc(VNULL, 1, sizeof(VclistCell));
    VASSERT( thee != VNULL);
    VASSERT( VclistCell_ctor2(thee, natoms) );

    return thee;
}

VPUBLIC int VclistCell_ctor2(VclistCell *thee, int natoms) {

    if (thee == VNULL) {
        Vnm_print(2, "VclistCell_ctor2:  NULL thee!\n");
        return 0;
    }

    thee->natoms = natoms;
    if (thee->natoms > 0) {
        thee->atoms = Vmem_malloc(VNULL, natoms, sizeof(Vatom *));
        if (thee->atoms == VNULL) {
            Vnm_print(2, 
          "VclistCell_ctor2:  unable to allocate space for %d atom pointers!\n",
              natoms);
            return 0;
        }
    }

    return 1;

}

VPUBLIC void VclistCell_dtor(VclistCell **thee) {
    
    if ((*thee) != VNULL) {
        VclistCell_dtor2(*thee);
        Vmem_free(VNULL, 1, sizeof(VclistCell), (void **)thee);
        (*thee) = VNULL;
    }

}

/** Main (stub) destructor */
VPUBLIC void VclistCell_dtor2(VclistCell *thee) {

    if (thee->natoms > 0) {
        Vmem_free(VNULL, thee->natoms, sizeof(Vatom *), 
                (void **)&(thee->atoms));
    }

}

