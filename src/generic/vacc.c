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

VPRIVATE int ivdwAccExclus(Vacc *thee, double center[3], double radius, int atomID);

VPUBLIC Vacc* Vacc_ctor(Valist *alist, double max_radius, int nx,
    int ny, int nz, int nsphere) {


    Vacc *thee = VNULL;

    /* Set up the structure */
    thee = Vmem_malloc(VNULL, 1, sizeof(Vacc) );
    VASSERT( thee != VNULL);
    VASSERT( Vacc_ctor2(thee, alist, max_radius, nx, ny, nz, nsphere));
    return thee;
}

VPUBLIC Vacc* Vacc_ctorFocus(Valist *alist, double max_radius, int nx,
    int ny, int nz, int nsphere,double x_min, double y_min, 
	double z_min, double x_max, double y_max, double z_max) {


    Vacc *thee = VNULL;

    /* Set up the structure */
    thee = Vmem_malloc(VNULL, 1, sizeof(Vacc) );
    VASSERT( thee != VNULL);
    VASSERT( Vacc_ctor2Focus(thee, alist, max_radius, nx, ny, nz, nsphere,
						x_min, y_min, z_min, x_max, y_max, z_max));

    return thee;
}

/** Get the dimensions of the molecule stored in thee->alist */
VPRIVATE void Vacc_getMolDims(
        Vacc *thee,
        double *x_min, /** Set to min molecule x */
        double *x_max, /** Set to max molecule x */
        double *y_min, /** Set to min molecule y */
        double *y_max, /** Set to max molecule y */
        double *z_min, /** Set to min molecule z */
        double *z_max, /** Set to max molecule z */
        double *r_max /** Set to max atom radius */
        ) {

    int i;
    double x, y, z;
    Valist *alist;
    Vatom *atom;

    alist = thee->alist;

    *x_max = -VLARGE;
    *y_max = -VLARGE;
    *z_max = -VLARGE;
    *x_min = -VLARGE;
    *y_min = -VLARGE;
    *z_min = VLARGE;
    *r_max = -1.0;
    for (i=0; i<Valist_getNumberAtoms(alist); i++) {
        atom = Valist_getAtom(alist, i);
        x = (Vatom_getPosition(atom))[0];
        y = (Vatom_getPosition(atom))[1];
        z = (Vatom_getPosition(atom))[2];
        if (x < *x_min) *x_min = x;
        if (y < *y_min) *y_min = y;
        if (z < *z_min) *z_min = z;
        if (x > *x_max) *x_max = x;
        if (y > *y_max) *y_max = y;
        if (z > *z_max) *z_max = z;
        if (Vatom_getRadius(atom) > *r_max) *r_max = Vatom_getRadius(atom);
    }
}

/** Setup Vacc grid */
VPRIVATE int Vacc_setupGrid(Vacc *thee) {
    /* Inflation factor ~ sqrt(2)*/
    #define VACC_INFLATE 1.42

    double x_max, x_min, y_max, y_min, z_max, z_min, r_max;
    double xlen, ylen, zlen;

    /* Get molecule dimensions */
    Vacc_getMolDims(thee, &x_max, &x_min, &y_max, &y_min, &z_max, 
            &z_min, &r_max);
 
    /* Set up grid spacings */
    xlen = x_max-x_min + 2.0*VACC_INFLATE*(r_max+thee->max_radius);
    ylen = y_max-y_min + 2.0*VACC_INFLATE*(r_max+thee->max_radius);
    zlen = z_max-z_min + 2.0*VACC_INFLATE*(r_max+thee->max_radius);
    thee->hx = xlen/((double)(thee->nx - 1));
    thee->hy = ylen/((double)(thee->ny - 1));
    thee->hzed = zlen/((double)(thee->nz - 1));
    Vnm_print(0, "Vacc_ctor2:  Grid lengths = (%g, %g, %g)\n",
            xlen, ylen, zlen);
 
    /* Setup lower corner */
    (thee->grid_lower_corner)[0] = 
        x_min - VACC_INFLATE*(r_max+thee->max_radius);
    (thee->grid_lower_corner)[1] = 
        y_min - VACC_INFLATE*(r_max+thee->max_radius);
    (thee->grid_lower_corner)[2] = 
        z_min - VACC_INFLATE*(r_max+thee->max_radius);
    Vnm_print(0, "Vacc_ctor2:  Grid lower corner = (%g, %g, %g)\n",
      (thee->grid_lower_corner)[0], (thee->grid_lower_corner)[1],
      (thee->grid_lower_corner)[2]);

    return 1;
}

/** Check and store parameters passed to constructor */
VPRIVATE int Vacc_storeParms(Vacc *thee,
        Valist *alist,
        double max_radius,
        int nx, int ny, int nz,
        int nsphere) {

    if (alist == VNULL) {
        Vnm_print(2, "Vacc_ctor2:  Got NULL Valist!\n");
        return 0;
    } else thee->alist = alist;
    if (nx < 3) {
        Vnm_print(2, "Vacc_ctor2:  nx (%d) must be greater than 2!\n", nx);
        return 0;
    } else thee->nx = nx;
    if (ny < 3) {
        Vnm_print(2, "Vacc_ctor2:  ny (%d) must be greater than 2!\n", ny);
        return 0;
    } else thee->ny = ny;
    if (nz < 3) {
        Vnm_print(2, "Vacc_ctor2:  nz (%d) must be greater than 2!\n", nz);
        return 0;
    } else thee->nz = nz;
    thee->n = nx*ny*nz;
    Vnm_print(0, "Vacc_ctor2:  Using %d x %d x %d hash table\n", nx, ny, nz);

    thee->max_radius = max_radius;
    Vnm_print(0, "Vacc_ctor2:  Using %g max radius\n", max_radius);

    thee->nsphere = nsphere;
    Vnm_print(0, "Vacc_ctor2:  Using %d-point probe sphere\n", nsphere);

    return 1;
}

/** Allocate (and clear) space for storage */
VPRIVATE int Vacc_allocate(Vacc *thee) {

    int i;

    thee->natoms = Vmem_malloc(thee->vmem, thee->n, sizeof(int));
    if (thee->natoms == VNULL) {
        Vnm_print(2, "Vacc_ctor2:  Failed to allocate %d (int)s for thee->natoms!\n",
                thee->n);
        return 0;
    }
    for (i=0; i<thee->n; i++) (thee->natoms)[i] = 0;

    thee->atomFlags = Vmem_malloc(thee->vmem, thee->n, sizeof(int));
    if (thee->atomFlags == VNULL) {
        Vnm_print(2, "Vacc_ctor2:  Failed to allocate %d (int)s for thee->atomFlags!\n", 
                thee->n);
        return 0;
    }
    for (i=0; i<thee->n; i++) (thee->atomFlags)[i] = 0;

    thee->atomIDs = Vmem_malloc(thee->vmem, thee->n, sizeof(int *));
    if (thee->atomIDs == VNULL) {
        Vnm_print(2, "Vacc_ctor2:  Failed to allocate %d (int *)s for thee->atomIDs!\n",
                thee->n);
        return 0;
    }
    for (i=0; i<thee->n; i++) (thee->atomIDs)[i] = VNULL;

    thee->area = Vmem_malloc(thee->vmem, Valist_getNumberAtoms(thee->alist), 
        sizeof(double));
    if (thee->area == VNULL) {
        Vnm_print(2, "Vacc_ctor2:  Failed to allocate %d (double)s for thee->area!\n", 
                Valist_getNumberAtoms(thee->alist));
        return 0;
    }
    for (i=0; i<Valist_getNumberAtoms(thee->alist); i++) thee->area[i] = 0;

    return 1;
}

/** Calculate the gridpoints an atom spans */
VPRIVATE void Vacc_gridSpan(Vacc *thee, 
        Vatom *atom, /** Atom */
        int *i_min, /** Set to min x-grid index */
        int *i_max, /** Set to max x-grid index */
        int *j_min, /** Set to min y-grid index */
        int *j_max, /** Set to max y-grid index */
        int *k_min, /** Set to min z-grid index */
        int *k_max  /** Set to max z-grid index */
        ) {

    double *coord;
    double x, y, z, rtot;

    /* Get the position in the grid's frame of reference */
    coord = Vatom_getPosition(atom);
    x = coord[0] - (thee->grid_lower_corner)[0];
    y = coord[1] - (thee->grid_lower_corner)[1];
    z = coord[2] - (thee->grid_lower_corner)[2];

    /* Get the range the atom radius + probe radius spans */
    rtot = Vatom_getRadius(atom) + thee->max_radius;
    
    /* Calculate the range of grid points the inflated atom spans in the x 
     * direction. */
    *i_max = (int)(ceil((x + rtot)/(thee->hx) ));
    *i_max = VMIN2(*i_max, thee->nx-1);
    *i_min = (int)(floor( (x - rtot)/(thee->hx) ));
    *i_min = VMAX2(*i_min, 0);
    *j_max = (int)( ceil( (y + rtot)/(thee->hy) ));
    *j_max = VMIN2(*j_max, thee->ny-1);
    *j_min = (int)(floor( (y - rtot)/(thee->hy) ));
    *j_min = VMAX2(*j_min, 0);
    *k_max = (int)( ceil( (z + rtot)/(thee->hzed) ));
    *k_max = VMIN2(*k_max, thee->nz-1);
    *k_min = (int)(floor( (z - rtot)/(thee->hzed) ));
    *k_min = VMAX2(*k_min, 0);

#if 0
    Vnm_print(0, "VACC DEBUG: %d <= i <= %d\n", i_min, i_max);
    Vnm_print(0, "VACC DEBUG: %d <= j <= %d\n", j_min, j_max);
    Vnm_print(0, "VACC DEBUG: %d <= k <= %d\n", k_min, k_max);
#endif

}

/** Get index for grid point in array */
VPRIVATE int Vacc_gridIndex(Vacc *thee,
        int i, /** Grid point x-index */
        int j, /** Grid point y-index */
        int k /** Grid point z-index */
        ) {

    return (thee->nz)*(thee->ny)*i + (thee->nz)*j + k;

}

/** Assign atoms to grid points */
VPRIVATE int Vacc_assignAtoms(Vacc *thee) {

    int iatom, i, j, k, ui, inext;
    int i_max, i_min, j_max, j_min, k_max, k_min;
    int totatoms;
    Vatom *atom;


    /* Find out how many atoms are associated with each grid point */
    totatoms = 0;
    for (iatom=0; iatom<Valist_getNumberAtoms(thee->alist); iatom++) { 

        /* Get grid span for atom */
        atom = Valist_getAtom(thee->alist, iatom);
        Vacc_gridSpan(thee, atom, 
                &i_min, &i_max, &j_min, &j_max, &k_min, &k_max); 

        /* Now find and assign the grid points */
        for ( i = i_min; i <= i_max; i++) {
            for ( j = j_min; j <= j_max; j++) {
                for ( k = k_min; k <= k_max; k++) {
                    /* Get index to array */
                    ui = Vacc_gridIndex(thee, i, j, k);
                    /* Increment number of atoms for this grid point */
                    (thee->natoms[ui])++;
                    totatoms += thee->natoms[ui];
                } 
            } 
        } 
    } 
    Vnm_print(0, "Vacc_ctor2:  Have %d atom entries\n", totatoms);

    /* Allocate the space to store the pointers to the atoms */
    for (ui=0; ui<thee->n; ui++) {
        if ((thee->natoms)[ui] > 0) {
            thee->atomIDs[ui] = Vmem_malloc(thee->vmem, thee->natoms[ui],
              sizeof(int));
            if (thee->atomIDs[ui] == VNULL) {
                Vnm_print(2, "Failed to allocate space for %d (int)s for thee->atomIDs[%d]!\n",
                        thee->natoms[ui], ui);
                return 0;
            }
        }
        /* Clear the counter for later use */
        thee->natoms[ui] = 0;
    }
 
    /* Assign the atoms to grid points */
    for (iatom=0; iatom<Valist_getNumberAtoms(thee->alist); iatom++) {

        /* Get grid span for atom */
        atom = Valist_getAtom(thee->alist, iatom);
        Vacc_gridSpan(thee, atom, 
                &i_min, &i_max, &j_min, &j_max, &k_min, &k_max); 

        /* Now find and assign the grid points */
        for (i = i_min; i <= i_max; i++) {
            for (j = j_min; j <= j_max; j++) {
                for (k = k_min; k <= k_max; k++) {
                    /* Get index to array */
                    ui = Vacc_gridIndex(thee, i, j, k);
                    /* Index of next available array location */
                    inext = thee->natoms[ui];
                    thee->atomIDs[ui][inext] = iatom;
                    /* Increment number of atoms */
                    (thee->natoms[ui])++;
                }
            }
        }
    } 

    return 1;
}

VPUBLIC int Vacc_ctor2(Vacc *thee, Valist *alist, double max_radius,
    int nx, int ny, int nz, int nsphere) {

    /* Check and store parameters */
    if (!Vacc_storeParms(thee, alist, max_radius, nx, ny, nz, nsphere)) {
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
    thee->sphere = Vacc_sphere(thee, &(thee->nsphere));
    if (thee->sphere == VNULL) {
        Vnm_print(2, "Vacc_ctor2:  probe sphere setup failed!\n");
        return 0;
    }
 
    /* Allocate space */
    if (!Vacc_allocate(thee)) {
        Vnm_print(2, "Vacc_ctor2:  memory allocation failed!\n");
        return 0;
    }

    /* Setup the grid */
    if (!Vacc_setupGrid(thee)) {
        Vnm_print(2, "Vacc_ctor2:  grid setup failed!\n");
        return 0;
    }

    /* Assign atoms to grid points */
    if (!Vacc_assignAtoms(thee)) {
        Vnm_print(2, "Vacc_ctor2:  atom assignment failed!\n");
        return 0;
    }



    return 1;
}

/* ///////////////////////////////////////////////////////////////////////////
// Routine:  Vacc_ctor2Focus
//
// Author:   Nathan Baker, Todd Dolinsky
/////////////////////////////////////////////////////////////////////////// */
VPUBLIC int Vacc_ctor2Focus(Vacc *thee, Valist *alist, double max_radius,
    int nx, int ny, int nz, int nsphere, double x_min, double y_min, 
	double z_min, double x_max, double y_max, double z_max) {

    /* Grid variables */
    int i;
    double x, y, z, *coord;
    int ii, jj, kk, totatoms;
    int i_min, j_min, k_min;
    int i_max, j_max, k_max;
    /* Natural grid coordinate (array position) */
    int ui;
    /* Atom radius */
    double rmax;
    double rtot, rtot2;
    Vatom *atom;

    VASSERT(alist != VNULL);
    thee->alist = alist;

    /* Set up memory management object */
    thee->vmem = Vmem_ctor("APBS::VACC");

    /* Set up grid dimensions */
    thee->nx = nx;
    thee->ny = ny;
    thee->nz = nz;
    thee->n = nx*ny*nz;
    if ((nx < 3) || (ny < 3) || (nz < 3)) {
        Vnm_print(2, "Vacc_ctorFocus:  nx, ny, nz must be greater than 2!\n");
        return 0;
    }
    Vnm_print(0, "Vacc_ctorFocus:  Using %d x %d x %d hash table\n", nx, ny, nz);
    Vnm_print(0, "Vacc_ctorFocus:  Using %g max radius\n", max_radius);
 
    /* Set up probe information */
    thee->nsphere = nsphere;
    thee->max_radius = max_radius;
    Vnm_print(0, "Vacc_ctorFocus:  Constructing sphere...\n");
    thee->sphere = Vacc_sphere(thee, &(thee->nsphere));
    VASSERT(thee->sphere != VNULL);
 
    /* Allocate space */
    thee->natoms = Vmem_malloc(thee->vmem, thee->n, sizeof(int));
    VASSERT(thee->natoms != VNULL);
    for (i=0; i<thee->n; i++) (thee->natoms)[i] = 0;
    thee->atomFlags = Vmem_malloc(thee->vmem, thee->n, sizeof(int));
    VASSERT(thee->atomFlags != VNULL);
    for (i=0; i<thee->n; i++) (thee->atomFlags)[i] = 0;
    thee->atomIDs = Vmem_malloc(thee->vmem, thee->n, sizeof(int *));
    VASSERT(thee->atomIDs != VNULL);
    for (i=0; i<thee->n; i++) (thee->atomIDs)[i] = VNULL;
    thee->area = Vmem_malloc(thee->vmem, Valist_getNumberAtoms(alist), 
        sizeof(double));
    VASSERT(thee->area != VNULL);
    for (i=0; i<Valist_getNumberAtoms(alist); i++) thee->area[i] = 0;

    /* Find dimensions of protein and atoms */
    rmax = -1.0;
    for (i=0; i<Valist_getNumberAtoms(alist); i++) {
        atom = Valist_getAtom(alist, i);
        if (Vatom_getRadius(atom) > rmax) rmax = Vatom_getRadius(atom);
    }

    /* Set up grid spacings, 2.84 > 2*sqrt(2) */
    thee->hx = (x_max - x_min + 2.84*(rmax+thee->max_radius))/(thee->nx - 1);
    thee->hy = (y_max - y_min + 2.84*(rmax+thee->max_radius))/(thee->ny - 1);
    thee->hzed = (z_max - z_min + 2.84*(rmax+thee->max_radius))/(thee->nz - 1);
 
    /* Inflate the grid a bit 1.42 > sqrt(2) */
    (thee->grid_lower_corner)[0] = x_min-1.42*(rmax+thee->max_radius);
    (thee->grid_lower_corner)[1] = y_min-1.42*(rmax+thee->max_radius);
    (thee->grid_lower_corner)[2] = z_min-1.42*(rmax+thee->max_radius);

    Vnm_print(0, "Vacc_ctorFocus:  Grid lower corner = (%g, %g, %g)\n",
      (thee->grid_lower_corner)[0], (thee->grid_lower_corner)[1],
      (thee->grid_lower_corner)[2]);
    Vnm_print(0, "Vacc_ctorFocus:  Grid lengths = (%g, %g, %g)\n",
      thee->hx*(thee->nx - 1), thee->hy*(thee->ny - 1), 
      thee->hzed*(thee->nz - 1));

    /* Find out how many atoms are associated with each grid point */
    totatoms = 0;
    for (i=0;i<Valist_getNumberAtoms(alist);i++) { 
        atom = Valist_getAtom(alist, i);
        /* Get the position in the grid's frame of reference */
        coord = Vatom_getPosition(atom);
        x = coord[0] - (thee->grid_lower_corner)[0];
        y = coord[1] - (thee->grid_lower_corner)[1];
        z = coord[2] - (thee->grid_lower_corner)[2];

        /* Get the range the atom radius + probe radius spans */
        rtot = Vatom_getRadius(atom) + thee->max_radius;
    
        /* Calculate the range of grid points the inflated atom spans in the x 
         * direction. */
        i_max = (int)( ceil( (x + rtot)/(thee->hx) ));
        i_max = VMIN2(i_max, nx-1);
        i_min = (int)(floor( (x - rtot)/(thee->hx) ));
        i_min = VMAX2(i_min, 0);
        j_max = (int)( ceil( (y + rtot)/(thee->hy) ));
        j_max = VMIN2(j_max, ny-1);
        j_min = (int)(floor( (y - rtot)/(thee->hy) ));
        j_min = VMAX2(j_min, 0);
        k_max = (int)( ceil( (z + rtot)/(thee->hzed) ));
        k_max = VMIN2(k_max, nz-1);
        k_min = (int)(floor( (z - rtot)/(thee->hzed) ));
        k_min = VMAX2(k_min, 0);

#if 0
        Vnm_print(0, "VACC DEBUG: %d <= i <= %d\n", i_min, i_max);
        Vnm_print(0, "VACC DEBUG: %d <= j <= %d\n", j_min, j_max);
        Vnm_print(0, "VACC DEBUG: %d <= k <= %d\n", k_min, k_max);
#endif

        /* Now find and assign the grid points */
        for ( ii = i_min; ii <= i_max; ii++) {
            for ( jj = j_min; jj <= j_max; jj++) {
                for ( kk = k_min; kk <= k_max; kk++) {
                    ui = (thee->nz)*(thee->ny)*ii + (thee->nz)*jj + kk;
                    (thee->natoms[ui])++;
                    totatoms += thee->natoms[ui];
                } 
            } 
        } 
    } /* for i =0:natoms */
    Vnm_print(0, "Vacc_ctorFocus:  Have %d atom entries\n", totatoms);

    /* Allocate the space to store the pointers to the atoms */
    for (i=0; i<thee->n; i++) {
        if ((thee->natoms)[i] > 0) {
            thee->atomIDs[i] = Vmem_malloc(thee->vmem, thee->natoms[i],
              sizeof(int));
            VASSERT(thee->atomIDs[i] != VNULL);
        }
        /* Clear the counter for later use */
        thee->natoms[i] = 0;
    }
 
    /* Assign the atoms to grid points */
    for (i=0; i<Valist_getNumberAtoms(alist); i++) {
        atom = Valist_getAtom(alist, i);
        /* Get the position in the grid's frame of reference */
        x = (Vatom_getPosition(atom))[0] - (thee->grid_lower_corner)[0];
        y = (Vatom_getPosition(atom))[1] - (thee->grid_lower_corner)[1];
        z = (Vatom_getPosition(atom))[2] - (thee->grid_lower_corner)[2];

        /* Get the range the atom radius + probe radius spans */
        rtot = Vatom_getRadius(atom) + thee->max_radius;
        rtot2 = VSQR(rtot);

        /* Now find and assign the grid points */
        i_max = (int)( ceil( (x + rtot)/(thee->hx) ));
        i_max = VMIN2(i_max, nx-1);
        i_min = (int)(floor( (x - rtot)/(thee->hx) ));
        i_min = VMAX2(i_min, 0);
        j_max = (int)( ceil( (y + rtot)/(thee->hy) ));
        j_max = VMIN2(j_max, ny-1);
        j_min = (int)(floor( (y - rtot)/(thee->hy) ));
        j_min = VMAX2(j_min, 0);
        k_max = (int)( ceil( (z + rtot)/(thee->hzed) ));
        k_max = VMIN2(k_max, nz-1);
        k_min = (int)(floor( (z - rtot)/(thee->hzed) ));
        k_min = VMAX2(k_min, 0);
        /* Now find and assign the grid points */
        for ( ii = i_min; ii <= i_max; ii++) {
            for ( jj = j_min; jj <= j_max; jj++) {
                for ( kk = k_min; kk <= k_max; kk++) {

                    ui = (thee->nz)*(thee->ny)*ii + (thee->nz)*jj + kk;
                    thee->atomIDs[ui][thee->natoms[ui]] = i;
                    (thee->natoms[ui])++;
                }
            }
        }
    } /* for i =0:natoms */

    return 1;
}


/* ///////////////////////////////////////////////////////////////////////////
// Routine:  Vacc_dtor
//
// Author:   Nathan Baker
/////////////////////////////////////////////////////////////////////////// */
VPUBLIC void Vacc_dtor(Vacc **thee) {
    
    if ((*thee) != VNULL) {
        Vacc_dtor2(*thee);
        Vmem_free(VNULL, 1, sizeof(Vacc), (void **)thee);
        (*thee) = VNULL;
    }

}

/* ///////////////////////////////////////////////////////////////////////////
// Routine:  Vacc_dtor2
//
// Author:   Nathan Baker
/////////////////////////////////////////////////////////////////////////// */
VPUBLIC void Vacc_dtor2(Vacc *thee) {

    int i;
    for (i=0; i<thee->n; i++) {
	if (thee->natoms[i] > 0)  Vmem_free(thee->vmem, (thee->natoms)[i],
          sizeof(int), (void **)&(thee->atomIDs[i]));
    }
    Vmem_free(thee->vmem, thee->n, sizeof(int *), (void **)&(thee->atomIDs));
    Vmem_free(thee->vmem, thee->n, sizeof(int), (void **)&(thee->natoms));
    Vmem_free(thee->vmem, thee->n, sizeof(int), (void **)&(thee->atomFlags));
    for (i=0; i<thee->nsphere; i++) 
      Vmem_free(thee->vmem, 3, sizeof(double), (void **)&(thee->sphere[i]));
    Vmem_free(thee->vmem, thee->nsphere, sizeof(double *), 
      (void **)&(thee->sphere));
    Vmem_free(thee->vmem, Valist_getNumberAtoms(thee->alist),
      sizeof(double), (void **)&(thee->area));

    Vmem_dtor(&(thee->vmem));
}

/* ///////////////////////////////////////////////////////////////////////////
// Routine:  Vacc_vdwAcc
//
// Author:   Nathan Baker
/////////////////////////////////////////////////////////////////////////// */
VPUBLIC double Vacc_vdwAcc(Vacc *thee, double center[3]) {

    int centeri, centerj, centerk;  /* Grid-based coordinates */
    int ui;                         /* Natural array coordinates */
    int iatom;                      /* Counters */
    double dist;
    Vatom *atom;
    double *apos;

    /* Convert to grid based coordinates */
    centeri = (int)( (center[0] - (thee->grid_lower_corner)[0])/thee->hx);
    centerj = (int)( (center[1] - (thee->grid_lower_corner)[1])/thee->hy);
    centerk = (int)( (center[2] - (thee->grid_lower_corner)[2])/thee->hzed);
   
    /* Check to make sure we're on the grid; if not, we're definitely 
     * accessible */ 
    if ((centeri < 0) || (centeri >= thee->nx) || 
        (centerj < 0) || (centerj >= thee->ny) || 
        (centerk < 0) || (centerk >= thee->nz)) return 1.0;
   
    /* If we're still here, then we need to check each atom until we find an
     * overlap at which point we can determine that the point is not 
     * accessible */
    ui = (thee->nz)*(thee->ny)*centeri + (thee->nz)*centerj + centerk;
    for (iatom=0;iatom<(thee->natoms)[ui];iatom++) {
        atom = Valist_getAtom(thee->alist, thee->atomIDs[ui][iatom]);
        apos = Vatom_getPosition(atom);
        dist = VSQR(center[0]-apos[0]) + VSQR(center[1]-apos[1])
               + VSQR(center[2]-apos[2]);
        if (dist < VSQR(Vatom_getRadius(atom))) return 0.0;
    }

    /* If we're still here, then the point is accessible */
    return 1.0;
}

/* ///////////////////////////////////////////////////////////////////////////
// Routine:  Vacc_ivdwAcc
//
// Author:   Nathan Baker
/////////////////////////////////////////////////////////////////////////// */
VPUBLIC double Vacc_ivdwAcc(Vacc *thee, double center[3], double radius) {

    return (double)ivdwAccExclus(thee, center, radius, -1);

}

/* ///////////////////////////////////////////////////////////////////////////
// Routine:  ivdwAccExclus
//
// Purpose:  Determines if a point is within the union of the spheres centered
//           at the atomic centers with radii equal to the sum of their van
//           der Waals radii and the probe radius.  Does not include
//           contributions from the specified atom.
//
// Args:     center => point to be tested
//           radius => radius to inflate by
//           atomID  => atom to ignore (-1 if none)
//           
//           Returns 1 if accessible (outside the molecular volume).
//
// Author:   Nathan Baker
/////////////////////////////////////////////////////////////////////////// */
VPRIVATE int ivdwAccExclus(Vacc *thee, double center[3], double radius, 
  int atomID) {

    int centeri, centerj, centerk;  /* Grid-based coordinates */
    int ui;                         /* Natural array coordinates */
    int iatom;                      /* Counters */
    double dist, *apos;
    Vatom *atom;

    /* We can only test probes with radii less than the max specified */
    VASSERT(thee != VNULL);
    if (radius > thee->max_radius) {
        Vnm_print(2, "Vacc_ivdwAcc: got radius (%g) bigger than max radius (%g)\n", 
          radius, thee->max_radius);
         VASSERT(0);
    }

    /* Convert to grid based coordinates */
    centeri = (int)( (center[0] - (thee->grid_lower_corner)[0])/thee->hx);
    centerj = (int)( (center[1] - (thee->grid_lower_corner)[1])/thee->hy);
    centerk = (int)( (center[2] - (thee->grid_lower_corner)[2])/thee->hzed);

    /* Check to make sure we're on the grid; if not, we're definitely
     * accessible */
    if ((centeri < 0) || (centeri >= thee->nx) || \
        (centerj < 0) || (centerj >= thee->ny) || \
        (centerk < 0) || (centerk >= thee->nz)) return 1;

    /* If we're still here, then we need to check each atom until we find an
     * overlap at which point we can determine that the point is not
     * accessible */
    ui = (thee->nz)*(thee->ny)*centeri + (thee->nz)*centerj + centerk;
    for (iatom=0;iatom<(thee->natoms)[ui];iatom++) {
        if (thee->atomIDs[ui][iatom] != atomID) {
            atom = Valist_getAtom(thee->alist, thee->atomIDs[ui][iatom]);
            apos = Vatom_getPosition(atom);
            dist = VSQR(apos[0]-center[0]) + VSQR(apos[1]-center[1])
                   + VSQR(apos[2]-center[2]);
            /* Only atoms with non-zero radii can contribute to solvent 
             * inaccessibility */
            if (Vatom_getRadius(atom) > 0.0) {
                if (dist < VSQR(Vatom_getRadius(atom)+radius)) return 0;
            }
        }
    }

    /* If we're still here, then the point is accessible */
    return 1;
}

/* ///////////////////////////////////////////////////////////////////////////
// Routine:  Vacc_splineAccGradAtom
//
// Author:   Nathan Baker
/////////////////////////////////////////////////////////////////////////// */
VPUBLIC void Vacc_splineAccGradAtom(Vacc *thee, double center[3], double win,
  double infrad, int atomID, double *grad) {

    int centeri, centerj, centerk;  /* Grid-based coordinates */
    double dist, *apos, arad, sm, sm2, w2i, w3i, mygrad;
    Vatom *atom;
    double mychi = 1.0;           /* Char. func. value for given atom */


    VASSERT(thee != NULL);

    /* Inverse squared window parameter */
    w2i = 1.0/(win*win);
    w3i = 1.0/(win*win*win);

    /* The grad is zero by default */
    grad[0] = 0.0;
    grad[1] = 0.0;
    grad[2] = 0.0;

    /* *** CALCULATE THE CHARACTERISTIC FUNCTION VALUE FOR THIS ATOM AND THE
     * *** MAGNITUDE OF THE FORCE *** */
    atom = Valist_getAtom(thee->alist, atomID);
    apos = Vatom_getPosition(atom);
    /* Zero-radius atoms don't contribute */
    if (Vatom_getRadius(atom) > 0.0) {
        arad = Vatom_getRadius(atom) + infrad;
        dist = VSQRT(VSQR(apos[0]-center[0]) + VSQR(apos[1]-center[1])
          + VSQR(apos[2]-center[2]));
        /* If we're inside an atom, the entire characteristic function
         * will be zero and the grad will be zero, so we can stop */
        if (dist <= (arad - win)) return;
        /* Likewise, if we're outside the smoothing window, the characteristic
         * function is unity and the grad will be zero, so we can stop */
        else if (dist >= (arad + win)) return;
        /* If we're inside the smoothing window */
        else {
            sm = dist - arad + win;
            sm2 = VSQR(sm);
            mychi = 0.75*sm2*w2i -0.25*sm*sm2*w3i;
            mygrad = 1.5*sm*w2i - 0.75*sm2*w3i;
        }
        /* Now assemble the grad vector */
        VASSERT(mychi > 0.0);
        grad[0] = -(mygrad/mychi)*((center[0] - apos[0])/dist);
        grad[1] = -(mygrad/mychi)*((center[1] - apos[1])/dist);
        grad[2] = -(mygrad/mychi)*((center[2] - apos[2])/dist);
    } 


    
}


/* ///////////////////////////////////////////////////////////////////////////
// Routine:  Vacc_splineAccAtom
//
// Author:   Nathan Baker
/////////////////////////////////////////////////////////////////////////// */
VPUBLIC double Vacc_splineAccAtom(Vacc *thee, double center[3], double win, 
  double infrad, int atomID) {

    int centeri, centerj, centerk;  
    double dist, *apos, arad, sm, sm2, w2i, w3i, value, stot, sctot;
    Vatom *atom;


    VASSERT(thee != NULL);

    /* Inverse squared window parameter */
    w2i = 1.0/(win*win);
    w3i = 1.0/(win*win*win);

    atom = Valist_getAtom(thee->alist, atomID);
    apos = Vatom_getPosition(atom);
    /* Zero-radius atoms don't contribute */
    if (Vatom_getRadius(atom) > 1.0) {
        arad = Vatom_getRadius(atom) + infrad;
        stot = arad + win;
        sctot = VMAX2(0, (arad - win));
        dist = VSQRT(VSQR(apos[0]-center[0]) + VSQR(apos[1]-center[1])
          + VSQR(apos[2]-center[2]));
        /* If we're inside an atom, the entire characteristic function
         * will be zero */
        if (dist <= sctot) {
            value = 0.0;
            return value;
        /* We're outside the smoothing window */
        } else if (dist >= stot) {
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

/* ///////////////////////////////////////////////////////////////////////////
// Routine:  Vacc_splineAcc
//
// Author:   Nathan Baker
/////////////////////////////////////////////////////////////////////////// */
VPUBLIC double Vacc_splineAcc(Vacc *thee, double center[3], double win, 
  double infrad) {

    int centeri, centerj, centerk, ui, iatom, atomID;      
    double value = 1.0;            


    VASSERT(thee != NULL);

    if (thee->max_radius < (win + infrad)) {
        Vnm_print(2, "Vacc_splineAcc:  ERROR -- Vacc constructed with max_radius=%g;\n",
          thee->max_radius);
        Vnm_print(2, "Vacc_splineAcc:  ERROR -- Insufficient for window=%g, inflation radius=%g\n", 
          win, infrad);
        VASSERT(0);
    }

    /* Convert to grid based coordinates */
    centeri = (int)( (center[0] - (thee->grid_lower_corner)[0])/thee->hx);
    centerj = (int)( (center[1] - (thee->grid_lower_corner)[1])/thee->hy);
    centerk = (int)( (center[2] - (thee->grid_lower_corner)[2])/thee->hzed);

    /* Check to make sure we're on the grid; if not, then our characteristic
     * function is definitely unity */
    if ((centeri < 0) || (centeri >= thee->nx) || \
        (centerj < 0) || (centerj >= thee->ny) || \
        (centerk < 0) || (centerk >= thee->nz)) {
        return 1;
    }

    /* If we're still here, then we need to check each atom until we find an
     * overlap at which point we can determine that the point is not
     * accessible */
    ui = (thee->nz)*(thee->ny)*centeri + (thee->nz)*centerj + centerk;
    /* First, reset the list of atom flags */
    for (iatom=0;iatom<(thee->natoms)[ui];iatom++) 
      thee->atomFlags[thee->atomIDs[ui][iatom]] = 0;
    /* Now loop through the atoms assembling the characteristic function */
    for (iatom=0;iatom<(thee->natoms)[ui];iatom++) {
        /* Check to see if we've counted this atom already */
        if (!(thee->atomFlags[thee->atomIDs[ui][iatom]])) {
            thee->atomFlags[thee->atomIDs[ui][iatom]] = 1;

            atomID = thee->atomIDs[ui][iatom];
            value *= Vacc_splineAccAtom(thee, center, win, infrad, atomID);
            
            if (value < VSMALL) return value;
        } 
    }
 
    return value;
}

/* ///////////////////////////////////////////////////////////////////////////
// Routine:  Vacc_molAcc
//
// Author:   Nathan Baker
/////////////////////////////////////////////////////////////////////////// */
VPUBLIC double Vacc_molAcc(Vacc *thee, double center[3], double radius) {

    int ipt;
    double vec[3];

    /* ******* CHECK IF OUTSIDE ATOM+PROBE RADIUS SURFACE ***** */
    if (Vacc_ivdwAcc(thee, center, radius) == 1.0) return 1;

    /* ******* CHECK IF INSIDE ATOM RADIUS SURFACE ***** */
    if (Vacc_vdwAcc(thee, center) == 0.0) return 0;

    /* ******* CHECK IF OUTSIDE MOLECULAR SURFACE ***** */
    /* Let S be the sphere of radius radius centered at the point we are
     * testing.  We are outside the molecular surface if there is a point on
     * the surface of S that is outside the atom+probe radius surface */
    VASSERT(thee->sphere != VNULL);
    for (ipt=0; ipt<thee->nsphere; ipt++) {
        vec[0] = radius*thee->sphere[ipt][0] + center[0];
        vec[1] = radius*thee->sphere[ipt][1] + center[1];
        vec[2] = radius*thee->sphere[ipt][2] + center[2];
        if (Vacc_ivdwAcc(thee, vec, radius)) return 1.0;
    }

    /* If all else failed, we are not inside the molecular surface */
    return 0.0;
}

/* ///////////////////////////////////////////////////////////////////////////
// Routine:  Vacc_fastMolAcc
//
// Author:   Nathan Baker
/////////////////////////////////////////////////////////////////////////// */
VPUBLIC double Vacc_fastMolAcc(Vacc *thee, double center[3], double radius) {

    int ipt;
    double vec[3];

    /* ******* CHECK IF OUTSIDE MOLECULAR SURFACE ***** */
    /* Let S be the sphere of radius radius centered at the point we are
     * testing.  We are outside the molecular surface if there is a point on
     * the surface of S that is outside the atom+probe radius surface */
    VASSERT(thee->sphere != VNULL);
    for (ipt=0; ipt<thee->nsphere; ipt++) {
        vec[0] = radius*thee->sphere[ipt][0] + center[0];
        vec[1] = radius*thee->sphere[ipt][1] + center[1];
        vec[2] = radius*thee->sphere[ipt][2] + center[2];
        if (Vacc_ivdwAcc(thee, vec, radius)) return 1.0;
    }

    /* If all else failed, we are not inside the molecular surface */
    return 0.0;
}


#if defined(HAVE_MC_H)
/* ///////////////////////////////////////////////////////////////////////////
// Routine:  Vacc_writeGMV
//
// Purpose:  Write out the chosen accessibility data at each vertex.  The the
//           appropriate isosurface routine, this would generate a
//           representation of the molecular surface as ``seen" by the PBE
//           solver.
//
// Arguments: radius   Radius of sphere to test
//            meth     Plot accessibility for molecular surface (meth=0),
//                     inflated van der Waals (meth=1), or van der Waals
//                     (meth=2)
//            gm       Gem object with mesh data
//            iodev    Device (usually "FILE")
//            iofmt    Format (usually "ASC")
//            iohost   Host   (usually "localhost")
//            iofile   Filename
//
// Author:   Nathan Baker
/////////////////////////////////////////////////////////////////////////// */
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

/* ///////////////////////////////////////////////////////////////////////////
// Routine:  Vacc_sphere
//
// Author:   Nathan Baker (original FORTRAN routine from UHBD by Michael
//           Gilson)
/////////////////////////////////////////////////////////////////////////// */
VPUBLIC double** Vacc_sphere(Vacc *thee, int *npts) {

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

/* ///////////////////////////////////////////////////////////////////////////
// Routine:  Vacc_totalSASA
//
// Author:   Nathan Baker (original FORTRAN routine from UHBD by Brock Luty)
/////////////////////////////////////////////////////////////////////////// */
VPUBLIC double Vacc_totalSASA(Vacc *thee, double radius) { 

    int i;
    double area = 0.0;

    for (i=0; i<Valist_getNumberAtoms(thee->alist); i++) {
        thee->area[i] = Vacc_atomSASA(thee, radius, i);
        area += (thee->area[i]);
    }

    return area;

}

/* ///////////////////////////////////////////////////////////////////////////
// Routine:  Vacc_atomSASA
//
// Purpose:  Calculates the contribution to the PROBE-CENTERED
//           solvent-accessible area from this atom
//
// Args:     radius  The radius of the solvent probe in Angstroms
//           iatom   Index of the atom in thee->alist
//
// Author:   Nathan Baker (original FORTRAN routine from UHBD by Brock Luty)
/////////////////////////////////////////////////////////////////////////// */
VPUBLIC double Vacc_atomSASA(Vacc *thee, double srad, int iatom) { 

    int ipt;
    double area = 0.0;
    double *tPos, tRad, vec[3];
    Vatom *thisAtom;

    /* Get the atom information */
    thisAtom = Valist_getAtom(thee->alist, iatom);
    tPos = Vatom_getPosition(thisAtom);
    tRad = Vatom_getRadius(thisAtom);

    for (ipt=0; ipt<thee->nsphere; ipt++) {
        vec[0] = (tRad+srad)*thee->sphere[ipt][0] + tPos[0];
        vec[1] = (tRad+srad)*thee->sphere[ipt][1] + tPos[1];
        vec[2] = (tRad+srad)*thee->sphere[ipt][2] + tPos[2];
        if (ivdwAccExclus(thee, vec, srad, iatom)) area += 1.0;
    }

    /* We will return UHBD's asas2: probe-centered solvent-accessible surface
     * area */
    area = area/((double)(thee->nsphere))*4.0*VPI*(tRad+srad)*(tRad+srad);

    return area;

}
