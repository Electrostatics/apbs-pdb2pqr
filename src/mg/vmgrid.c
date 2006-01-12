/**
 *  @file    vmgrid.c
 *  @author  Nathan Baker
 *  @brief   Class Vmgrid methods
 *  @ingroup Vmgrid
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
#include "apbs/vmgrid.h"

VEMBED(rcsid="$Id$")

/* ///////////////////////////////////////////////////////////////////////////
// Routine:  Vmgrid_ctor
// Author:   Nathan Baker
/////////////////////////////////////////////////////////////////////////// */
VPUBLIC Vmgrid* Vmgrid_ctor() {

    Vmgrid *thee = VNULL;

    thee = Vmem_malloc(VNULL, 1, sizeof(Vmgrid));
    VASSERT(thee != VNULL);
    VASSERT(Vmgrid_ctor2(thee));

    return thee;
}

/* ///////////////////////////////////////////////////////////////////////////
// Routine:  Vmgrid_ctor2
// Author:   Nathan Baker
/////////////////////////////////////////////////////////////////////////// */
VPUBLIC int Vmgrid_ctor2(Vmgrid *thee) {

    int i;

    if (thee == VNULL) return 0;

    thee->ngrids = 0;
    for (i=0; i<VMGRIDMAX; i++) thee->grids[i] = VNULL;
    
    return 1;
}

/* ///////////////////////////////////////////////////////////////////////////
// Routine:  Vmgrid_dtor
// Author:   Nathan Baker
/////////////////////////////////////////////////////////////////////////// */
VPUBLIC void Vmgrid_dtor(Vmgrid **thee) {

    if ((*thee) != VNULL) {
        Vmgrid_dtor2(*thee);
        Vmem_free(VNULL, 1, sizeof(Vmgrid), (void **)thee);
        (*thee) = VNULL;
    }
}

/* ///////////////////////////////////////////////////////////////////////////
// Routine:  Vmgrid_dtor2
// Author:   Nathan Baker
/////////////////////////////////////////////////////////////////////////// */
VPUBLIC void Vmgrid_dtor2(Vmgrid *thee) { ; }

/* ///////////////////////////////////////////////////////////////////////////
// Routine:  Vmgrid_value
// Author:   Nathan Baker
/////////////////////////////////////////////////////////////////////////// */
VPUBLIC int Vmgrid_value(Vmgrid *thee, double pt[3], double *value) {

    int i, rc;
    double tvalue;
  
    VASSERT(thee != VNULL);

    for (i=0; i<thee->ngrids; i++) {
        rc = Vgrid_value(thee->grids[i], pt, &tvalue);
        if (rc) {
            *value = tvalue;
            return 1;
        }
    }

    Vnm_print(2, "Vmgrid_value:  Point (%g, %g, %g) not found in \
hiearchy!\n", pt[0], pt[1], pt[2]);

    return 0;
}

/* ///////////////////////////////////////////////////////////////////////////
// Routine:  Vmgrid_curvature
//
//   Notes:  cflag=0 ==> Reduced Maximal Curvature
//           cflag=1 ==> Mean Curvature (Laplace)
//           cflag=2 ==> Gauss Curvature
//           cflag=3 ==> True Maximal Curvature
//
// Authors:  Nathan Baker
/////////////////////////////////////////////////////////////////////////// */
VPUBLIC int Vmgrid_curvature(Vmgrid *thee, double pt[3], int cflag, 
  double *value) {

    int i, rc;
    double tvalue;

    VASSERT(thee != VNULL);

    for (i=0; i<thee->ngrids; i++) {
        rc = Vgrid_curvature(thee->grids[i], pt, cflag, &tvalue);
        if (rc) {
            *value = tvalue;
            return 1;
        }
    }

    Vnm_print(2, "Vmgrid_curvature:  Point (%g, %g, %g) not found in \
hiearchy!\n", pt[0], pt[1], pt[2]);

    return 0;


}

/* ///////////////////////////////////////////////////////////////////////////
// Routine:  Vmgrid_gradient
//
// Authors:  Nathan Baker
/////////////////////////////////////////////////////////////////////////// */
VPUBLIC int Vmgrid_gradient(Vmgrid *thee, double pt[3], double grad[3]) {

    int i, j, rc;
    double tgrad[3];

    VASSERT(thee != VNULL);

    for (i=0; i<thee->ngrids; i++) {
        rc = Vgrid_gradient(thee->grids[i], pt, tgrad);
        if (rc) {
            for (j=0; j<3; j++) grad[j] = tgrad[j];
            return 1;
        }
    }

    Vnm_print(2, "Vmgrid_gradient:  Point (%g, %g, %g) not found in \
hiearchy!\n", pt[0], pt[1], pt[2]);

    return 0;


}

/* ///////////////////////////////////////////////////////////////////////////
// Routine:  Vmgrid_addGrid
//
// Authors:  Nathan Baker
/////////////////////////////////////////////////////////////////////////// */
VPUBLIC int Vmgrid_addGrid(Vmgrid *thee, Vgrid *grid) {

    int i, j, rc;
    double tgrad[3];

    VASSERT(thee != VNULL);

    if (grid == VNULL) {
        Vnm_print(2, "Vmgrid_addGrid:  Not adding VNULL grid!\n");
        return 0;
    }

    if (thee->ngrids >= VMGRIDMAX) {
        Vnm_print(2, "Vmgrid_addGrid:  Too many grids in hierarchy (max = \
%d)!\n", VMGRIDMAX);
        Vnm_print(2, "Vmgrid_addGrid:  Not adding grid!\n");
        return 0;
    }

    thee->grids[thee->ngrids] = grid;
    (thee->ngrids)++;

    return 1;

}
