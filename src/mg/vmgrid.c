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
#include "apbs/vmgrid.h"
#include "apbs/vstring.h"

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
