/**
 *  @file    femparm.c
 *  @ingroup FEMparm
 *  @author  Nathan Baker
 *  @brief   Class FEMparm methods
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
 * Copyright (c) 2003.  Washington University in St. Louis.
 * All Rights Reserved.
 * Portions Copyright (c) 1999-2003.  The Regents of the University of
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
#include "apbs/femparm.h"

/* ///////////////////////////////////////////////////////////////////////////
// Class FEMparm: Inlineable methods
/////////////////////////////////////////////////////////////////////////// */
#if !defined(VINLINE_MGPARM)

#endif /* if !defined(VINLINE_MGPARM) */

/* ///////////////////////////////////////////////////////////////////////////
// Class FEMparm: Non-inlineable methods
/////////////////////////////////////////////////////////////////////////// */

/* ///////////////////////////////////////////////////////////////////////////
// Routine:  FEMparm_ctor
//
// Author: Nathan Baker
/////////////////////////////////////////////////////////////////////////// */
VPUBLIC FEMparm* FEMparm_ctor() {

    /* Set up the structure */
    FEMparm *thee = VNULL;
    thee = Vmem_malloc(VNULL, 1, sizeof(FEMparm));
    VASSERT( thee != VNULL);
    VASSERT( FEMparm_ctor2(thee) );

    return thee;
}

/* ///////////////////////////////////////////////////////////////////////////
// Routine:  FEMparm_ctor2
//
// Author:   Nathan Baker
/////////////////////////////////////////////////////////////////////////// */
VPUBLIC int FEMparm_ctor2(FEMparm *thee) {

    if (thee == VNULL) return 0;

    thee->parsed = 0;

    return 1; 
}

/* ///////////////////////////////////////////////////////////////////////////
// Routine:  FEMparm_dtor
//
// Author:   Nathan Baker
/////////////////////////////////////////////////////////////////////////// */
VPUBLIC void FEMparm_dtor(FEMparm **thee) {
    if ((*thee) != VNULL) {
        FEMparm_dtor2(*thee);
        Vmem_free(VNULL, 1, sizeof(FEMparm), (void **)thee);
        (*thee) = VNULL;
    }
}

/* ///////////////////////////////////////////////////////////////////////////
// Routine:  FEMparm_dtor2
//
// Author:   Nathan Baker
/////////////////////////////////////////////////////////////////////////// */
VPUBLIC void FEMparm_dtor2(FEMparm *thee) { ; }

/* ///////////////////////////////////////////////////////////////////////////
// Routine:  FEMparm_check
//
// Author:   Nathan Baker
/////////////////////////////////////////////////////////////////////////// */
VPUBLIC int FEMparm_check(FEMparm *thee) { 

    /* Check to see if we were even filled... */
    if (!thee->parsed) {
        Vnm_print(2, "FEMparm_check:  not filled!\n");
        return 0;
    }

    return 1;

}
