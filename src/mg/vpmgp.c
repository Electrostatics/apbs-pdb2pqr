/**
 *  @file    vpmgp.c
 *  @author  Nathan Baker
 *  @brief   Class Vpmgp methods
 *  @ingroup Vpmgp
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
#include "apbs/vpmgp.h"

VEMBED(rcsid="$Id$")

/* ///////////////////////////////////////////////////////////////////////////
// Class Vpmgp: Inlineable methods
/////////////////////////////////////////////////////////////////////////// */
#if !defined(VINLINE_VACC)
#endif /* if !defined(VINLINE_VACC) */

/* ///////////////////////////////////////////////////////////////////////////
// Class Vpmgp: Non-inlineable methods
/////////////////////////////////////////////////////////////////////////// */

/* ///////////////////////////////////////////////////////////////////////////
// Routine:  Vpmgp_ctor
//
// Author:   Nathan Baker
/////////////////////////////////////////////////////////////////////////// */
VPUBLIC Vpmgp* Vpmgp_ctor(int nx, int ny, int nz, int nlev, double hx, 
  double hy, double hzed, int nonlin) {

    Vpmgp *thee = VNULL;

    /* Set up the structure */
    thee = Vmem_malloc(VNULL, 1, sizeof(Vpmgp) );
    VASSERT( thee != VNULL);
    VASSERT(Vpmgp_ctor2(thee, nx, ny, nz, nlev, hx, hy, hzed, nonlin));

    return thee;
}

/* ///////////////////////////////////////////////////////////////////////////
// Routine:  Vpmgp_ctor2
//
// Author:   Nathan Baker
/////////////////////////////////////////////////////////////////////////// */
VPUBLIC int Vpmgp_ctor2(Vpmgp *thee, int nx, int ny, int nz, int nlev,
  double hx, double hy, double hzed, int nonlin) {

    /* Specified parameters */
    thee->nx = nx;
    thee->ny = ny;
    thee->nz = nz;
    thee->hx = hx;
    thee->hy = hy;
    thee->hzed = hzed;
    thee->xlen = ((double)(nx-1))*hx;
    thee->ylen = ((double)(ny-1))*hy;
    thee->zlen = ((double)(nz-1))*hzed;
    thee->nlev = nlev; 
    thee->nonlin = nonlin;
    if (nonlin == 0) thee->ipkey = -1;
    else thee->ipkey = 0;


    /* Default parameters */
    thee->errtol = 1.0e-6;   /* Here are a few comments.  Mike had this set to
                              * 1e-9; convential wisdom sets this at 1e-6 for
                              * the PBE; Ray Luo sets this at 1e-3 for his
                              * accelerated PBE (for dynamics, etc.) */
    thee->itmax = 200;
    thee->istop = 1;
    thee->iinfo = 1;         /* I'd recommend either 1 (for debugging LPBE) or 
                              * 2 (for debugging NPBE), higher values give 
                              * too much output */
    thee->bcfl = BCFL_SDH;
    thee->key = 0;
    thee->iperf = 0;
    if (thee->nonlin == 1) { 
        Vnm_print(0, "Vpmp_ctor2:  Using meth = 1, mgcoar = 2, mgsolv = 0\n");
        thee->meth = 1;
        thee->mgcoar = 2;
        thee->mgsolv = 0;
    } else {                 
        Vnm_print(0, "Vpmp_ctor2:  Using meth = 0, mgcoar = 2, mgsolv = 0\n");
#if 0                               /* Fastest convergence, but lots of mem */
        thee->meth = 0;
        thee->mgcoar = 2;
        thee->mgsolv = 0;
#else                               /* Most rigorous (good for testing) */
        thee->meth = 2;
        thee->mgcoar = 2;
        thee->mgsolv = 1;
#endif
    }
    thee->mgkey = 0;
    thee->nu1 = 2;
    thee->nu2 = 2;
    thee->mgsmoo = 1;
    thee->mgprol = 0;
    thee->mgdisc = 0;
    thee->omegal = 8.0e-1;
    thee->omegan = 9.0e-1;
    thee->ipcon = 3;
    thee->irite = 8;
    thee->xcent = 0.0;
    thee->ycent = 0.0;
    thee->zcent = 0.0;

    return 1;
}

/* ///////////////////////////////////////////////////////////////////////////
// Routine:  Vpmgp_dtor
//
// Author:   Nathan Baker
/////////////////////////////////////////////////////////////////////////// */
VPUBLIC void Vpmgp_dtor(Vpmgp **thee) {
    
    if ((*thee) != VNULL) {
        Vpmgp_dtor2(*thee);
        Vmem_free(VNULL, 1, sizeof(Vpmgp), (void **)thee);
        (*thee) = VNULL;
    }

}

/* ///////////////////////////////////////////////////////////////////////////
// Routine:  Vpmgp_dtor2
//
// Author:   Nathan Baker
/////////////////////////////////////////////////////////////////////////// */
VPUBLIC void Vpmgp_dtor2(Vpmgp *thee) { ; }
