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
 * Copyright (c) 2002-2009, Washington University in St. Louis.
 * Portions Copyright (c) 2002-2009.  Nathan A. Baker
 * Portions Copyright (c) 1999-2002.  The Regents of the University of California.
 * Portions Copyright (c) 1995.  Michael Holst
 *
 * All rights reserved.
 * 
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are met: 
 *
 * -  Redistributions of source code must retain the above copyright notice, this
 * list of conditions and the following disclaimer.  
 * 
 * - Redistributions in binary form must reproduce the above copyright notice,
 * this list of conditions and the following disclaimer in the documentation
 * and/or other materials provided with the distribution.
 * 
 * - Neither the name of Washington University in St. Louis nor the names of its
 * contributors may be used to endorse or promote products derived from this
 * software without specific prior written permission.
 * 
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
 * "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
 * LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
 * A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR
 * CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
 * EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
 * PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
 * PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
 * LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
 * NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
 * SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 *
 * @endverbatim
 */

#include "apbscfg.h"
#include "apbs/vpmgp.h"
#include "apbs/mgparm.h"

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
VPUBLIC Vpmgp* Vpmgp_ctor(MGparm *mgparm) {

    Vpmgp *thee = VNULL;

    /* Set up the structure */
    thee = Vmem_malloc(VNULL, 1, sizeof(Vpmgp) );
    VASSERT( thee != VNULL);
    VASSERT(Vpmgp_ctor2(thee,mgparm));

    return thee;
}

/* ///////////////////////////////////////////////////////////////////////////
// Routine:  Vpmgp_ctor2
//
// Author:   Nathan Baker
/////////////////////////////////////////////////////////////////////////// */
VPUBLIC int Vpmgp_ctor2(Vpmgp *thee,MGparm *mgparm) {
	
	/* Specified parameters */
    thee->nx = mgparm->dime[0];
    thee->ny = mgparm->dime[1];
    thee->nz = mgparm->dime[2];
    thee->hx = mgparm->grid[0];
    thee->hy = mgparm->grid[1];
    thee->hzed = mgparm->grid[2];
    thee->xlen = ((double)(mgparm->dime[0]-1))*mgparm->grid[0];
    thee->ylen = ((double)(mgparm->dime[1]-1))*mgparm->grid[1];
    thee->zlen = ((double)(mgparm->dime[2]-1))*mgparm->grid[2];
    thee->nlev = mgparm->nlev;
	
    thee->nonlin = mgparm->nonlintype;
	thee->meth = mgparm->method;
	
#ifdef DEBUG_MAC_OSX_OCL
#include "mach_chud.h"
	if(kOpenCLAvailable)
		thee->meth = 4;
#endif
	
    if (thee->nonlin == NONLIN_LPBE) thee->ipkey = IPKEY_LPBE; /* LPBE case */
	else if(thee->nonlin == NONLIN_SMPBE) thee->ipkey = IPKEY_SMPBE; /* SMPBE case */
    else thee->ipkey = IPKEY_NPBE; /* NPBE standard case */
	
    /* Default parameters */
    if (mgparm->seterrtol) { /* If errtol is set by the user in APBS input file, */
                             /* then use this custom-defined errtol */
        thee->errtol = mgparm->errtol; 
        Vnm_print(1, "  Error tolerance (errtol) is now set to user-defined \
value: %g \n", thee->errtol);
        Vnm_print(0, "Error tolerance (errtol) is now set to user-defined \
value: %g \n", thee->errtol);
    } else thee->errtol = 1.0e-6;   /* Here are a few comments.  Mike had this set to
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
	thee->mgcoar = 2;
	thee->mgkey = 0;
	thee->nu1 = 2;
	thee->nu2 = 2;
	thee->mgprol = 0;
	thee->mgdisc = 0;
	thee->omegal = 19.4e-1;
	thee->omegan = 9.0e-1;
	thee->ipcon = 3;
	thee->irite = 8;
	thee->xcent = 0.0;
	thee->ycent = 0.0;
	thee->zcent = 0.0;
	
	/* Default value for all APBS runs */
	thee->mgsmoo = 1;
    if (thee->nonlin == NONLIN_NPBE || thee->nonlin == NONLIN_SMPBE) { 
		/* SMPBE Added - SMPBE needs to mimic NPBE */
		Vnm_print(0, "Vpmp_ctor2:  Using meth = 1, mgsolv = 0\n");
        thee->mgsolv = 0;
	}else{        
		/* Most rigorous (good for testing) */
		Vnm_print(0, "Vpmp_ctor2:  Using meth = 2, mgsolv = 1\n");
        thee->mgsolv = 1;
    }

	/* TEMPORARY USEAQUA */
	/* If we are using aqua, our solution method is either VSOL_CGMGAqua or VSOL_NewtonAqua
	 * so we need to temporarily override the mgsolve value and set it to 0
	 */
	if(mgparm->useAqua == 1) thee->mgsolv = 0;
	
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
