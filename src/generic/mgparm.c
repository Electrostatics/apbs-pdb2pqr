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
// File:     mgparm.c
//
// Purpose:  Class MGparm: methods. 
//
// Author:   Nathan Baker
/////////////////////////////////////////////////////////////////////////// */

#include "apbscfg.h"
#include "apbs/mgparm.h"

/* ///////////////////////////////////////////////////////////////////////////
// Class MGparm: Inlineable methods
/////////////////////////////////////////////////////////////////////////// */
#if !defined(VINLINE_MGPARM)

#endif /* if !defined(VINLINE_MGPARM) */

/* ///////////////////////////////////////////////////////////////////////////
// Class MGparm: Non-inlineable methods
/////////////////////////////////////////////////////////////////////////// */

/* ///////////////////////////////////////////////////////////////////////////
// Routine:  MGparm_ctor
//
// Author: Nathan Baker
/////////////////////////////////////////////////////////////////////////// */
VPUBLIC MGparm* MGparm_ctor() {

    /* Set up the structure */
    MGparm *thee = VNULL;
    thee = Vmem_malloc(VNULL, 1, sizeof(MGparm));
    VASSERT( thee != VNULL);
    VASSERT( MGparm_ctor2(thee) );

    return thee;
}

/* ///////////////////////////////////////////////////////////////////////////
// Routine:  MGparm_ctor2
//
// Purpose:  Construct the MGparm object
//
// Notes:    Constructor broken into two parts for FORTRAN users.
//
// Returns:  1 if sucessful, 0 otherwise
//
// Author:   Nathan Baker
/////////////////////////////////////////////////////////////////////////// */
VPUBLIC int MGparm_ctor2(MGparm *thee) {

    int i;

    if (thee == VNULL) return 0;

    thee->nlev = -1;
    for (i=0; i<3; i++) thee->dime[i] = -1;

    thee->parsed = 0;

    thee->setdime = 0;
    thee->setnlev = 0;
    thee->setgrid = 0;
    thee->setglen = 0;
    thee->setgcent = 0;  
    thee->setmolid = 0;
    thee->setnonlin = 0;
    thee->setbcfl = 0;
    thee->setnion = 0;
    for (i=0; i<MAXION; i++) thee->setion[i] = 0;
    thee->setpdie = 0;
    thee->setsdie = 0;
    thee->setsrfm = 0;
    thee->setsrad = 0;
    thee->setswin = 0; 
    thee->settemp = 0;
    thee->setgamma = 0;
    thee->setcalcenergy = 0;      
    thee->setcalcforce = 0;       
    thee->setwritepot = 0; 
    thee->setwriteacc = 0; 
    thee->nion = 0;
    thee->swin = 0;
    thee->srad = 1.4;

    return 1; 
}

/* ///////////////////////////////////////////////////////////////////////////
// Routine:  MGparm_dtor
//
// Author:   Nathan Baker
/////////////////////////////////////////////////////////////////////////// */
VPUBLIC void MGparm_dtor(MGparm **thee) {
    if ((*thee) != VNULL) {
        MGparm_dtor2(*thee);
        Vmem_free(VNULL, 1, sizeof(MGparm), (void **)thee);
        (*thee) = VNULL;
    }
}

/* ///////////////////////////////////////////////////////////////////////////
// Routine:  MGparm_dtor2
//
// Purpose:  Destroy the atom object
//
// Author:   Nathan Baker
/////////////////////////////////////////////////////////////////////////// */
VPUBLIC void MGparm_dtor2(MGparm *thee) { ; }

/* ///////////////////////////////////////////////////////////////////////////
// Routine:  MGparm_check
//
// Purpose:  Check the parameter settings for internal consistency
//
// Returns:  1 if OK, 0 otherwise
//
// Author:   Nathan Baker
/////////////////////////////////////////////////////////////////////////// */
VPUBLIC int MGparm_check(MGparm *thee) { 

    int i, tdime[3], ti, tnlev[3], nlev;

    /* Check to see if we were even filled... */
    if (!thee->parsed) {
        Vnm_print(2, "MGparm_check:  not filled!\n");
        return 0;
    }
    /* Check to see that everything was set */
    if (!thee->setdime) {
        Vnm_print(2, "NOsh:  DIME not set!\n");
        return 0;
    }
    if (!thee->setnlev) {
        Vnm_print(2, "NOsh:  NLEV not set!\n");
        return 0;
    }
    if ((!thee->setgrid) && (!thee->setglen)) {
        Vnm_print(2, "NOsh:  Neither GRID nor GLEN set!\n");
        return 0;
    }
    if ((thee->setgrid) && (thee->setglen)) {
        Vnm_print(2, "NOsh:  Both GRID and GLEN set!\n");
        return 0;
    }
    if (!thee->setgcent) {
        Vnm_print(2, "NOsh:  GCENT not set!\n");
        return 0;
    }
    if (!thee->setmolid) {
        Vnm_print(2, "NOsh:  MOL not set!\n");
        return 0;
    }
    if (!thee->setnonlin) {
        Vnm_print(2, "NOsh:  LPBE or NPBE not set!\n");
        return 0;
    }
    if (!thee->setbcfl) {
        Vnm_print(2, "NOsh:  BCFL not set!\n");
        return 0;
    }
    if (!thee->setnion) {
        thee->setnion = 1;
        thee->nion = 0;
    } 
    for (i=0; i<thee->nion; i++) {
        if (!thee->setion[i]) {
            Vnm_print(2, "NOsh: ION #%d not set!\n",i);
            return 0;
        }
    }
    if (!thee->setpdie) {
        Vnm_print(2, "NOsh: PDIE not set!\n");
        return 0;
    }
    if (!thee->setsdie) {
        Vnm_print(2, "NOsh: SDIE not set!\n");
        return 0;
    }
    if (!thee->setsrfm) {
        Vnm_print(2, "NOsh: SRFM not set!\n");
        return 0;
    }
    if (((thee->srfm==0) || (thee->srfm==1)) && (!thee->setsrad)) {
        Vnm_print(2, "NOsh: SRAD not set!\n");
        return 0;
    }
    if ((thee->srfm==2) && (!thee->setswin)) {
        Vnm_print(2, "NOsh: SWIN not set!\n");
        return 0;
    }
    if (!thee->settemp) {
        Vnm_print(2, "NOsh: TEMP not set!\n");
        return 0;
    }
    if (!thee->setgamma) {
        Vnm_print(2, "NOsh: GAMMA not set!\n");
        return 0;
    }
    if (!thee->setcalcenergy) thee->calcenergy = 0;
    if (!thee->setcalcforce) thee->calcforce = 0;
    if (!thee->setwritepot) thee->writepot = 0;
    if (!thee->setwriteacc) thee->writeacc = 0;

    /* We're going to perform a sanity check on nlev and dime, but only in the
     * case where nlev and dime have been set to non-negative values.  This
     * allows us to ignore phony settings as established temporarily during
     * automatic focusing. */
    if ((thee->nlev >= 0) && (thee->dime[0] >= 0) &&
        (thee->dime[1] >= 0) && (thee->dime[2] >= 0)) {
	/* Calculate the actual number of grid points and nlev to satisfy the
	 * formula:  n = c * 2^(l+1) + 1, where n is the number of grid points,
	 * c is an integer, and l is the number of levels */
        for (i=0; i<3; i++) {
            tdime[i] = thee->dime[i];
            ti = tdime[i] - 1;
            tnlev[i] = 0;
            /* Find the maximum number of times this dimension can be divided by
             * two */
            while (1) {
                if (VODD(ti)) break;
                (tnlev[i])++;
                ti = (int)ceil(0.5*ti);
            }
            (tnlev[i])--;
            /* We'd like to have at least 4 levels in the multigrid hierarchy.
	     * This means that the dimension needs to be c*32 + 1, where c is
	     * an integer. */
            if (tnlev[i] < 4) {
                Vnm_print(2, "NOsh:  Bad dime[%d]  = %d (%d nlev)!\n",
                  i, tdime[i], tnlev[i]);
                ti = (int)(tdime[i]/32.0);
                if (ti < 1) ti = 1;
                tdime[i] = ti*32 + 1;
                tnlev[i] = 4;
                Vnm_print(2, "NOsh:  Reset dime[%d] to %d and nlev to 4.\n", 
                  i, tdime[i]);
            }
        }
	/* The actual number of levels we'll be using is the smallest number of
         * possible levels in any dimensions */
        nlev = VMIN2(tnlev[0], tnlev[1]);
        nlev = VMIN2(nlev, tnlev[2]);
        /* Set the number of levels and dimensions */
        Vnm_print(0, "NOsh:  nlev = %d, dime = (%d, %d, %d)\n", nlev, tdime[0],
          tdime[1], tdime[2]);
        thee->nlev = nlev;
        for (i=0; i<3; i++) thee->dime[i] = tdime[i];
    }

    return 1;
}

/* ///////////////////////////////////////////////////////////////////////////
// Routine:  MGparm_copy
//  
// Purpose:  Copy parm into thee
//
// Args:     parm    object to copy into thee
//
// Author:   Nathan Baker
/////////////////////////////////////////////////////////////////////////// */
VPUBLIC void MGparm_copy(MGparm *thee, MGparm *parm) {

    int i;

    VASSERT(thee != VNULL);
    VASSERT(parm != VNULL);

    for (i=0; i<3; i++) thee->dime[i] = parm->dime[i];
    thee->setdime = parm->setdime;
    thee->nlev = parm->nlev;
    thee->setnlev = parm->setnlev;
    for (i=0; i<3; i++) thee->grid[i] = parm->grid[i];
    thee->setgrid = parm->setgrid;
    for (i=0; i<3; i++) thee->glen[i] = parm->glen[i];
    thee->setglen = parm->setglen;
    thee->cmeth = parm->cmeth;
    for (i=0; i<3; i++) thee->center[i] = parm->center[i];
    thee->setgcent = parm->setgcent;
    thee->centmol = parm->centmol;
    thee->molid = parm->molid;
    thee->setmolid = parm->setmolid;
    thee->nonlin = parm->nonlin; 
    thee->setnonlin = parm->setnonlin;
    thee->bcfl = parm->bcfl;
    thee->setbcfl = parm->setbcfl;
    thee->nion = parm->nion;
    thee->setnion = parm->setnion;
    for (i=0; i<MAXION; i++) {
        thee->ionq[i] = parm->ionq[i];
        thee->ionc[i] = parm->ionc[i];
        thee->ionr[i] = parm->ionr[i];
        thee->setion[i] = parm->setion[i];
    };
    thee->pdie = parm->pdie;
    thee->setpdie = parm->setpdie;
    thee->sdie = parm->sdie;
    thee->setsdie = parm->setsdie;
    thee->srfm = parm->srfm;
    thee->setsrfm = parm->setsrfm;
    thee->srad = parm->srad;
    thee->setsrad = parm->setsrad;
    thee->swin = parm->swin;
    thee->setswin = parm->setswin;
    thee->temp = parm->temp;
    thee->settemp = parm->settemp;
    thee->gamma = parm->gamma;
    thee->setgamma = parm->setgamma;
    thee->calcenergy = parm->calcenergy;
    thee->setcalcenergy = parm->setcalcenergy;
    thee->calcforce = parm->calcforce;
    thee->setcalcforce = parm->setcalcforce;
    thee->writepot = parm->writepot;
    thee->setwritepot = parm->setwritepot;
    for (i=0; i<VMAX_ARGLEN; i++) thee->writepotstem[i] = parm->writepotstem[i];
    thee->writepotfmt = parm->writepotfmt;
    thee->writeacc = parm->writeacc;
    thee->setwriteacc = parm->setwriteacc;
    for (i=0; i<VMAX_ARGLEN; i++) thee->writeaccstem[i] = parm->writeaccstem[i];
    thee->writeaccfmt = parm->writeaccfmt;
    thee->parsed = parm->parsed;

}
