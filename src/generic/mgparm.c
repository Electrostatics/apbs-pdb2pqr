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

    int i;

    /* Check to see if we were even filled... */
    if (!thee->parsed) {
        Vnm_print(2, "MGparm_check:  not filled!\n");
        return 0;
    }
    /* Check to see that everything was set */
    if (!thee->setdime) {
        Vnm_print(2, "NOsh: DIME not set!\n");
        return 0;
    }
    if (!thee->setnlev) {
        Vnm_print(2, "NOsh: NLEV not set!\n");
        return 0;
    }
    if ((!thee->setgrid) && (!thee->setglen)) {
        Vnm_print(2, "NOsh: Neither GRID nor GLEN set!\n");
        return 0;
    }
    if ((thee->setgrid) && (thee->setglen)) {
        Vnm_print(2, "NOsh: Both GRID and GLEN set!\n");
        return 0;
    }
    if (!thee->setgcent) {
        Vnm_print(2, "NOsh: GCENT not set!\n");
        return 0;
    }
    if (!thee->setmolid) {
        Vnm_print(2, "NOsh: MOL not set!\n");
        return 0;
    }
    if (!thee->setnonlin) {
        Vnm_print(2, "NOsh: LPBE or NPBE not set!\n");
        return 0;
    }
    if (!thee->setbcfl) {
        Vnm_print(2, "NOsh: BCFL not set!\n");
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

    return 0;
}
