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
/// Copyright � 1999. The Regents of the University of California (Regents).
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
// File:     femparm.c
//
// Purpose:  Class FEMparm: methods. 
//
// Author:   Nathan Baker
/////////////////////////////////////////////////////////////////////////// */

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
