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
// File:     vatom.c
//
// Purpose:  Class Vatom: methods.
//
// Author:   Nathan Baker
/////////////////////////////////////////////////////////////////////////// */

#include "apbscfg.h"
#include "apbs/vatom.h"

VEMBED(rcsid="$Id$")

/* ///////////////////////////////////////////////////////////////////////////
// Class Vatom: Inlineable methods
/////////////////////////////////////////////////////////////////////////// */
#if !defined(VINLINE_VATOM)
/* ///////////////////////////////////////////////////////////////////////////
// Routine:  Vatom_getPosition
//
// Purpose:  Return atom position (as 3-array)
//
// Author:   Nathan Baker
/////////////////////////////////////////////////////////////////////////// */
VPUBLIC double *Vatom_getPosition(Vatom *thee) { 

   VASSERT(thee != VNULL);
   return thee->position; 

}

/* ///////////////////////////////////////////////////////////////////////////
// Routine:  Vatom_setRadius
//
// Purpose:  Set atom radius 
//
// Author:   Nathan Baker
/////////////////////////////////////////////////////////////////////////// */
VPUBLIC void Vatom_setRadius(Vatom *thee, double radius) { 

   VASSERT(thee != VNULL);
   thee->radius = radius; 

}

/* ///////////////////////////////////////////////////////////////////////////
// Routine:  Vatom_getRadius
//
// Purpose:  Return atom radius 
//
// Author:   Nathan Baker
/////////////////////////////////////////////////////////////////////////// */
VPUBLIC double Vatom_getRadius(Vatom *thee) { 

   VASSERT(thee != VNULL);
   return thee->radius; 

}

/* ///////////////////////////////////////////////////////////////////////////
// Routine:  Vatom_setCharge
//
// Purpose:  Set atom charge 
//
// Author:   Nathan Baker
/////////////////////////////////////////////////////////////////////////// */
VPUBLIC void Vatom_setCharge(Vatom *thee, double charge) { 

   VASSERT(thee != VNULL);
   thee->charge = charge; 

}

/* ///////////////////////////////////////////////////////////////////////////
// Routine:  Vatom_getCharge
//
// Purpose:  Return atom charge 
//
// Author:   Nathan Baker
/////////////////////////////////////////////////////////////////////////// */
VPUBLIC double Vatom_getCharge(Vatom *thee) { 

   VASSERT(thee != VNULL);
   return thee->charge; 

}

/* ///////////////////////////////////////////////////////////////////////////
// Routine:  Vatom_memChk
//
// Purpose:  Return the number of bytes used by this object
//
// Author:   Nathan Baker
/////////////////////////////////////////////////////////////////////////// */
VPUBLIC int Vatom_memChk(Vatom *thee) { return sizeof(Vatom); }

#endif /* if !defined(VINLINE_VATOM) */

/* ///////////////////////////////////////////////////////////////////////////
// Class Vatom: Non-inlineable methods
/////////////////////////////////////////////////////////////////////////// */

/* ///////////////////////////////////////////////////////////////////////////
// Routine:  Vatom_ctor
//
// Purpose:  Construct the atom object
//
// Author:   Nathan Baker
/////////////////////////////////////////////////////////////////////////// */
VPUBLIC Vatom* Vatom_ctor() {

    /* Set up the structure */
    Vatom *thee = VNULL;
    thee = (Vatom *)Vmem_malloc( VNULL, 1, sizeof(Vatom) );
    VASSERT( thee != VNULL);
    VASSERT( Vatom_ctor2(thee));

    return thee;
}

/* ///////////////////////////////////////////////////////////////////////////
// Routine:  Vatom_ctor2
//
// Purpose:  Construct the atom object
//
// Notes:    Constructor broken into two parts for FORTRAN users.
//
// Author:   Nathan Baker
/////////////////////////////////////////////////////////////////////////// */
VPUBLIC int Vatom_ctor2(Vatom *thee) { 
    thee->partID = -1;
    return 1; 
}

/* ///////////////////////////////////////////////////////////////////////////
// Routine:  Vatom_dtor
//
// Purpose:  Destroy the atom object
//
// Author:   Nathan Baker
/////////////////////////////////////////////////////////////////////////// */
VPUBLIC void Vatom_dtor(Vatom **thee) {
    if ((*thee) != VNULL) {
        Vatom_dtor2(*thee);
        Vmem_free(VNULL, 1, sizeof(Vatom), (void **)thee);
        (*thee) = VNULL;
    }
}

/* ///////////////////////////////////////////////////////////////////////////
// Routine:  Vatom_dtor2
//
// Purpose:  Destroy the atom object
//
// Author:   Nathan Baker
/////////////////////////////////////////////////////////////////////////// */
VPUBLIC void Vatom_dtor2(Vatom *thee) { ; }

/* ///////////////////////////////////////////////////////////////////////////
// Routine:  Vatom_setPosition
//
// Purpose:  Set atom position (as 3-array)
//
// Author:   Nathan Baker
/////////////////////////////////////////////////////////////////////////// */
VPUBLIC void Vatom_setPosition(Vatom *thee, double position[3]) { 

   VASSERT(thee != VNULL);
   (thee->position)[0] = position[0]; 
   (thee->position)[1] = position[1]; 
   (thee->position)[2] = position[2]; 

}

