/* ///////////////////////////////////////////////////////////////////////////
/// MC = < Manifold Code >
/// Copyright (C) 1993  Michael Holst
/// Copyright (C) 1999  Contributing authors listed in the code documentation
///                     retain the copyrights to their code contributions.
///
/// The program MC is copyrighted software, and there are restrictions on its
/// use, modification, and distribution; please refer to the document COPYING
/// which accompanies this software.
///
/// This program is distributed in the hope that it will be useful,
/// but WITHOUT ANY WARRANTY; without even the implied warranty of
/// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
/// See the document COPYING for more details.
///
/// You should have received a copy of the COPYING file with this program;
/// if not, write to the author Michael Holst at:  mholst@math.ucsd.edu
///
/// rcsid="$Id$"
/// /////////////////////////////////////////////////////////////////////// */

/* ///////////////////////////////////////////////////////////////////////////
// File:     vatom.c
//
// Purpose:  Class Vatom: methods.
//
// Author:   Nathan Baker
/////////////////////////////////////////////////////////////////////////// */

#include "apbs/vatom.h"

VEMBED(rcsid="$Id$")

/* ///////////////////////////////////////////////////////////////////////////
// Class Vatom: Inlineable methods
/////////////////////////////////////////////////////////////////////////// */
#if !defined(VINLINE_VATOM)
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
    thee = (Vatom *)Vram_ctor( 1, sizeof(Vatom) );
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
VPUBLIC int Vatom_ctor2(Vatom *thee) { return 1; }

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
        Vram_dtor((Vram **)thee, 1, sizeof(Vatom));
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
// Purpose:  Set atom position (as Vec3 vector)
//
// Author:   Nathan Baker
/////////////////////////////////////////////////////////////////////////// */
VPUBLIC void Vatom_setPosition(Vatom *thee, Vec3 position) { 

   VASSERT(thee != VNULL);
   (thee->position)[0] = position[0]; 
   (thee->position)[1] = position[1]; 
   (thee->position)[2] = position[2]; 

}

/* ///////////////////////////////////////////////////////////////////////////
// Routine:  Vatom_getPosition
//
// Purpose:  Return atom position (as Vec3 vector)
//
// Author:   Nathan Baker
/////////////////////////////////////////////////////////////////////////// */
VPUBLIC double * Vatom_getPosition(Vatom *thee) { 

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
