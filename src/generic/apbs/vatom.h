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
// File:     vatom.h    < vatom.c >
//
// Purpose:  
//    Class Vatom:  
//      An atom class for interfacing MC with PDB files and  
//      other things in an efficient manner.  Basically nice
//      temporary storage for stuff read from a PDB while  
//      assigning charges and stuff from a parameter file.
//
// Author:   Nathan Baker
/////////////////////////////////////////////////////////////////////////// */

#ifndef _VATOM_H_
#define _VATOM_H_

#include "mc/vhal.h"
#include "mc/vram.h"
#include "mc/vec3.h"

/* ///////////////////////////////////////////////////////////////////////////
// Class Vatom: Parameters and datatypes
/////////////////////////////////////////////////////////////////////////// */

/* ///////////////////////////////////////////////////////////////////////////
// Structure Vatom: Definition
/////////////////////////////////////////////////////////////////////////// */

typedef struct Vatom {

    Vec3 position;     /* Atomic position */
    double radius;     /* Atomic radius   */
    double charge;     /* Atomic charge   */

} Vatom;

/* ///////////////////////////////////////////////////////////////////////////
// Class Vatom: Inlineable methods (vatom.c)
/////////////////////////////////////////////////////////////////////////// */

#if !defined(VINLINE_VATOM)
#else /* if defined(VINLINE_VATOM) */
#endif /* if !defined(VINLINE_VATOM) */

/* ///////////////////////////////////////////////////////////////////////////
// Class Vatom: Non-Inlineable methods (vatom.c)
/////////////////////////////////////////////////////////////////////////// */

VEXTERNC Vatom * Vatom_ctor();
VEXTERNC int     Vatom_ctor2(Vatom *thee);
VEXTERNC void    Vatom_dtor(Vatom **thee);
VEXTERNC void    Vatom_dtor2(Vatom *thee);

VEXTERNC void   Vatom_setPosition(Vatom *thee, Vec3 position);
VEXTERNC double* Vatom_getPosition(Vatom *thee);
VEXTERNC void     Vatom_setRadius(Vatom *thee, double radius);
VEXTERNC double   Vatom_getRadius(Vatom *thee);
VEXTERNC void   Vatom_setCharge(Vatom *thee, double charge);
VEXTERNC double   Vatom_getCharge(Vatom *thee);

#endif /* ifndef _VATOM_H_ */
