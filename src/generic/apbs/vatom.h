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
#include "mc/vmem.h"
#include "mc/vec3.h"

#include "apbs/vhal.h"

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
    VEXTERNC double* Vatom_getPosition(Vatom *thee);
    VEXTERNC void    Vatom_setRadius(Vatom *thee, double radius);
    VEXTERNC double  Vatom_getRadius(Vatom *thee);
    VEXTERNC void    Vatom_setCharge(Vatom *thee, double charge);
    VEXTERNC double  Vatom_getCharge(Vatom *thee);
    VEXTERNC int    Vatom_memChk(Vatom *thee);
#else /* if defined(VINLINE_VATOM) */
#   define Vatom_getPosition(thee) ((thee)->position)
#   define Vatom_setRadius(thee, radius) ((thee)->radius = radius)
#   define Vatom_getRadius(thee) ((thee)->radius)
#   define Vatom_setCharge(thee, charge) ((thee)->charge = charge)
#   define Vatom_getCharge(thee) ((thee)->charge)
#   define Vatom_memChk(thee) (sizeof(Vatom))
#endif /* if !defined(VINLINE_VATOM) */

/* ///////////////////////////////////////////////////////////////////////////
// Class Vatom: Non-Inlineable methods (vatom.c)
/////////////////////////////////////////////////////////////////////////// */

VEXTERNC Vatom * Vatom_ctor();
VEXTERNC int     Vatom_ctor2(Vatom *thee);
VEXTERNC void    Vatom_dtor(Vatom **thee);
VEXTERNC void    Vatom_dtor2(Vatom *thee);

VEXTERNC void   Vatom_setPosition(Vatom *thee, Vec3 position);

#endif /* ifndef _VATOM_H_ */
